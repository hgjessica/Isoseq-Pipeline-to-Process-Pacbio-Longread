"""Shared helpers for the Iso-Seq step scripts (config, running commands, sample selection)."""
import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path


class Ctx:
    def __init__(self, cfg, dry_run):
        self.cfg = cfg
        self.dry = dry_run
        self.work = Path(cfg["work_dir"])
        self.threads = str(cfg.get("threads", 8))

    def path(self, *parts):
        return self.work.joinpath(*parts)


def load_config(path):
    with open(path) as fh:
        cfg = json.load(fh)
    for key in ("work_dir", "samples"):
        if key not in cfg:
            sys.exit(f"config is missing required key: {key}")
    return cfg


def run(ctx, cmd, **kw):
    """Run a command (list form, no shell). Returns the exit code."""
    print("$ " + " ".join(shlex.quote(str(c)) for c in cmd), flush=True)
    if ctx.dry:
        return 0
    return subprocess.run([str(c) for c in cmd], **kw).returncode


def need(ctx, *files):
    """Abort with a clear message if an input is missing (skipped in dry run)."""
    for f in files:
        if not ctx.dry and not Path(f).exists():
            sys.exit(f"[ERROR] missing: {f}")


def count_reads(ctx, bam):
    if ctx.dry:
        return "?"
    out = subprocess.run(["samtools", "view", "-c", str(bam)],
                         capture_output=True, text=True)
    return out.stdout.strip() or "?"


def sm_tag(bam):
    hdr = subprocess.run(["samtools", "view", "-H", str(bam)],
                         capture_output=True, text=True).stdout
    for line in hdr.splitlines():
        if line.startswith("@RG"):
            for field in line.split("\t"):
                if field.startswith("SM:"):
                    return field[3:]
    return ""


def kept_samples(cfg):
    return [s for s in cfg["samples"] if s.get("keep", True)]


def pick_sample(cfg, name):
    """--sample NAME, or the SLURM array task id as an index into kept samples."""
    kept = kept_samples(cfg)
    if name:
        for s in kept:
            if s["name"] == name:
                return s
        sys.exit(f"unknown (or keep=false) sample: {name}")
    task = os.environ.get("SLURM_ARRAY_TASK_ID")
    if task is None:
        sys.exit("give --sample NAME (or run as a SLURM array task)")
    return kept[int(task)]


def fl_bam(ctx, sample):
    return Path(sample["fl_bam"]) if "fl_bam" in sample else ctx.path("fl", f"{sample['name']}.fl.bam")


def task_tmp(ctx, tag):
    """Per-task scratch directory so concurrent array tasks cannot collide."""
    base = Path(ctx.cfg.get("scratch_dir") or ctx.path(".tmp"))
    d = base / f"{tag}_{os.environ.get('SLURM_ARRAY_TASK_ID', 'local')}"
    if not ctx.dry:
        d.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(d)
    return d


def ensure_fai(ctx):
    fa = ctx.cfg["reference_fasta"]
    need(ctx, fa)
    if ctx.dry or Path(fa + ".fai").exists():
        return
    if run(ctx, ["samtools", "faidx", fa]) != 0:
        sys.exit("[ERROR] samtools faidx failed")


def base_parser(doc):
    """Parser with the options every step script shares."""
    ap = argparse.ArgumentParser(description=doc.strip().split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="\n".join(doc.strip().split("\n\n")[1:]))
    ap.add_argument("--config", required=True, help="JSON config (see config.example.json)")
    ap.add_argument("--dry-run", action="store_true", help="print commands, do not run")
    return ap


def add_submit(sub, choices):
    """Add the `submit` subcommand (SLURM array, one task per kept sample)."""
    p = sub.add_parser("submit", help="submit a per-sample step as a SLURM array")
    p.add_argument("step_name", choices=sorted(choices))
    p.add_argument("step_args", nargs=argparse.REMAINDER,
                   help="extra args passed to the step, e.g. --input-stage flnc")
    return p


def submit_array(ctx, args, script):
    """Write and submit a SLURM array script (one task per kept sample) for a step."""
    sl = ctx.cfg.get("slurm", {})
    n = len(kept_samples(ctx.cfg))
    conc = sl.get("array_concurrency", 3)
    extra = " ".join(shlex.quote(a) for a in args.step_args)
    log_dir = ctx.path("slurm_logs")
    script = f"""#!/bin/bash
#SBATCH -p {sl.get('partition', 'PARTITION_NAME')}
#SBATCH -N 1
#SBATCH -n {ctx.threads}
#SBATCH --mem={sl.get('mem', '64G')}
#SBATCH -t {sl.get('time', '48:00:00')}
#SBATCH --export=ALL
#SBATCH -J {args.step_name}
#SBATCH --array=0-{n - 1}%{conc}
#SBATCH -o {log_dir}/{args.step_name}_%A_%a.out
#SBATCH -e {log_dir}/{args.step_name}_%A_%a.err

source {sl.get('conda_sh', 'Path_to_conda.sh')}
conda activate {sl.get('conda_env', 'isoseq')}
python3 {Path(script).resolve()} --config {Path(args.config).resolve()} {args.step_name} {extra}
"""
    sbatch_file = ctx.path("slurm", f"{args.step_name}.sbatch")
    print(script)
    if ctx.dry:
        return
    log_dir.mkdir(parents=True, exist_ok=True)
    sbatch_file.parent.mkdir(parents=True, exist_ok=True)
    sbatch_file.write_text(script)
    subprocess.run(["sbatch", str(sbatch_file)])
