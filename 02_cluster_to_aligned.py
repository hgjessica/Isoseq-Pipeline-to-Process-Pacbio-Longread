#!/usr/bin/env python3
"""Part 2: FLNC reads -> clustered transcripts -> aligned BAM.

    python 02_cluster_to_aligned.py --config config.json cluster --sample sample_1
    python 02_cluster_to_aligned.py --config config.json fix-sm [--dir clustered] [--apply]
    python 02_cluster_to_aligned.py --config config.json align --sample sample_1 [--input-stage clustered]
    python 02_cluster_to_aligned.py --config config.json submit cluster
    python 02_cluster_to_aligned.py --config config.json submit align --input-stage flnc

Steps: cluster (isoseq cluster2) -> fix-sm (set SM:/LB: to the sample name, optional) ->
align (pbmm2 --preset ISOSEQ). --input-stage chooses what to align: clustered (default),
flnc (read-level depth, needed for collapse) or fl (before refine, keeps the polyA tail).
Output: clustered/<sample>.clustered.bam, aligned/<sample>.aligned.bam (+ .bai)
"""
import re
import subprocess
import sys
from pathlib import Path

from common import *  # noqa: F401,F403  (helpers shared by all steps)
from common import submit_array

ALIGN_STAGES = {
    # stage: (input dir, input suffix, output dir, output suffix)
    "clustered": ("clustered", ".clustered.bam", "aligned", ".aligned.bam"),
    "flnc": ("flnc", ".flnc.bam", "aligned_flnc", ".flnc.aligned.bam"),
    "fl": ("fl", ".fl.bam", "aligned_fl", ".fl.aligned.bam"),
}


def step_cluster(ctx, args):
    """FLNC -> clustered transcripts (isoseq cluster2)."""
    s = pick_sample(ctx.cfg, args.sample)
    task_tmp(ctx, "cluster")
    src = ctx.path("flnc", f"{s['name']}.flnc.bam")
    out = ctx.path("clustered", f"{s['name']}.clustered.bam")
    need(ctx, src)
    if not ctx.dry:
        out.parent.mkdir(exist_ok=True)
    rc = run(ctx, ["isoseq", "cluster2", "-j", ctx.threads, "--log-level", "INFO", src, out])
    if rc:
        sys.exit(f"[ERROR] cluster2 failed for {s['name']}")
    print(f"[OK] {out}  transcripts: {count_reads(ctx, out)}")


def step_fix_sm(ctx, args):
    """Set SM:/LB: in each BAM header to the file's sample name (BAMs are rewritten)."""
    d = ctx.path(args.dir)
    if not d.is_dir():
        sys.exit(f"no such directory: {d}")
    print("MODE:", "apply" if args.apply else "dry run (pass --apply)")
    for bam in sorted(d.glob("*.bam")):
        sample = bam.name.split(".")[0]
        cur = sm_tag(bam)
        if cur == sample:
            print(f"  {bam.name}: SM already {sample}")
            continue
        print(f"  {bam.name}: SM {cur or '?'} -> {sample}")
        if not args.apply:
            continue
        hdr = subprocess.run(["samtools", "view", "-H", str(bam)],
                             capture_output=True, text=True).stdout
        lines = []
        for line in hdr.splitlines():
            if line.startswith("@RG"):
                line = re.sub(r"\tSM:[^\t]*", f"\tSM:{sample}", line)
                line = re.sub(r"\tLB:[^\t]*", f"\tLB:{sample}", line)
            lines.append(line)
        hdr_file = d / f".hdr_{sample}.sam"
        tmp = Path(str(bam) + ".reheader.tmp")
        hdr_file.write_text("\n".join(lines) + "\n")
        with open(tmp, "wb") as out:
            rc = subprocess.run(["samtools", "reheader", str(hdr_file), str(bam)],
                                stdout=out, stderr=subprocess.DEVNULL).returncode
        hdr_file.unlink()
        if rc != 0 or tmp.stat().st_size == 0:
            tmp.unlink(missing_ok=True)
            print("      FAILED - original left untouched")
            continue
        tmp.replace(bam)
        pbi, bai = Path(str(bam) + ".pbi"), Path(str(bam) + ".bai")
        if pbi.exists():
            pbi.unlink()
            subprocess.run(["pbindex", str(bam)])
        if bai.exists():
            bai.unlink()
            subprocess.run(["samtools", "index", "-@", "4", str(bam)])


def step_align(ctx, args):
    """Align to the reference genome with pbmm2 (--input-stage clustered|flnc|fl)."""
    s = pick_sample(ctx.cfg, args.sample)
    task_tmp(ctx, "align")
    in_dir, in_sfx, out_dir, out_sfx = ALIGN_STAGES[args.input_stage]
    src = fl_bam(ctx, s) if args.input_stage == "fl" else ctx.path(in_dir, s["name"] + in_sfx)
    out = ctx.path(out_dir, s["name"] + out_sfx)
    ensure_fai(ctx)
    need(ctx, src)
    if not ctx.dry:
        out.parent.mkdir(exist_ok=True)
    rc = run(ctx, ["pbmm2", "align", "--preset", "ISOSEQ", "--sort", "-j", ctx.threads,
                   "-J", "8", "--log-level", "INFO", ctx.cfg["reference_fasta"], src, out])
    if rc:
        sys.exit(f"[ERROR] pbmm2 failed for {s['name']}")
    if not ctx.dry and not Path(str(out) + ".bai").exists():
        run(ctx, ["samtools", "index", "-@", "8", out])
    print(f"[OK] {out}  records: {count_reads(ctx, out)}")


def main():
    ap = base_parser(__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("cluster", help="FLNC -> clustered transcripts")
    p.add_argument("--sample", help="sample name (default: SLURM array task id)")
    p.set_defaults(fn=step_cluster)

    p = sub.add_parser("fix-sm", help="set BAM SM:/LB: to the sample name")
    p.add_argument("--dir", default="clustered")
    p.add_argument("--apply", action="store_true")
    p.set_defaults(fn=step_fix_sm)

    p = sub.add_parser("align", help="align with pbmm2")
    p.add_argument("--sample", help="sample name (default: SLURM array task id)")
    p.add_argument("--input-stage", choices=sorted(ALIGN_STAGES), default="clustered")
    p.set_defaults(fn=step_align)

    add_submit(sub, ["cluster", "align"]).set_defaults(fn=lambda c, a: submit_array(c, a, __file__))

    args = ap.parse_args()
    args.fn(Ctx(load_config(args.config), args.dry_run), args)


if __name__ == "__main__":
    main()
