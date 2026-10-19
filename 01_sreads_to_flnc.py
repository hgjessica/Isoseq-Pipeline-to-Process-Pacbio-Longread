#!/usr/bin/env python3
"""Part 1: S-reads -> FL reads -> FLNC reads.

    python 01_sreads_to_flnc.py --config config.json lima
    python 01_sreads_to_flnc.py --config config.json rename [--apply] [--drop-others]
    python 01_sreads_to_flnc.py --config config.json refine --sample sample_1
    python 01_sreads_to_flnc.py --config config.json submit refine

Steps: lima (per cell, demultiplex by barcode) -> rename (barcode -> sample name) ->
refine (primer + polyA trimming, per sample). Skip lima/rename if you already have
demultiplexed FL reads (set "fl_bam" for the sample in the config).
Output: fl/<sample>.fl.bam, flnc/<sample>.flnc.bam
"""
import re
import subprocess
import sys
from pathlib import Path

from common import *  # noqa: F401,F403  (helpers shared by all steps)
from common import submit_array


def step_lima(ctx, args):
    """S-reads -> FL reads, one lima run per cell (demultiplexes by barcode)."""
    primers = ctx.cfg["primers_fasta"]
    need(ctx, primers)
    if not ctx.dry:
        ctx.path("fl").mkdir(exist_ok=True)
    for cell in ctx.cfg["cells"]:
        seg = ctx.path(cell, "segmented.bam")
        tar = ctx.path(f"{cell}.tar")
        if not seg.exists() and tar.exists():
            print(f"[EXTRACT] {tar}")
            if run(ctx, ["tar", "-xf", tar, "-C", ctx.work]) != 0:
                print(f"[ERROR] extract failed for {cell}")
                continue
        if not seg.exists() and not ctx.dry:
            print(f"[SKIP] neither {seg} nor {tar} found")
            continue
        rc = run(ctx, ["lima", "--isoseq", "--peek-guess", "-j", ctx.threads,
                       "--log-level", "INFO", seg, primers,
                       ctx.path("fl", f"{cell}.fl.bam")])
        print(f"[{'OK' if rc == 0 else 'ERROR'}] cell {cell}")
    print("Next: rename --apply (barcode -> sample name)")


def step_rename(ctx, args):
    """Rename lima's barcode-named outputs to <sample>.fl.bam using the config."""
    fl_dir = ctx.path("fl")
    if not fl_dir.is_dir():
        sys.exit(f"no such directory: {fl_dir}")
    by_bc = {s["barcode"]: s["name"] for s in ctx.cfg["samples"] if "barcode" in s}

    if args.apply:
        # lima writes *.build files while indexing; renaming mid-write breaks it.
        try:
            jobs = subprocess.run(["squeue", "-u", os.environ.get("USER", ""), "-h", "-o", "%j"],
                                  capture_output=True, text=True).stdout
            if "lima" in jobs:
                sys.exit("REFUSING: a lima job is still in the queue. Wait, then rerun.")
        except FileNotFoundError:
            pass  # no SLURM here
        if list(fl_dir.glob("*.build")):
            sys.exit("REFUSING: fl/ contains *.build files (index build in progress).")
    print("MODE:", "apply" if args.apply else "dry run (pass --apply)")

    pat = re.compile(r"(IsoSeqX_bc\d+_5p)--IsoSeqX_3p(.*)$")
    n = skipped = 0
    for f in sorted(fl_dir.iterdir()):
        if f.name.endswith(".build"):
            continue
        m = pat.search(f.name)
        if not m:
            continue
        sample = by_bc.get(m.group(1))
        if not sample:
            print(f"  SKIP {f.name} (no mapping for {m.group(1)})")
            skipped += 1
            continue
        new = fl_dir / f"{sample}.fl{m.group(2)}"
        if new.exists():
            print(f"  SKIP {f.name} -> {new.name} (target exists)")
            skipped += 1
            continue
        print(f"  {f.name} -> {new.name}")
        if args.apply:
            f.rename(new)
        n += 1
    print(f"{n} files {'renamed' if args.apply else 'would be renamed'}, {skipped} skipped")

    if args.drop_others:
        print("\nsamples with keep=false:")
        for s in ctx.cfg["samples"]:
            if s.get("keep", True):
                continue
            for f in fl_dir.glob(f"{s['name']}.fl*"):
                print(f"  {'deleting' if args.apply else 'would delete'} {f.name}")
                if args.apply:
                    f.unlink()


def step_refine(ctx, args):
    """FL -> FLNC (primer + polyA trimming)."""
    s = pick_sample(ctx.cfg, args.sample)
    task_tmp(ctx, "refine")
    src = fl_bam(ctx, s)
    out = ctx.path("flnc", f"{s['name']}.flnc.bam")
    need(ctx, src, ctx.cfg["primers_fasta"])
    if not ctx.dry:
        out.parent.mkdir(exist_ok=True)
    rc = run(ctx, ["isoseq", "refine", "-j", ctx.threads, "--require-polya",
                   src, ctx.cfg["primers_fasta"], out])
    if rc:
        sys.exit(f"[ERROR] refine failed for {s['name']}")
    print(f"[OK] {out}  FLNC reads: {count_reads(ctx, out)}")


def main():
    ap = base_parser(__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("lima", help="S-reads -> FL reads, per cell").set_defaults(fn=step_lima)

    p = sub.add_parser("rename", help="barcode -> sample name")
    p.add_argument("--apply", action="store_true")
    p.add_argument("--drop-others", action="store_true",
                   help="delete demultiplexed files of samples with keep=false")
    p.set_defaults(fn=step_rename)

    p = sub.add_parser("refine", help="FL -> FLNC")
    p.add_argument("--sample", help="sample name (default: SLURM array task id)")
    p.set_defaults(fn=step_refine)

    add_submit(sub, ["refine"]).set_defaults(fn=lambda c, a: submit_array(c, a, __file__))

    args = ap.parse_args()
    args.fn(Ctx(load_config(args.config), args.dry_run), args)


if __name__ == "__main__":
    main()
