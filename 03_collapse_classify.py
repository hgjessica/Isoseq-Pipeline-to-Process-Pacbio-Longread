#!/usr/bin/env python3
"""Part 3 (optional): aligned FLNC reads -> collapsed isoforms -> pigeon classification.

    python 03_collapse_classify.py --config config.json pigeon-prepare     # once
    python 03_collapse_classify.py --config config.json collapse --sample sample_1
    python 03_collapse_classify.py --config config.json submit collapse

Needs aligned FLNC from part 2 (align --input-stage flnc) and flnc/<sample>.flnc.bam.
Output: pigeon/<sample>/
"""
import re
import subprocess
import sys
from pathlib import Path

from common import *  # noqa: F401,F403  (helpers shared by all steps)
from common import submit_array


def sorted_gtf(ctx):
    gtf = ctx.cfg["annotation_gtf"]
    return re.sub(r"\.gtf$", ".sorted.gtf", gtf)


def step_pigeon_prepare(ctx, args):
    """One-off: sort + index the annotation for pigeon. Run BEFORE collapse."""
    gtf, fa = ctx.cfg["annotation_gtf"], ctx.cfg["reference_fasta"]
    need(ctx, gtf, fa)
    ensure_fai(ctx)
    srt = Path(sorted_gtf(ctx))
    if srt.exists() and Path(str(srt) + ".pgi").exists():
        print(f"[HAVE] {srt} and .pgi already exist - nothing to do")
        return
    if run(ctx, ["pigeon", "prepare", gtf, fa]) != 0:
        sys.exit("[ERROR] pigeon prepare failed")


def step_collapse(ctx, args):
    """Collapse aligned FLNC into isoforms, then pigeon classify/filter/report."""
    s = pick_sample(ctx.cfg, args.sample)
    name = s["name"]
    task_tmp(ctx, "pigeon")
    aln = ctx.path("aligned_flnc", f"{name}.flnc.aligned.bam")
    flnc = ctx.path("flnc", f"{name}.flnc.bam")
    gtf, fa = sorted_gtf(ctx), ctx.cfg["reference_fasta"]
    out = ctx.path("pigeon", name)
    need(ctx, aln, flnc, gtf, gtf + ".pgi", fa + ".fai")
    if not ctx.dry:
        out.mkdir(parents=True, exist_ok=True)

    gff = out / f"{name}.collapsed.gff"
    steps = [("[1/5] isoseq collapse", ["isoseq", "collapse", "-j", ctx.threads,
                                          "--log-level", "INFO", aln, flnc, gff]),
             ("[2/5] pigeon prepare", ["pigeon", "prepare", gff])]
    for label, cmd in steps:
        print(label)
        if run(ctx, cmd) != 0:
            sys.exit(f"[ERROR] {label} failed")

    sorted_gff = out / f"{name}.collapsed.sorted.gff"
    if not sorted_gff.exists():
        sorted_gff = gff
    abund = out / f"{name}.collapsed.abundance.txt"
    print("[3/5] pigeon classify")
    cmd = ["pigeon", "classify", "-j", ctx.threads, "--log-level", "INFO"]
    if abund.exists():
        cmd += ["--fl", abund]
    cmd += ["-o", out / name, sorted_gff, gtf, fa]
    if run(ctx, cmd) != 0:
        sys.exit("[ERROR] pigeon classify failed")

    cls = out / f"{name}_classification.txt"
    if not cls.exists() and not ctx.dry:
        found = sorted(out.glob("*_classification.txt"))
        cls = found[0] if found else cls
    print("[4/5] pigeon filter")
    if run(ctx, ["pigeon", "filter", cls]) != 0:
        print("[WARN] pigeon filter returned non-zero")

    print("[5/5] pigeon report")
    filt = sorted(out.glob("*filtered_lite_classification.txt")) if not ctx.dry else []
    if filt:
        if run(ctx, ["pigeon", "report", filt[0], out / f"{name}_saturation.txt"]) != 0:
            print("[WARN] pigeon report failed")
    print(f"[OK] {out}")


def main():
    ap = base_parser(__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("pigeon-prepare", help="one-off annotation prep").set_defaults(fn=step_pigeon_prepare)

    p = sub.add_parser("collapse", help="collapse + pigeon classify/filter/report")
    p.add_argument("--sample", help="sample name (default: SLURM array task id)")
    p.set_defaults(fn=step_collapse)

    add_submit(sub, ["collapse"]).set_defaults(fn=lambda c, a: submit_array(c, a, __file__))

    args = ap.parse_args()
    args.fn(Ctx(load_config(args.config), args.dry_run), args)


if __name__ == "__main__":
    main()
