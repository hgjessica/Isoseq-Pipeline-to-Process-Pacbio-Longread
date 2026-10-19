# Isoseq-Pipeline-to-Process-Pacbio-Longread
Python scripts that run the PacBio Iso-Seq workflow one step at a time, locally or as SLURM job
arrays, split into two main parts (plus an optional third):

| Script | Covers |
|---|---|
| `01_sreads_to_flnc.py` | `lima` → `rename` → `refine` (S-reads → FL → FLNC) |
| `02_cluster_to_aligned.py` | `cluster` → `fix-sm` → `align` (FLNC → clustered → aligned BAM) |
| `03_collapse_classify.py` | `pigeon-prepare` → `collapse` (optional isoform classification) |

`common.py` holds the shared helpers and must stay in the same folder. All paths, sample names and cluster settings live in
one JSON config, so the code itself contains nothing specific to a project.

```
S-reads ──lima──▶ FL reads ──refine──▶ FLNC ──cluster2──▶ clustered ──pbmm2──▶ aligned (IGV)
                                         └──pbmm2──▶ aligned FLNC ──collapse + pigeon──▶ classified isoforms
```

## 1. Requirements

- Python ≥ 3.8 (standard library only)
- Tools on `PATH`: `lima`, `isoseq`, `pbmm2`, `pigeon`, `samtools` (and `pbindex` for `fix-sm`)
- Optionally SLURM (`sbatch`, `squeue`) for the `submit` command

Create the tool environment once from the included file:

```bash
conda env create -f environment.yml
conda activate isoseq
```

You also need:
- a reference genome FASTA (indexed automatically with `samtools faidx` if the `.fai` is missing)
- a gene annotation GTF
- the Iso-Seq cDNA primer FASTA used for your library prep

## 2. Configure

Copy `config.example.json` to `config.json` and fill in every `Path_to_...` value.

| Key | Meaning |
|---|---|
| `work_dir` | Directory holding the cell folders; all outputs are written here |
| `reference_fasta`, `annotation_gtf`, `primers_fasta` | Inputs described above |
| `threads` | Threads per job |
| `scratch_dir` | Optional. Per-task temp dir (use a disk with plenty of free space) |
| `slurm` | Partition, memory, time limit, array concurrency, conda setup (used by `submit`) |
| `cells` | Cell folders in `work_dir`, each with `segmented.bam` (or a `<cell>.tar` that contains it). Only needed for the `lima` step |
| `samples` | One entry per sample, see below |

Sample entries:

```json
{"name": "sample_1", "cell": "cell_1", "barcode": "IsoSeqX_bc01_5p", "keep": true}
```

- `barcode` maps lima's per-barcode output file to the sample name.
- `keep: false` marks samples on the same cell that are not yours. They are skipped by every
  per-sample step, and `rename --drop-others` can delete their demultiplexed files.
- Already-demultiplexed FL data (no `lima` step needed): give `"fl_bam": "<path>"` instead of
  `cell`/`barcode`.

Keep `config.json` out of version control if it contains real names or paths.

## 3. Run

Add `--dry-run` before the subcommand to print the commands without running anything.
Per-sample steps take `--sample NAME`; under SLURM they use the array task id instead.

```bash
P1="python 01_sreads_to_flnc.py --config config.json"
P2="python 02_cluster_to_aligned.py --config config.json"
P3="python 03_collapse_classify.py --config config.json"
```

**Part 1 — S-reads → FLNC** (`01_sreads_to_flnc.py`)

| Step | Command | Output |
|---|---|---|
| lima | `$P1 lima` | `fl/<cell>.fl.IsoSeqX_bcNN_5p--IsoSeqX_3p.bam` |
| rename (barcode → sample) | `$P1 rename` (check), then `$P1 rename --apply` | `fl/<sample>.fl.bam` |
| refine | `$P1 submit refine` | `flnc/<sample>.flnc.bam` |

**Part 2 — cluster → aligned** (`02_cluster_to_aligned.py`)

| Step | Command | Output |
|---|---|---|
| cluster | `$P2 submit cluster` | `clustered/<sample>.clustered.bam` |
| fix-sm (if the SM label is wrong) | `$P2 fix-sm` (check), then `$P2 fix-sm --apply` | rewritten in place |
| align clustered transcripts | `$P2 submit align` | `aligned/<sample>.aligned.bam` + `.bai` |
| align FLNC reads (for part 3) | `$P2 submit align --input-stage flnc` | `aligned_flnc/<sample>.flnc.aligned.bam` |

**Part 3 (optional) — collapse and classify** (`03_collapse_classify.py`)

| Step | Command | Output |
|---|---|---|
| prepare annotation (once) | `$P3 pigeon-prepare` | sorted GTF + `.pgi` |
| collapse + classify | `$P3 submit collapse` | `pigeon/<sample>/` |

Notes:
- `submit STEP [args]` (available in each script for its per-sample steps) writes `slurm/STEP.sbatch` (one array task per kept sample, limited to
  `array_concurrency` at once), prints it, and calls `sbatch`. Run a step directly (without
  `submit`) to test one sample first.
- Run `lima` and `rename` on a login node or as one job; they work per cell, not per sample.
- Align from a different stage with `--input-stage fl` (aligns FL reads before refine, so the
  polyA tail shows as a soft clip) or `--input-stage flnc`.
- `align` and `collapse` need their input stage to exist. Wait for the previous array to finish
  before submitting the next step.

## 4. Things that commonly go wrong

- **A cell is a pool of several samples.** The barcode lives in the read sequence, not a BAM
  tag, so the tar or cell name is *not* the sample name. Confirm the barcode → sample table with
  your sequencing provider and check lima's `*.lima.guess`/summary files before trusting it.
- **Wrong sample label inside the BAM.** The instrument may stamp one `SM:` value on the whole
  cell, and lima does not rewrite it when splitting by barcode. Run `fix-sm` on `clustered/`
  (small files) rather than `flnc/` (very large) before aligning, because downstream tools read
  `SM:`, not the filename.
- **Rename safety.** `rename --apply` refuses to run while a job named `lima` is in the queue or
  `fl/` contains `*.build` files (an index being written). Wait rather than bypassing it.
- **Full node scratch.** If node-local `/tmp` is full, set `scratch_dir` to a larger disk.
  Each task uses its own subdirectory so concurrent tasks do not collide.
- **Pigeon reference.** Run `pigeon-prepare` once before `collapse`. Running it inside the array
  would have every task writing the same files.
- **Why align FLNC and not only clustered.** Clustering collapses ~20 reads into one record and
  drops singletons, which are often the rare isoforms. Collapse needs aligned FLNC to give
  per-isoform read counts.
