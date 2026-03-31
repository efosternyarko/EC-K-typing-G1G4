#!/usr/bin/env python3
"""
kpsm_screen.py — Group 2/3 pre-screen for E. coli G1/G4 K-locus typing.

RATIONALE
---------
Group 1 and 4 capsule loci (wzy-dependent, cps) are structurally distinct
from Group 2/3 capsule loci (ABC-transporter, kps). However, Group 2/3
strains carry a wzy-dependent O-antigen locus that can spuriously match
the G1/G4 Kaptive database, producing false-positive or low-confidence
G1/G4 calls (the "wzy-interference artefact").

kpsM encodes an inner membrane component of the Group 2/3 ABC-transporter
capsule export system. It is present in all Group 2/3 strains and absent
from Group 1/4 strains. kpsM_3 is the Group 3-specific paralogue. Both
are included as query sequences so that the screen captures Group 2 and
Group 3 strains. Detection of either gene (>=90% identity, >=80% coverage)
identifies a Group 2/3 organism for which G1/G4 typing should be suppressed.

USAGE
-----
    # Screen a single assembly:
    python kpsm_screen.py --assembly genome.fa --kpsm-ref DB/kpsM_reference.fasta

    # Screen a directory of assemblies and run Kaptive on kpsM-negative ones:
    python kpsm_screen.py \
        --assembly-dir assemblies/ \
        --kpsm-ref DB/kpsM_reference.fasta \
        --kaptive-db DB/EC-K-typing_group1and4_v1.2.gbk \
        --output-dir results/ \
        --run-kaptive

OUTPUT
------
    kpsm_screen_results.tsv  — per-assembly kpsM status (columns below)
    kaptive_output.tsv       — Kaptive results for kpsM-negative assemblies only

TSV COLUMNS
-----------
    assembly        Path to assembly FASTA
    kpsM_hit        True / False
    kpsM_pident     % identity of best kpsM hit (NA if no hit)
    kpsM_qcov       % query coverage of best hit (NA if no hit)
    group           Group2_3 / G1_G4_candidate
    kaptive_run     True / False (was Kaptive run?)
"""

import argparse
import subprocess
import sys
from pathlib import Path

# ── Thresholds ────────────────────────────────────────────────────────────────
KPSM_PIDENT_THRESH = 90.0   # % nucleotide identity for a kpsM hit
KPSM_QCOV_THRESH   = 80.0   # % query coverage for a kpsM hit


def minimap2_kpsm(assembly: Path, kpsm_ref: Path):
    """
    Map kpsM reference sequences against the assembly using minimap2.
    Returns (hit: bool, pident: float|None, qcov: float|None).

    Parses PAF output. Identity = matching_bases / alignment_block_len.
    Query coverage = aligned_query_bases / query_len.
    """
    result = subprocess.run(
        ["minimap2", "-c", "--secondary=no", str(assembly), str(kpsm_ref)],
        capture_output=True, text=True
    )

    best_pident, best_qcov = None, None
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 11:
            continue
        query_len    = int(parts[1])
        query_start  = int(parts[2])
        query_end    = int(parts[3])
        n_matches    = int(parts[9])
        aln_block    = int(parts[10])

        if aln_block == 0:
            continue

        pident = 100.0 * n_matches / aln_block
        qcov   = 100.0 * (query_end - query_start) / query_len

        if best_pident is None or pident > best_pident:
            best_pident = pident
            best_qcov   = qcov

    if best_pident is None:
        return False, None, None

    hit = (best_pident >= KPSM_PIDENT_THRESH) and (best_qcov >= KPSM_QCOV_THRESH)
    return hit, best_pident, best_qcov


def run_kaptive(assemblies: list, kaptive_db: Path, output_dir: Path):
    """Run Kaptive on a list of assembly paths."""
    kaptive_out = output_dir / "kaptive_output.tsv"
    cmd = [
        "kaptive", "assembly",
        str(kaptive_db),
        *[str(a) for a in assemblies],
        "-o", str(kaptive_out),
    ]
    print(f"  Running Kaptive on {len(assemblies)} assembly/assemblies...")
    subprocess.run(cmd, check=True)
    return kaptive_out


def screen_assemblies(assemblies: list, kpsm_ref: Path):
    """Screen a list of assemblies for kpsM. Returns list of result dicts."""
    results = []
    for i, asm in enumerate(assemblies, 1):
        asm = Path(asm)
        print(f"  [{i}/{len(assemblies)}] {asm.name} ...", end=" ", flush=True)
        hit, pident, qcov = minimap2_kpsm(asm, kpsm_ref)

        if hit:
            group = "Group2_3"
            print(f"kpsM+ (pident={pident:.1f}%, qcov={qcov:.1f}%) -> Group 2/3, G1/G4 suppressed")
        elif pident is not None:
            group = "G1_G4_candidate"
            print(f"kpsM- (best pident={pident:.1f}%) -> G1/G4 candidate")
        else:
            group = "G1_G4_candidate"
            print("no kpsM hit -> G1/G4 candidate")

        results.append({
            "assembly":    str(asm),
            "kpsM_hit":    hit,
            "kpsM_pident": f"{pident:.1f}" if pident is not None else "NA",
            "kpsM_qcov":   f"{qcov:.1f}"   if qcov   is not None else "NA",
            "group":        group,
        })
    return results


def write_tsv(results: list, path: Path):
    """Write screening results to TSV."""
    cols = ["assembly", "kpsM_hit", "kpsM_pident", "kpsM_qcov", "group", "kaptive_run"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in results:
            kaptive_run = r.get("kaptive_run", False)
            f.write("\t".join([
                r["assembly"],
                str(r["kpsM_hit"]),
                r["kpsM_pident"],
                r["kpsM_qcov"],
                r["group"],
                str(kaptive_run),
            ]) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Screen assemblies for kpsM (Group 2/3 marker) before G1/G4 Kaptive typing."
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--assembly", type=Path, help="Single assembly FASTA.")
    src.add_argument("--assembly-dir", type=Path,
                     help="Directory of assembly FASTAs (*.fa, *.fasta, *.fna).")
    src.add_argument("--assembly-list", type=Path,
                     help="Text file listing one assembly path per line.")

    parser.add_argument("--kpsm-ref", type=Path,
                        default=Path(__file__).parent.parent / "DB" / "kpsM_reference.fasta",
                        help="kpsM/kpsM_3 reference FASTA (default: DB/kpsM_reference.fasta). "
                             "Contains both kpsM (Group 2) and kpsM_3 (Group 3) references.")
    parser.add_argument("--kaptive-db", type=Path,
                        help="Kaptive G1/G4 database (.gbk). Required with --run-kaptive.")
    parser.add_argument("--output-dir", type=Path, default=Path("kpsm_screen_results"),
                        help="Output directory.")
    parser.add_argument("--run-kaptive", action="store_true",
                        help="Run Kaptive on kpsM-negative assemblies after screening.")
    parser.add_argument("--pident", type=float, default=KPSM_PIDENT_THRESH,
                        help=f"kpsM hit identity threshold (default: {KPSM_PIDENT_THRESH}).")
    parser.add_argument("--qcov", type=float, default=KPSM_QCOV_THRESH,
                        help=f"kpsM hit coverage threshold (default: {KPSM_QCOV_THRESH}).")

    args = parser.parse_args()

    global KPSM_PIDENT_THRESH, KPSM_QCOV_THRESH
    KPSM_PIDENT_THRESH = args.pident
    KPSM_QCOV_THRESH   = args.qcov

    # Collect assemblies
    if args.assembly:
        assemblies = [args.assembly]
    elif args.assembly_dir:
        assemblies = sorted(
            list(args.assembly_dir.glob("*.fa")) +
            list(args.assembly_dir.glob("*.fasta")) +
            list(args.assembly_dir.glob("*.fna"))
        )
        if not assemblies:
            sys.exit(f"No FASTA files found in {args.assembly_dir}")
    else:
        assemblies = [Path(l.strip()) for l in args.assembly_list.read_text().splitlines() if l.strip()]

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nkpsM pre-screen: {len(assemblies)} assembly/assemblies")
    print(f"  kpsM reference : {args.kpsm_ref}")
    print(f"  Thresholds     : pident >= {KPSM_PIDENT_THRESH}%, qcov >= {KPSM_QCOV_THRESH}%\n")

    results = screen_assemblies(assemblies, args.kpsm_ref)

    # Optionally run Kaptive on negatives
    candidates = [r["assembly"] for r in results if r["group"] == "G1_G4_candidate"]
    if args.run_kaptive:
        if not args.kaptive_db:
            sys.exit("--kaptive-db required when using --run-kaptive")
        if candidates:
            run_kaptive([Path(a) for a in candidates], args.kaptive_db, args.output_dir)
        else:
            print("No G1/G4 candidates — Kaptive not run.")
        for r in results:
            r["kaptive_run"] = r["group"] == "G1_G4_candidate"

    # Write results
    tsv_out = args.output_dir / "kpsm_screen_results.tsv"
    write_tsv(results, tsv_out)

    # Summary
    n_g23  = sum(1 for r in results if r["group"] == "Group2_3")
    n_cand = sum(1 for r in results if r["group"] == "G1_G4_candidate")
    print(f"\nSummary")
    print(f"  Group 2/3 (kpsM+, G1/G4 suppressed) : {n_g23}")
    print(f"  G1/G4 candidates (kpsM-)             : {n_cand}")
    print(f"  Results written to                   : {tsv_out}")


if __name__ == "__main__":
    main()
