"""Reproduce aggregate host counts from an explicitly supplied public GenBank file.

This script never downloads data. Compound/fuzzy/partial/exception CDSs are excluded
with counts, rather than flattened. Starts and terminal stops are not counted as
elongation codons. Bundled profiles contain aggregate counts, not genomic sequences.
"""

import argparse
import hashlib
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.SeqFeature import ExactPosition, SimpleLocation

from dna_design.profiles import CODONS, CodonProfile


def derive(path: Path, identifier: str, organism: str, strain: str, date: str) -> CodonProfile:
    record = SeqIO.read(path, "genbank")
    if record.annotations.get("molecule_type") != "DNA":
        raise ValueError("Codon profile derivation requires an explicit DNA record")
    counts: Counter[str] = Counter({codon: 0 for codon in CODONS})
    exclusions: Counter[str] = Counter()
    included = 0
    for feature in record.features:
        if feature.type != "CDS":
            continue
        qualifiers = feature.qualifiers
        if {"pseudo", "pseudogene", "exception", "transl_except", "ribosomal_slippage"} & set(
            qualifiers
        ):
            exclusions["exception_or_pseudogene"] += 1
            continue
        location = feature.location
        if (
            not isinstance(location, SimpleLocation)
            or not isinstance(location.start, ExactPosition)
            or not isinstance(location.end, ExactPosition)
            or location.ref is not None
            or location.strand not in (1, -1)
        ):
            exclusions["unsupported_location"] += 1
            continue
        if qualifiers.get("transl_table", ["1"]) != ["11"] or qualifiers.get(
            "codon_start", ["1"]
        ) != ["1"]:
            exclusions["other_code_or_frame"] += 1
            continue
        dna = feature.extract(record.seq)
        try:
            protein = str(dna.translate(table=11, cds=True))
        except (ValueError, TranslationError):
            exclusions["invalid_complete_cds"] += 1
            continue
        if str(dna) != str(dna).upper() or set(str(dna)) - set("ACGT"):
            exclusions["ambiguous_dna"] += 1
            continue
        translation = qualifiers.get("translation")
        if translation is not None and translation != [protein]:
            exclusions["translation_disagreement"] += 1
            continue
        counts.update(str(dna[i : i + 3]) for i in range(3, len(dna) - 3, 3))
        included += 1
    return CodonProfile(
        id=identifier,
        version=f"{record.id}-genome-cds-v1-{date}",
        organism=organism,
        strain=strain,
        reference_kind="genome_wide",
        counts=dict(counts),
        source_url=f"https://www.ncbi.nlm.nih.gov/nuccore/{record.id}",
        source_accession=record.id,
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        accessed_date=date,
        derivation="derive_codon_profile.py v1: exact complete table-11 CDSs; exclude start/stop",
        included_cds=included,
        excluded_cds=dict(exclusions),
        limitations=(
            "Genome-wide reference, not a high-expression gene set or expression predictor",
            "Excluded exceptional/partial/compound CDSs are documented in excluded_cds",
            "Annotation snapshot is bound to the downloaded record checksum",
            "NCBI public data policy applies; third-party rights are not independently adjudicated",
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--id", required=True)
    parser.add_argument("--organism", required=True)
    parser.add_argument("--strain", required=True)
    parser.add_argument("--accessed-date", required=True)
    args = parser.parse_args()
    profile = derive(args.input, args.id, args.organism, args.strain, args.accessed_date)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(profile.model_dump_json(indent=2) + "\n")
    print(profile.id, profile.included_cds, profile.excluded_cds, profile.fingerprint)
