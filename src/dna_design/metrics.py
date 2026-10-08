"""Transparent local sequence and codon statistics, separate from expression claims."""

import math
from collections import Counter
from typing import Annotated

from Bio.Data import CodonTable
from pydantic import Field

from .models import CDSMode, Construct, FrozenModel
from .profiles import CodonProfile, relative_codon_weights
from .translation import extract_feature, translate_feature


def gc_fraction(sequence: str) -> float:
    if not sequence or set(sequence) - set("ACGT"):
        raise ValueError("GC statistics require nonempty concrete uppercase DNA")
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def local_gc(sequence: str, window: int) -> tuple[tuple[int, float], ...]:
    gc_fraction(sequence)
    if type(window) is not int or window < 1:
        raise ValueError("GC window must be a positive integer")
    return tuple(
        (start, gc_fraction(sequence[start : start + window]))
        for start in range(len(sequence) - window + 1)
    )


def sequence_entropy(sequence: str) -> float:
    gc_fraction(sequence)
    return -sum(
        (n / len(sequence)) * math.log2(n / len(sequence)) for n in Counter(sequence).values()
    )


class CDSCodonMetrics(FrozenModel):
    feature_id: str
    codon_counts: dict[str, int]
    adaptation_score: Annotated[float, Field(ge=0, le=1)]
    cai: Annotated[float, Field(ge=0, le=1)] | None
    rare_fraction: Annotated[float, Field(ge=0, le=1)] | None
    rare_codon_indices: tuple[int, ...]
    rare_clusters: tuple[tuple[int, ...], ...]
    warnings: tuple[str, ...]


def codon_metrics(
    reference: Construct,
    sequence: str,
    profile: CodonProfile,
    *,
    rare_threshold: float | None = None,
) -> tuple[CDSCodonMetrics, ...]:
    if rare_threshold is not None and (
        not math.isfinite(rare_threshold) or not 0 <= rare_threshold <= 1
    ):
        raise ValueError("Rare threshold must be explicitly supplied in [0,1]")
    weights = relative_codon_weights(profile)
    table = CodonTable.unambiguous_dna_by_id[profile.translation_table]
    synonym_counts = Counter(table.forward_table.values())
    result: list[CDSCodonMetrics] = []
    for feature in reference.features:
        if not feature.is_cds:
            continue
        if feature.translation_table != profile.translation_table:
            raise ValueError(f"CDS {feature.id} genetic code differs from the codon reference")
        translate_feature(sequence, feature)
        dna = extract_feature(sequence, feature)[feature.frame_offset :]
        codons = [(i // 3, dna[i : i + 3]) for i in range(0, len(dna), 3)]
        if feature.cds_mode == CDSMode.COMPLETE:
            codons = codons[1:-1]
        codons = [(i, c) for i, c in codons if c not in table.stop_codons]
        values = [weights[c] for _, c in codons if synonym_counts[table.forward_table[c]] > 1]
        score = (
            1.0
            if not values
            else 0.0
            if 0.0 in values
            else math.exp(sum(map(math.log, values)) / len(values))
        )
        rare = tuple(
            i for i, c in codons if rare_threshold is not None and weights[c] < rare_threshold
        )
        runs: list[list[int]] = []
        for index in rare:
            if runs and runs[-1][-1] == index - 1:
                runs[-1].append(index)
            else:
                runs.append([index])
        warnings = (
            ()
            if profile.reference_kind == "high_expression"
            else ("Genome/custom reference adaptation is not high-expression-reference CAI",)
        )
        result.append(
            CDSCodonMetrics(
                feature_id=feature.id,
                codon_counts=dict(Counter(c for _, c in codons)),
                adaptation_score=score,
                cai=score if profile.reference_kind == "high_expression" else None,
                rare_fraction=len(rare) / len(codons)
                if rare_threshold is not None and codons
                else 0.0
                if rare_threshold is not None
                else None,
                rare_codon_indices=rare,
                rare_clusters=tuple(tuple(run) for run in runs if len(run) > 1),
                warnings=warnings,
            )
        )
    return tuple(result)
