"""Overlapping, strand-aware IUPAC recognition matching, including circular origins."""

import re

from Bio import Restriction
from Bio.Data import IUPACData
from Bio.Restriction.Restriction import AllEnzymes
from Bio.Seq import Seq

from .models import (
    FrozenModel,
    MotifConstraint,
    MotifMode,
    MotifStrand,
    Region,
    Strand,
    Topology,
)


class MotifOccurrence(FrozenModel):
    start: int
    length: int
    strand: Strand

    def positions(self, sequence_length: int) -> tuple[int, ...]:
        return tuple((self.start + offset) % sequence_length for offset in range(self.length))


def _pattern(motif: str) -> str:
    return "".join(f"[{IUPACData.ambiguous_dna_values[symbol]}]" for symbol in motif)


def find_motif(
    sequence: str, constraint: MotifConstraint, topology: Topology = Topology.LINEAR
) -> tuple[MotifOccurrence, ...]:
    """Occurrence start is always a 0-based position on the forward genomic axis."""
    length = len(constraint.motif)
    if not sequence or length > len(sequence):
        return ()
    search = sequence
    if topology == Topology.CIRCULAR:
        search += sequence[: length - 1]
    patterns: list[tuple[str, Strand]] = []
    if constraint.strand != MotifStrand.REVERSE:
        patterns.append((constraint.motif, Strand.FORWARD))
    if constraint.strand != MotifStrand.FORWARD:
        patterns.append((str(Seq(constraint.motif).reverse_complement()), Strand.REVERSE))
    matches: list[MotifOccurrence] = []
    for motif, strand in patterns:
        for match in re.finditer(f"(?=({_pattern(motif)}))", search):
            start = match.start()
            if start >= len(sequence):
                continue
            if constraint.position is not None and start != constraint.position:
                continue
            occurrence = MotifOccurrence(start=start, length=length, strand=strand)
            if constraint.region and not all(
                constraint.region.start <= p < constraint.region.end
                for p in occurrence.positions(len(sequence))
            ):
                continue
            matches.append(occurrence)
    return tuple(sorted(matches, key=lambda item: (item.start, item.strand)))


def restriction_constraint(
    enzyme: str,
    mode: MotifMode,
    *,
    id: str | None = None,
    region: Region | None = None,
    position: int | None = None,
    strand: MotifStrand = MotifStrand.BOTH,
) -> MotifConstraint:
    """Resolve a recognition motif, not restriction cuts or methylation behavior."""
    enzyme_class = getattr(Restriction, enzyme, None)
    if enzyme_class is None or enzyme_class not in AllEnzymes:
        raise ValueError(f"Unknown restriction enzyme: {enzyme}")
    return MotifConstraint(
        id=id or f"{mode.value}:{enzyme}",
        motif=str(enzyme_class.site),
        mode=mode,
        enzyme=enzyme,
        region=region,
        position=position,
        strand=strand,
    )
