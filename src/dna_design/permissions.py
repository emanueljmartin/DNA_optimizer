"""Compose per-base controls. Permission alone never proves synonymy."""

from .models import Construct, FrozenModel, MotifMode, Mutability, Strand
from .motifs import find_motif


class BasePermission(FrozenModel):
    position: int
    state: Mutability
    feature_ids: tuple[str, ...]
    cds_ids: tuple[str, ...]
    reasons: tuple[str, ...]


def effective_permissions(reference: Construct) -> tuple[BasePermission, ...]:
    preserved: dict[int, list[str]] = {}
    for constraint in reference.constraints:
        if constraint.mode == MotifMode.PRESERVE:
            for occurrence in find_motif(reference.sequence, constraint, reference.topology):
                for position in occurrence.positions(len(reference.sequence)):
                    preserved.setdefault(position, []).append(f"PRESERVE motif {constraint.id}")
    result: list[BasePermission] = []
    for position in range(len(reference.sequence)):
        features = [f for f in reference.features if f.start <= position < f.end]
        cds = [f for f in features if f.is_cds]
        locks: list[str] = list(preserved.get(position, ()))
        if not features:
            locks.append("unannotated DNA")
        for feature in features:
            if feature.mutability == Mutability.LOCKED:
                locks.append(f"LOCKED feature {feature.id}")
            if feature.mutability == Mutability.SYNONYMOUS_ONLY and not cds:
                locks.append(f"no CDS interpretation for {feature.id}")
            if feature.is_cds and feature.frame_offset:
                offset = (
                    feature.end - 1 - position
                    if feature.strand == Strand.REVERSE
                    else position - feature.start
                )
                if offset < feature.frame_offset:
                    locks.append(f"untranslated frame-offset base in {feature.id}")
            if feature.rbs_mask is not None:
                offset = (
                    feature.end - 1 - position
                    if feature.strand == Strand.REVERSE
                    else position - feature.start
                )
                if feature.rbs_mask[offset] != "N":
                    locks.append(f"fixed RBS mask base in {feature.id}")
        if locks:
            state = Mutability.LOCKED
        elif cds:
            state = Mutability.SYNONYMOUS_ONLY
        elif any(f.mutability == Mutability.SYNONYMOUS_ONLY for f in features):
            state = Mutability.LOCKED
        elif any(f.mutability == Mutability.DESIGNABLE_RBS for f in features):
            state = Mutability.DESIGNABLE_RBS
        else:
            state = Mutability.FREE_MUTATION
        result.append(
            BasePermission(
                position=position,
                state=state,
                feature_ids=tuple(f.id for f in features),
                cds_ids=tuple(f.id for f in cds),
                reasons=tuple(locks),
            )
        )
    return tuple(result)
