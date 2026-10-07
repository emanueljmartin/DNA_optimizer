"""Fresh CDS extraction/translation; never use stored translation as authority."""

from Bio.Data.CodonTable import TranslationError as BioTranslationError
from Bio.Seq import Seq

from .errors import TranslationError
from .models import CDSMode, Feature, Strand


def extract_feature(sequence: str, feature: Feature) -> str:
    if feature.end > len(sequence):
        raise TranslationError(f"Feature {feature.id} is outside the sequence")
    dna = sequence[feature.start : feature.end]
    return str(Seq(dna).reverse_complement()) if feature.strand == Strand.REVERSE else dna


def translate_feature(sequence: str, feature: Feature) -> str:
    if not feature.is_cds:
        raise TranslationError(f"Feature {feature.id} is not a CDS")
    dna = extract_feature(sequence, feature)[feature.frame_offset :]
    if set(dna) - set("ACGT"):
        raise TranslationError("Translation requires concrete uppercase DNA")
    try:
        protein = str(
            Seq(dna).translate(
                table=feature.translation_table, cds=feature.cds_mode == CDSMode.COMPLETE
            )
        )
    except (BioTranslationError, ValueError) as exc:
        raise TranslationError(f"CDS {feature.id}: {exc}") from exc
    if feature.cds_mode == CDSMode.FRAGMENT and "*" in protein[:-1]:
        raise TranslationError(f"CDS {feature.id}: internal stop in coding fragment")
    return protein
