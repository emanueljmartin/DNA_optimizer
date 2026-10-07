"""Local single-record import with explicit representational limits."""

import warnings
from io import StringIO

from Bio import BiopythonParserWarning, SeqIO
from Bio.SeqFeature import CompoundLocation, ExactPosition
from pydantic import ValidationError

from .errors import ParseError, TranslationError, UnsupportedLocationError
from .models import CDSMode, Construct, Feature, Strand, Topology
from .translation import translate_feature


def _normalize(text: str) -> str:
    dna = "".join(text.split()).upper()
    if not dna or set(dna) - set("ACGT"):
        raise ParseError("Expected nonempty DNA containing only A/C/G/T (case is normalized)")
    return dna


def parse_raw(
    text: str, *, name: str = "untitled", topology: Topology = Topology.LINEAR
) -> Construct:
    try:
        return Construct(name=name, sequence=_normalize(text), topology=topology)
    except ValidationError as exc:
        raise ParseError(str(exc)) from exc


def parse_fasta(text: str, *, topology: Topology = Topology.LINEAR) -> Construct:
    lines = text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ParseError("FASTA must begin with a header")
    if any(line.startswith(">") for line in lines[1:]):
        raise ParseError("Exactly one FASTA record is required; multiple ORFs use annotations")
    name = lines[0][1:].strip()
    if not name:
        raise ParseError("FASTA header cannot be empty")
    try:
        return Construct(
            name=name,
            sequence=_normalize("\n".join(lines[1:])),
            topology=topology,
            input_format="fasta",
            original_record_text=text,
        )
    except ValidationError as exc:
        raise ParseError(str(exc)) from exc


def _single_qualifier(qualifiers: dict[str, tuple[str, ...]], key: str, default: str) -> str:
    values = qualifiers.get(key, (default,))
    if len(values) != 1:
        raise ParseError(f"GenBank qualifier {key} must have exactly one value")
    return values[0]


def parse_genbank(text: str) -> Construct:
    """Reject all features we cannot faithfully interpret; never flatten joins."""
    if not text.lstrip().startswith("LOCUS ") or not text.rstrip().endswith("//"):
        raise ParseError("GenBank must begin with LOCUS and end with //; surrounding data rejected")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            records = list(SeqIO.parse(StringIO(text), "genbank"))
    except (ValueError, BiopythonParserWarning) as exc:
        raise ParseError(f"Malformed GenBank: {exc}") from exc
    if len(records) != 1:
        raise ParseError("Exactly one GenBank record is required")
    record = records[0]
    try:
        dna = _normalize(str(record.seq))
    except Exception as exc:
        raise ParseError(f"GenBank sequence is unavailable or invalid: {exc}") from exc
    features: list[Feature] = []
    for index, biofeature in enumerate(record.features):
        location = biofeature.location
        if location is None or isinstance(location, CompoundLocation):
            raise UnsupportedLocationError(
                f"Feature {index} ({biofeature.type}): missing/compound location unsupported"
            )
        if (
            not isinstance(location.start, ExactPosition)
            or not isinstance(location.end, ExactPosition)
            or location.ref is not None
            or location.ref_db is not None
        ):
            raise UnsupportedLocationError(
                f"Feature {index} ({biofeature.type}): fuzzy/remote location unsupported"
            )
        qualifiers = {key: tuple(values) for key, values in biofeature.qualifiers.items()}
        unsupported = {"transl_except", "exception", "ribosomal_slippage", "pseudo", "pseudogene"}
        if biofeature.type == "CDS" and unsupported.intersection(qualifiers):
            raise ParseError(f"Feature {index}: unsupported CDS exception/pseudogene qualifier")
        # Source annotation is retained as metadata and in the original text; treating
        # its whole-record interval as a functional lock would prevent all CDS edits.
        if biofeature.type == "source":
            continue
        name = next(
            (qualifiers[key][0] for key in ("label", "gene", "locus_tag") if qualifiers.get(key)),
            f"{biofeature.type}_{index}",
        )
        try:
            offset = int(_single_qualifier(qualifiers, "codon_start", "1")) - 1
            table = int(_single_qualifier(qualifiers, "transl_table", "11"))
            feature = Feature(
                id=f"gb:{index}",
                name=name,
                feature_type=biofeature.type,
                start=int(location.start),
                end=int(location.end),
                strand=Strand(location.strand or 0),
                frame_offset=offset,
                translation_table=table,
                cds_mode=CDSMode.FRAGMENT if offset else CDSMode.COMPLETE,
                description="; ".join(qualifiers.get("note", ())) or None,
                qualifiers=qualifiers,
            )
            if feature.is_cds:
                translated = translate_feature(dna, feature)
                if "translation" in qualifiers:
                    annotated = "".join(_single_qualifier(qualifiers, "translation", "").split())
                    if annotated != translated.removesuffix("*"):
                        raise ParseError(
                            f"Feature {index}: annotated translation disagrees with DNA"
                        )
        except (ValidationError, ValueError, TranslationError) as exc:
            raise ParseError(f"Feature {index}: {exc}") from exc
        features.append(feature)
    # SeqIO Reference objects are not JSON values. Preserve the complete original
    # text and an explicit textual reference representation rather than drop them.
    annotations = dict(record.annotations)
    if "references" in annotations:
        annotations["references"] = [str(ref) for ref in annotations["references"]]
    annotations["source_features"] = [
        {
            "location": str(f.location),
            "qualifiers": f.qualifiers,
        }
        for f in record.features
        if f.type == "source"
    ]
    try:
        return Construct(
            name=record.name,
            sequence=dna,
            topology=record.annotations.get("topology", "linear"),
            features=tuple(features),
            input_format="genbank",
            record_annotations=annotations,
            original_record_text=text,
        )
    except ValidationError as exc:
        raise ParseError(f"GenBank controls invalid: {exc}") from exc
