"""Typed controls for Phase 1. Canonical DNA is concrete uppercase A/C/G/T."""

from enum import IntEnum, StrEnum
from typing import Annotated, Self

from Bio.Data import CodonTable, IUPACData
from pydantic import BaseModel, ConfigDict, Field, JsonValue, field_validator, model_validator

Coordinate = Annotated[int, Field(ge=0, strict=True)]
Identifier = Annotated[str, Field(min_length=1)]


class FrozenModel(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


class Topology(StrEnum):
    LINEAR = "linear"
    CIRCULAR = "circular"


class Strand(IntEnum):
    REVERSE = -1
    UNKNOWN = 0
    FORWARD = 1


class Mutability(StrEnum):
    LOCKED = "LOCKED"
    SYNONYMOUS_ONLY = "SYNONYMOUS_ONLY"
    FREE_MUTATION = "FREE_MUTATION"
    DESIGNABLE_RBS = "DESIGNABLE_RBS"


class CDSMode(StrEnum):
    COMPLETE = "complete"
    FRAGMENT = "fragment"


class MotifMode(StrEnum):
    PRESERVE = "PRESERVE"
    AVOID = "AVOID"
    REQUIRED = "REQUIRED"


class MotifStrand(StrEnum):
    FORWARD = "forward"
    REVERSE = "reverse"
    BOTH = "both"


class Region(FrozenModel):
    start: Coordinate
    end: Coordinate

    @model_validator(mode="after")
    def nonempty(self) -> Self:
        if self.end <= self.start:
            raise ValueError("Intervals must be nonempty 0-based half-open")
        return self


class Feature(Region):
    id: Identifier
    name: Identifier
    feature_type: Identifier
    strand: Strand = Strand.UNKNOWN
    mutability: Mutability | None = None
    description: str | None = None
    translation_table: Annotated[int, Field(strict=True)] = 11
    frame_offset: Annotated[int, Field(ge=0, le=2, strict=True)] = 0
    cds_mode: CDSMode = CDSMode.COMPLETE
    rbs_mask: str | None = None
    qualifiers: dict[str, tuple[str, ...]] = Field(default_factory=dict)

    @model_validator(mode="after")
    def biological_controls(self) -> Self:
        if self.mutability is None:
            default = Mutability.SYNONYMOUS_ONLY if self.is_cds else Mutability.LOCKED
            object.__setattr__(self, "mutability", default)
        if self.is_cds:
            if self.strand == Strand.UNKNOWN:
                raise ValueError("CDS strand must be explicitly +1 or -1")
            if self.translation_table not in CodonTable.unambiguous_dna_by_id:
                raise ValueError("Unknown NCBI genetic code")
            table = CodonTable.unambiguous_dna_by_id[self.translation_table]
            if any(codon in table.forward_table for codon in table.stop_codons):
                raise ValueError("Context-dependent stop/sense genetic codes are unsupported")
            if self.cds_mode == CDSMode.COMPLETE and self.frame_offset:
                raise ValueError("A complete CDS cannot have a leading frame offset")
            coding_length = self.end - self.start - self.frame_offset
            if coding_length < 3 or coding_length % 3:
                raise ValueError("CDS coding length after frame offset must be a multiple of 3")
        elif self.frame_offset or self.cds_mode != CDSMode.COMPLETE:
            raise ValueError("Frame offset and CDS mode apply only to CDS features")
        if self.mutability == Mutability.DESIGNABLE_RBS:
            if self.rbs_mask is None:
                raise ValueError("DESIGNABLE_RBS requires an explicit N/ACGT mask")
            if len(self.rbs_mask) != self.end - self.start:
                raise ValueError("RBS mask length must equal its feature interval")
            if set(self.rbs_mask) - set("ACGTN"):
                raise ValueError("RBS mask uses N for mutable and A/C/G/T for fixed bases")
            if self.strand == Strand.UNKNOWN:
                raise ValueError("A designable RBS requires an explicit strand")
        elif self.rbs_mask is not None:
            raise ValueError("RBS masks require DESIGNABLE_RBS mutability")
        return self

    @property
    def is_cds(self) -> bool:
        return self.feature_type == "CDS"


class MotifConstraint(FrozenModel):
    id: Identifier
    mode: MotifMode
    motif: str
    strand: MotifStrand = MotifStrand.BOTH
    region: Region | None = None
    position: Coordinate | None = None
    enzyme: str | None = None
    description: str | None = None

    @field_validator("motif")
    @classmethod
    def valid_motif(cls, value: str) -> str:
        if not value or set(value) - set(IUPACData.ambiguous_dna_values):
            raise ValueError("Motifs must be nonempty uppercase IUPAC DNA")
        return value


class VersionedProfile(FrozenModel):
    id: Identifier
    version: Identifier
    source: str | None = None


class Construct(FrozenModel):
    name: Identifier = "untitled"
    sequence: str
    topology: Topology = Topology.LINEAR
    source_sequence: str | None = None
    features: tuple[Feature, ...] = ()
    constraints: tuple[MotifConstraint, ...] = ()
    optimization_history: tuple[dict[str, str], ...] = ()
    host_profile: VersionedProfile | None = None
    synthesis_profile: VersionedProfile | None = None
    rbs_settings: dict[str, str | int | float | bool] = Field(default_factory=dict)
    fragmentation_settings: dict[str, str | int | float | bool] = Field(default_factory=dict)
    random_seed: Annotated[int, Field(strict=True)] | None = None
    software_version: str = "0.1.0"
    model_versions: dict[str, str] = Field(default_factory=dict)
    input_format: str = "raw"
    record_annotations: dict[str, JsonValue] = Field(default_factory=dict)
    original_record_text: str | None = None
    import_warnings: tuple[str, ...] = ()

    @field_validator("sequence", "source_sequence")
    @classmethod
    def concrete_dna(cls, value: str | None) -> str | None:
        if value is not None and (not value or set(value) - set("ACGT")):
            raise ValueError("DNA must be nonempty concrete uppercase A/C/G/T")
        return value

    @model_validator(mode="after")
    def controls_in_bounds(self) -> Self:
        if self.source_sequence is None:
            object.__setattr__(self, "source_sequence", self.sequence)
        if len({f.id for f in self.features}) != len(self.features):
            raise ValueError("Feature IDs must be unique")
        if len({c.id for c in self.constraints}) != len(self.constraints):
            raise ValueError("Constraint IDs must be unique")
        for feature in self.features:
            if feature.end > len(self.sequence):
                raise ValueError(f"Feature {feature.id} exceeds construct length")
        for constraint in self.constraints:
            if len(constraint.motif) > len(self.sequence):
                raise ValueError("A motif cannot exceed one complete construct length")
            if constraint.region and constraint.region.end > len(self.sequence):
                raise ValueError("Motif region exceeds construct length")
            if constraint.position is not None and constraint.position >= len(self.sequence):
                raise ValueError("Motif position exceeds construct length")
        return self
