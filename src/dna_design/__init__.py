"""DNA Design Platform: Phase 1 public domain interfaces."""

from .coordinates import internal_to_ui, ui_to_internal
from .errors import (
    DNAError,
    InvalidDesignError,
    ParseError,
    TranslationError,
    UnsupportedLocationError,
)
from .models import (
    CDSMode,
    Construct,
    Feature,
    MotifConstraint,
    MotifMode,
    MotifStrand,
    Mutability,
    Region,
    Strand,
    Topology,
    VersionedProfile,
)
from .motifs import MotifOccurrence, find_motif, restriction_constraint
from .parsing import parse_fasta, parse_genbank, parse_raw
from .permissions import BasePermission, effective_permissions
from .translation import extract_feature, translate_feature
from .validation import (
    CheckStatus,
    Finding,
    ValidationReport,
    check_feasibility,
    validate_candidate,
)

__version__ = "0.1.0"

__all__ = [
    "BasePermission",
    "CDSMode",
    "CheckStatus",
    "Construct",
    "DNAError",
    "Feature",
    "Finding",
    "InvalidDesignError",
    "MotifConstraint",
    "MotifMode",
    "MotifOccurrence",
    "MotifStrand",
    "Mutability",
    "ParseError",
    "Region",
    "Strand",
    "Topology",
    "TranslationError",
    "UnsupportedLocationError",
    "ValidationReport",
    "VersionedProfile",
    "check_feasibility",
    "effective_permissions",
    "extract_feature",
    "find_motif",
    "internal_to_ui",
    "parse_fasta",
    "parse_genbank",
    "parse_raw",
    "restriction_constraint",
    "translate_feature",
    "ui_to_internal",
    "validate_candidate",
]
