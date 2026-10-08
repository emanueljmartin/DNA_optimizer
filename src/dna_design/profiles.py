"""Local, provenance-backed codon references; no network access at runtime."""

import hashlib
import json
from importlib.resources import files
from itertools import product
from typing import Annotated, Literal, Self

from Bio.Data import CodonTable
from pydantic import Field, model_validator

from .models import FrozenModel, Identifier

CODONS = tuple("".join(bases) for bases in product("ACGT", repeat=3))


class CodonProfile(FrozenModel):
    id: Identifier
    version: Identifier
    organism: Identifier
    strain: Identifier
    translation_table: Annotated[int, Field(strict=True)] = 11
    reference_kind: Literal["genome_wide", "high_expression", "custom"]
    counts: dict[str, Annotated[int, Field(ge=0, strict=True)]]
    source_url: Identifier
    source_accession: Identifier
    source_sha256: Annotated[str, Field(pattern=r"^[0-9a-f]{64}$")]
    accessed_date: Annotated[str, Field(pattern=r"^\d{4}-\d{2}-\d{2}$")]
    derivation: Identifier
    included_cds: Annotated[int, Field(ge=1, strict=True)]
    excluded_cds: dict[str, Annotated[int, Field(ge=0, strict=True)]]
    limitations: tuple[str, ...] = ()

    @model_validator(mode="after")
    def valid_counts(self) -> Self:
        if set(self.counts) != set(CODONS):
            raise ValueError("Codon reference must explicitly contain all 64 DNA codons")
        table = CodonTable.unambiguous_dna_by_id.get(self.translation_table)
        if table is None or any(c in table.forward_table for c in table.stop_codons):
            raise ValueError("Unknown or context-dependent genetic code")
        if not sum(self.counts.values()):
            raise ValueError("Codon reference must contain observed codons")
        for amino_acid in set(table.forward_table.values()):
            if not any(self.counts[c] for c, aa in table.forward_table.items() if aa == amino_acid):
                raise ValueError(f"No reference codons for amino acid {amino_acid}")
        return self

    @property
    def fingerprint(self) -> str:
        payload = json.dumps(self.model_dump(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode()).hexdigest()


def load_codon_profile(host: str) -> CodonProfile:
    """Load a packaged profile carrying its source checksum; no network access."""
    aliases = {"ecoli": "ecoli-mg1655", "pao1": "pseudomonas-pao1"}
    identifier = aliases.get(host, host)
    if identifier not in {"ecoli-mg1655", "pseudomonas-pao1"}:
        raise ValueError(f"Unknown built-in host: {host}")
    resource = files("dna_design").joinpath("data", f"{identifier}.json")
    return CodonProfile.model_validate_json(resource.read_text())


def relative_codon_weights(profile: CodonProfile) -> dict[str, float]:
    """Observed count / maximum synonymous count. Zero counts stay zero."""
    table = CodonTable.unambiguous_dna_by_id[profile.translation_table]
    maxima = {
        aa: max(profile.counts[c] for c, amino in table.forward_table.items() if amino == aa)
        for aa in set(table.forward_table.values())
    }
    return {c: profile.counts[c] / maxima[aa] for c, aa in table.forward_table.items()}
