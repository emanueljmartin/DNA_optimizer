"""Amino-acid input integrity and seeded initialization, not expression benchmarks."""

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from dna_design import CDSMode, translate_feature, validate_candidate
from dna_design.profiles import CodonProfile, load_codon_profile
from dna_design.reverse_translation import reverse_translate


@pytest.mark.parametrize("host", ["ecoli", "pao1"])
def test_reverse_translation_preserves_input_and_records_provenance(host: str) -> None:
    profile = load_codon_profile(host)
    result = reverse_translate("M ACDEFGHIKLMNPQRSTVWY*", profile, seed=2)
    assert translate_feature(result.sequence, result.features[0]) == "MACDEFGHIKLMNPQRSTVWY"
    assert result == reverse_translate("MACDEFGHIKLMNPQRSTVWY", profile, seed=2)
    assert result.random_seed == 2
    assert result.host_profile.version == profile.version
    assert result.model_versions["codon_profile_sha256"] == profile.fingerprint
    assert validate_candidate(result, result.sequence).valid


@pytest.mark.parametrize("protein", ["", "MXX", "M*E", "MUO", "-", "EAA"])
def test_unsupported_complete_proteins_are_rejected(protein: str) -> None:
    with pytest.raises(ValueError):
        reverse_translate(protein, load_codon_profile("ecoli"), seed=1)


def test_fragments_and_alternative_code_are_explicit() -> None:
    data = load_codon_profile("ecoli").model_dump()
    data["translation_table"] = 4
    profile = CodonProfile.model_validate(data)
    for protein in ("WEE", "WEE*"):
        fragment = reverse_translate(protein, profile, seed=1, cds_mode=CDSMode.FRAGMENT)
        assert translate_feature(fragment.sequence, fragment.features[0]) == protein
    complete = reverse_translate("MWEE", profile, seed=1)
    assert translate_feature(complete.sequence, complete.features[0]) == "MWEE"
    with pytest.raises(ValueError, match="seed"):
        reverse_translate("ME", profile, seed=True)


@settings(max_examples=20, deadline=None)
@given(st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=40), st.integers(0, 65535))
def test_seeded_reverse_translation_invariant(suffix: str, seed: int) -> None:
    protein = "M" + suffix
    construct = reverse_translate(protein, load_codon_profile("ecoli"), seed=seed)
    assert translate_feature(construct.sequence, construct.features[0]) == protein
    assert construct == reverse_translate(protein, load_codon_profile("ecoli"), seed=seed)
