"""Seeded amino-acid input conversion; independently verify the generated CDS."""

import random

from Bio.Data import CodonTable

from .models import CDSMode, Construct, Feature, Strand, VersionedProfile
from .profiles import CodonProfile
from .translation import translate_feature
from .validation import validate_candidate


def reverse_translate(
    protein: str,
    profile: CodonProfile,
    *,
    seed: int,
    cds_mode: CDSMode = CDSMode.COMPLETE,
    name: str = "reverse-translated protein",
) -> Construct:
    """Sample observed synonymous codons, then establish the new DNA reference.

    Complete proteins must begin with M; an optional terminal * is interpreted as
    termination and a stop is appended when absent. Fragment mode preserves an
    explicitly supplied terminal *. Internal stops/ambiguous amino acids fail.
    No optimization or expression claim is made by this initialization step.
    """
    if type(seed) is not int:
        raise ValueError("Reverse translation requires an integer random seed")
    profile = CodonProfile.model_validate(profile.model_dump())
    cds_mode = CDSMode(cds_mode)
    protein = "".join(protein.split()).upper()
    has_stop = protein.endswith("*")
    sense = protein[:-1] if has_stop else protein
    if not sense or set(sense) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError(
            "Protein must contain the 20 canonical amino acids, without internal stops"
        )
    if cds_mode == CDSMode.COMPLETE and not sense.startswith("M"):
        raise ValueError(
            "Complete protein must begin with M; choose explicit fragment mode otherwise"
        )
    table = CodonTable.unambiguous_dna_by_id[profile.translation_table]
    rng = random.Random(seed)
    codons: list[str] = []
    for index, amino in enumerate(sense):
        if index == 0 and cds_mode == CDSMode.COMPLETE:
            if "ATG" not in table.start_codons:
                raise ValueError("This genetic code has no supported ATG initiation codon")
            codons.append("ATG")
            continue
        choices = sorted(c for c, aa in table.forward_table.items() if aa == amino)
        codons.append(rng.choices(choices, weights=[profile.counts[c] for c in choices], k=1)[0])
    if cds_mode == CDSMode.COMPLETE or has_stop:
        # Stops were deliberately excluded from the reference's elongation counts.
        # Uniform seeded selection is an engineering policy, not inferred stop usage.
        codons.append(rng.choice(sorted(table.stop_codons)))
    sequence = "".join(codons)
    feature = Feature(
        id="cds",
        name=name,
        feature_type="CDS",
        start=0,
        end=len(sequence),
        strand=Strand.FORWARD,
        translation_table=profile.translation_table,
        cds_mode=cds_mode,
    )
    expected = sense + ("*" if has_stop and cds_mode == CDSMode.FRAGMENT else "")
    if translate_feature(sequence, feature) != expected:
        raise ValueError("Generated DNA failed independent amino-acid verification")
    construct = Construct(
        name=name,
        sequence=sequence,
        features=(feature,),
        random_seed=seed,
        host_profile=VersionedProfile(
            id=profile.id, version=profile.version, source=profile.source_url
        ),
        model_versions={"codon_profile_sha256": profile.fingerprint},
        input_format="amino_acid",
        import_warnings=("Seeded codon sampling initializes DNA; it does not predict expression",),
    )
    validate_candidate(construct, sequence).require_valid()
    return construct
