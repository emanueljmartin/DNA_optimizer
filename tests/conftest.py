from io import StringIO

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from dna_design import Construct, Feature, Strand


@pytest.fixture
def cds() -> Feature:
    return Feature(id="cds", name="synthetic CDS", feature_type="CDS", start=0, end=9, strand=1)


@pytest.fixture
def coding_construct(cds: Feature) -> Construct:
    return Construct(name="synthetic", sequence="ATGGAATAA", features=(cds,))


def genbank_text(
    *,
    sequence: str = "ATGGAATAA",
    strand: int = 1,
    location: SimpleLocation | None = None,
    qualifiers: dict | None = None,
    extra_features: tuple[SeqFeature, ...] = (),
    topology: str = "linear",
) -> str:
    dna = str(Seq(sequence).reverse_complement()) if strand == Strand.REVERSE else sequence
    record = SeqRecord(Seq(dna), id="synthetic.1", name="synthetic", description="Public test DNA")
    record.annotations = {"molecule_type": "DNA", "topology": topology}
    record.features = [
        SeqFeature(
            SimpleLocation(0, len(dna)), type="source", qualifiers={"organism": ["synthetic"]}
        ),
        SeqFeature(
            location or SimpleLocation(0, len(dna), strand=strand),
            type="CDS",
            qualifiers=qualifiers
            or {"gene": ["test"], "transl_table": ["11"], "translation": ["ME"]},
        ),
        *extra_features,
    ]
    output = StringIO()
    SeqIO.write(record, output, "genbank")
    return output.getvalue()
