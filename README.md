# DNA Design Platform

A local, constraint-preserving DNA sequence-design library. **Phase 1 only**:
annotated constructs, biological mutation permissions, imports and independent final
validation. No codon optimizer, expression predictor, RBS model, server, or UI yet.

Start future sessions with [AGENTS.md](AGENTS.md), the complete
[product specification](docs/PRODUCT_SPEC.md), [architecture](docs/architecture.md),
[development plan](docs/development_plan.md), [decisions](docs/DECISIONS.md), and
[scientific provenance](docs/scientific_provenance.md). The full initial request is
preserved in [docs/prompts/INITIAL_CODEX_BOOTSTRAP.md](docs/prompts/INITIAL_CODEX_BOOTSTRAP.md).

## Setup and checks

Python >=3.12; development currently uses 3.12.3. Python 3.13 is preferred when
available and supported by the pinned dependencies. From the repository root:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-dev.lock
.venv/bin/python -m pip install --no-build-isolation -e '.[dev]'
.venv/bin/python -m pytest
.venv/bin/python -m compileall -q src tests
.venv/bin/ruff check .
.venv/bin/mypy src
```

If a Codespace interpreter lacks ensurepip, use `virtualenv --python python3.12 .venv`
from an interpreter with virtualenv installed. No administrator access to the client
machine is needed. Never commit .venv or private biological sequences.

## Example

```python
from dna_design import Construct, Feature, validate_candidate

cds = Feature(
    id="cds", name="synthetic", feature_type="CDS", start=0, end=9, strand=1, translation_table=11
)
original = Construct(sequence="ATGGAATAA", features=(cds,))
report = validate_candidate(original, "ATGGAGTAA")
report.require_valid()  # synonymous; protein and all original controls preserved
print(report.model_dump_json(indent=2))
```

Coordinates are 0-based half-open internally and 1-based inclusive in user-facing
input; call explicit conversion functions. Unannotated bases and functional non-coding
features default to LOCKED. FREE_MUTATION needs an explicit annotation. All overlapping
CDSs retain their proteins, even under a free-mutation annotation. Candidate annotations
cannot redefine original controls. Sequence data stay local; no external services run.

`parse_raw`, `parse_fasta`, and `parse_genbank` import one construct. Models serialize
to/from JSON via Pydantic. Complete CDSs require valid start/terminal stop codons;
fragments require explicit mode. Both strands and explicit genetic-code IDs work.
IUPAC motifs and circular-origin motif occurrences are supported; recognition-site
compatibility uses Biopython's restriction database.

Limits: equal-length substitutions, concrete A/C/G/T DNA, exact simple feature intervals.
Compound/fuzzy/remote GenBank locations and recoding exceptions are rejected explicitly.
`check_feasibility` detects provable conflicts; absence of conflicts does not prove a
solution exists. Unsupported features are never silently flattened. See the architecture
for source-annotation handling and provenance details.

Phase 2 will add sourced host profiles, codon metrics and seeded multiobjective
optimization after explicit authorization. Do not automatically begin it.
