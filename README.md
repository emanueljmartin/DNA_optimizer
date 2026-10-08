# DNA Design Platform

A local DNA design workbench with a browser interface, codon and synthesis design,
thermodynamic RBS prediction/libraries, Gibson fragmentation, reports and validated
exports. Original proteins, protected bases, motifs and masks remain authoritative.
Phases 1–8 have an implemented engineering checkpoint; experimental accuracy has
not been established. The reliability study and data choices will be planned with
the user after this build. Executed checks are in the development plan.

Start future sessions with [AGENTS.md](AGENTS.md), the complete
[product specification](docs/PRODUCT_SPEC.md), [architecture](docs/architecture.md),
[development plan](docs/development_plan.md), [decisions](docs/DECISIONS.md), and
[scientific provenance](docs/scientific_provenance.md). The full initial request is
preserved in [docs/prompts/INITIAL_CODEX_BOOTSTRAP.md](docs/prompts/INITIAL_CODEX_BOOTSTRAP.md).

## Setup and checks

Python >=3.12; the rebuilt Codespace currently uses 3.13.5. From the repository root:

```bash
bash scripts/setup_app.sh
.venv/bin/dna-design serve
```

Open **http://127.0.0.1:8000**. In Codespaces, use the private forwarded port 8000.
The installer preserves the working `.venv`, installs declared/pinned local API and
RBS dependencies, and builds the browser with Node >=24 and `npm ci`. It does not
publish or deploy. Rebuilt devcontainers run this installer automatically.
The packaged browser is included in locally built wheels; source edits require
`cd frontend && npm run build` to refresh it. Bundled React notices are included.
Use `bash scripts/build_package.sh` for a clean local source archive/wheel; it builds
in a fresh temporary directory so cached setuptools assets cannot enter the package.

Use **Load synthetic demo** to explore the workflow. Import your construct, curate
features in **Sequence**, run **Design** or **RBS**, plan **Assembly**, then inspect
**Reports**. Annotation edits explicitly establish new controls and clear prior
designs. **Save project** downloads a draft containing original/final DNA, settings,
history and locks; save it before closing or refreshing the browser. Final FASTA,
GenBank, JSON and ordering CSV exports run fresh validation and reject failures.
Reopening a project restores each stage's recorded controls, including separate
codon/RBS search budgets and complete Tm solution assumptions.
JSON retains the full project; GenBank retains typed feature IDs/mutability/masks
through documented custom qualifiers, but does not replace the JSON project history.

The solver exposes seeds, alternatives, preference weights and advanced parameters.
Presets weight codon adaptation and selected synthesis risks; unconfigured synthesis
criteria remain unassessed. Generic/public profiles are advisory. Custom codon JSON
must include counts and source provenance; high-expression CAI requires an explicitly
sourced high-expression reference. Other-host RBS prediction requires sourced anti-SD
parameters and explicit uncalibrated opt-in. Do not interpret TIR as measured proteins.

Run engineering checks from the root after activating the environment:

```bash
source .venv/bin/activate
which python
python --version
pytest -q
python -m compileall -q src tests
ruff check .
ruff format --check .
mypy src
python -m pip check
npm --prefix frontend run build
npm --prefix frontend run format:check
npm --prefix frontend test
```

For local thermodynamic RBS prediction, install the optional backend into the same
environment and rerun the suite:

```bash
python -m pip install -r requirements-rbs.lock
python -m pip install --no-build-isolation -e '.[dev,rbs]'
pytest -q
```

OSTIR 1.1.3 and ViennaRNA 2.7.2 run locally through the Python RNA API; no prediction
service or RNAfold CLI is required. Without this extra, prediction fails explicitly
and three real-backend Python integration tests skip; core and test-double tests still run.
Dependency licenses and model limitations are recorded in
[scientific provenance](docs/scientific_provenance.md#local-rbs-model-and-license-boundary).

Chromium browser tests require `cd frontend && npm exec -- playwright install --with-deps chromium`.
This installs test tools in the development container; normal application use needs
only your existing browser. A stale Yarn apt key blocked the stock dependency helper
in this Codespace; browser libraries were installed from signed Debian sources, with
no trust bypass or source-list changes. CI installs browser dependencies separately.

The Reliability page evaluates local CSV predictions using an explicit provenance
manifest and disjoint calibration/held-out groups. It does not download experimental
datasets or generate accuracy evidence automatically. It shows per-host/condition
metrics and evidence warnings and can download the full report. The
[proposed reliability study](docs/reliability_plan.md) lists the biological questions
and data decisions to review next. See
[benchmark workflow](docs/benchmark_workflow.md) for the schema, public source index,
units and leakage checks. CLI commands are `dna-design validate`, `dna-design export`
and `dna-design benchmark`; output commands refuse to overwrite existing files.

Setup creates an absent `.venv` with `python3 -m venv --copies .venv` and reuses a
working environment. An unusable environment stops setup before dependency installation;
move it aside before explicitly recreating it. After activation, `which python` must
point to `/workspaces/DNA_optimizer/.venv/bin/python` in this workspace. Never commit
.venv or private biological sequences. The devcontainer uses a prebuilt Python 3.13
image and Node 24; no Dockerfile or user override is needed. The running `vscode` user
and configuration were verified after the user's rebuild. A further container rebuild
and remote CI were not performed during this verification.

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

## Codon and synthesis design

```python
from dna_design import (
    SearchConfig,
    SynthesisProfile,
    load_codon_profile,
    optimize,
    repair_synthesis,
    reverse_translate,
    validate_candidate,
)

profile = load_codon_profile("ecoli")  # "pao1" is also bundled; no download
original = reverse_translate("M" + "E" * 20, profile, seed=7)
config = SearchConfig(seed=7, iterations=200, alternatives=2)
codon = optimize(original, config, profile)
if not codon.designs:
    raise RuntimeError(codon.validation.model_dump_json())
repair = repair_synthesis(
    original,
    SynthesisProfile(homopolymer_limit=5),
    config,
    codon_profile=profile,
    initial_sequence=codon.designs[0].sequence,
)
final = repair.search.designs[0]
validate_candidate(original, final.sequence).require_valid()
print(repair.after[0].model_dump_json(indent=2))
```

The example's homopolymer limit is an explicitly selected engineering policy.
Generic policies invent no thresholds; unconfigured criteria remain unassessed.
`twist-public` and `genscript-public` are incomplete snapshots of cited public advice,
with no vendor-acceptance claims. A valid biological design can retain synthesis risks;
inspect per-finding results. `require_policy()` also rejects unassessed criteria.

Built-in host data are genome-wide counts, so their score is labeled codon adaptation.
CAI is returned only for explicitly sourced high-expression reference profiles.
Zero-count codons retain weight zero; rare thresholds and GC targets are explicit.
Search supports weighted losses, hill climbing or simulated annealing, distinct
alternatives, and custom scorers. Search exhaustion does not prove infeasibility.
The default first-30-nt CDS protection comes from the product specification. Start
triplets stay fixed; additional designed locks and overlapping proteins always apply.

## RBS design and later codon changes

`RBSContext` identifies original RBS/CDS features, an upstream interval adjacent to
the oriented CDS start, and the downstream context length. `predict_rbs` returns
relative model TIR, free-energy components, spacing, a labeled unbound RNA structure,
versions and warnings. `design_rbs` respects the original N/ACGT mask and optionally
redesigns synonymous bases in the first requested CDS interval. `design_rbs_library`
selects distinct RBS sequences for explicit targets; inspect residual log-fold errors.

Always pass the **same original construct** to every stage. After `design_rbs`, carry
the matching `initiation_locks[index]` into any later codon or synthesis search:

```python
from dna_design import SearchConfig, optimize


def continue_codon_design(original, rbs_result, profile):
    index = rbs_result.selected_design_index
    return optimize(
        original,
        SearchConfig(seed=8, additional_locks=(rbs_result.initiation_locks[index],)),
        profile,
        initial_sequence=rbs_result.search.designs[index].sequence,
    )
```

These locks cover the entire supplied upstream and downstream prediction context,
including any larger joint design region. Never replace the original annotations
with a designed construct to establish new biological controls. E. coli uses OSTIR's
published calibration; this project has no held-out expression validation. PAO1 RBS
transfer requires explicitly sourced host parameters and `allow_uncalibrated=True`.
No default PAO1 calibration or experimental accuracy is claimed.
