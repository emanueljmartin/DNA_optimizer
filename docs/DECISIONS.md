# Architectural and scientific decisions

## 2026-10-08 — Resume each stage faithfully and bind reports to evaluated inputs

Decision: recover the latest codon, RBS and assembly controls from existing project
request history, retain unexposed request restrictions/budgets on rerun, give codon
and RBS independent search controls, and expose every Tm solution parameter. Loading
another reference resets controls to defined application defaults. Reason: the
last-operation-only restoration lost earlier stage settings; an asymmetric Tm
solution was silently made symmetric and omitted salt parameters fell back to
defaults. Permanent browser regressions reproduce both failures and verify requests
sent to the backend. Host changes also clear restored codon references and previous
anti-SD parameters: a regression reproduced PAO1 selection leaving E. coli counts
active, and now checks the chosen profile recorded in the saved result. No biological
default, domain gate or saved schema was changed.

Decision: render backend benchmark metrics with evidence/independence warnings and
allow exact JSON download. Display results only for their exact evaluated CSV and
manifest text, including delayed responses. Reason: a scientist should be able to
review results without interpreting raw JSON, and edited inputs must not carry an
old result. Synthetic ratio/ordinal and delayed-response regressions check this;
they make no experimental accuracy claim. A proposed reliability study is documented
in reliability_plan.md; scientific thresholds, datasets and splits remain for review.

## 2026-10-07 — Phase boundaries and persistent specification

Decision: preserve complete B–AC requirements and full bootstrap; implement only
Phase 1 as a domain library. Reason: biological correctness precedes optimization.
Alternatives: combined optimizer/UI implementation. Implications: no expression or
RBS predictions and no host codon tables yet. Reference: PRODUCT_SPEC.md AA–AC.

## 2026-10-07 — Original construct is the validation authority

Decision: frozen typed controls, independently re-evaluated reference translations,
candidate metadata equality checks, equal-length substitutions only. Reason: a
candidate cannot redefine protected DNA/proteins. Alternatives: validate candidate
annotations alone; cache translations. Implications: indels and coordinate remapping
are deferred. Reference: PRODUCT_SPEC.md AB, X.

## 2026-10-07 — Conservative composed mutability

Decision: LOCKED dominates, unannotated bases lock, all overlapping CDSs constrain
protein identity; non-CDS SYNONYMOUS_ONLY without a CDS locks. Explicit masks govern
RBS designability. Reason: absence of annotation is not evidence of non-function.
Alternatives: last annotation wins, blanket free non-coding mutation. Implications:
GenBank gene annotations can lock CDSs until explicitly curated. Reference: E.

## 2026-10-07 — Translation and faithful import

Decision: Biopython/NCBI genetic codes, table 11 default for manually annotated
bacterial CDSs, INSDC default table 1 for imported CDSs without transl_table, explicit
fragment mode, reject unsupported exception/fuzzy/compound/remote locations and
context-dependent genetic codes. Reason: no faithful simple-interval representation
for these features. Alternatives: flatten joins or guess recoding. Implications:
some real-world records fail with explicit limitations. References:
[NCBI codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi),
[Biopython sequence tutorial](https://biopython.org/docs/latest/Tutorial/chapter_seq_objects.html).

Regression review corrected the initial GenBank omission handling to the documented
table 1 default. A dedicated test rejects a bacterial-only start when no alternate
table is specified; the public manual annotation default remains bacterial table 11.

## 2026-10-07 — Motifs and feasibility

Decision: literal nucleotide preservation at reference occurrences; both strands
and circular-origin matching; ERROR for hard AVOID violations; prove only explicit
conflicts without claiming general solver completeness. Reason: degenerate matches
must not weaken preservation. Alternatives: preserve motif match only, call every
failed candidate infeasible. Implications: unresolved synonymous repair feasibility
remains unknown; no optimization algorithm introduced.

## 2026-10-07 — Environment and dependency reproducibility

Decision: local Python 3.12 environment (3.13 not installed), exact dependency pins,
CI on 3.12/3.13; lightweight devcontainer after Phase 1 verification. Alternatives:
use system Python 3.14 immediately or download another interpreter. Implications:
supported Python floor remains 3.12; no client-machine modifications.

## 2026-10-07 — Acceptance checkpoint and feasibility boundary

Decision: mark Phase 1 complete after 123 passing tests, compilation/import checks,
Ruff and strict typing. Reason: executable evidence covers the specified invariants
and adversarial reference checks. Alternatives: stop at scaffolding, or implement
Phase 2 before validating controls. Implications: Phase 2 stays unstarted; CI and a
container rebuild remain unexecuted and are explicitly distinguished from local checks.
Feasibility additionally proves impossible REQUIRED motifs when no scoped window can
match under locked positions, even after relaxing all unlocked bases to arbitrary DNA.
Absence of a proof remains unknown, never a guarantee of biological feasibility.
References: development_plan.md checkpoint table; testing_strategy.md; Phase 1 tests.

## 2026-10-08 — Resume existing Phase 1 and preserve virtual environments

Decision: reverify the existing core on the rebuilt Codespace's Python 3.13.5, install
the already-declared development dependencies into the working `.venv`, and avoid
recreating it. Setup creates absent environments with `--copies`, checks the required
Python version and virtual-environment prefix, and stops with recovery instructions
for unusable existing paths. Reason: the initial collection failure came from an
external pytest selected when the local test tools were missing; repeating environment
creation can unnecessarily replace an existing interpreter. Alternatives: rebuild the
application, overwrite every `.venv`, or change biological behavior to address setup.
Implications: seven permanent setup regressions supplement the existing 123 scientific
tests; all 130 pass locally. Application semantics and Phase 1 representation limits
remain unchanged. Phase 2 requires separate authorization; remote CI and another
container rebuild remain unexecuted. Reference: development_plan.md resumed verification.

## 2026-10-08 — Explicit authorization for Phases 2–4 and DNA-only imports

Decision: honor the user's request to finish Phase 1 checks and continue through
codon optimization, synthesis and local RBS functionality, keeping later phases out
of scope. Phase 1 additionally rejects non-DNA GenBank molecule types even when their
letters fit A/C/G/T. Reason: molecule identity cannot be inferred from the alphabet;
unsupported representations must not be reinterpreted. Three regressions reproduce
RNA, mRNA and protein records. Reference:
[NCBI molecule-type documentation](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/#MoleculeTypeB).

## 2026-10-08 — Codon references and bounded seeded search

Decision: bundle reproducible aggregate genome-wide elongation counts from NCBI
NC_000913.3 and NC_002516.2 with dates, checksums, exclusions and profile fingerprints.
Keep genome-wide adaptation separate from high-expression-reference CAI; accept
explicit sourced custom/high-expression profiles. Zero counts remain zero, rare/GC
settings are explicit, and no high-expression dataset is inferred. Reason: host
composition is not a measured expression model. Sources/transforms are recorded in
scientific_provenance.md; both snapshots were rederived exactly from fresh downloads.

Decision: a configurable seeded hill-climbing/annealing interface scores only valid
designs and independently rechecks all original proteins, overlaps, locks and motifs.
Only unresolved REQUIRED/AVOID states can be traversed during motif repair. Report
bounded-search exhaustion separately from Phase 1's proven conflicts. Distinct returned
alternatives include objective components and per-base changes. Protect the first
30 nt by default and independently gate configured mutation boundaries. Generate
synonyms from valid starting DNA while retaining original proteins, which avoids
incorrectly freezing a partially locked codon that was redesigned in an earlier stage.
Reason: proposal generators and score caches cannot be the final biological authority.

Decision: amino-acid input initializes a new reference with seed-controlled sampling
of observed counts, explicit complete/fragment semantics, ATG initiation for complete
M-starting proteins and uniform stop selection. Fresh translation verifies the input
protein. This is initialization, not a greedy codon-replacement optimizer or an
expression claim. Reference: PRODUCT_SPEC.md H–K, AA; Sharp and Li (1987).

## 2026-10-08 — Transparent synthesis policies and residual risk

Decision: generic profiles invent no thresholds. Callers select explicit GC,
homopolymer, repeat and entropy policies; public GenScript/Twist snapshots carry
sources, dates, product scope and omitted criteria. Provide per-finding before/after
reports and same-solver repair with original-reference validation. A biologically
valid result can retain risks; policy certification fails on residual/unassessed
criteria. Reason: public guidance is incomplete and proprietary acceptance cannot be
reconstructed. Alternatives: opaque synthetic difficulty scores or invented vendor
limits. Implications: vendor acceptance and repeat melting-temperature modeling
remain unclaimed. References: scientific_provenance.md, PRODUCT_SPEC.md O.

## 2026-10-08 — Optional local thermodynamic RBS model and complete-context locks

Decision: integrate pinned OSTIR 1.1.3/ViennaRNA 2.7.2 behind a replaceable local
protocol, with lazy imports, explicit missing-backend errors, energy/spacing/version
outputs and a labeled unbound full-context structure. Published E. coli parameters
are consumed from the upstream model; PAO1/custom hosts require sourced anti-SD and
explicit uncalibrated opt-in. No project experimental accuracy or absolute-expression
claim is made. OSTIR's obsolete CLI-presence import warning is ignored only by its
exact text after checking the supported versions; execution uses the real RNA Python
API and is covered by real integration tests.

Decision: masked target search, optional synonymous joint CDS changes and target
libraries share the original-reference validation gate. Libraries require distinct
RBS sequences and retain their search history/selected index. Final predictions run
afresh. Later codon/synthesis optimization carries an AdditionalLock covering every
supplied upstream/downstream prediction base and any larger joint region. Regression
review caught that a 30-nt lock could allow changes in a 35-nt prediction context;
both-strand tests now verify complete-context protection. Reason: keeping only the
RBS motif or a shorter CDS interval does not preserve the model prediction.

Dependency licenses were reviewed from installed distributions and primary sources:
OSTIR GPLv3 and ViennaRNA's attribution/redistribution/commercial-product conditions.
Keep dependencies optional and do not vendor their code/binaries. This authorizes
local research execution; commercial redistribution needs a separate license decision.
References: scientific_provenance.md; PRODUCT_SPEC.md L; OSTIR (2021).

## 2026-10-08 — Full-tool continuation and assembly

The user explicitly authorized Phases 5–8 and deferred the experimental reliability
study/data choices until after the full build. Engineering tests continue throughout.
No automatic commit, push, publication or deployment is authorized.

Decision: deterministic dynamic programming starts at the geometric minimum fragment
count and covers all eligible boundaries, considering overlap uniqueness, composition,
GC/Tm preferences, protected DNA, synthesis policy and balance. Circular closure adds
one shared overlap to the ordering geometry. The fixed origin is never rotated silently.
Tm uses Biopython DNA_NN3/Allawi–SantaLucia under displayed solution assumptions, not
invented reaction concentrations. Independent reconstruction and fresh geometry/metric
checks gate ordering export; an ambiguous diagnostic plan cannot be exported as valid.
Reason: ordering length includes overlap and physical reconstruction is a hard invariant.

## 2026-10-08 — Local browser, project snapshots and exports

Decision: add optional pinned FastAPI/Uvicorn and React/TypeScript/Vite around the
existing domain library. All science remains in Python. One background job worker,
bounded queue/retention, cooperative progress/cancellation and local files make job
state explicit. Use same-origin API calls, local-client headers, no telemetry exporter
or external assets/prediction calls, and no sequence access logs. The launcher defaults
to loopback; no hosted deployment is performed. Interactive size/budget limits are
engineering resource limits rather than weakened biological constraints.

Draft JSON can save unfinished work with its validation findings. Every final format
freshly rechecks original proteins, locks, motifs, masks, carried context locks and any
assembly. JSON records source/final DNA, settings/history, complete selected codon
profile, fingerprints and software versions; GenBank preserves metadata and explicit
validated custom qualifiers for feature IDs, mutability, masks and fragment CDSs.
Annotation edits preserve notes/qualifiers, establish controls from original DNA and
clear earlier design state. Combined codon/synthesis weights are transparent preferences.
Reason: convenient saving must not imply an invalid design passed export validation.

Decision: retain built browser assets for local wheel use without Node at runtime;
keep editable frontend source and npm lock for reproducible rebuilds. Include bundled
React MIT notices. Current Starlette requires httpx2 for warning-free test clients;
pin it and its resolved dependencies. The installer preserves the existing Python
interpreter, installs declared extras, runs npm ci/build and does not deploy. A stale
Yarn apt key affected browser test dependencies; use signed Debian sources without
changing trust configuration. Application use itself needs no test Chromium install.

## 2026-10-08 — Auditable evaluation and bounded performance

Decision: implement a versioned local mapped-CSV loader, primary-source candidate index,
mandatory provenance/license/checksum/units/split manifest, tied-rank/linear/log metrics,
and fold errors only for a defensible ratio scale. Calibration is explicitly separate
and fitted per host/condition using positive calibration data. Reject sequence/group
and declared model-training leakage; unknown upstream overlap remains unverified.
Published ordinal Boël scores cannot receive fold-error interpretation. Do not bundle
experimental records, invent source mappings or claim performance from synthetic data.
This is evaluator infrastructure; the experimental study remains a post-build decision.

Decision: keep only the best requested distinct valid search alternatives and cache
at most 512 oriented RBS contexts. Preserve proposal trajectories/random draws and
fresh final predictions/validation. Test equivalence to retaining all proposals and
measure actual synthetic multi-kb performance with tracemalloc. Reason: retaining
every full construct grows memory with the search budget; context reuse is sufficient
for prediction caching. Timings are engineering measurements, not accuracy benchmarks.

Packaging regression: cached setuptools `build/lib` retained an obsolete browser
bundle after a source rebuild. The supported package script now builds a source archive
then a wheel in a fresh temporary directory. A permanent test requires packaged asset
files to match exactly the index manifest; source tests/fixtures and dependency notices
are explicitly included. Existing workspace/build artifacts are preserved.
