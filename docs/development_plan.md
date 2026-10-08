# Development milestones

The user authorized the remaining Phases 5–8 on 2026-10-08: build the full local
tool, then jointly plan independent reliability tests and data. Implement sequential
milestones, retaining biological gates and running engineering tests at every stage.
Complete requirements are in PRODUCT_SPEC.md; executable evidence is recorded below.

| Milestone | Scope / implementation targets | Dependencies | Tests / definition of done | Intentionally deferred |
| --- | --- | --- | --- | --- |
| 1: core sequence engine | Models, coordinates, strands/codes, mutability, locks, motifs, parsers, independent validation, proven conflicts | Persistent documentation; Python >=3.12; Pydantic/Biopython | Complete pytest including properties, import/compilation, configured lint/types; all AC invariants covered; limitations documented | Optimizers, host usage, server, UI, fragment design |
| 2: codon optimization | Versioned host profiles, provenance-backed usage/CAI, rare codons, GC, modular multiobjective solver, seeded proposals, alternatives | Accepted Phase 1; explicit new authorization; sourced PAO1/MG1655 data | Protein/lock/motif integrity, repeatable seeds, objective behavior, genuine alternative diversity; fresh final validation | Vendor rules, RBS claims, UI |
| 3: synthesis | Transparent generic GC/repeat/homopolymer/complexity scanner and repair; versioned public/custom profiles | Phase 2; defensible profile sources | Inserted defects detected and improved without protein/lock/motif regression; per-finding report | Proprietary acceptance claims |
| 4: RBS | Backend abstraction, host settings, predict/design, masks, joint CDS, libraries | Phases 1–3; model/license review; host validation | Fixed masks, synonymous CDS, target distance behavior; model provenance and PAO1 limits | Unvalidated absolute-expression claims |
| 5: Gibson | Minimum fragment count, breakpoints, overlap metrics, exact reconstruction and ordering export | Validated construct; public Tm/quality model | Ordered lengths include overlap; exact complete reconstruction; uniqueness/repeat diagnostics | Other assembly methods |
| 6: browser UI/API | FastAPI; React/TypeScript; input/viewer/annotations/settings/RBS/fragments/reports/help/exports | Stable domain contracts; validated scientific modules | API/schema, UI interactions, responsiveness, validation gate on export; local-default privacy | Frontend biological calculators |
| 7: benchmarks | Versioned public loaders, held-out quantitative reports | Licensed datasets; model versions; explicit split policy | Reproducible metrics; no sequence-equality criterion; calibration/validation disjoint | Fabricated performance values |
| 8: polish | Profiling/caching, progress/cancellation, errors, packaging/installers/usability | Functional tested previous phases | Performance budgets, reproducible setup, regression suite, docs continuity | Unrequested scientific scope |

## Historical Phase 1 checkpoint — 2026-10-07

**Phase 1 complete, 2026-10-07.** Bootstrap documents were created before application
functionality. The implementation is a tested local domain library; Phase 2 has not
started and requires explicit authorization.

Executed in `/workspaces/DNA_optimizer` on branch `phase-1-core`, using local Python
3.12.3 (system default 3.14.2; Python 3.13 not installed). Startup confirmed a clean
repository, Node 24.21.0 and npm 11.19.0. The supplied expected path
`/workspaces/dna-design-platform` differs from this actual, user-selected workspace;
work stayed in the provided repository. No client/host PC changes were made.

| Executed command / check | Evidence |
| --- | --- |
| `.venv/bin/python -m pytest` | 123 passed, 0 failed, 0 skipped; includes Hypothesis properties |
| `.venv/bin/python -m compileall -q src tests` | Exit 0 |
| `.venv/bin/ruff check .` | All checks passed |
| `.venv/bin/ruff format --check .` | 25 files already formatted |
| `.venv/bin/mypy src` | No issues in 9 source files |
| `.venv/bin/python -m pip check` | No broken requirements |
| Import and integration sanity | All domain modules imported; synonymous validation of a 3,006-nt synthetic CDS passed |
| Specification preservation | Full prompt equals supplied text after newline normalization; entire B–AC block is present verbatim in PRODUCT_SPEC.md |
| Container JSON / shell syntax | JSON parsed; setup script checked with bash syntax validation |
| `git status`, `git diff --stat`, `git diff --check` | Reviewed before the checkpoint; edits/new files remain for user review; no merge/push performed |

Pinned runtime/dev/build dependencies are in pyproject.toml and requirements-dev.lock.
CI config covers Python 3.12 and 3.13. A lightweight devcontainer specifies Python 3.13
and Node 24. CI execution and an actual container rebuild were not performed in this
session; local executable evidence is for Python 3.12.3 only. Container images/features
are version-tagged rather than immutable digest pins. Recreate .venv when switching
Python versions. No experimental expression/RBS benchmarks have been run.

Acceptance evidence maps to tests/test_models_coordinates.py (schema/coordinates/
round trips), test_translation_permissions.py (protein, strand, tables, overlap,
locks, masks and properties), test_motifs.py (modes, IUPAC, scopes, circular origin,
restriction names, conflicts), test_parsing_validation.py (imports, reference authority,
all constraints and structured reports), and test_regressions.py (import semantics,
literal multi-CDS fixture, serializer bypasses, required feasibility and enzyme labels).

Known representation limits: substitutions only; concrete DNA; exact single intervals;
no recoding/compound features; recognition motifs rather than enzyme cut simulation;
no general satisfiability solver. See architecture.md for the scientific implications.

## Historical resumed Phase 1 verification — 2026-10-08

Phase 1 remains complete; Phase 2 has not started. Work resumed in the existing
`/workspaces/DNA_optimizer` repository on `phase-1-core`, preserving all prior edits
and the existing commit. The implementation and 123 existing tests cover the Phase 1
acceptance requirements; no missing scientific component or failing application test
was identified. No scientific/application code was changed during this verification.

The rebuilt environment runs Python 3.13.5 (`/usr/local/bin/python3`). Its existing
`.venv` runs successfully, and activation selects
`/workspaces/DNA_optimizer/.venv/bin/python`. Runtime packages were already present,
but local development tools were missing: the initial `pytest -q` used an external
pytest and failed collection with `ModuleNotFoundError: Bio`. Installing
`requirements-dev.lock` and the editable `.[dev]` project fixed the environment.
The pip upgrade check found the installed version already current. The environment
was not deleted or recreated.

The devcontainer JSON parses successfully, uses a prebuilt Python image (so has no
Dockerfile), and contains no `remoteUser` or `containerUser`. `getent passwd vscode`
and `id` confirm a valid running user. Setup now creates absent environments with
`--copies`, checks Python >=3.12 and the intended virtual-environment prefix, preserves
existing environments, and stops before installation if an existing path is unusable.
Seven permanent regressions cover environment creation/reuse, unusable interpreters
and paths, installation failure, and guards under Python optimization. A real setup
rerun also preserved the existing interpreter's inode and modification timestamp.

| Executed command / check | Evidence |
| --- | --- |
| Activated `.venv`, `pytest -q` | 130 passed, 0 failed, 0 skipped; includes original Hypothesis invariants and seven setup regressions |
| `python -m compileall -q src tests` | Exit 0 |
| `ruff check .` / `ruff format --check .` | All checks passed; 26 files formatted |
| `mypy src` | No issues in 9 source files |
| `python -m pip check` | No broken requirements |
| Devcontainer JSON / `bash -n .devcontainer/setup.sh` | Passed |

No further container rebuild or remote CI run was performed. The sandbox command
runner could not create a user namespace; authorized checks ran outside that runner
without changing host security settings. No commit, push, reset, or merge was performed.
Before Phase 2: review this checkpoint, verify configured CI on 3.12/3.13, and obtain
explicit Phase 2 authorization. Host datasets still require documented provenance.

## Authorized continuation — 2026-10-08

The user requested completion of Phase 1 checks followed by the next few phases;
working scope is Phases 2–4. A follow-up import audit reproduced RNA/protein GenBank
records with A/C/G/T letters being silently treated as DNA. The parser now requires
an explicit DNA molecule type. Three permanent regressions were added; all 133 Phase 1
tests pass, as do Ruff, formatting, typing and compilation. The Codespace creation
log confirms successful container startup followed by the stale-interpreter setup
failure already repaired. This completes the additional Phase 1 check.

Implementation plan for the authorized continuation:

* Phase 2: derive versioned genome-wide codon counts from public reference records
  NC_000913.3 and NC_002516.2 with checksums and explicit exclusions; keep genome-wide
  adaptation distinct from high-expression-reference CAI. Build seeded, bounded local
  search with independent hard-constraint gates, configurable objectives and distinct
  alternatives. Report search exhaustion separately from proven infeasibility.
* Phase 3: transparent composition/repeat/complexity findings with explicit policies;
  no invented vendor thresholds. Repair through the same hard-constraint gate and
  validate against the original construct, including when starting from Phase 2 output.
* Phase 4: optional local OSTIR/ViennaRNA backend, dependency/license/source review,
  explicit host limitations, masked RBS design and optional joint initiation-CDS changes.
  Targets are relative model outputs; fresh final predictions and original-reference
  validation are mandatory. Carry explicit initiation locks into later codon redesign.

The entire original-reference validation suite remains required throughout. No sequence
data are sent to prediction services; only public source assets are downloaded.

## Historical checkpoint — Phases 1–4, 2026-10-08

Phase 1 was independently rechecked before continuation. Phases 2–4 now provide a
tested local domain checkpoint on the existing `phase-1-core` branch. No Phase 5
implementation, server, UI, commit, push or merge was started.

| Phase | Implemented evidence | Practical limits |
| --- | --- | --- |
| 1 | Original core invariants, parsers, genetic codes, motifs, composed mutability and fresh validation; explicit DNA-only GenBank guard | Concrete DNA, substitutions, exact single intervals; unsupported biology rejected |
| 2 | Versioned MG1655/PAO1 genome-wide counts; local reproducible derivation; custom/high-expression profiles; adaptation/conditional CAI, rare/GC metrics; seeded weighted solver, motif repair, alternatives, changes and original-reference fingerprints; amino-acid reverse translation | No fabricated high-expression built-in or expression accuracy; bounded local search is not a complete solver |
| 3 | Explicit versioned generic/custom/public policies; global/local GC and spread, homopolymers, direct/inverted repeats, entropy; seeded synonymous/free-region repair and per-finding before/after reports | Advisory public approximations; residual and unassessed risks remain visible; no proprietary acceptance claims |
| 4 | Optional local OSTIR/ViennaRNA protocol; sourced E. coli defaults; oriented predictions/energies/spacing/unbound structure; masks, relative targets, joint synonymous CDS design and distinct-RBS libraries; fresh final predictions and carried full-context locks | No held-out project validation or default PAO1 calibration; other hosts require sourced parameters and explicit uncalibrated opt-in |

Regression review found and covered invalid-reference serialization, missing-profile
rare objectives, partially locked starting codons, fresh search-region protection,
invalid-CDS derivation exclusions, and RBS locks that previously covered 30 nt while
prediction consumed 35 nt. Protection now covers the complete supplied prediction
context and larger joint regions. Cross-phase tests keep the same original construct
through all stages and subsequent codon changes, on both strands.

Local evidence uses the existing `.venv` Python 3.13.5 with the optional `rbs` extra.
The full suite currently has 190 passing tests, 0 failures and 0 skips. Three real
OSTIR/ViennaRNA integration cases run locally; core-only installations skip these
explicitly while executing mask/design protocol tests with a fictional test double.
Synthetic integrity tests are not experimental expression benchmarks. CI is configured
for Python 3.12/3.13 with core and optional extras; remote CI has not been run here.

Both packaged codon profiles were rederived from freshly downloaded public NCBI
records: full profile content and source checksums match, with 4,294 MG1655 and 5,570
PAO1 qualifying CDSs. Wheel packaging includes both aggregate JSON datasets; no genome
record or dependency binary is bundled in the platform wheel.

| Executed command / check | Current evidence |
| --- | --- |
| Activated `.venv`, `which python`, `python --version` | `/workspaces/DNA_optimizer/.venv/bin/python`, Python 3.13.5 |
| Activated `.venv`, `pytest -q` | 190 passed, 0 failed, 0 skipped, including all three real local RBS integration cases |
| Fresh copied-interpreter environment; development lock plus built platform wheel; pytest | 187 passed, 0 failed, 3 expected optional-backend skips; no editable source install |
| `python -m compileall -q src tests scripts` | Exit 0 |
| `ruff check .` / `ruff format --check .` | Passed; 39 files formatted, including README examples |
| `mypy src` | No issues in 15 source files |
| `python -m pip check`, both environments | No broken requirements |
| `python -m pip wheel --no-deps --no-build-isolation` / wheel inspection | Built successfully; includes both packaged profile JSON files |
| Public-source rederivation | Source checksums and complete profile contents match both built-ins |
| README Python examples | All executed successfully through the public API |
| Devcontainer JSON / `bash -n .devcontainer/setup.sh` | Passed; no invalid user override |
| `git diff --check`, `git status`, `git diff --stat` | Reviewed; changes remain unstaged/untracked for user review, preserving prior work |

The isolated core-only environment and wheel artifacts are outside the repository in
temporary directories. They do not replace `.venv`, bundle dependency binaries, or
change application biological controls. A further container rebuild and remote CI
run were not performed. Local Python 3.12 execution was not available in this resumed
environment; the configured matrix still requires remote verification on that version.

At that checkpoint, Phase 5 required new authorization; the user subsequently
authorized the full tool through Phase 8, as recorded below.

## Current full-tool checkpoint — Phases 1–8, 2026-10-08

Resumed the existing project and completed an engineering checkpoint through Phase 8
on `phase-1-core`. Original Phase 1–4 modules/tests were retained and extended, not
recreated. No automatic commit, push, merge, reset, publication or deployment occurred.
Experimental model accuracy has not been established; the user explicitly deferred
the reliability study/data/split choices until after the build.

| Phase | Completed implementation | Evidence / practical limits |
| --- | --- | --- |
| 1–4 | Existing biological core, sourced codon profiles/search, synthesis, local RBS and initiation locks retained | Full original invariant suite; real OSTIR/ViennaRNA integration; original reference remains authoritative |
| 5 | Minimum feasible Gibson tiling, boundary scoring, full overlap metrics, circular closure, fragment synthesis findings and ordering CSV | Exact reconstruction/properties; fresh DNA/geometry/Tm/metadata/uniqueness protection; no wet-lab yield claim |
| 6 | Local FastAPI + React/TypeScript interface; input/viewer/selection/annotations/motifs; seeded preferences/custom references; RBS/library; assembly; reports/help/drafts/final exports | API and Chromium workflows; backend-only science; JSON provenance/replay; final export revalidation |
| 7 | Versioned mapped public/local CSV loader, primary-source discovery index, required manifest/schema, tied-rank/linear/log metrics and defensible fold reports | Synthetic arithmetic/leakage tests only; exact public raw-column/license/split preparation and experimental evaluation remain for the study |
| 8 | Cooperative background progress/cancellation, bounded alternatives/context cache, installer, local package builder, browser notices, profiling, CI configuration and docs | Seed invariance, real cancellation, installer rerun, clean wheel/source archive and isolated install; CI/rebuild not executed remotely |

Regression coverage includes metadata-preserving annotation edits, oriented codon logs,
GenBank IDs/mutability/masks and reverse fragment interpretation, carried RBS starting
DNA in libraries, stale/forged fragment metrics, context-protection policies, CSV units
and split leakage, and stale compiled browser assets in cached wheel builds. The clean
package builder works from a fresh source archive, preserving previous workspace work.

| Executed check | Result |
| --- | --- |
| `.venv/bin/pytest -q` | **241 passed, 0 failed, 0 skipped** on Python 3.13.5, with real optional RBS/API dependencies |
| `npm run build` | Strict TypeScript and production Vite build passed; all assets local |
| `npm run format:check`, `npm test` | Formatting passed; **7 Chromium tests passed**, including real RBS prediction, restoration, cancellation and mobile use |
| Ruff checks/format, `mypy src` | Passed; 56 Python files formatted, 23 source files typed |
| Compilation, pip check, shell syntax, devcontainer JSON, `git diff --check` | Passed |
| `bash scripts/setup_app.sh` | Passed against the existing `.venv`; no environment deletion/recreation |
| Local package build / isolated installed wheel | Browser, notices, profiles and optional API status verified; without RBS: **238 passed, 3 intentional real-backend skips**; after installing pinned RBS extras: **241 passed, 0 failed, 0 skipped** |

Browser test dependencies initially failed because Chromium libraries were absent and
the existing Yarn apt source had an unavailable signing key. Required libraries were
installed inside the container using signed Debian sources, without bypassing signatures
or changing apt source/trust configuration. Current Starlette's test client requires
httpx2; the declared API extra/lock now pin it. The container has a valid `vscode` user,
no remote/container user override, and uses the full installer at post-create. Activated
Python resolves to `/workspaces/DNA_optimizer/.venv/bin/python`. No further rebuild ran.

Actual synthetic profiling on this machine (seed 7, no allocation tracing): 3,006-base
analysis took 0.105 s, 100 codon proposals took 3.206 s, and 6,000-base Gibson planning
(2,000-base ordered limit) took 0.838 s. A repeat during concurrent test/install load
took 0.300 / 8.280 / 1.201 s. Allocation tracing measured peaks of 3.599 / 4.331 / 17.101
MB of Python allocations and inflated timings to 0.772 / 32.457 / 7.819 s; native memory
is not included. Reproduce with `scripts/profile_local.py`, optionally `--memory`.
These are engineering measurements for these budgets, not biological benchmarks.

Next: review the scientist's workflow and select defensible public/user datasets and
independent host/condition/group splits, verify upstream calibration overlap, then run
the quantitative reliability study. Biological accuracy, absolute expression, proprietary
vendor acceptance, PAO1 RBS calibration and wet-lab assembly success remain unclaimed.
Implementation details, provenance and operational limits are in architecture.md,
scientific_provenance.md, testing_strategy.md and benchmark_workflow.md.
Known future scientific work includes sourced high-expression sets where justified,
held-out RBS/expression benchmarking with a separate calibration split, experimentally
defensible PAO1 RBS calibration, and a supported repeat-Tm model. Commercial packaging
requires a separate dependency license decision. None of these deferred items is
represented as completed or replaced by a fabricated parameter/benchmark.

## Current continuation checkpoint — workflow replay and study preparation, 2026-10-08

Continued the existing full-tool implementation. A browser audit reproduced two
scientific-control replay defects: last-operation-only restoration lost earlier
codon/RBS settings, and Tm restoration discarded unequal strand concentrations and
other salt assumptions. The browser now restores each stage from existing history,
keeps codon/RBS search budgets separate, resets controls when adopting another
reference, preserves unexposed restrictions/budgets on rerun, and exposes/restores
all Tm conditions. Saved schema and original biological gates remain unchanged.

Reliability now displays backend-generated counts, per-host/condition metrics,
evidence and independence warnings, and downloads the unchanged JSON report. A
report belongs to its exact CSV/manifest inputs; editing them hides earlier results,
including an earlier request still in flight. Ratio/ordinal arithmetic and stale
response regressions use clearly synthetic observations, never accuracy evidence.

The existing forward/reverse pipeline tests now run through unique-overlap Gibson
planning, exact reconstruction, draft restoration and all four final exports.
Seven new browser cases exercise multi-stage control replay, cross-project resets,
complete Tm replay, ratio/ordinal reports/downloads and delayed evaluation responses,
plus host changes clearing a restored codon reference and previous anti-SD parameters.
The host regression checks the profile recorded in the saved backend result.

Executed checks: **241 Python tests passed, 0 failed, 0 skipped**; **14 Chromium
tests passed**; strict TypeScript/production build, frontend formatting, Ruff
checks/formatting, mypy (23 source files), compilation, pip dependency checks, shell
syntax and Git whitespace checks passed. The existing Python 3.13.5 `.venv` and
devcontainer remain healthy; no environment recreation or further rebuild occurred.
The updated source archive/wheel built cleanly; a fresh temporary wheel installation
using the existing dependency environment served the current local browser assets
and passed all three selected packaging/both-strand pipeline checks. Neither the
repository environment nor previous temporary package artifacts were replaced.
Changes are saved in the workspace, **not committed or pushed**; `phase-1-core`
and the original commit remain intact.

[reliability_plan.md](reliability_plan.md) is a concrete proposal for biological
review: separate integrity, expression association/improvement, RBS target accuracy,
synthesis and assembly evidence; freeze provenance/context/groups/splits before
evaluation; agree acceptance criteria and uncertainty. The primary Boël source was
checked for host/assay applicability. No experimental dataset, split, threshold,
PAO1 calibration or performance result has been selected or claimed. Next: review
the intended first use/host/assay and approve a defensible data preparation protocol.
