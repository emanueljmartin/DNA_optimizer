# Testing strategy

Scientific invariants take priority over interface polish. Run `.venv/bin/python -m
pytest`, `.venv/bin/python -m compileall -q src tests`, `.venv/bin/ruff check .`, and
`.venv/bin/mypy src` before declaring a milestone complete. Also run Ruff formatting
checks, pip dependency checks and `git diff --check`.

Unit tests exercise concrete DNA validation, models, bounds, coordinate conversions,
both CDS strands, frame offsets, complete/fragment semantics, table-specific codons,
motif scopes/orientations, masks and composed permissions. Integration tests import
annotated GenBank then independently validate synonymous and invalid candidates.
Adversarial tests alter candidate metadata/source sequence to ensure it cannot change
the biological baseline. Infeasibility tests distinguish proven conflicts from an
invalid but potentially repairable candidate. Synthetic fixtures never contain
confidential/unpublished DNA.

Hypothesis checks coordinate round trips, DNA serialization/import round trips,
locked-region substitutions, and synonymous protein preservation across randomized
codon choices. Tests assert invariants rather than a unique synonymous sequence.
Every discovered bug receives a permanent regression test. Regression checks include
off-by-one bounds, same-position circular sites, overlapping genes, trailing whitespace,
and duplicate/ambiguous identifiers. Never weaken tests to obtain a passing suite.

Devcontainer setup regressions exercise copied-interpreter creation, reuse without
overwriting a working environment, recovery errors for unusable interpreters/files/
symlinks, and stopping after dependency-install failure. A real non-virtual-environment
interpreter is rejected even with Python optimization enabled. Installer commands are
otherwise simulated in temporary directories so these tests need no package downloads
and cannot alter the developer's `.venv`.

Phase 2 tests exercise sourced-profile schemas, reference exclusions, zero-count
weights, conditional CAI, rare/GC metrics, explicit objective parameters, seeded
alternative diversity, motif repair, overlap/lock gates, structured invalid references,
and budget exhaustion. Hypothesis exercises both CDS orientations and protein-input
reverse translation. Partial locks use valid starting DNA while preserving original
proteins. Fresh search validation rejects changes outside selected mutation regions
or inside configured initiation protection.

Phase 3 synthetic stress cases insert extreme global/local GC, homopolymers, direct
and inverted repeats and low entropy. Detection and seeded repair preserve original
motifs, CDSs and locked regions; immutable residual risk and unassessed criteria remain
visible. Public-guidance values are checked separately from heuristic execution.

Phase 4 uses an explicitly fictional counting test double to isolate masks, target
losses, joint CDS boundaries, original-validation gates, target-library diversity and
carried context locks. It never represents a biological model. Three optional real
OSTIR/ViennaRNA integration cases execute oriented prediction, report model outputs
and perform constrained target design. Core-only installations skip those three;
the optional-backend CI matrix must execute them. Full local tests use the real extra.
Cross-phase pipeline tests retain the SAME original construct throughout codon,
synthesis, RBS and subsequent codon design. Final prediction and validation bypass
caches, and both forward/reverse initiation contexts remain protected.

Phase 5 tests verify linear/circular ordered lengths, overlap geometry, minimum counts,
uniqueness, protected-boundary policies, exact reconstruction and fresh rejection of
forged DNA/coordinates/metrics. Hypothesis checks randomized tiling without missing or
duplicated bases. A public Biopython example checks the Tm model under explicit conditions.

API/project/CLI tests cover input defaults, UI coordinates, annotation metadata,
oriented codon logs, seed replay, alternatives, recorded profiles, job cancellation,
privacy/body limits, failed final exports, draft restoration, GenBank round trips and
refusal to overwrite output files. Chromium tests exercise real browser curation,
selection, local-only network requests, optimization, real OSTIR prediction, assembly,
downloads, input errors and mobile layout.
Browser regressions also check multi-stage history restoration, independent codon/RBS
budgets, full Tm solution replay, preservation of unexposed search restrictions, and
clearing controls when loading a different construct or selecting a different host.
The host case checks the exact scoring profile saved by the backend. Both-strand pipeline cases now
continue through unique-overlap Gibson reconstruction, draft reload and all exports.

Benchmark tests use explicitly synthetic observations to verify tied-rank correlation,
constant/zero handling, exact fold arithmetic, separate calibration, per-host/condition
reports, checksums, ordinal restrictions and sequence/group/training leakage rejection.
They establish evaluator arithmetic, never experimental model accuracy. Public datasets
require explicit source/units/column/split preparation and license provenance. No public
experimental benchmark has been executed here. The study/data choices are deferred to
the post-build discussion with the user.
Browser evaluation cases exercise actual backend synthetic ratio/ordinal reports,
warnings, unassessed metrics and exact report downloads; an intentionally delayed
response verifies that edited inputs cannot acquire an earlier evaluation's results.
The draft study and biological review questions are in reliability_plan.md.

Execution tests verify cancellation preserves source DNA and progress does not change
seeded results. Bounded alternative retention returns the same best designs as retaining
all proposals. `scripts/profile_local.py` records actual synthetic multi-kb timings and
Python allocations. No hard machine-dependent speed assertion is placed in unit tests.
CI covers Python 3.12/3.13 with API and optional RBS, plus Node 24 Chromium workflows;
local executable evidence does not imply remote CI has executed.
