# Testing strategy

Scientific invariants take priority over interface polish. Run `.venv/bin/python -m
pytest`, `.venv/bin/python -m compileall -q src tests`, `.venv/bin/ruff check .`, and
`.venv/bin/mypy src` before declaring Phase 1 complete.

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

Later phases add optimizer reproducibility/alternative diversity, synthesis inserted
defect-and-repair tests, RBS mask/target and protein integrity checks, Gibson exact
reconstruction/length tests, API tests, frontend interaction tests, and held-out public
benchmark reports. Training/calibration and independent validation remain separate;
no fabricated values or claimed experimental accuracy from unit tests.
