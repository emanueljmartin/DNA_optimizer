# DNA Design Platform: repository instructions

Read `docs/PRODUCT_SPEC.md`, `docs/architecture.md`, `docs/DECISIONS.md`, and
`docs/development_plan.md` before significant changes. PRODUCT_SPEC.md is
authoritative. The complete bootstrap is in docs/prompts/INITIAL_CODEX_BOOTSTRAP.md.

1. Never silently weaken hard biological constraints or trade them against scores.
2. LOCKED DNA must remain byte-identical to the original construct.
3. Independently verify every affected protein, including overlapping and reverse CDSs.
4. Unannotated DNA is not freely mutable by default. Functional non-coding features
   default to LOCKED. Compose overlapping rules; never use last annotation wins.
5. Internal coordinates are 0-based half-open; user coordinates are 1-based inclusive.
6. Do not invent scientific parameters, datasets, provenance, or benchmark results.
7. Keep calibration and independent benchmark validation distinct.
8. Every bug fix requires a permanent regression test; never weaken tests to pass.
9. Every stochastic algorithm requires a controllable seed.
10. Final designs require fresh independent validation against the ORIGINAL construct.
    Candidate annotations, constraints, or translations cannot redefine the reference.
11. Backend/domain logic is the single source of scientific truth. Do not duplicate
    biological calculations in frontend code.
12. Sequence data remain local by default. External transmission must be opt-in.
    Commit only synthetic/public test sequences, never confidential sequences.
13. Run the complete relevant tests before declaring a milestone complete.
14. Do not automatically start the next phase. This bootstrap authorizes Phase 1 only.
15. Update repository documentation and the decision log for consequential changes.
16. Work on the existing development branch; never force-push, rewrite shared history,
    or automatically merge to main. Exclude credentials, environments, and caches.

Use Python >=3.12, a local .venv, pytest, and property-based tests for invariants.
Keep parsing, translation, motif matching, permissions, and final validation modular.
Read docs/scientific_provenance.md before adding biological models or datasets.
Unsupported biological representations must raise explicit errors, never be flattened.
