# Development milestones

Only Phase 1 is authorized. Do not automatically start Phase 2. Complete requirements
are in PRODUCT_SPEC.md; executable evidence is recorded here when available.

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

## Current checkpoint

Bootstrap documents created before application functionality. Phase 1 implementation
and verification are in progress; do not mark complete until the exact commands and
results below are populated after execution.

Known representation limits: substitutions only; concrete DNA; exact single intervals;
no recoding/compound features; recognition motifs rather than enzyme cut simulation;
no general satisfiability solver. See architecture.md for the scientific implications.
