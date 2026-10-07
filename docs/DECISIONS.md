# Architectural and scientific decisions

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

Decision: Biopython/NCBI genetic codes, table 11 default for bacterial CDSs, explicit
fragment mode, reject unsupported exception/fuzzy/compound/remote locations and
context-dependent genetic codes. Reason: no faithful simple-interval representation
for these features. Alternatives: flatten joins or guess recoding. Implications:
some real-world records fail with explicit limitations. References:
[NCBI codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi),
[Biopython sequence tutorial](https://biopython.org/docs/latest/Tutorial/chapter_seq_objects.html).

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
