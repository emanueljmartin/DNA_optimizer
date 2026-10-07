# Architecture and Phase 1 boundaries

The platform designs annotated DNA while preserving hard biological constraints.
The authoritative complete requirements are in PRODUCT_SPEC.md. Phase 1 supplies
a local Python domain library; no optimizer, server, UI, or external sequence service
is authorized yet. Later FastAPI/Pydantic APIs and React/TypeScript views consume
this library rather than recalculate biology.

## Modules and interfaces

`src/dna_design/models.py`: frozen Pydantic models for constructs, features,
topology, strand, mutability, CDS interpretation, motifs, settings, and provenance.
Tuples prevent reordering/mutation of constraint and feature collections. Qualifiers
are retained for GenBank fidelity but do not override typed biological fields.

`coordinates.py`: explicit coordinate conversion with bounds checking.
`translation.py`: oriented feature extraction and table-aware CDS translation via
Biopython. Complete CDSs use biological start/stop validation; fragment mode must
be explicit and retains terminal stop symbols. Partial leading bases are locked.
`motifs.py`: IUPAC matching on both strands, positional scopes and circular origin
handling; restriction enzyme names resolve through Biopython's recognition database.
`permissions.py`: per-base composed substitution permissions. LOCKED dominates;
unannotated positions are locked; all CDS translations constrain allowed changes.
`parsing.py`: single-record raw/FASTA/GenBank import. Exact simple locations and
qualifiers are preserved. Unsupported locations and translation exceptions fail loudly.
`validation.py`: fresh candidate validation against reference sequence, annotations,
and constraints, with structured PASS/WARNING/ERROR results and provable conflicts.

Public interfaces: `parse_raw`, `parse_fasta`, `parse_genbank`, `translate_feature`,
`effective_permissions`, `validate_candidate`, `check_feasibility`. JSON round trips
use Pydantic's model serialization. The validator never trusts candidate source DNA,
annotations, embedded translations, or cached optimizer results.

## Constraint semantics

Phase 1 supports substitutions of equal-length, concrete A/C/G/T sequences. Indels
require future coordinate remapping. Import normalizes case/whitespace once; models
themselves accept canonical uppercase DNA only. Ambiguous input DNA is rejected;
IUPAC ambiguity is supported in motifs, not in design bases.

All overlapping restrictions compose. FREE_MUTATION permits substitutions only in
explicit designable annotations. A CDS still requires protein preservation even when
an overlapping feature is FREE_MUTATION. SYNONYMOUS_ONLY without a covering CDS
locks that position because no protein-preservation interpretation is available.
DESIGNABLE_RBS requires a matching explicit N/ACGT mask; fixed letters are immutable.
Functional features may be explicitly made designable; their defaults remain LOCKED.

PRESERVE binds reference occurrences (orientation and start) and their actual original
nucleotides, including degenerate motifs; REQUIRED demands a matching occurrence;
AVOID rejects all scoped occurrences. Both-strand matching is the default. Motif
regions contain the complete occurrence; position pins its forward genomic start.
Circular origin occurrences are supported for motifs up to one full construct turn.
Linear feature intervals cannot cross the origin; compound GenBank features are rejected.

Feasibility reports only proven contradictions: absent preserved sites, identical
required/forbidden recognition patterns with covering scopes, preserved forbidden
occurrences, and unavoidable AVOID occurrences entirely within immutable positions.
Other failures mean candidate invalidity or unknown feasibility, not proof that no
synonymous solution exists. No satisfiability solver is implemented in Phase 1.

## Ambiguities, risks, and decisions

* “AVOID where possible” versus never violating hard constraints: a forbidden motif
  causes ERROR. Proven immutable conflicts are marked infeasible; no silent relaxation.
* Complete versus partial CDS: default complete, requiring a start and terminal stop;
  fragments use explicit mode and frame offset. Fuzzy/compound/remote GenBank locations,
  transl_except, ribosomal slippage, and unrepresented exceptions are rejected.
* Alternative genetic codes: use NCBI IDs exposed by installed Biopython; codes with
  context-dependent dual stop/sense codons are unsupported rather than guessed.
* Annotated GenBank translations are checked against the nucleotide interpretation.
  Disagreement is an error, not permission to alter the protein reference.
* Generic `source` annotations are metadata, not functional locks on the entire genome.
  GenBank `gene` annotations are conservatively LOCKED: explicit user curation is
  needed to make their overlapping CDS designable. Imports prioritize fidelity.
* Python 3.12 is available; 3.13 is absent. Use 3.12 now and test 3.13 in CI.
* Host usage datasets, expression models, thermodynamic RBS parameters, synthesis
  thresholds and benchmark loaders remain deferred and require provenance.
* A score is not absolute expression; TIR will be relative model output, not measured
  proteins. PAO1 transferability must be separately validated.

## Validation and operational boundaries

Reference and candidate schema checks establish structure. Final validation re-extracts
and translates each reference CDS on both sequences, checks every locked position,
all motif constraints and mask literals, and rejects changed candidate control metadata.
Reference validity is checked before accepting any candidate. Invalid reports cannot
pass `require_valid()`, the boundary future exporters must call.

All computation stays local. No optimization cache, mutation acceptance score, or
candidate metadata supplies the biological baseline. Equal-length scanning is linear
in construct length times feature/motif count; multi-kb constructs are suitable for
Phase 1. Performance optimization follows profiling in later phases.
