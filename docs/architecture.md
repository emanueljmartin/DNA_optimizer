# Architecture and authorized phase boundaries

The platform designs annotated DNA while preserving hard biological constraints.
The authoritative complete requirements are in PRODUCT_SPEC.md. Phase 1 supplies
a local Python domain library. The user subsequently authorized Phases 2–4 and
then the full tool through Phases 5–8 on 2026-10-08. Codon, synthesis, local RBS,
Gibson, project/export and benchmark modules extend that core. FastAPI/Pydantic APIs
and React/TypeScript views consume the library; biological calculations stay in Python.

## Modules and interfaces

`src/dna_design/models.py`: frozen Pydantic models for constructs, features,
topology, strand, mutability, CDS interpretation, motifs, settings, and provenance.
Tuples prevent reordering/mutation of constraint and feature collections. Qualifiers
are retained for GenBank fidelity but do not override typed biological fields.

`coordinates.py`: explicit coordinate conversion with bounds checking.
`translation.py`: oriented feature extraction and table-aware CDS translation via
Biopython. Manually annotated bacterial CDSs default to table 11; GenBank imports
honor INSDC's omitted-transl_table default of 1. Complete CDSs use biological
start/stop validation; fragment mode must
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

`profiles.py` and `data/*.json`: extensible sourced codon references and two packaged
genome-wide profiles. `scripts/derive_codon_profile.py` reproduces aggregate counts
from local public GenBank snapshots with checksums and explicit exclusions.
`metrics.py`: global/local GC, Shannon base entropy, codon adaptation/conditional CAI,
explicit rare thresholds and clusters. `reverse_translation.py`: seeded weighted
initialization from canonical protein input, followed by independent translation.

`optimization.py`: configurable seeded hill climbing/simulated annealing, synonymous
and explicitly permitted non-coding proposals, original-reference gates, additional
locks, objective components, base-change logs and distinct alternatives. Motif repair
may traverse unresolved REQUIRED/AVOID states but never protein/lock/mask violations.
No invalid sequence is returned as a Design. Only the Phase 1 feasibility proofs
produce `infeasible`; insufficient exploration produces `budget_exhausted`. Invalid
input controls produce `invalid_reference`. `success` means the requested number of
valid designs was produced, not that an optimum or a biological target was reached.

`synthesis.py`: versioned generic/custom and public-guidance profiles, global/local
GC, GC spread, homopolymers, direct/inverted repeats and entropy findings; same-solver
repair with visible residual risk. Synthesis risk is a soft loss unless the caller
requires the selected policy via `require_policy()`; that gate fails on unresolved
or unassessed criteria. Public profiles are incomplete advisory approximations.

`rbs.py`: replaceable protocol, optional lazy OSTIR/ViennaRNA integration, oriented
contexts, relative TIR/energy/spacing/unbound structure reports, masked seeded target
search, optional joint CDS design and distinct-RBS target libraries. Missing/unsupported
backends fail explicitly. Final predictions bypass score caches. Library selection
retains the complete search result and a selected index, rather than rewriting history.

All sequential stages take the ORIGINAL construct plus an optional valid starting
DNA sequence. Proposal choices use that starting DNA, including partially locked
synonymous codons; original proteins/locked DNA remain the validation reference.
Fresh search validation separately enforces allowed-position boundaries, initiation
protection relative to the starting sequence and explicitly carried AdditionalLocks.
After RBS design, callers carry the corresponding initiation lock into later searches.
It covers the complete supplied prediction context and any larger joint design region,
so later codon changes cannot silently change the predicted initiation rate.

## Constraint semantics

`gibson.py` chooses the lowest feasible fragment count using deterministic dynamic
programming over all eligible boundaries. Ordered lengths include shared overlaps;
circular closure is accounted for separately. Unique full overlaps are mandatory
for final export. DNA_NN3 Tm uses explicit concentrations/salt assumptions, with
public Biopython defaults rather than an inferred Gibson reaction mix. Independent
validation rechecks DNA, tiling/order, overlaps, count/length limits, metadata and Tm.

`projects.py` preserves original controls, candidate DNA, stage settings/history,
profile contents/fingerprints, software versions, context locks and optional assembly.
The browser restores each stage's latest controls from full request history, so a
later assembly operation does not erase codon/RBS controls on reload. Codon and RBS
search budgets remain separate. Unexposed search restrictions and assembly budgets
are retained on rerun; all explicit Tm solution fields are shown and restored.
Selecting another host clears a restored/custom codon override and previous anti-SD
parameters, so the selected built-in host supplies the scoring reference and foreign
RBS settings still require explicit sourced parameters and opt-in.
Draft snapshots can retain unfinished work; final export separately enforces fresh
validation. FASTA, GenBank, project JSON and fragment CSV have different information
content. GenBank uses explicit validated custom qualifiers for platform feature IDs,
mutability, RBS masks and fragment-CDS interpretation; original source/literature
metadata is preserved. JSON carries cross-stage locks and complete settings/history.

`service.py` orchestrates domain operations and backend-only analysis/codon change
logs. Explicit annotation or motif curation creates a new baseline from the original
DNA and clears earlier design state. Sequential design operations retain the same
original reference and starting candidate. Combined codon/synthesis preferences are
soft engineering weights; they never weaken original biological gates.

`api.py` is an optional local FastAPI service. One background worker runs at most
four active/queued jobs, with in-memory progress, cooperative cancellation, opaque
job IDs and bounded retention (16 records / 30-minute completed-job expiry). Failed
or cancelled jobs leave the caller's project intact. Same-origin requests require a
local-client header; bodies have a 2 MB limit, constructs a 100,000-base interactive
limit and searches explicit iteration/alternative limits. No DNA access logs, CORS,
external prediction services, telemetry exporter, CDN assets or automatic persistence.
`cli.py` binds to loopback by default. The browser uses backend reports, coordinate
conversion for display, and a 700-base paged viewer; files stay in the user's browser
until explicitly saved or posted to the local process.

`benchmarks.py` provides versioned mapped public/local CSV loading and quantitative
held-out reports. Manifests record checksum, provenance/license, units, column mapping,
transforms, exclusions, group/split policy and model-training overlap. Calibration is
optional and fitted only to disjoint positive calibration rows per host/condition.
Ordinal scores receive rank/linear association metrics; log/fold errors require a
ratio scale and comparable units or explicit separate calibration. Replicate groups,
exact sequences and declared model-training overlaps are checked. No experimental
observations or claimed accuracy are bundled; public source conversion remains
explicit and study choices will be agreed with the user.
The browser renders backend metrics and evidence/independence warnings, and downloads
the unchanged report JSON. A report is displayed only for the exact CSV/manifest
text evaluated, preventing a delayed response from labeling edited study inputs.

`execution.py` supplies cooperative progress/cancellation without changing random
draws. Search retains only the best requested distinct alternatives; RBS scoring caches
up to 512 oriented contexts and bypasses the cache for final predictions. Reproducible
synthetic profiling lives in `scripts/profile_local.py`. Installer/build scripts and
optional dependency locks keep frontend packaging and local model execution explicit.

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
occurrences, unavoidable AVOID occurrences entirely within immutable positions,
and REQUIRED motifs with no compatible scoped window even when every unlocked
position is relaxed to arbitrary A/C/G/T.
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
* The rebuilt Codespace supplies Python 3.13.5 and a working local `.venv`. The original
  checkpoint used 3.12.3; configured CI covers 3.12/3.13. Setup reuses valid environments,
  creates absent ones with copied interpreters, and rejects unusable existing paths.
* Host codon counts, synthesis guidance snapshots and pinned RBS backend parameters
  have documented provenance. High-expression built-ins, tAI/codon-pair models,
  experimentally calibrated PAO1 RBS parameters remain deferred. Public experimental
  study execution and exact source-column/split selections remain pending user review.
* A codon score is not absolute expression. TIR is a relative model output, not measured
  proteins. No held-out experimental validation has been performed in this project;
  PAO1 transfer requires explicit uncalibrated opt-in and sourced parameters.

## Validation and operational boundaries

Reference and candidate schema checks establish structure. Final validation re-extracts
and translates each reference CDS on both sequences, checks every locked position,
all motif constraints and mask literals, and rejects changed candidate control metadata.
Reference validity is checked before accepting any candidate. Invalid reports cannot
pass `require_valid()`, which all final exporters call.

All computation stays local. No optimization cache, mutation acceptance score, or
candidate metadata supplies the biological baseline. Equal-length scanning is linear
in construct length times feature/motif count; multi-kb constructs are suitable for
the current domain implementation. Search additionally repeats independent validation
per proposal; thermodynamic design is bounded by a caller-selected budget. General
solver completeness remains outside scope. Recorded synthetic timings describe this
machine and budget only; they do not establish biological reliability.
