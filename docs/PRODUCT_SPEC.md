# DNA Design Platform — authoritative product specification

The complete product requirements below are preserved from sections B–AC of the initial bootstrap. Only Phase 1 is currently authorized. Implementation status and limitations live in architecture.md and development_plan.md; they do not weaken these requirements.

# B. PERMANENT ROLE

Act as a senior:

- computational synthetic biology engineer
- bioinformatics developer
- scientific software architect
- optimization engineer
- backend developer
- frontend developer
- test engineer

The objective is not to create a toy codon optimizer.

The objective is to build a scientifically defensible, extensible DNA sequence-design platform for synthetic biology.

Initial biological focus:

1. *Pseudomonas aeruginosa* PAO1
2. *Escherichia coli* K-12/MG1655

The architecture must later support additional:

- organisms
- codon-reference datasets
- synthesis vendors
- RBS models
- cloning methods
- optimization algorithms
- sequence-design objectives

---

# C. HIGH-LEVEL PRODUCT

The application must eventually support the following conceptual workflow:

```text
Input construct
    ↓
feature annotation
    ↓
constraint definition
    ↓
expression/codon optimization
    ↓
RBS / translation-initiation optimization if requested
    ↓
synthesis optimization
    ↓
independent biological validation
    ↓
Gibson fragmentation if required
    ↓
export
```

Every step must understand which sequence regions it is permitted to alter.

No mutation may be accepted simply because it improves an optimization score if it violates a higher-priority biological or user-defined constraint.

---

# D. HARD CONSTRAINTS VERSUS SOFT OBJECTIVES

This distinction is fundamental.

## Hard constraints

Hard constraints may NEVER be silently violated.

Examples include:

- protected nucleotide regions must remain exactly unchanged;
- CDSs must encode exactly the same amino-acid sequence unless explicit protein modification was requested;
- overlapping CDSs must preserve all encoded proteins;
- required restriction sites/motifs must remain present;
- forbidden sites/motifs must be absent when a biologically valid synonymous solution exists;
- fixed positions in RBS design masks must remain unchanged;
- Gibson fragments must not exceed the specified maximum ordered length;
- user-defined immutable regions must remain byte-identical.

If hard constraints conflict, report the design as infeasible and explain the conflict.

Do not trade hard constraints against optimization scores.

## Soft objectives

Soft objectives may trade off against one another.

Examples include:

- host codon usage
- CAI
- rare-codon frequency
- rare-codon clustering
- global GC
- local/sliding-window GC
- homopolymers
- direct repeats
- inverted repeats
- low sequence complexity
- RNA structure
- synthesis difficulty
- codon-pair preferences
- sequence diversity
- RBS target TIR
- problematic internal sequence motifs

Never present the combined optimization score as a direct prediction of absolute protein expression.

---

# E. CORE MUTABILITY MODEL

Every nucleotide must effectively be governed by an explicit mutability state.

Implement at least:

## LOCKED

Nucleotide sequence must remain exactly unchanged.

Use for:

- protected regions
- cloning overhangs
- functional motifs
- user-protected sequences
- default functional non-coding features

## SYNONYMOUS_ONLY

DNA may change only if all affected CDS translations remain exactly unchanged.

This must be enforced biologically, not merely by looking at individual nucleotide positions.

The implementation must account for:

- strand
- reading frame
- translation table
- overlapping CDSs

## FREE_MUTATION

Arbitrary nucleotide substitutions may be permitted only when the sequence has explicitly been declared non-functional/designable.

Never assume unannotated non-coding DNA is freely mutable.

## DESIGNABLE_RBS

Sequence is available specifically for RBS optimization according to an explicit design mask and context.

Promoters, operators, terminators, RBSs and other functional non-coding features must default to LOCKED unless explicitly made designable.

When annotations overlap, do not use a naïve "last annotation wins" model.

The strictest applicable biological rules must compose.

---

# F. COORDINATES

Use exactly:

```text
Internal:
0-based, half-open coordinates

User interface:
1-based, inclusive coordinates
```

Never silently mix these conventions.

Provide explicit conversion functions and tests.

---

# G. CONSTRUCT DATA MODEL

A construct must eventually contain at least:

- nucleotide sequence
- construct name
- topology
- source sequence
- annotations/features
- ORFs/CDSs
- strand information
- optimization history
- active constraints
- organism/host profile
- synthesis profile
- RBS settings
- fragmentation settings
- random seed
- software version
- model versions

Each annotation should include at minimum:

- ID
- name
- feature type
- start
- end
- strand
- mutability
- optional description

Useful feature types include:

- CDS
- promoter
- RBS
- 5' UTR
- terminator
- operator
- linker
- tag
- cloning overhang
- restriction site
- spacer
- protected region
- mutable non-coding region
- custom feature

---

# H. INPUT FORMATS

Eventually support:

- raw DNA
- FASTA
- GenBank
- amino-acid sequence
- constructs containing multiple ORFs

When GenBank is imported, preserve annotations whenever they can be represented faithfully.

If a complex/compound GenBank location cannot yet be represented correctly, explicitly report that limitation.

Do not silently flatten biological annotations incorrectly.

Amino-acid reverse translation belongs to the codon-optimization phase.

---

# I. HOST PROFILES

Initial built-in hosts:

1. *Pseudomonas aeruginosa* PAO1
2. *Escherichia coli* K-12/MG1655

Host profiles must be extensible.

They should eventually support:

- genome-wide codon usage
- highly expressed gene codon usage where scientifically justified
- custom codon tables
- genetic code
- future tRNA/tAI models
- relevant RBS model parameters
- version information
- provenance/references

Do NOT invent high-expression codon tables.

Every built-in biological dataset must have documented provenance.

---

# J. EXPRESSION / CODON OPTIMIZATION

Do NOT implement codon optimization as:

"replace every codon with the most frequent synonymous codon."

Build a multi-objective optimization framework.

Eventually include:

- host codon usage
- CAI
- rare codons
- rare-codon clusters
- global GC
- local GC
- homopolymers
- direct repeats
- inverted repeats where practical
- sequence complexity
- duplicated sequence
- optional codon-pair bias
- 5' CDS RNA structure
- internal Shine-Dalgarno-like sequence penalties where scientifically justified
- synthesis constraints
- user-defined motifs

Provide presets eventually:

- Balanced
- Expression priority
- Synthesis priority
- Custom

Detailed objective weights should live in advanced settings.

The optimizer must:

1. parse hard constraints;
2. generate a valid starting sequence;
3. score it;
4. propose synonymous/local mutations;
5. reject hard-constraint violations;
6. accept/reject according to a configurable search algorithm;
7. run until convergence/budget;
8. independently validate the final result from scratch.

Potential algorithms include:

- stochastic local search
- simulated annealing
- evolutionary search
- iterative motif repair

Keep the solver interface modular.

Optimization must be reproducible from a random seed.

Eventually support:

- best design
- N alternative designs

Alternative designs must be genuinely different while satisfying all hard constraints.

---

# K. N-TERMINAL CODING REGION

Treat the first coding region downstream of the start codon specially.

Default conceptual length:

```text
30 nt
```

Eventually support modes such as:

- lock first X nt
- conservative optimization
- joint RBS/CDS optimization
- normal codon optimization

Do not assume codon frequency alone determines expression.

RNA structure around the translation-initiation region matters.

---

# L. RBS MODULE

RBS functionality is scientifically important.

Eventually support:

## Predict

Inputs:

- host
- upstream sequence
- RBS sequence
- downstream CDS context
- start codon

Outputs:

- predicted relative translation-initiation rate
- thermodynamic components when available
- relevant RNA structure
- spacing information
- warnings
- model version

Never claim the predicted TIR is an experimentally measured absolute number of proteins.

## Design

Allow:

- target relative TIR
- constant upstream sequence
- mutable RBS pattern
- downstream CDS
- host
- optional synonymous redesign of first X nt of CDS

Example mask:

```text
NNNNNNNAGGAGGNNNNNNN
```

Meaning:

```text
N = mutable
A/C/G/T = fixed
```

`AGGAGG` must remain fixed.

## RBS library

Eventually allow generating designs around targets such as:

```text
250
500
1000
2500
5000
10000
25000
```

## Joint RBS/CDS optimization

Allow synonymous redesign of the first X nt of CDS.

Protein sequence must remain identical.

If full-CDS codon optimization occurs afterward, it must either:

- preserve the optimized initiation region, or
- include TIR preservation explicitly as a constraint/objective.

Do not implement a simplistic Shine-Dalgarno heuristic and label it as a rigorous thermodynamic RBS model.

Use an abstraction layer so model backends can be substituted.

ViennaRNA/OSTIR or another defensible implementation may later be used where scientifically and legally appropriate.

Document limitations for *P. aeruginosa* if its model parameters are less validated than *E. coli*.

---

# M. RESTRICTION SITES AND MOTIFS

Implement conceptually distinct motif constraints:

## PRESERVE

An existing site/motif must remain at its position.

## AVOID

Occurrences should be removed where biologically possible, and new occurrences must not be created.

## REQUIRED

At least one occurrence must exist, optionally in a specified region/position.

Support:

- restriction-enzyme names
- arbitrary nucleotide motifs
- multiple motifs
- future IUPAC/degenerate motifs

For sites inside coding DNA, removal should use synonymous mutation when possible.

If removal is impossible because of hard constraints, explain why.

Circular constructs must eventually handle motifs spanning the sequence origin correctly.

---

# N. PROTECTED REGIONS

Eventually allow protecting regions through:

- graphical sequence selection
- coordinates
- motif search and occurrence selection

Protected regions must remain exactly nucleotide-identical.

---

# O. SYNTHESIS OPTIMIZATION

Create a synthesis-quality model independent from expression optimization.

Eventually support profiles:

- Generic
- Twist-public
- GenScript-public
- Custom

Never claim to reproduce proprietary synthesis algorithms.

Vendor profiles may include only publicly documented information.

Each vendor profile must include:

- profile version
- source
- date/version accessed
- explicit approximation statement

If a vendor threshold is not public, do not invent one.

Generic synthesis analysis should eventually inspect:

- global GC
- local GC
- homopolymers
- direct repeats
- repeated blocks
- tandem repeats
- inverted-repeat/hairpin risk where practical
- low complexity
- within-construct sequence similarity
- user-forbidden motifs
- problematic sequence structures

Mutation rules:

```text
coding DNA → synonymous changes
FREE_MUTATION non-coding DNA → arbitrary permitted changes
LOCKED DNA → never change
```

After synthesis optimization, all biological constraints must be reevaluated.

Do not collapse synthesis assessment into one unexplained score.

Use transparent findings such as:

```text
Global GC                     PASS
Local GC                      PASS
Homopolymers                  PASS
Direct repeats                WARNING
Forbidden motifs              PASS
Protected sequence integrity  PASS
Protein sequence integrity    PASS
```

---

# P. GIBSON FRAGMENTATION

Eventually split long constructs into synthesis-ready Gibson fragments.

Default conceptual overlap:

```text
30 nt
```

Maximum ordered fragment length includes overlap.

For total sequence length L, maximum ordered fragment length M and overlap O, the minimum feasible number of fragments must satisfy:

```text
n * M - (n - 1) * O >= L
```

Determine minimum feasible n before breakpoint optimization.

Breakpoint scoring should eventually consider:

- overlap GC
- overlap Tm
- uniqueness
- repeats
- composition
- secondary structure where practical
- protected regions
- synthesis quality

Advanced settings should eventually include:

- overlap length
- GC range
- preferred Tm
- breakpoint search window
- minimum fragment length

Every fragmentation result must validate:

- every fragment <= maximum ordered length
- overlaps are correct
- adjacent fragments reconstruct the complete sequence
- no bases are missing
- no unintended duplication exists
- fragment order is unambiguous

In-silico reconstruction must exactly reproduce the final construct.

---

# Q. USER INTERFACE

Eventually build a browser-based UI.

Backend:

- Python >=3.12
- FastAPI
- Pydantic
- Biopython where appropriate

Frontend:

- React
- TypeScript
- Vite or similar

Eventually provide an interactive annotated sequence viewer.

The user should be able to:

- click features
- select bases
- create annotations
- lock/unlock regions
- assign mutability
- mark CDSs
- define RBSs
- mark freely mutable spacers
- protect restriction sites
- search motifs
- view translation
- highlight mutations

Provide both:

- graphical feature view
- raw sequence view

For long sequences, keep the UI responsive.

Normal mode should expose only major settings.

Advanced mode should expose technical parameters.

Do not duplicate scientific calculations in the frontend.

The backend is the single source of truth.

---

# R. CONTEXTUAL HELP AND REFERENCES

Eventually every meaningful scientific metric/parameter should have contextual help.

Examples:

- CAI
- Target TIR
- Local GC
- Protected region
- Codon reference
- overlap length
- rare-codon threshold

Help should explain:

1. definition
2. what it influences
3. why changing it matters
4. important limitations
5. scientific reference where appropriate

Manage help text centrally.

Every major scientific model needs provenance.

Important literature includes at least:

Salis HM, Mirsky EA, Voigt CA. 2009.
Automated design of synthetic ribosome binding sites to control protein expression.
Nature Biotechnology.
DOI: 10.1038/nbt.1568

Kudla et al. 2009.
Coding-sequence determinants of gene expression in Escherichia coli.
Science.
DOI: 10.1126/science.1170160

Boël et al. 2016.
Codon influence on protein expression in E. coli correlates with mRNA levels.
Nature.
DOI: 10.1038/nature16509

Farasat et al. 2014.
DOI: 10.15252/msb.20134955

Johns et al. 2018.
Metagenomic mining of regulatory elements enables programmable species-selective gene expression.
Nature Methods.
DOI: 10.1038/nmeth.4633

Do not use references decoratively.

Tie citations to the algorithms/assumptions they support.

---

# S. OUTPUTS

Eventually provide:

## Final optimized construct

Full DNA sequence.

## Original versus optimized

Show:

- changed nucleotides
- changed codons
- unchanged protected regions
- original translation
- optimized translation

## Metrics

At minimum eventually:

- CAI/selected codon score
- codon statistics
- GC
- local GC
- rare codons/clusters
- repeat warnings
- synthesis warnings
- motif checks
- RBS/TIR metrics

## Mutation/change log

For each mutation where feasible:

- position
- original nucleotide/codon
- optimized nucleotide/codon
- amino acid
- reason/category

## Gibson fragments

Report:

- fragment number
- full ordered sequence
- unique sequence length
- overlap sequence
- total ordered length
- GC
- overlap metrics

---

# T. EXPORTS

Eventually support:

- FASTA
- GenBank
- CSV
- JSON project file

GenBank should preserve annotations.

JSON projects must contain enough information to reproduce a design:

- settings
- constraints
- feature annotations
- seed
- model versions
- source sequence
- final sequence

CSV should be useful for fragment ordering.

---

# U. BENCHMARKING

Do NOT validate the optimizer by requiring it to generate exactly the same synonymous DNA sequence as another optimizer.

Many synonymous sequences may be valid.

Evaluate:

- experimental ranking
- trends
- expression differences
- target TIR behavior
- constraint satisfaction
- synthesis-risk detection

Build a benchmark framework.

Potential public datasets include:

## Kudla 2009

Use synonymous GFP sequences to evaluate whether scoring distinguishes high versus low experimental expression.

## Boël 2016

Use heterologous expression data for codon/expression evaluation.

## Salis 2009

RBS benchmark/calibration.

## Farasat 2014

RBS validation.

## Johns 2018

Especially important because it includes *Pseudomonas aeruginosa* measurements.

Useful RBS benchmark statistics:

- Spearman correlation
- Pearson correlation on appropriate transformed values
- median fold error
- fraction within 2-fold
- fraction within 5-fold
- fraction within 10-fold

For values spanning orders of magnitude, use appropriate log-space analyses.

Training/calibration datasets and held-out validation datasets must remain explicitly separate.

Do not train on a dataset and present evaluation against that same dataset as independent validation.

Do not fabricate benchmark values.

---

# V. SYNTHESIS STRESS TESTS

Eventually construct synthetic test sequences deliberately containing:

- long homopolymers
- extreme local GC
- repeated blocks
- long direct repeats
- inverted repeats
- low complexity

Verify that scanners identify the inserted defect.

Then verify repair:

```text
problematic sequence
    ↓
synthesis optimizer
    ↓
repaired sequence
```

Confirm:

- intended defect improves
- protein remains unchanged for CDS
- protected regions remain unchanged
- no new forbidden motifs appear

---

# W. AUTOMATED TESTING

Scientific correctness is more important than UI polish.

Use:

- unit tests
- integration tests
- property-based tests
- regression tests
- API tests later
- frontend tests later
- benchmark tests later

Important invariants include:

## Sequence

- synonymous optimization preserves protein
- locked sequence remains byte-identical
- required motifs remain
- removable forbidden motifs can be removed synonymously

## Coordinates

- feature coordinates remain correct
- UI/internal conversion is correct

## RBS

Eventually:

- fixed mask positions never change
- mutable N positions may change
- CDS changes remain synonymous
- target distance improves relative to input

## Gibson

Eventually:

- reconstruction reproduces full sequence exactly
- fragment lengths satisfy limits
- overlaps are correct

Every bug fix must gain a permanent regression test.

---

# X. FINAL VALIDATION

Every final design must be independently validated from scratch.

Never trust only the optimizer's internal state.

Eventually check:

- valid DNA alphabet
- correct length
- CDS translation integrity
- protected sequence integrity
- required motifs/sites
- forbidden motifs/sites
- RBS fixed positions
- synthesis constraints
- Gibson reconstruction if used

If final validation fails, do not export the design as valid.

---

# Y. PRIVACY

Sequence data must remain local to the application by default.

Do not upload DNA to external web services silently.

Optional online functionality must be opt-in and clearly disclose what is transmitted.

During development, use synthetic/public sequences for committed tests.

Do not commit confidential/unpublished biological sequences unless explicitly instructed.

---

# Z. PERFORMANCE

The application must eventually remain practical for multi-kb constructs.

Use:

- caching where appropriate
- incremental calculations where appropriate
- progress reporting
- cancellation for long jobs
- deterministic seeds

Avoid needless recalculation.

---

# AA. DEVELOPMENT PHASES

Development must proceed in phases.

## Phase 1 — Core sequence engine

Implement:

- construct model
- feature annotations
- coordinate system
- CDS translation
- mutation permissions
- protected regions
- restriction/motif constraints
- raw DNA/FASTA/GenBank parsing
- validation framework
- extensive tests

NO codon optimizer before these invariants are reliable.

## Phase 2 — Codon optimization

Later implement:

- host profiles
- codon usage
- CAI
- synonymous mutation engine
- GC metrics
- rare codons
- multi-objective scoring
- reproducible optimizer
- alternative designs
- integrity tests

## Phase 3 — Synthesis optimization

Later implement:

- generic synthesis scanner
- repeats
- homopolymers
- local/global GC
- complexity
- freely mutable non-coding regions
- versioned synthesis profiles

## Phase 4 — RBS prediction/design

Later implement:

- RBS model abstraction
- host configuration
- prediction
- masks
- target optimization
- joint CDS redesign
- RBS libraries

## Phase 5 — Gibson fragmentation

Later implement:

- minimum fragment count
- breakpoint optimization
- overlap scoring
- reconstruction
- export

## Phase 6 — UI

Later build:

- sequence input
- sequence viewer
- annotations
- protected-region selection
- optimization settings
- advanced settings
- RBS interface
- Gibson interface
- reports
- exports
- help

## Phase 7 — Benchmarks

Implement public benchmark loaders and quantitative reports.

## Phase 8 — Polish

Improve:

- performance
- errors
- installers/packaging
- documentation
- usability

Do not collapse these into one uncontrolled implementation task.

---

# AB. PHASE 1 ACCEPTANCE REQUIREMENTS

Your immediate implementation target is Phase 1 ONLY.

Implement:

- central Pydantic/data models
- construct model
- feature annotations
- topology
- strand handling
- four mutability states
- coordinate conversion
- raw DNA parsing
- FASTA parsing
- GenBank parsing
- CDS translation
- reverse-strand translation
- genetic-code support
- protected regions
- motif constraints
- PRESERVE / AVOID / REQUIRED
- arbitrary nucleotide motifs
- IUPAC motifs where practical
- restriction-enzyme compatibility
- final validation
- infeasibility reporting where provable

Important design requirement:

Final validation must compare a candidate against the ORIGINAL/reference construct.

A candidate must never be allowed to redefine what "protected sequence" or "original protein" means.

---

# AC. REQUIRED PHASE 1 TESTS

At minimum test:

## Sequence validation

- legal DNA accepted
- illegal DNA rejected

## CDS preservation

- synonymous mutation passes
- nonsynonymous mutation fails
- reverse-strand CDS preserved
- alternative translation tables respected
- overlapping CDS invariants preserved

## Locked sequence

- locked DNA remains byte-identical
- candidate mutation inside lock fails

## Coordinates

- 1-based inclusive UI → 0-based half-open internal
- reverse conversion
- boundary cases
- off-by-one regression cases

## Annotation overlap

- LOCKED dominates
- FREE_MUTATION cannot weaken an overlapping functional annotation
- overlapping CDSs preserve all encoded proteins
- unannotated DNA is conservatively non-mutable

## Motifs

- AVOID
- REQUIRED
- PRESERVE
- multiple motifs
- IUPAC matching
- circular origin-spanning motif where implemented
- forbidden motif trapped inside immutable DNA reports infeasibility
- contradictory required + forbidden motif reports conflict

## Parsing

- raw DNA
- FASTA
- GenBank
- forward feature
- reverse feature
- annotation preservation
- unsupported compound locations explicitly reported

## Final validation

- original sequence is the preservation reference
- candidate cannot redefine original translation
- candidate cannot redefine protected DNA
- simultaneous hard constraints
- structured PASS/WARNING/ERROR reporting

## Property-based testing

Use it where helpful for:

- coordinate conversions
- translation invariants
- lock preservation
- sequence round trips

---

