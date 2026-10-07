You are taking long-term ownership of a new scientific software project called:

# DNA Design Platform

You are operating inside the GitHub repository for this project, in a GitHub Codespace.

Your first responsibility is to make this repository self-contained so that future Codex sessions can continue development without relying on this conversation.

Do not merely acknowledge these instructions. Execute them.

---

# A. BOOTSTRAP THE PROJECT KNOWLEDGE FIRST

Before implementing application functionality, create the following permanent repository files:

```text
docs/PRODUCT_SPEC.md
docs/architecture.md
docs/development_plan.md
docs/scientific_provenance.md
docs/testing_strategy.md
docs/DECISIONS.md
AGENTS.md
README.md
```

These files must become the persistent source of truth for future development.

Also create:

```text
docs/prompts/
```

and save this complete bootstrap instruction as:

```text
docs/prompts/INITIAL_CODEX_BOOTSTRAP.md
```

Future Codex sessions must be able to understand the project simply by reading:

```text
AGENTS.md
docs/PRODUCT_SPEC.md
docs/architecture.md
docs/development_plan.md
docs/DECISIONS.md
```

---

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

# AD. AGENTS.md CONTENT

Create an `AGENTS.md` that future Codex sessions must obey.

At minimum it must say:

1. Read `docs/PRODUCT_SPEC.md` before major architectural changes.
2. PRODUCT_SPEC.md is authoritative.
3. Never silently weaken hard biological constraints.
4. LOCKED DNA is byte-identical.
5. Protein preservation must be independently verified.
6. Unannotated DNA is not freely mutable by default.
7. Internal coordinates are 0-based half-open.
8. User coordinates are 1-based inclusive.
9. Do not invent scientific parameters or provenance.
10. Calibration and independent benchmark validation must remain distinct.
11. Every bug fix requires a regression test.
12. Every stochastic algorithm requires a controllable seed.
13. Final designs require independent validation.
14. Backend/domain logic is the single source of scientific truth.
15. Do not duplicate biological calculations in frontend code.
16. Sequence data remain local by default.
17. Run tests before declaring a milestone complete.
18. Do not weaken tests to obtain a passing suite.
19. Do not automatically start the next development phase.
20. Keep repository documentation updated as architectural decisions change.

Also instruct future agents to read:

```text
docs/DECISIONS.md
docs/development_plan.md
```

before significant changes.

---

# AE. DECISION LOG

Create:

```text
docs/DECISIONS.md
```

Maintain a lightweight architectural/scientific decision log.

For consequential decisions record:

```text
Date
Decision
Reason
Alternatives considered
Scientific/technical implications
References if applicable
```

Future Codex sessions should update this when they make meaningful architecture/scientific decisions.

---

# AF. DEVELOPMENT PLAN

Create:

```text
docs/development_plan.md
```

Turn the phase structure above into a practical milestone plan.

Each milestone must have:

- scope
- dependencies
- implementation targets
- tests
- definition of done
- intentionally deferred work

Never mark something complete without executable evidence where testing is possible.

---

# AG. SCIENTIFIC PROVENANCE FILE

Create:

```text
docs/scientific_provenance.md
```

Use this as the index of:

- biological datasets
- parameters
- equations
- scientific models
- vendor public rules
- benchmark datasets
- DOI/reference
- version/access date
- where each is used in code

Do not populate unsupported biological values simply to fill the table.

Use `TBD — requires provenance` where appropriate.

---

# AH. DEVELOPMENT ENVIRONMENT

This repository is being developed in GitHub Codespaces.

Use the Codespace terminal for deterministic computation.

Confirm at startup:

```bash
pwd
git status
git branch --show-current
python3 --version
node --version
npm --version
```

The workspace should be under:

```text
/workspaces/dna-design-platform
```

Set up a project-local Python virtual environment if appropriate.

Do not modify the university host PC.

Do not require administrator access to the client machine.

Eventually create:

```text
.devcontainer/devcontainer.json
```

so the development environment becomes reproducible.

The devcontainer should eventually pin/document:

- Python version
- Node version
- required system dependencies
- useful VS Code extensions
- setup commands

Do not over-engineer the devcontainer before Phase 1 is working.

---

# AI. GIT RULES

Use Git carefully.

Do not:

- rewrite shared history
- force push
- commit credentials
- commit `.venv`
- commit `node_modules`
- commit caches
- commit giant generated benchmark outputs
- commit confidential biological sequences

Work on the existing development branch.

Before suggesting a checkpoint:

```bash
git status
git diff --stat
```

and run the complete relevant test suite.

At the end of Phase 1, suggest a commit message but do not merge to `main` automatically.

---

# AJ. WORK STYLE

Work autonomously.

Do NOT stop after:
- creating a scaffold,
- writing one model,
- or creating a handful of tests,

to ask whether you should continue.

Continue until the current milestone is a meaningful tested checkpoint.

For routine engineering decisions, make a defensible choice and document it.

Only ask for user input when there is a genuinely consequential scientific choice that cannot be resolved responsibly from:

- this specification,
- established software-engineering practice,
- or reliable scientific evidence.

Do not invent an answer merely to avoid asking.

Do not start Phase 2 automatically.

---

# AK. IMMEDIATE EXECUTION PLAN

Execute the following now:

## Step 1 — environment inspection

Run:

```bash
pwd
git status
git branch --show-current
python3 --version
node --version
npm --version
```

Confirm that you are operating in the expected Codespace repository.

## Step 2 — persistent project knowledge

Create:

```text
docs/PRODUCT_SPEC.md
AGENTS.md
docs/architecture.md
docs/development_plan.md
docs/scientific_provenance.md
docs/testing_strategy.md
docs/DECISIONS.md
docs/prompts/INITIAL_CODEX_BOOTSTRAP.md
```

Populate them with the relevant permanent information from this prompt.

`docs/PRODUCT_SPEC.md` must contain the complete product requirements from sections B through AC above in a well-organized form.

Do not replace the requirements with a short summary.

`docs/prompts/INITIAL_CODEX_BOOTSTRAP.md` should preserve this full bootstrap prompt for reference.

## Step 3 — architecture

Critically review the product requirements.

Document:

- contradictions
- ambiguities
- risks
- unresolved scientific assumptions
- dependency concerns
- final architecture
- key interfaces
- validation strategy
- testing strategy

## Step 4 — repository scaffold

Create the maintainable project structure.

Create appropriate:

```text
pyproject.toml
.gitignore
.gitattributes
README.md
```

Use Python >=3.12.

Prefer Python 3.13 for this Codespace if dependency compatibility permits.

## Step 5 — Phase 1 implementation

Implement the Phase 1 core sequence engine.

Do not implement the Phase 2 optimizer.

## Step 6 — tests

Write and execute Phase 1 tests.

Use pytest.

Use property-based tests where useful.

Every discovered bug must receive a regression test.

## Step 7 — final verification

Before declaring the milestone complete, run:

- full pytest suite
- Python import/compilation sanity check
- relevant static/type checks if configured
- `git status`
- `git diff --stat`

Inspect failures and fix them.

## Step 8 — final report

Only after the checkpoint is genuinely tested, report:

### Architecture
What was created.

### Implementation
What modules/features were implemented.

### Scientific invariants
What is actually enforced.

### Assumptions
What assumptions were made.

### Tests
Exact commands and:

```text
passed:
failed:
skipped:
```

### Unresolved issues
Anything still scientifically or technically uncertain.

### Git
Current branch and working-tree status.

### Suggested commit
Suggested commit message.

### Next milestone
State what Phase 2 would contain, but DO NOT begin it.

---

# AL. CONTINUITY REQUIREMENT

A critical success criterion for this bootstrap is:

After this task finishes, a future Codex session must NOT need this chat.

It should be able to start by reading:

```text
AGENTS.md
docs/PRODUCT_SPEC.md
docs/architecture.md
docs/development_plan.md
docs/DECISIONS.md
docs/scientific_provenance.md
```

and understand:

- what the application is;
- why it exists;
- biological invariants;
- software architecture;
- current implementation state;
- completed work;
- next milestone;
- known unresolved scientific issues.

Build the repository documentation accordingly.

Begin now.

Do not merely describe what you intend to do.

Inspect the Codespace and execute the bootstrap and Phase 1 implementation.