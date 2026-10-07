# Scientific provenance index

Sources reviewed on 2026-10-07. Code paths below are the Phase 1 domain library.
Installed dependency versions are pinned in pyproject.toml. No host usage dataset,
vendor threshold, expression model or RBS prediction dataset is currently bundled.

| Item | Source / version | Use and limits | Code |
| --- | --- | --- | --- |
| Genetic code IDs, start/stop codons | [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi); installed Biopython CodonTable | Bacterial default table 11; per-CDS explicit override; dual stop/sense codes rejected | translation.py, models.py |
| Complete-CDS translation | [Biopython sequence tutorial](https://biopython.org/docs/latest/Tutorial/chapter_seq_objects.html) | Alternative starts translate as M in complete CDS; strict start/stop validation; fragments explicitly distinct | translation.py |
| GenBank annotation representation | [Biopython annotations tutorial](https://biopython.org/docs/latest/Tutorial/chapter_seq_annot.html) | Exact simple locations, strand, codon_start and qualifiers; unrepresentable annotations rejected | parsing.py |
| Restriction recognition sequences | [Biopython restriction documentation](https://biopython.org/docs/latest/api/Bio.Restriction.html) | Recognition compatibility only, not cut chemistry or digest planning; bundled dependency database | motifs.py |
| IUPAC DNA symbols | Biopython IUPACData.ambiguous_dna_values | Used only for motif interpretation; input DNA remains concrete | motifs.py |
| Coordinates and immutability | PRODUCT_SPEC.md E–F, AB | Engineering contracts, not empirical biological parameters | coordinates.py, permissions.py, validation.py |
| Initial hosts: PAO1, MG1655 | TBD — requires provenance | Profiles/codon tables deferred to Phase 2; do not invent high-expression references | Not implemented |
| CAI / rare codon thresholds | TBD — requires provenance | Reference set and zero-frequency policy required before implementation | Not implemented |
| RBS thermodynamic parameters, host calibration | TBD — requires provenance | Backend/license and PAO1 validity review required | Not implemented |
| Generic synthesis thresholds | TBD — requires provenance | Transparent configurable heuristics, not measured vendor acceptance | Not implemented |
| Twist-public / GenScript-public | TBD — requires provenance | Public sources/access dates/approximation statement required; proprietary rules unavailable | Not implemented |
| 30 nt initiation region, 30 nt overlap | User specification K, P | Conceptual defaults, not scientifically validated universal optima | Deferred |

## Planned literature and benchmark index

The following references were specified by the user. They are queued for source and
data-access review, not presented as validation of any implemented expression model.
Before use, record exact datasets, versions/checksums, licenses, transforms, exclusions,
calibration split, held-out split and the code/model version. Numerical benchmark
results are currently absent.

| Reference | Proposed use | Current status |
| --- | --- | --- |
| Salis, Mirsky, Voigt (2009), [10.1038/nbt.1568](https://doi.org/10.1038/nbt.1568) | Thermodynamic RBS abstraction and held-out RBS validation | Deferred; requires full source/data review |
| Kudla et al. (2009), [10.1126/science.1170160](https://doi.org/10.1126/science.1170160) | Synonymous GFP ranking; initiation structure assumptions | Deferred; no accuracy claims |
| Boël et al. (2016), [10.1038/nature16509](https://doi.org/10.1038/nature16509) | Codon/mRNA/expression evaluation | Deferred; dataset provenance required |
| Farasat et al. (2014), [10.15252/msb.20134955](https://doi.org/10.15252/msb.20134955) | RBS calibration/validation separation | Deferred |
| Johns et al. (2018), [10.1038/nmeth.4633](https://doi.org/10.1038/nmeth.4633) | Species-selective regulatory models, PAO1 assessment | Deferred; host-specific validity unresolved |

Benchmark metrics will include Spearman, appropriately transformed Pearson, median
fold error and fractions within 2/5/10-fold. Units, zero-value handling, and log-space
transform assumptions must be stated; calibration data cannot be reused as independent
validation. Literature is tied to its proposed model role rather than added decoratively.
