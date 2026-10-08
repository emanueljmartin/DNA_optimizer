# Scientific provenance index

Sources reviewed on 2026-10-07 and 2026-10-08. Code paths cover Phases 1–8.
Dependency versions are pinned in pyproject.toml, development/RBS/API lock files and
frontend/package-lock.json.
Only aggregate public host codon counts are bundled; no experimental expression
benchmark dataset or private sequence is bundled. All prediction/design runs are local.

| Item | Source / version | Use and limits | Code |
| --- | --- | --- | --- |
| Genetic code IDs, start/stop codons | [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi); Biopython 1.88 CodonTable | Manual bacterial CDS default 11; imported GenBank omitted transl_table defaults to 1; per-CDS override; dual stop/sense codes rejected | translation.py, models.py, parsing.py |
| Complete-CDS translation | [Biopython sequence tutorial](https://biopython.org/docs/latest/Tutorial/chapter_seq_objects.html) | Alternative starts translate as M in complete CDS; strict start/stop validation; fragments explicitly distinct | translation.py |
| GenBank annotation representation | [Biopython annotations tutorial](https://biopython.org/docs/latest/Tutorial/chapter_seq_annot.html) | Exact simple locations, strand, codon_start and qualifiers; unrepresentable annotations rejected | parsing.py |
| Restriction recognition sequences | [Biopython restriction documentation](https://biopython.org/docs/latest/api/Bio.Restriction.html) | Recognition compatibility only, not cut chemistry or digest planning; bundled dependency database | motifs.py |
| IUPAC DNA symbols | Biopython IUPACData.ambiguous_dna_values | Used only for motif interpretation; input DNA remains concrete | motifs.py |
| Coordinates and immutability | PRODUCT_SPEC.md E–F, AB | Engineering contracts, not empirical biological parameters | coordinates.py, permissions.py, validation.py |
| GenBank default genetic code | [NCBI default transl_table convention](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/genetic-codes/) | Omission means table 1; do not infer code 11 from the application's initial host focus | parsing.py |
| Initial hosts: PAO1, MG1655 | NCBI NC_002516.2 / NC_000913.3; snapshots/checksums below | Aggregate complete-CDS elongation counts; no fabricated high-expression set | profiles.py, data/*.json, scripts/derive_codon_profile.py |
| CAI / codon adaptation | [Sharp and Li (1987)](https://doi.org/10.1093/nar/15.3.1281) | Relative synonymous weights; conditional high-expression CAI; explicit zero policy below | metrics.py, profiles.py |
| Rare-codon threshold | Explicit caller-selected value | Relative weight below threshold; consecutive indices form clusters; no universal biological cutoff | metrics.py, optimization.py |
| RBS parameters / host calibration | [OSTIR paper](https://doi.org/10.21105/joss.03362); OSTIR 1.1.3 / ViennaRNA 2.7.2 | Published E. coli calibration; transfer and license limits below | rbs.py |
| Generic synthesis thresholds | Explicit user-selected engineering policies | Defaults leave numeric criteria unassessed; no empirical acceptance predictor | synthesis.py |
| Twist-public / GenScript-public | Public sources below, accessed 2026-10-08 | Product-specific incomplete guidance snapshots, not proprietary algorithms | synthesis.py |
| 30 nt initiation region | User specification K | Configurable protection default, not a universal optimum; RBS locks extend to full context | optimization.py, rbs.py |
| 30 nt assembly overlap | User specification P | Explicit configurable engineering overlap; ordered lengths include shared bases | gibson.py |
| Overlap Tm | [Biopython MeltingTemp API](https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html), Biopython 1.88 | DNA_NN3 (Allawi–SantaLucia 1997), displayed strand/salt assumptions; not measured Gibson yield | gibson.py |
| Benchmark metrics | Explicit evaluator definitions in benchmark_workflow.md | Average tied-rank Spearman, linear/log Pearson, positive ratio-scale fold error; explicit calibration/units/splits | benchmarks.py |

## Public codon-reference snapshots

Downloaded on 2026-10-08 from NCBI EFetch using the versioned accessions below
(`db=nuccore`, `rettype=gbwithparts`, `retmode=text`). These downloads contain only
public records; runtime profile loading never downloads data. NCBI imposes no use
restrictions on its molecular data but does not adjudicate all submitters' third-party
rights: [NCBI policy](https://www.ncbi.nlm.nih.gov/home/about/policies/).

| Profile | Source record | Included CDSs | Excluded CDSs | Source file SHA-256 |
| --- | --- | --- | --- | --- |
| ecoli-mg1655 | [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3) | 4,294 | 24 exceptional/pseudogene | `313752bfe10d4aaf04abdb541a1309bce9a6b66f059e75fe215898e736a5a8a5` |
| pseudomonas-pao1 | [NC_002516.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_002516.2) | 5,570 | 2 exceptional/pseudogene | `8b4388ea984a0d7df4850e556aebc17acf77f2982acf2f192088df01befb4a6c` |

Derivation v1 accepts exact single-interval, explicit-strand, complete table-11 CDSs
with codon_start=1 and canonical DNA. It excludes pseudo/pseudogene, exception,
transl_except, ribosomal slippage, unsupported/partial/remote locations, other codes
or frames, invalid complete translations and nucleotide/annotated-protein disagreement.
Every exclusion is counted by reason. Internal elongation codons are counted; the
initiator and terminal stop are omitted. Overlapping qualifying CDSs each contribute
their coding observations. No highly expressed gene set is inferred from annotation.

Reproduce from a local downloaded source snapshot:

```bash
python scripts/derive_codon_profile.py /path/NC_000913.3.gb /tmp/ecoli-profile.json \
  --id ecoli-mg1655 --organism 'Escherichia coli' --strain 'K-12/MG1655' \
  --accessed-date 2026-10-08
python scripts/derive_codon_profile.py /path/NC_002516.2.gb /tmp/pao1-profile.json \
  --id pseudomonas-pao1 --organism 'Pseudomonas aeruginosa' --strain PAO1 \
  --accessed-date 2026-10-08
```

Exact reproduction requires the source file checksum to match. Accession sequence
versions do not guarantee immutability of annotation snapshots. The packaged JSON
records contain the checksum, date, version, transform and exclusions, and expose
a SHA-256 fingerprint of their entire canonical content. Only these small aggregate
JSON files are committed; downloaded genome records remain outside the repository.

## Codon-score and initialization policies

For each sense codon, weight is its observed count divided by the maximum count for
its synonymous amino acid. Counts of zero remain zero, with no fabricated pseudocount.
The gene score is the geometric mean, omitting single-synonym amino acids under its
selected table, and complete-CDS start/stop. A contributing zero yields score zero;
an empty contributing set yields one with no evidence of expression improvement.

The result is named **codon adaptation** for genome-wide/custom references. The `cai`
field is populated only when the explicit sourced profile is marked `high_expression`;
users are responsible for supplying a defensible reference set. No high-expression
profile is bundled. Rare thresholds use these relative weights and must be explicit.
GC targets, objective weights and annealing temperature are engineering preferences,
not fitted expression parameters. All stochastic methods record a caller-supplied seed.

Amino-acid reverse translation samples observed synonymous counts to initialize a
new DNA reference. Complete input must start with M; ATG initiation and a terminal
stop are explicit policies, without adding amino acids. Stops are sampled uniformly
because elongation profiles deliberately omit stop usage. Fragment mode must be
explicit. Every generated CDS is freshly translated against the input protein.

## Synthesis policies and limits

Generic/custom profiles require explicit numerical policies for global/local GC,
GC spread, homopolymers, direct/inverted repeats and Shannon base entropy. Entropy
measures mononucleotide composition, not a full empirical sequence-complexity model.
Repeat scans use exact words at the selected minimum length; inverted scans report
one non-overlapping pair per word/reverse-complement family. Circular scans include
origin-spanning windows/runs bounded by one construct turn. Unconfigured or overly
large windows remain visibly unassessed. Findings show positions, observed values
and selected-policy penalties. Their sum is a transparent engineering loss rather
than a vendor acceptance probability. Repair can retain immutable or residual risks.

`genscript-public`, snapshot 2026-10-08, uses public GC warning boundaries 20–80%
and repeat length 20 bp from [GenScript gene synthesis guidance](https://www.genscript.com/gene_synthesis.html).
The page also names homopolymers/palindromes without numerical cutoffs; none are invented.

`twist-public`, snapshot 2026-10-08, represents **clonal gene/vector design guidance**:
25–65% GC, 20-bp repeats and a maximum 50-percentage-point spread among 50-bp GC
windows. Source: [Twist vector-design technical note, revision 4](https://www.twistbioscience.com/content/dam/twistbioscience/resources/2023-06/DOC-001081_TechNote_TwistTipVectorDesign-REV4-singles.pdf).
The guidance's alternative repeat melting-temperature condition is omitted because
no repeat Tm model is implemented. This is not a universal Twist product policy or
a reconstruction of proprietary scoring. Vendor rules can change after these snapshots.

## Local RBS model and license boundary

Backend: [OSTIR 1.1.3](https://github.com/barricklab/ostir/tree/v1.1.3), using the
[ViennaRNA 2.7.2](https://github.com/ViennaRNA/ViennaRNA/tree/v2.7.2) Python API.
Credit: Cameron T. Roots, Alexandra Lukasiewicz and Jeffrey E. Barrick,
[OSTIR: open source translation initiation rate prediction (2021)](https://doi.org/10.21105/joss.03362);
ViennaRNA authors and the Institute for Theoretical Chemistry, University of Vienna.
The backend is optional and lazily imported; its absence never triggers a fabricated
Shine-Dalgarno substitute. It accepts table-11 complete CDS contexts and ATG/GTG/TTG.

The upstream [factory source](https://github.com/barricklab/ostir/blob/v1.1.3/ostir/ostir_factory.py)
documents calibration on Salis2009, RNA2004 folding settings, 37°C and 35-nt cutoff.
The E. coli anti-SD default `ACCTCCTTA` comes from its
[public API source](https://github.com/barricklab/ostir/blob/v1.1.3/ostir/ostir.py).
Scientific constants are consumed from upstream, without copying or refitting them.
Relative TIR, energy terms and spacing are returned with model versions. The displayed
unbound full-context MFE structure is separately labeled and does not purport to be
the ribosome-bound conformation or OSTIR's cropped mRNA energy component.

No project calibration or held-out experimental validation was performed. Upstream
calibration is not independent validation here. Synthetic tests establish execution,
mask/strand/protein integrity, deterministic search and target-loss behavior only.
There is no built-in PAO1 RBS calibration. Other hosts require a sourced nine-base
anti-SD, an uncalibrated label and explicit opt-in; no experimentally verified host
transfer or absolute protein-expression accuracy is claimed. Context truncation is
reported rather than hidden. Designed initiation locks cover all supplied context.

Licenses were inspected in the installed distributions and upstream sources:
[OSTIR GPLv3](https://github.com/barricklab/ostir/blob/master/LICENSE) includes
copyleft/source obligations for redistribution; [ViennaRNA COPYING](https://github.com/ViennaRNA/ViennaRNA/blob/v2.7.2/COPYING)
permits research/educational/commercial use and modification with attribution, limits
fee-based redistribution, and directs commercial-product incorporation to the authors.
Local use/testing proceeds without redistributing their binaries or source in this
repository. Optional dependency status does not remove their licensing conditions.
Commercial packaging/distribution needs a separate license compatibility decision;
this checkpoint does not grant permission or declare a platform distribution license.

## Planned literature and benchmark index

The following references were specified by the user. Primary paper/supplement discovery
links and a versioned mapped-CSV loader are implemented. Exact raw column/assay mappings,
licenses and experimental splits still require study preparation with the user; these
references are not presented as validation of the implemented expression model.
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

The evaluator now provides Spearman, linear/ratio-scale log Pearson, median fold error
and fractions within 2/5/10-fold when defensible. Units, zero handling, transforms and
calibration are explicit. Boël's published ordinal expression scores support rankings,
not fold error; regulatory activity in Johns is not isolated initiation rate. Calibration
data cannot be reused as independent validation; upstream calibration overlap remains
separate from file splits. No experimental numerical results are reported. See
[benchmark workflow](benchmark_workflow.md) for requirements and source discovery.

## Assembly solution assumptions and browser assets

Tm uses Biopython DNA_NN3 nearest-neighbor DNA/DNA estimates, self-complementarity
detection, and explicit concentrations. Defaults are each strand 25 nM, Na 50 mM,
K/Tris/Mg/dNTP 0 mM, salt correction 5, as documented in the API. Salt method 7 is an
explicit alternative. These are public software defaults, not inferred Gibson reaction
conditions. The documented example produces 60.32 °C. Overlap scoring includes exact
word uniqueness/orientation, GC/Tm preferences, mononucleotide entropy, protected-base
fraction and selected synthesis risks. It does not model shorter off-target homologies
or predict measured assembly efficiency. Reconstruction is an engineering invariant.

React 19.3.0, React DOM 19.3.0 and their resolved scheduler are bundled with their
MIT notices in frontend/public/THIRD_PARTY_NOTICES.txt and the built web package.
No external fonts, CDNs, images, telemetry exporter or prediction service is used.
Optional scientific dependencies retain their separate licensing conditions above;
no project distribution license or commercial redistribution decision is inferred.
