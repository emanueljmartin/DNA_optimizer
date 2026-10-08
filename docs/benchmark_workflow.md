# Local experimental evaluation workflow

The framework is implemented; no experimental accuracy study has been run. We will
agree on datasets, host/construct conditions, exclusions, training overlap and held-out
splits with the user before making reliability claims. Synthetic tests establish
arithmetic and integrity only. Do not require identical synonymous DNA from another
optimizer as evidence of correctness.

Use local files. The app never downloads benchmark DNA or sends it to a service.
The evaluator consumes **already generated predictions** with their explicit model
version; it does not label codon adaptation as calibrated expression. Generate those
predictions with the selected local domain model and retain its full configuration.

## Versioned mapped loader

`public-mapped-csv-v1` accepts public source tables after explicit local preparation,
or user data in the same schema. Native paper spreadsheets/PDFs are not silently
interpreted. Map verified source columns in the manifest; record every conversion,
sequence normalization, measurement aggregation and exclusion. Do not infer units,
license, host, context or training labels from a paper title.

Required CSV columns (names may be explicitly mapped):

| Column | Meaning |
| --- | --- |
| id | Unique observation ID; replicates retain separate IDs |
| sequence | Exact canonical A/C/G/T model-input sequence; no inferred flanking context |
| group | Related-construct/replicate leakage group, as defined in the manifest |
| host | Actual assayed host/strain label |
| condition | Assay/normalization condition; evaluations are separated by host and condition |
| split | `calibration`, `heldout` or `excluded`, chosen before evaluation |
| measured | Finite nonnegative observed value; its meaning/units must be declared |
| predicted | Finite nonnegative output from the declared model/version |

The manifest's schema is generated from the backend:
[benchmark_manifest.schema.json](benchmark_manifest.schema.json). It requires dataset
ID/version, source, license and license source, access date, exact SHA-256 of the UTF-8
CSV text, experimental/synthetic evidence label, measured and predicted definitions/
units, model version, split policy and group definition. Optional fields specify column
mapping, transforms, exclusions, ratio/ordinal scale, comparable units, calibration and
known training fingerprints/groups. Use the checksum of the exact text you submit;
newline/encoding changes alter it. UTF-8 BOMs and implicit sequence normalization are
not performed by the loader.

Source candidates from the specification are exposed by `/api/benchmarks/catalog`:

| Candidate | Primary source | Role and caution |
| --- | --- | --- |
| Kudla 2009 | [Author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC3902468/) | Synonymous GFP ranking; review exact measurements and sequence context |
| Boël 2016 | [Paper and supplementary CSV link](https://pmc.ncbi.nlm.nih.gov/articles/PMC5054687/) | Heterologous expression ranking; published scores are ordinal, not fold changes |
| Salis 2009 | [Paper](https://doi.org/10.1038/nbt.1568) | RBS calibration/held-out subsets; OSTIR documents upstream use for calibration |
| Farasat 2014 | [Paper and supplements](https://link.springer.com/article/10.15252/msb.20134955) | Host/context-specific RBS evaluation; raw table meanings require review |
| Johns 2018 | [Paper](https://doi.org/10.1038/nmeth.4633) | Species-selective regulatory activity; includes transcription, not isolated TIR |

The registry gives discovery/proposed roles, not approved experimental splits or a
blanket data license. Exact raw records/headers and reuse rights have not been validated
for this project. No paper tables are bundled or fabricated. Each study's manifest must
record the exact prepared data checksum and verified license/source.

## Metrics and separation

Spearman uses average tied ranks; linear Pearson describes association. Constant or
insufficient pairs return unassessed correlation. Ratio-scale positive pairs additionally
support natural-log Pearson. Zero pairs remain in rankings and are explicitly excluded
from log/fold metrics, with no pseudocount. Negative/nonfinite input is rejected.

Median symmetric fold error and fractions within 2/5/10-fold require ratio-scale
measurements and defensible comparable units or explicitly separate calibration.
Multiplicative calibration fits the mean log measured/predicted difference using only
at least two positive calibration records per host/condition; it never fits held-out
records. Ordinal Boël scores cannot be assigned ratio-scale fold errors. Reports retain
counts, exclusions, calibration factors, definitions, source manifest and warnings.

Calibration versus held-out **group and exact-sequence** overlap is rejected, including
replicates. Declared model-training sequences/groups cannot appear in held-out records.
Group IDs must capture related variants/proteins as scientifically appropriate; exact
sequence checks alone cannot detect all homology leakage. The file split is separate
from upstream model-training/calibration overlap. Unknown upstream overlap is explicitly
unverified; a documented disjoint claim requires its verification source and remains a
user declaration. Salis2009/OSTIR reuse always carries the upstream-calibration caveat.

## Run and preserve the result

Use **Reliability** to load a local CSV plus manifest and evaluate, or:

```bash
.venv/bin/dna-design benchmark /path/data.csv /path/manifest.json --output /path/report.json
```

The CLI refuses to overwrite an existing report. The browser displays evidence and
independence warnings, record counts and per-host/condition ranking, log and eligible
fold metrics; the full backend JSON/provenance is available to inspect and download.
Changing either input hides the earlier report, including a response still in flight.
Keep the exact inputs, manifest, report and prediction settings together. Do not
commit private/unpublished sequences. The proposed [reliability plan](reliability_plan.md)
sets out the biological decisions to review before the study. A later reliability
study should cover protein/lock integrity, measured expression ranking, RBS target
error, synthesis findings and assembly reconstruction separately, with appropriately
independent experimental evidence where biological accuracy is claimed.
