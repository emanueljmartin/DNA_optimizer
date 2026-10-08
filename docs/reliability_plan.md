# Reliability study proposal

Status: prepared for biological review on 2026-10-08. This is a study proposal,
not an executed benchmark or an experimental accuracy claim. The engineering
checkpoint is recorded in [development_plan.md](development_plan.md). The local
data and evaluation contract is in [benchmark_workflow.md](benchmark_workflow.md).

## Questions and evidence

| Question | Data and comparison | What establishes success | Current evidence |
| --- | --- | --- | --- |
| Does redesign preserve the construct's biological rules? | Annotated synthetic cases: both strands, overlaps, alternative codes, locked non-coding DNA, RBS masks, all motif modes, circular origins; adversarial candidate metadata | Every returned/exported design independently satisfies the original controls; any integrity failure blocks release | Automated examples, property tests, independent validation and rejection tests |
| Do consecutive stages preserve earlier decisions? | Codon → synthesis → RBS → subsequent codon → Gibson → draft reload → final exports; replay stored controls, including masks and thermodynamic conditions | Same original proteins/locked bases, carried initiation context, exact fragment reconstruction, exact saved settings sent on rerun | Forward/reverse pipeline and browser regressions; fictional scorer isolates control integrity, real local backend tests separately exercise execution |
| Is the codon score associated with measured expression? | Public or user observations with actual host/strain and assay labels; compare a frozen score with observed expression per condition | Predeclared ranking/association criterion and uncertainty, reviewed before evaluation | Evaluator arithmetic tested; no experimental association measured |
| Does optimization improve expression? | Paired original/designed constructs for multiple proteins under matched experimental conditions; appropriate controls | Predeclared biological improvement and acceptable variability; independent observations | No paired experimental study; a higher adaptation score alone does not answer this |
| Does RBS prediction rank initiation rates and reach requested targets? | Full, verified upstream/CDS model inputs with measured outcomes; held-out constructs separate from model calibration | Predeclared rank/target-error criteria, with fold errors only for defensible ratio-scale units or separate calibration | Real local model execution and synthetic target-search integrity; no project experimental accuracy result |
| Are synthesis findings useful? | Defined defects plus independently recorded vendor outcomes for the chosen product/policy | Correct defect detection/repair with preserved biology; vendor outcome criterion agreed separately | Synthetic stress tests; public guidance is advisory and incomplete |
| Will fragments assemble as intended? | Exact in-silico reconstruction plus separately documented experimental assemblies | Exact software reconstruction and an agreed experimental success measure | Reconstruction, length/overlap/Tm checks; no wet-lab yield result |

No synonymous-sequence identity requirement is imposed when comparing another
optimizer. Many different DNA sequences can satisfy the same biological rules.
Engineering invariants require exact correctness; expression and assembly performance
require observations and uncertainty appropriate to the chosen question.

## Prepare each dataset before scoring

1. Record the exact source/version, reuse license, raw files and checksums. Review
   actual columns and units, strain, promoter/vector, CDS boundaries, flanking
   sequence, assay and replicates. Keep missing-context cases explicit; never
   invent upstream DNA or infer the assayed strain from our codon-profile choice.
2. Define the prediction being tested. Genome-wide adaptation is a score, not an
   absolute expression prediction. RBS output is a relative model TIR. Keep
   fluorescence, soluble protein, transcript abundance and isolated initiation
   claims distinct according to what the experiment measured.
3. Inventory upstream training/calibration records for the exact model version.
   Exclude known overlap from an independent validation claim. If upstream
   overlap cannot be established, label independence unverified.
4. Choose leakage groups before evaluation: at least technical/biological
   replicates and related sequence variants; protein/family grouping where the
   question requires generalization to unseen proteins. The loader checks exact
   sequences and declared groups; it cannot infer biological relatedness.
5. Freeze calibration/development and held-out groups, exclusion reasons,
   normalization and model settings. Determine sample sizes/replicates and an
   uncertainty method from the intended effect and assay variability; none is
   inferred here. Calibrate only on the calibration subset, by actual host and
   condition. No tuning on held-out outcomes.
6. Generate predictions locally with retained model version, complete context,
   settings and seed where relevant. The present evaluator reads supplied
   predictions; it does not reconstruct missing experimental context or run a
   new biological model from a paper table automatically.
7. Run the frozen evaluation and preserve inputs, manifest, exclusions, prediction
   settings, quantitative report and original/final projects. Report failures,
   zero observations and unassessed metrics. If model/settings change after
   seeing results, version the study and disclose that the same held-out set has
   been inspected.

Use local study files outside Git for unpublished/private sequences. Only sourced
public or clearly synthetic fixtures belong in the repository.

## Public candidates and host applicability

The specification's candidate index and primary links are in
[benchmark_workflow.md](benchmark_workflow.md). They are candidates, not selected
study datasets. Kudla/Boël concern expression; Salis/Farasat concern RBS-related
evaluation; Johns requires separating broader regulatory activity from initiation.
Raw sequence context, measurement meanings, reuse rights and upstream overlap still
require study-specific review.

One applicability issue was verified in the primary Boël paper: its large expression
dataset uses BL21(DE3) with pMGK and a T7/pET expression system, with ordinal 0–5 gel
scores. Thus it can inform a scoped ranking study, but cannot establish MG1655 or
PAO1 performance or ratio-scale fold accuracy. These are study-design implications
of the reported conditions, not an executed evaluation.
[Boël et al., Methods and large-scale analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC5054687/).

Keep PAO1 evaluation separate. The tool has a sourced PAO1 codon profile, but no
project-validated PAO1 RBS calibration. Supplying an anti-SD and choosing uncalibrated
transfer provides a computation, not evidence of accuracy in that host.

## Biological decisions to agree before experimental evaluation

* First intended use: expression ranking, paired optimization benefit, RBS target
  design, synthesis ordering, or assembly; a separate criterion for each claim.
* Actual first host/strain, construct context and assay; which public observations
  and/or local paired data are relevant to that use.
* Acceptable ranking/target error or improvement, and the uncertainty required for
  a practical decision. Existing 2/5/10-fold summaries are descriptive, not preset
  acceptance thresholds.
* Related-construct grouping, calibration/held-out allocation, replicate handling
  and independent upstream-training verification.

Next review should decide these biological questions. Implementation details and
automated verification remain the agent's responsibility. No experimental study
split, performance target or host calibration has been chosen on the user's behalf.
