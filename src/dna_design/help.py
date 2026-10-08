"""Central scientific help consumed by the browser without duplicated calculations."""

HELP = {
    "codon_adaptation": {
        "title": "Codon reference and CAI",
        "definition": "Relative synonymous usage compared with a versioned "
        "reference. CAI is reported only for an explicit "
        "high-expression reference.",
        "effect": "Weights prefer codons represented in the selected host reference.",
        "limitations": "Genome usage is not an experimentally calibrated "
        "expression predictor; zero counts stay zero.",
        "reference": "https://doi.org/10.1093/nar/15.3.1281",
    },
    "rare_threshold": {
        "title": "Rare-codon threshold",
        "definition": "Codons with relative synonymous weight below the selected threshold.",
        "effect": "Affects reported rare frequency and the optional rare-codon loss.",
        "limitations": "No universal threshold is supplied; values are explicit "
        "engineering choices.",
        "reference": "docs/scientific_provenance.md",
    },
    "gc": {
        "title": "Global and local GC",
        "definition": "Fraction of concrete DNA bases that are G or C, globally or in a "
        "specified sliding window.",
        "effect": "Selected ranges or targets guide composition preferences.",
        "limitations": "GC alone does not predict expression or synthesis acceptance; vendor "
        "guidance is product-specific.",
        "reference": "docs/scientific_provenance.md",
    },
    "locked": {
        "title": "Protected DNA",
        "definition": "LOCKED bases must be identical to the original DNA. Unannotated "
        "bases default to locked.",
        "effect": "All overlapping rules compose and all original CDS proteins must remain "
        "identical.",
        "limitations": "Changing annotations starts a new reference and clears previous "
        "design results.",
        "reference": "docs/PRODUCT_SPEC.md",
    },
    "tir": {
        "title": "Target relative TIR",
        "definition": "Local OSTIR model output estimating relative translation initiation, "
        "with energy components.",
        "effect": "RBS search minimizes absolute log-fold distance to a positive target.",
        "limitations": "Not measured absolute protein expression. E. coli calibration is "
        "upstream; PAO1 transfer requires explicit sourced parameters and "
        "opt-in.",
        "reference": "https://doi.org/10.21105/joss.03362",
    },
    "initiation": {
        "title": "Initiation context protection",
        "definition": "Default codon search protects the first 30 coding bases. RBS "
        "design carries a lock over its complete supplied prediction "
        "context.",
        "effect": "Later searches automatically preserve that designed context in this app.",
        "limitations": "The 30-nt default is a product setting, not a universal "
        "biological optimum.",
        "reference": "docs/PRODUCT_SPEC.md",
    },
    "synthesis": {
        "title": "Synthesis findings",
        "definition": "Explicit GC, repeat, homopolymer and entropy policies with "
        "per-finding observations.",
        "effect": "Repair optimizes selected risks while preserving original biological "
        "constraints.",
        "limitations": "Unconfigured criteria remain unassessed. Public vendor profiles "
        "are incomplete approximations; residual risks can remain.",
        "reference": "docs/scientific_provenance.md",
    },
    "overlap": {
        "title": "Gibson overlaps and ordered length",
        "definition": "Adjacent ordered fragments share exact overlaps. The maximum "
        "length includes those overlaps.",
        "effect": "The planner starts at the geometric minimum count and chooses valid "
        "boundaries; circular plans include a closing overlap.",
        "limitations": "Exact reconstruction and unique full overlaps do not prove "
        "wet-lab assembly yield or eliminate shorter off-target "
        "homologies.",
        "reference": "docs/PRODUCT_SPEC.md",
    },
    "tm": {
        "title": "Overlap melting temperature",
        "definition": "Biopython DNA_NN3 DNA/DNA nearest-neighbor estimate under displayed "
        "strand/salt concentrations.",
        "effect": "An optional preferred temperature influences boundary scoring.",
        "limitations": "Default solution values come from Biopython, not a measured Gibson "
        "reaction mix; conditions must match the intended experiment.",
        "reference": "https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html",
    },
    "validation": {
        "title": "Independent final validation",
        "definition": "Rechecks the final DNA against original proteins, locks, "
        "motifs, masks, carried initiation locks and any fragment plan.",
        "effect": "Failed final validation blocks all final-design exports; draft "
        "snapshots remain available.",
        "limitations": "These integrity checks are distinct from experimental "
        "expression or assembly reliability validation.",
        "reference": "docs/testing_strategy.md",
    },
}
