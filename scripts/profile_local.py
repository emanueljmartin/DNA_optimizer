"""Reproducible synthetic engineering timings; no expression/assembly accuracy claim."""

import argparse
import json
import platform
import random
import time
import tracemalloc

from dna_design import Construct, Feature, SearchConfig, load_codon_profile, optimize
from dna_design.gibson import GibsonSettings, plan_gibson
from dna_design.projects import DesignProject
from dna_design.service import AnalysisRequest, analyze


def measure(name, function, memory):
    if memory:
        tracemalloc.start()
    start = time.perf_counter()
    function()
    elapsed = time.perf_counter() - start
    result = {
        "operation": name,
        "seconds": round(elapsed, 3),
    }
    if memory:
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        result["peak_python_mb"] = round(peak / 1_000_000, 3)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--memory", action="store_true", help="Trace Python allocations; adds timing overhead"
    )
    memory = parser.parse_args().memory
    sequence = "ATG" + "GAA" * 1000 + "TAA"
    original = Construct(
        sequence=sequence,
        features=(
            Feature(
                id="cds",
                name="Synthetic CDS",
                feature_type="CDS",
                start=0,
                end=len(sequence),
                strand=1,
            ),
        ),
    )
    profile = load_codon_profile("ecoli")
    rng = random.Random(7)
    assembly_dna = "".join(rng.choices("ACGT", k=6000))
    assembly = Construct(sequence=assembly_dna)
    measurements = [
        measure(
            "3006-base analysis",
            lambda: analyze(AnalysisRequest(project=DesignProject(original=original))),
            memory,
        ),
        measure(
            "3006-base codon search, 100 proposals, seed 7",
            lambda: optimize(
                original, SearchConfig(seed=7, iterations=100, alternatives=3), profile
            ),
            memory,
        ),
        measure(
            "6000-base Gibson, max 2000 bp",
            lambda: plan_gibson(
                assembly, assembly_dna, GibsonSettings(maximum_ordered_length=2000)
            ),
            memory,
        ),
    ]
    print(
        json.dumps(
            {
                "kind": "synthetic engineering measurements",
                "python": platform.python_version(),
                "instrumentation": "tracemalloc; Python allocations only; includes overhead"
                if memory
                else "Wall clock without allocation tracing",
                "measurements": measurements,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
