"""Progress/cancellation leave deterministic science unchanged."""

import pytest

from dna_design.execution import ExecutionControl, OperationCancelled
from dna_design.gibson import GibsonSettings, plan_gibson
from dna_design.optimization import SearchConfig, optimize

from .test_optimization import synthetic_profile


def test_progress_seed_and_bounded_alternatives_have_identical_results(coding_construct):
    reference = coding_construct
    profile = synthetic_profile()
    calls = []
    config = SearchConfig(seed=5, iterations=100, alternatives=2, n_terminal_lock_nt=0)
    original = optimize(reference, config, profile)
    progress = optimize(
        reference,
        config,
        profile,
        execution=ExecutionControl(progress=lambda *args: calls.append(args)),
    )
    all_results = optimize(reference, config.model_copy(update={"alternatives": 1000}), profile)
    assert progress == original
    assert original.designs == all_results.designs[:2]
    assert calls[0][0] == 0 and calls[-1][0] == 100
    assert len(calls) >= 100


def test_cancel_interrupts_optimization_and_assembly_without_changing_original(coding_construct):
    reference = coding_construct
    saved = reference.model_dump_json()
    control = ExecutionControl(cancelled=lambda: True)
    with pytest.raises(OperationCancelled):
        optimize(reference, SearchConfig(seed=5), synthetic_profile(), execution=control)
    with pytest.raises(OperationCancelled):
        plan_gibson(
            reference,
            reference.sequence,
            GibsonSettings(maximum_ordered_length=100),
            execution=control,
        )
    assert reference.model_dump_json() == saved
