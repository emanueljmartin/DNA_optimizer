"""Cooperative, deterministic progress and cancellation hooks for local jobs."""

from collections.abc import Callable
from dataclasses import dataclass

from .errors import DNAError


class OperationCancelled(DNAError):
    pass


@dataclass(frozen=True)
class ExecutionControl:
    progress: Callable[[int, int, str], None] | None = None
    cancelled: Callable[[], bool] | None = None

    def checkpoint(self, completed: int, total: int, stage: str) -> None:
        if self.cancelled and self.cancelled():
            raise OperationCancelled("Operation cancelled; original project was preserved")
        if self.progress:
            self.progress(completed, total, stage)
