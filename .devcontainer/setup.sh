#!/usr/bin/env bash
set -euo pipefail
if [[ ! -e .venv && ! -L .venv ]]; then
    python3 -m venv --copies .venv
fi
if ! .venv/bin/python -c '
import sys
from pathlib import Path
if sys.version_info < (3, 12):
    raise SystemExit("Python >=3.12 required")
if sys.prefix == sys.base_prefix:
    raise SystemExit("Interpreter is not in a virtual environment")
if Path(sys.prefix).resolve() != Path(".venv").resolve():
    raise SystemExit("Wrong virtual environment")
'; then
    printf '%s\n' \
        'Existing .venv is unusable or requires Python >=3.12.' \
        'Move the old environment aside, then recreate it with python3 -m venv --copies .venv.' >&2
    exit 1
fi
.venv/bin/python -m pip install -r requirements-dev.lock
.venv/bin/python -m pip install --no-build-isolation -e '.[dev]'
