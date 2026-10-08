#!/usr/bin/env bash
# Build from a fresh source archive to avoid setuptools' stale build/lib assets.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
(
    cd frontend
    npm ci
    npm run build
)
package_stage="$(mktemp -d /tmp/dna-design-package.XXXXXXXX)"
export DNA_DESIGN_PACKAGE_STAGE="$package_stage"
.venv/bin/python - <<'PY'
import os
from setuptools.build_meta import build_sdist
build_sdist(os.environ["DNA_DESIGN_PACKAGE_STAGE"])
PY
.venv/bin/python -m pip wheel --no-deps --no-build-isolation "$package_stage"/*.tar.gz -w "$package_stage"
printf '%s\n' "Local source archive and wheel: $package_stage"
