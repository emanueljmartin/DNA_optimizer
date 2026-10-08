#!/usr/bin/env bash
# Local full-tool installer; preserves a healthy .venv and never deploys.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
bash .devcontainer/setup.sh
.venv/bin/python -m pip install -r requirements-app.lock -r requirements-rbs.lock
.venv/bin/python -m pip install --no-build-isolation -e '.[dev,app,rbs]'
if ! command -v node >/dev/null || ! command -v npm >/dev/null; then
    printf '%s\n' 'Install Node 24 and npm to build the browser.' >&2
    exit 1
fi
node -e 'if (Number(process.versions.node.split(".")[0]) < 24) { console.error("Node >=24 required"); process.exit(1); }'
(
    cd frontend
    npm ci
    npm run build
)
.venv/bin/python -m pip check
printf '%s\n' 'Ready. Run .venv/bin/dna-design serve and open http://127.0.0.1:8000.'
