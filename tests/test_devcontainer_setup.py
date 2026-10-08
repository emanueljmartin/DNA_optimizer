"""Setup must preserve existing environments and fail before installing into broken ones."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

SETUP = Path(__file__).resolve().parents[1] / ".devcontainer" / "setup.sh"


def run_setup(
    tmp_path: Path, *, existing: bool, usable: bool = True, install_fails: bool = False
) -> tuple[subprocess.CompletedProcess[str], Path, Path, Path]:
    """Fake Python commands exercise setup without downloads or touching the real .venv."""
    interpreter = tmp_path / "test-python"
    interpreter.write_text(
        "#!/usr/bin/env bash\n"
        'if [[ "${1:-}" == "-c" ]]; then\n'
        f"    exit {0 if usable else 1}\n"
        "fi\n"
        'printf "%s\\n" "$*" >> "$SETUP_TEST_INSTALL_LOG"\n'
        f"exit {17 if install_fails else 0}\n"
    )
    interpreter.chmod(0o755)
    commands = tmp_path / "commands"
    commands.mkdir()
    python3 = commands / "python3"
    python3.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        'printf "%s\\n" "$*" >> "$SETUP_TEST_CREATE_LOG"\n'
        "mkdir -p .venv/bin\n"
        'cp "$SETUP_TEST_INTERPRETER" .venv/bin/python\n'
    )
    python3.chmod(0o755)
    venv_python = tmp_path / ".venv" / "bin" / "python"
    if existing:
        venv_python.parent.mkdir(parents=True)
        venv_python.write_bytes(interpreter.read_bytes())
        venv_python.chmod(0o755)
        (tmp_path / ".venv" / "keep-existing").write_text("existing environment")
    create_log = tmp_path / "create.log"
    install_log = tmp_path / "install.log"
    result = subprocess.run(
        ["bash", str(SETUP)],
        cwd=tmp_path,
        env={
            **os.environ,
            "PATH": f"{commands}{os.pathsep}{os.environ['PATH']}",
            "SETUP_TEST_INTERPRETER": str(interpreter),
            "SETUP_TEST_CREATE_LOG": str(create_log),
            "SETUP_TEST_INSTALL_LOG": str(install_log),
        },
        capture_output=True,
        text=True,
        check=False,
    )
    return result, create_log, install_log, venv_python


def test_setup_creates_copied_environment_and_installs_declared_dependencies(
    tmp_path: Path,
) -> None:
    result, create_log, install_log, venv_python = run_setup(tmp_path, existing=False)
    assert result.returncode == 0, result.stderr
    assert create_log.read_text().splitlines() == ["-m venv --copies .venv"]
    assert venv_python.is_file()
    assert install_log.read_text().splitlines() == [
        "-m pip install -r requirements-dev.lock",
        "-m pip install --no-build-isolation -e .[dev]",
    ]


def test_setup_reuses_existing_environment_without_recreating_it(tmp_path: Path) -> None:
    result, create_log, install_log, venv_python = run_setup(tmp_path, existing=True)
    assert result.returncode == 0, result.stderr
    assert not create_log.exists()
    assert venv_python.read_bytes() == (tmp_path / "test-python").read_bytes()
    assert (tmp_path / ".venv" / "keep-existing").read_text() == "existing environment"
    assert install_log.exists()


def test_setup_preserves_broken_environment_and_reports_recovery(tmp_path: Path) -> None:
    result, create_log, install_log, venv_python = run_setup(tmp_path, existing=True, usable=False)
    assert result.returncode != 0
    assert "Move the old environment aside" in result.stderr
    assert not create_log.exists()
    assert not install_log.exists()
    assert venv_python.read_bytes() == (tmp_path / "test-python").read_bytes()
    assert (tmp_path / ".venv" / "keep-existing").exists()


def test_setup_stops_when_dependency_installation_fails(tmp_path: Path) -> None:
    result, _, install_log, _ = run_setup(tmp_path, existing=True, install_fails=True)
    assert result.returncode == 17
    assert install_log.read_text().splitlines() == ["-m pip install -r requirements-dev.lock"]


@pytest.mark.parametrize("entry_type", ["file", "dangling-symlink"])
def test_setup_does_not_overwrite_unusable_environment_path(
    tmp_path: Path, entry_type: str
) -> None:
    environment = tmp_path / ".venv"
    if entry_type == "file":
        environment.write_text("preserve this file")
    else:
        environment.symlink_to(tmp_path / "missing-environment")
    result, create_log, install_log, _ = run_setup(tmp_path, existing=False)
    assert result.returncode != 0
    assert "Move the old environment aside" in result.stderr
    assert not create_log.exists()
    assert not install_log.exists()
    if entry_type == "file":
        assert environment.read_text() == "preserve this file"
    else:
        assert environment.is_symlink()


def test_setup_rejects_non_virtualenv_interpreter_even_with_python_optimization(
    tmp_path: Path,
) -> None:
    interpreter = tmp_path / ".venv" / "bin" / "python"
    interpreter.parent.mkdir(parents=True)
    interpreter.symlink_to(sys.executable)
    result = subprocess.run(
        ["bash", str(SETUP)],
        cwd=tmp_path,
        env={**os.environ, "PYTHONOPTIMIZE": "1"},
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "Interpreter is not in a virtual environment" in result.stderr
    assert "pip" not in result.stderr
    assert interpreter.is_symlink()
