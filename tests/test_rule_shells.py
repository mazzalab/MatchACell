"""Static checks on the shell: blocks of the Snakemake rules."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RULES = sorted((ROOT / "workflow/rules").glob("*.smk"))
SHELL = re.compile(r'shell:\s*r?"""(.*?)"""', re.S)
ONLY_PARAMS = re.compile(r"^(\{params\.\w+\}\s*)+$")


@pytest.mark.parametrize("rule_file", RULES, ids=lambda p: p.name)
def test_shell_does_not_end_on_an_optional_flag_line(rule_file):
    # Snakemake strips trailing whitespace from the formatted command: when the
    # last line holds only flags that render empty, the command ends in "\",
    # which bash passes to the script as a literal argument.
    for block in SHELL.findall(rule_file.read_text()):
        lines = [ln.strip() for ln in block.strip().splitlines()]
        if len(lines) >= 2:
            assert not (ONLY_PARAMS.match(lines[-1]) and lines[-2].endswith("\\")), \
                f"{rule_file.name}: put optional flags on the last argument's line"
