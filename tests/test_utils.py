import sys

import pytest
from pytest_mock import MockerFixture

from sequana_pipelines.lora.src.utils import (
    base_version_parser,
    ccs_version_parser,
    checkm_version_parser,
    get_apptainer_version,
    get_tools_versions,
    run_tools_versions,
    run_version,
    seqtk_version_parser,
    strip_version_parser,
    unicycler_version_parser,
)


@pytest.mark.asyncio
async def test_run_version():
    _, version = await run_version("python", "python --version", False, {})
    assert version == f"{sys.version.split()[0]}"


@pytest.mark.asyncio
async def test_run_version_tool_not_found():
    _, version = await run_version("( ͡° ͜ʖ ͡°)", "lora_is_pretty_good --version", False, {})
    assert version == "Tool not found"


@pytest.mark.asyncio
async def test_run_version_not_found(mocker: MockerFixture):
    async def mock_create_subprocess_shell(cmd, stdout, stderr):
        class MockProcess:
            """Mock Process object return by create_subprocess_shell"""

            def __init__(self):
                self.returncode = 1

            async def communicate(self):
                return "( ͡° ͜ʖ ͡°)", "( ͡° ʖ̯ ͡°)"

        return MockProcess()

    mocker.patch("sequana_pipelines.lora.src.utils.asyncio.create_subprocess_shell", mock_create_subprocess_shell)
    _, version = await run_version("( ͡° ͜ʖ ͡°)", "lora_is_pretty_good --version", False, {})
    assert version == "No version found"


@pytest.mark.asyncio
async def test_run_tools_versions():
    used_tools = {"( ͡° ͜ʖ ͡°)": "lora_is_pretty_good --version", "python": "python --version"}
    results = await run_tools_versions(used_tools, {})
    assert results[0][1] in ["No version found", "Tool not found"]
    assert results[1][1] == f"{sys.version.split()[0]}"


@pytest.mark.asyncio
async def test_run_get_apptainer_version(mocker: MockerFixture):
    mocker.patch("sequana_pipelines.lora.src.utils.shutil.which", return_value="ola")
    used_tools = {"( ͡° ͜ʖ ͡°)": "lora_is_pretty_good --version"}
    results = await run_tools_versions(used_tools, {"( ͡° ͜ʖ ͡°)": "/path/to/lenny_418.0.img"})
    assert results[0][1] == "418.0"


def test_get_tools_versions(mocker: MockerFixture):
    # create fake lora config
    config = {
        "assembler": "flye",
        "ccs": {"do": False},
        "blast": {"do": True},
        "busco": {"do": True},
        "canu_correction": {"do": False},
        "circlator": {"do": False},
        "prokka": {"do": False},
        "sequana_coverage": {"do": True},
    }

    async def mock_run_version(tool, cmd, has_apptainer, apptainers):
        return tool, "1337"

    mocker.patch("sequana_pipelines.lora.src.utils.run_version", mock_run_version)
    results = get_tools_versions(config)
    expected_keys = {
        "flye",
        "blast",
        "busco",
        "minimap2",
        "sequana_coverage",
        "seqkit",
        "samtools",
        "quast",
        "seqtk",
        "multiqc",
    }
    for tool, _ in results:
        assert tool in expected_keys


def test_get_tools_versions_with_ccs(mocker: MockerFixture):
    """ccs.do=True exercises the dict-of-subtools branch in _iter_optional_tools."""
    config = {
        "assembler": "flye",
        "ccs": {"do": True},
        "blast": {"do": False},
        "busco": {"do": False},
        "canu_correction": {"do": False},
        "circlator": {"do": False},
        "prokka": {"do": False},
        "sequana_coverage": {"do": False},
    }

    async def mock_run_version(tool, cmd, has_apptainer, apptainers):
        return tool, "1337"

    mocker.patch("sequana_pipelines.lora.src.utils.run_version", mock_run_version)
    results = get_tools_versions(config)
    tool_names = {tool for tool, _ in results}
    assert "ccs" in tool_names
    assert "pbindex" in tool_names


# ---------------------------------------------------------------------------
# Version parsers
# ---------------------------------------------------------------------------


def test_base_version_parser():
    assert base_version_parser("tool 1.2.3\nsome extra line") == "1.2.3"


def test_strip_version_parser():
    assert strip_version_parser("  2.9.1  \n") == "2.9.1"


def test_seqtk_version_parser():
    # seqtk prints "Usage" on line 0 and "Version: 1.3-r106" on line 1
    output = "Usage: seqtk <command>\nVersion: 1.3-r106\n"
    assert seqtk_version_parser(output) == "1.3-r106"


def test_unicycler_version_parser():
    assert unicycler_version_parser("Unicycler 0.5.0\nmore info") == "0.5.0"


def test_ccs_version_parser():
    assert ccs_version_parser("ccs 6.4.0 (commit ...)\nmore") == "6.4.0"


def test_checkm_version_parser():
    # checkm prints a header line, then "checkm version X.Y.Z" on line 1
    output = "----------- CheckM -----------\ncheckm version 1.2.2\n"
    assert checkm_version_parser(output) == "1.2.2"


# ---------------------------------------------------------------------------
# get_apptainer_version
# ---------------------------------------------------------------------------


def test_get_apptainer_version_direct():
    result = get_apptainer_version("flye", {"flye": "/path/to/flye_2.9.img"})
    assert result == "2.9"


def test_get_apptainer_version_correspondance():
    # canu_correction maps to the "canu" apptainer key
    result = get_apptainer_version("canu_correction", {"canu": "/path/to/canu_2.2.img"})
    assert result == "2.2"


def test_get_apptainer_version_missing():
    result = get_apptainer_version("unknown_tool", {})
    assert result == "No version found"


# ---------------------------------------------------------------------------
# run_version edge cases
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_run_version_returncode_none(mocker: MockerFixture):
    async def mock_create_subprocess_shell(cmd, stdout, stderr):
        class MockProcess:
            def __init__(self):
                self.returncode = None

            async def communicate(self):
                return b"", b""

        return MockProcess()

    mocker.patch("sequana_pipelines.lora.src.utils.asyncio.create_subprocess_shell", mock_create_subprocess_shell)
    _, version = await run_version("sometool", "sometool --version", False, {})
    assert version == "Tool not found"


@pytest.mark.asyncio
async def test_run_version_index_error(mocker: MockerFixture):
    """Empty stdout triggers IndexError in base_version_parser → 'No version found'."""

    async def mock_create_subprocess_shell(cmd, stdout, stderr):
        class MockProcess:
            def __init__(self):
                self.returncode = 0

            async def communicate(self):
                return b"", b""

        return MockProcess()

    mocker.patch("sequana_pipelines.lora.src.utils.asyncio.create_subprocess_shell", mock_create_subprocess_shell)
    # "samtools" is not in VERSION_PARSER so it falls back to base_version_parser;
    # empty stdout causes IndexError → "No version found"
    _, version = await run_version("samtools", "samtools --version", False, {})
    assert version == "No version found"


@pytest.mark.asyncio
async def test_run_version_apptainer_fallback(mocker: MockerFixture):
    """returncode 127 with apptainer present → reads version from apptainer path."""
    mocker.patch("sequana_pipelines.lora.src.utils.shutil.which", return_value="/usr/bin/apptainer")

    async def mock_create_subprocess_shell(cmd, stdout, stderr):
        class MockProcess:
            def __init__(self):
                self.returncode = 127

            async def communicate(self):
                return b"", b"command not found"

        return MockProcess()

    mocker.patch("sequana_pipelines.lora.src.utils.asyncio.create_subprocess_shell", mock_create_subprocess_shell)
    _, version = await run_version("flye", "flye --version", True, {"flye": "/path/to/flye_2.9.img"})
    assert version == "2.9"
