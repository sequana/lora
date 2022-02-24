import pytest
import sys
from pytest_mock import MockerFixture

from sequana_pipelines.lora.src.utils import get_tools_versions, run_tools_versions, run_version


@pytest.mark.asyncio
async def test_run_version():
    tool, version = await run_version("python", "python --version")
    assert version == f"{sys.version.split()[0]}"


@pytest.mark.asyncio
async def test_run_version_tool_not_found():
    tool, version = await run_version("( ͡° ͜ʖ ͡°)", "lora_is_pretty_good --version")
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
    tool, version = await run_version("( ͡° ͜ʖ ͡°)", "lora_is_pretty_good --version")
    assert version == "No version found"


@pytest.mark.asyncio
async def test_run_tools_versions():
    used_tools = {"( ͡° ͜ʖ ͡°)": "lora_is_pretty_good --version", "python": "python --version"}
    results = await run_tools_versions(used_tools)
    assert results[0][1] == "Tool not found"
    assert results[1][1] == f"{sys.version.split()[0]}"


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

    async def mock_run_version(tool, cmd):
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
        "multiqc"
    }
    for tool, _ in results:
        assert tool in expected_keys
