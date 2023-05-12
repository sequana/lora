import os

import pytest

from sequana_pipelines.lora.src.report import create_reports, create_summary, get_quast_information

from . import test_dir


@pytest.fixture(scope="module")
def config():
    yield {
        "assembler": "flye",
        "ccs": {"do": False},
        "blast": {"do": True},
        "busco": {"do": True},
        "canu_correction": {"do": False},
        "circlator": {"do": False},
        "prokka": {"do": False},
        "sequana_coverage": {"do": True},
    }


def test_lora_create_lora_summary(mocker, tmpdir, config):
    """Test LORA summary report creation"""

    async def mock_run_version(tool, cmd, has_apptainer, apptainers):
        return tool, "1337"

    mocker.patch("sequana_pipelines.lora.src.utils.run_version", mock_run_version)

    lora_dir = test_dir / "resources"
    summary = tmpdir.join("summary.html")
    samples = ["toto", "tata"]
    create_summary(summary, "lora.html", get_quast_information(samples, lora_dir), config, lora_dir)
    assert os.path.exists(summary)


def test_lora_create_lora_report(mocker, tmpdir, config):
    """Test LORA report creation"""

    def mock_create_summary(summary_name, lora_name, quast_info, config, lora_dir):
        pass

    mocker.patch("sequana_pipelines.lora.src.report.create_summary", mock_create_summary)
    lora_dir = test_dir / "resources"
    lora_report = tmpdir.join("lora.html")
    samples = ["toto", "tata"]
    create_reports("summary.html", lora_report, samples, config, lora_dir=lora_dir)
    assert os.path.exists(lora_report)


def test_lora_create_both_reports(mocker, tmpdir, config):
    """Test LORA summary report creation"""

    async def mock_run_version(tool, cmd, has_apptainer, apptainers):
        return tool, "1337"

    mocker.patch("sequana_pipelines.lora.src.utils.run_version", mock_run_version)

    lora_dir = test_dir / "resources"
    summary = tmpdir.join("summary.html")
    lora_report = tmpdir.join("lora.html")
    samples = ["toto", "tata"]
    create_reports(summary, lora_report, samples, config, lora_dir)
    assert os.path.exists(summary)
    assert os.path.exists(lora_report)
