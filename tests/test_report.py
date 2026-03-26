import os
from collections import defaultdict

import pytest

from sequana_pipelines.lora.src.report import (
    _iter_value_to_float,
    create_reports,
    create_summary,
    fill_blast_result,
    fill_sequana_coverage,
    get_busco_information,
    get_quast_information,
)

from . import test_dir


@pytest.fixture(scope="module")
def config():
    yield {
        "assembler": "flye",
        "ccs": {"do": False},
        "blast": {"do": True},
        "busco": {"do": True},
        "checkm": {"do": False},
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


def test_lora_create_reports_no_blast_no_busco(mocker, tmpdir, config):
    """Reports with blast and busco disabled should still succeed."""
    cfg = {**config, "blast": {"do": False}, "busco": {"do": False}}

    async def mock_run_version(tool, cmd, has_apptainer, apptainers):
        return tool, "1337"

    mocker.patch("sequana_pipelines.lora.src.utils.run_version", mock_run_version)

    lora_dir = test_dir / "resources"
    summary = tmpdir.join("summary.html")
    lora_report = tmpdir.join("lora.html")
    samples = ["toto", "tata"]
    create_reports(summary, lora_report, samples, cfg, lora_dir)
    assert os.path.exists(summary)
    assert os.path.exists(lora_report)


# ---------------------------------------------------------------------------
# Unit tests for individual report helpers
# ---------------------------------------------------------------------------


def test_iter_value_to_float_integer():
    assert list(_iter_value_to_float(["42"])) == [42]


def test_iter_value_to_float_float():
    result = list(_iter_value_to_float(["3.14159"]))
    assert result == ["3.142"]


def test_iter_value_to_float_string_passthrough():
    assert list(_iter_value_to_float(["N/A"])) == ["N/A"]


def test_iter_value_to_float_mixed():
    result = list(_iter_value_to_float(["10", "2.5", "na"]))
    assert result[0] == 10
    assert result[1] == "2.500"
    assert result[2] == "na"


def test_get_quast_information():
    lora_dir = test_dir / "resources"
    samples = ["tata", "toto"]
    results = get_quast_information(samples, lora_dir)
    assert set(results.keys()) == {"tata", "toto"}
    # QUAST_KEY has 5 entries → 5 values per sample
    for sample, values in results.items():
        assert len(values) == 5
    # Mapped (%) must be capped at 100.0
    from sequana_pipelines.lora.src.enums import QUAST_KEY

    idx = QUAST_KEY.index("Mapped (%)")
    for values in results.values():
        assert float(values[idx]) <= 100.0


def test_get_busco_information():
    lora_dir = test_dir / "resources"
    samples = ["tata", "toto"]
    results = list(get_busco_information(samples, lora_dir))
    assert len(results) == 2
    for sample, values in results:
        assert sample in samples
        # BUSCO_KEY has 4 entries
        assert len(values) == 4


def test_fill_blast_result():
    analysis = defaultdict(lambda: defaultdict(dict))
    lora_dir = test_dir / "resources"
    samples = ["tata", "toto"]
    fill_blast_result(analysis, samples, lora_dir)
    # tata.tsv has hits for tig00000001
    assert "tig00000001" in analysis["tata"]
    assert "blast" in analysis["tata"]["tig00000001"]
    assert isinstance(analysis["tata"]["tig00000001"]["blast"], list)


def test_fill_sequana_coverage():
    analysis = defaultdict(lambda: defaultdict(dict))
    lora_dir = test_dir / "resources"
    samples = ["tata", "toto"]
    fill_sequana_coverage(analysis, samples, lora_dir)
    for sample in samples:
        assert "tig00000001" in analysis[sample]
        cov = analysis[sample]["tig00000001"]["coverage"]
        assert "length" in cov
        assert "DOC" in cov
        assert "BOC" in cov
        assert "cov_image" in analysis[sample]["tig00000001"]
