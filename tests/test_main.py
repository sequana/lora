import subprocess

import pytest
import yaml
from click.testing import CliRunner

import sequana_pipelines.lora.main as main

from . import test_dir


def test_standalone_subprocess(tmpdir):

    input_dir = test_dir
    cmd = [
        "test",
        "--input-directory",
        str(input_dir),
        "--working-directory",
        str(tmpdir),
        "--force",
        "--assembler",
        "flye",
        "--genome-size",
        "1m",
    ]
    subprocess.call(cmd)


def test_standalone_script(tmpdir):
    input_dir = test_dir
    runner = CliRunner()

    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(input_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--data-type",
            "pacbio-corr",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
        ],
    )

    assert results.exit_code == 0


def test_standalone_script_nanopore(tmpdir):
    input_dir = test_dir
    runner = CliRunner()

    args = [
        "--input-directory",
        str(input_dir),
        "--working-directory",
        str(tmpdir),
        "--force",
        "--data-type",
        "nano-hq",
        # "--mode",
        # "eukaryota",
        "--assembler",
        "flye",
        "--genome-size",
        "1m",
    ]
    results = runner.invoke(main.main, args)
    assert results.exit_code == 0


def test_busco_lineage_not_overridden_by_mode(tmpdir):
    """--busco-lineage must take precedence over the default set by --mode."""
    runner = CliRunner()

    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--mode",
            "bacteria",
            "--busco-lineage",
            "streptomyces",
        ],
    )
    assert results.exit_code == 0

    import yaml

    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert (
        "streptomyces" in config["busco"]["lineage"]
    ), f"Expected streptomyces lineage but got: {config['busco']['lineage']}"


def test_version():
    cmd = ["sequana_lora", "--version"]
    subprocess.call(cmd)


def test_help():
    cmd = ["sequana_lora", "--help"]
    subprocess.call(cmd)


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("data_type", ["nano-raw", "nano-corr", "pacbio-raw", "pacbio-hifi"])
def test_data_type_variants(tmpdir, data_type):
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--data-type",
            data_type,
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
        ],
    )
    assert results.exit_code == 0


# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------


def test_mode_bacteria_sets_prokka_and_coverage(mocker, tmpdir):
    mocker.patch(
        "sequana_pipelines.lora.main.get_busco_lineages_and_urls",
        return_value={"bacteria": "http://example.com/bacteria_odb12.tar.gz"},
    )
    mocker.patch("sequana_pipelines.lora.main.download_and_extract_tar_gz")
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--mode",
            "bacteria",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert config["prokka"]["do"] is True
    assert config["sequana_coverage"]["do"] is True


def test_mode_eukaryota_sets_quast_options(mocker, tmpdir):
    mocker.patch(
        "sequana_pipelines.lora.main.get_busco_lineages_and_urls",
        return_value={"eukaryota": "http://example.com/eukaryota_odb12.tar.gz"},
    )
    mocker.patch("sequana_pipelines.lora.main.download_and_extract_tar_gz")
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--mode",
            "eukaryota",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert "-e" in config["quast"]["options"]


# ---------------------------------------------------------------------------
# Individual flags
# ---------------------------------------------------------------------------


def test_do_prokka(tmpdir):
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--do-prokka",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert config["prokka"]["do"] is True


def test_do_circlator(tmpdir):
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--do-circlator",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert config["circlator"]["do"] is True


def test_canu_without_circlator_warns(tmpdir):
    """Canu without --do-circlator should succeed but log a warning."""
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "canu",
            "--genome-size",
            "1m",
            "--data-type",
            "pacbio-corr",
        ],
    )
    assert results.exit_code == 0


def test_do_coverage(tmpdir):
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--do-coverage",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert config["sequana_coverage"]["do"] is True


def test_do_correction(tmpdir):
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--do-correction",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert config["canu_correction"]["do"] is True


def test_do_ccs(tmpdir):
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "pacbio-hifi",
            "--pacbio-build-ccs",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert config["ccs"]["do"] is True


def test_do_polypolish(tmpdir):
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--do-polypolish",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert config["polypolish"]["do"] is True


# ---------------------------------------------------------------------------
# BLAST options
# ---------------------------------------------------------------------------


def test_blastdb(tmp_path, tmpdir):
    blast_dir = tmp_path / "blastdb"
    blast_dir.mkdir()
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--blastdb",
            str(blast_dir),
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert config["blast"]["do"] is True


def test_blast_email(tmpdir):
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--blast-email",
            "test@example.com",
            "--blast-remote-db",
            "nt",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert config["blast"]["do"] is True
    assert config["blast"]["remote"] is True
    assert config["blast"]["email"] == "test@example.com"
    assert config["blast"]["remote_db"] == "nt"


# ---------------------------------------------------------------------------
# Genome size
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("genome_size", ["5000000", "100k", "2g"])
def test_genome_size_variants(tmpdir, genome_size):
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            genome_size,
            "--data-type",
            "nano-hq",
        ],
    )
    assert results.exit_code == 0


def test_invalid_genome_size(tmpdir):
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "5x",
            "--data-type",
            "nano-hq",
        ],
    )
    assert results.exit_code != 0


# ---------------------------------------------------------------------------
# BUSCO options
# ---------------------------------------------------------------------------


def test_busco_print_lineages(mocker):
    mocker.patch(
        "sequana_pipelines.lora.main.get_busco_lineages_and_urls",
        return_value={
            "bacteria": "http://example.com/bacteria_odb12.tar.gz",
            "eukaryota": "http://example.com/eukaryota_odb12.tar.gz",
        },
    )
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--busco-print-lineages",
        ],
    )
    assert results.exit_code == 0
    assert "bacteria" in results.output


def test_invalid_busco_lineage(mocker, tmpdir):
    mocker.patch(
        "sequana_pipelines.lora.main.get_busco_lineages_and_urls",
        return_value={"bacteria": "http://example.com/bacteria_odb12.tar.gz"},
    )
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--busco-lineage",
            "totally_invalid_xyz123",
        ],
    )
    assert results.exit_code != 0


def test_busco_download_retries_then_succeeds(mocker, tmpdir):
    """Download fails twice then succeeds on the third attempt."""
    mocker.patch(
        "sequana_pipelines.lora.main.get_busco_lineages_and_urls",
        return_value={"bacteria": "http://example.com/bacteria_odb12.tar.gz"},
    )
    mocker.patch("sequana_pipelines.lora.main.time.sleep")
    call_count = {"n": 0}

    def flaky_download(url, dest):
        call_count["n"] += 1
        if call_count["n"] < 3:
            raise OSError("network error")

    mocker.patch("sequana_pipelines.lora.main.download_and_extract_tar_gz", side_effect=flaky_download)
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--busco-lineage",
            "bacteria",
        ],
    )
    assert results.exit_code == 0
    assert call_count["n"] == 3


def test_busco_download_fails_all_retries(mocker, tmpdir):
    """All three download attempts fail → exit code != 0."""
    mocker.patch(
        "sequana_pipelines.lora.main.get_busco_lineages_and_urls",
        return_value={"bacteria": "http://example.com/bacteria_odb12.tar.gz"},
    )
    mocker.patch("sequana_pipelines.lora.main.time.sleep")
    mocker.patch(
        "sequana_pipelines.lora.main.download_and_extract_tar_gz",
        side_effect=OSError("network error"),
    )
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--busco-lineage",
            "bacteria",
        ],
    )
    assert results.exit_code != 0


def test_busco_lineage_reuses_existing_dir(mocker, tmpdir):
    """When the busco lineage dir already exists locally, skip download."""
    lineage_url = "http://example.com/bacteria_odb12.tar.gz"
    mocker.patch(
        "sequana_pipelines.lora.main.get_busco_lineages_and_urls",
        return_value={"bacteria": lineage_url},
    )
    mock_download = mocker.patch("sequana_pipelines.lora.main.download_and_extract_tar_gz")
    # pre-create the expected busco download dir so it appears already present

    busco_dir = tmpdir.mkdir("busco_downloads").mkdir("bacteria_odb12")
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--busco-lineage",
            "bacteria",
        ],
    )
    assert results.exit_code == 0
    mock_download.assert_not_called()
    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert "busco_downloads/bacteria_odb12" in config["busco"]["lineage"]


# ---------------------------------------------------------------------------
# CheckM options
# ---------------------------------------------------------------------------


def test_checkm_valid_name(tmpdir):
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--checkm-name",
            "Bacillus",
            "--checkm-rank",
            "genus",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((tmpdir / "config.yaml").read_text("utf-8"))
    assert config["checkm"]["do"] is True
    assert config["checkm"]["taxon_name"] == "Bacillus"


def test_checkm_invalid_name(tmpdir):
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--checkm-name",
            "InvalidXYZ123Organism",
            "--checkm-rank",
            "genus",
        ],
    )
    assert results.exit_code != 0


# ---------------------------------------------------------------------------
# ChoiceOrDir validation
# ---------------------------------------------------------------------------


def test_choice_or_dir_invalid_directory(tmpdir):
    """A directory path without dataset.cfg should be rejected by ChoiceOrDir."""
    runner = CliRunner()
    results = runner.invoke(
        main.main,
        [
            "--input-directory",
            str(test_dir),
            "--working-directory",
            str(tmpdir / "output"),
            "--force",
            "--assembler",
            "flye",
            "--genome-size",
            "1m",
            "--data-type",
            "nano-hq",
            "--busco-lineage",
            str(test_dir),  # real dir, no dataset.cfg
        ],
    )
    assert results.exit_code != 0
