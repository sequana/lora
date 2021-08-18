import os

from sequana_pipelines.lora import create_report
from . import test_dir


def test_lora_create_report(tmpdir):
    """Test LORA report creation"""
    lora_dir = os.path.join(test_dir, 'resources')
    summary = tmpdir.join('summary.html')
    samples = ['toto']
    create_report(summary, samples, lora_dir=lora_dir)
    assert os.path.exists(summary)
