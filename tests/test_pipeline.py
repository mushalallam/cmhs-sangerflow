import json
from pathlib import Path

import pytest

from sangerflow.models import ReadData, SampleInput
from sangerflow.pipeline import RunConfig, _auto_orient, run_pipeline


def prepared_read(path: Path) -> ReadData:
    return ReadData(
        name=path.stem,
        path=path.resolve(),
        sequence="AACCGG",
        qualities=[35] * 6,
        trim_start=0,
        trim_end=6,
        trimmed_sequence="AACCGG",
        trimmed_qualities=[35] * 6,
    )


def test_pipeline_writes_auditable_outputs(tmp_path: Path, monkeypatch):
    forward = tmp_path / "sample_F.ab1"
    reverse = tmp_path / "sample_R.ab1"
    forward.write_bytes(b"forward")
    reverse.write_bytes(b"reverse")
    reference = tmp_path / "reference.fasta"
    reference.write_text(">ref\nTTTTAACCGGTTTT\n", encoding="utf-8")

    def fake_prepare(path, config, *, reverse=False):
        return prepared_read(path)

    monkeypatch.setattr("sangerflow.pipeline._prepare_read", fake_prepare)
    output = tmp_path / "results"
    progress = []
    metadata = run_pipeline(
        [SampleInput("sample", forward, reverse)],
        reference,
        output,
        RunConfig(min_read_length=1, min_overlap=1),
        progress_callback=lambda *values: progress.append(values),
    )

    assert metadata["passed_samples"] == 1
    assert metadata["failed_samples"] == 0
    assert (output / "reports/report.html").is_file()
    assert (output / "reports/variants.vcf").is_file()
    assert (output / "consensus/all_consensus.fasta").is_file()
    assert (output / "traces/sample.forward.quality.svg").is_file()
    assert (output / "traces/batch_qc.svg").is_file()
    loaded = json.loads((output / "reports/run_metadata.json").read_text())
    assert loaded["input_sha256"][str(forward.resolve())]
    assert progress == [(0, 1, "sample", "RUNNING"), (1, 1, "sample", "PASS")]


def test_pipeline_records_sample_failure(tmp_path: Path, monkeypatch):
    bad = tmp_path / "bad_F.ab1"
    bad.write_bytes(b"not an ABI")
    reference = tmp_path / "reference.fasta"
    reference.write_text(">ref\nAACCGG\n", encoding="utf-8")

    def fail_prepare(path, config, *, reverse=False):
        raise ValueError("invalid chromatogram")

    monkeypatch.setattr("sangerflow.pipeline._prepare_read", fail_prepare)
    output = tmp_path / "results"
    metadata = run_pipeline(
        [SampleInput("bad", bad, None)], reference, output, RunConfig(min_read_length=1)
    )
    assert metadata["failed_samples"] == 1
    assert "invalid chromatogram" in (output / "reports/samples.tsv").read_text()


def test_pipeline_refuses_nonempty_output(tmp_path: Path):
    reference = tmp_path / "reference.fasta"
    reference.write_text(">ref\nAACCGG\n", encoding="utf-8")
    output = tmp_path / "results"
    output.mkdir()
    (output / "keep.txt").write_text("do not overwrite", encoding="utf-8")
    with pytest.raises(ValueError, match="not empty"):
        run_pipeline([], reference, output, RunConfig())
    assert (output / "keep.txt").read_text() == "do not overwrite"


def test_pipeline_rejects_normalized_name_collision(tmp_path: Path):
    reference = tmp_path / "reference.fasta"
    reference.write_text(">ref\nAACCGG\n", encoding="utf-8")
    samples = [SampleInput("a b"), SampleInput("a_b")]
    with pytest.raises(ValueError, match="not unique"):
        run_pipeline(samples, reference, tmp_path / "results", RunConfig())


def test_pipeline_rejects_low_identity_reference(tmp_path: Path, monkeypatch):
    forward = tmp_path / "sample_F.ab1"
    forward.write_bytes(b"forward")
    reference = tmp_path / "reference.fasta"
    reference.write_text(">ref\nTTTTTT\n", encoding="utf-8")
    monkeypatch.setattr(
        "sangerflow.pipeline._prepare_read", lambda path, config: prepared_read(path)
    )
    output = tmp_path / "results"
    metadata = run_pipeline(
        [SampleInput("sample", forward, None)], reference, output, RunConfig(min_read_length=1)
    )
    assert metadata["failed_samples"] == 1
    assert "verify sample" in (output / "reports/samples.tsv").read_text()


def test_auto_orientation_corrects_read_against_reference(tmp_path: Path):
    value = ReadData(
        "read",
        tmp_path / "read.ab1",
        "AATCGG",
        [35] * 6,
        trim_start=0,
        trim_end=6,
        trimmed_sequence="AATCGG",
        trimmed_qualities=[35] * 6,
    )
    _auto_orient(value, "CCGATT", RunConfig(min_read_length=1))
    assert value.trimmed_sequence == "CCGATT"
    assert value.orientation_corrected is True
