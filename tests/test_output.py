from pathlib import Path

from sangerflow.models import ReadData, Variant
from sangerflow.output import safe_name, write_fastq, write_trace_svg, write_vcf


def test_safe_name_removes_path_characters():
    assert safe_name("patient 1/../x") == "patient_1_.._x"


def test_fastq_encoding(tmp_path: Path):
    output = tmp_path / "read.fastq"
    write_fastq(output, "read", "AC", [0, 40])
    assert output.read_text().splitlines()[-1] == "!I"


def test_vcf_has_required_columns(tmp_path: Path):
    output = tmp_path / "calls.vcf"
    variant = Variant("sample", 2, "A", "G", "SNV", 30)
    write_vcf(output, "ref", 10, [variant])
    text = output.read_text()
    assert text.startswith("##fileformat=VCFv4.3")
    assert "ref\t2\t.\tA\tG\t30.0\tPASS" in text


def test_vcf_sorts_positions_and_sanitizes_sample(tmp_path: Path):
    output = tmp_path / "calls.vcf"
    variants = [
        Variant("sample one", 5, "A", "G", "SNV", 30),
        Variant("sample one", 2, "C", "T", "SNV", 30),
    ]
    write_vcf(output, "ref", 10, variants)
    records = [line for line in output.read_text().splitlines() if not line.startswith("#")]
    assert records[0].split("\t")[1] == "2"
    assert "SAMPLE=sample_one" in records[0]


def test_trace_svg(tmp_path: Path):
    read = ReadData(
        "trace",
        tmp_path / "trace.ab1",
        "AC",
        [30, 30],
        trace_channels={base: list(range(20)) for base in "ACGT"},
        peak_locations=[5, 15],
        trim_start=0,
        trim_end=2,
        trimmed_sequence="AC",
        trimmed_qualities=[30, 30],
    )
    output = tmp_path / "trace.svg"
    assert write_trace_svg(output, read)
    assert "<svg" in output.read_text()
    assert "green area retained" in output.read_text()
