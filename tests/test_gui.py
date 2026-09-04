import io
from pathlib import Path

from sangerflow.gui import GuiState, _config_from_payload, create_app


def test_gui_preview_and_synthetic_demo(tmp_path: Path):
    app = create_app(tmp_path / "outputs")
    state: GuiState = app.extensions["sangerflow_state"]
    client = app.test_client()
    prefix = f"/{state.token}"
    headers = {"X-SangerFlow-Token": state.token}

    page = client.get(prefix + "/")
    assert page.status_code == 200
    assert b"CMHS SangerFlow Pipeline" in page.data
    assert client.get(prefix + "/api/status").status_code == 403

    preview = client.post(
        prefix + "/api/preview",
        headers=headers,
        data={
            "reads": [
                (io.BytesIO(b"forward"), "patient-01_F.ab1"),
                (io.BytesIO(b"reverse"), "patient-01_R.abi"),
            ],
            "reference": (io.BytesIO(b">demo\nAACCGGTT\n"), "reference.fasta"),
        },
        content_type="multipart/form-data",
    )
    assert preview.status_code == 200
    rows = preview.get_json()["rows"]
    assert [(row["sample"], row["direction"]) for row in rows] == [
        ("patient-01", "forward"),
        ("patient-01", "reverse"),
    ]

    demo = client.post(prefix + "/api/demo", headers=headers)
    assert demo.status_code == 200
    result = demo.get_json()
    assert result["status"] == "COMPLETE"
    assert "No patient data" in result["message"]
    report = client.get(result["report"])
    assert report.status_code == 200
    assert b"synthetic demonstration" in report.data
    state.temporary.cleanup()


def test_gui_advanced_configuration():
    config = _config_from_payload(
        {
            "preset": "advanced",
            "settings": {
                "mask_below": "25",
                "min_overlap": "40",
                "min_reference_identity": "0.9",
            },
        }
    )
    assert config.mask_below == 25
    assert config.min_overlap == 40
    assert config.min_reference_identity == 0.9
