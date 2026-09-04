"""Local-only browser interface for CMHS SangerFlow Pipeline."""

# ruff: noqa: E501 -- the embedded HTML/CSS/JavaScript is kept readable in its native form.

from __future__ import annotations

import os
import platform
import secrets
import shutil
import subprocess
import tempfile
import threading
import webbrowser
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any
from uuid import uuid4

from flask import Flask, abort, jsonify, request, send_from_directory
from werkzeug.serving import make_server

from .demo import create_demo_results
from .models import SampleInput
from .pairing import infer_sample_and_direction
from .pipeline import RunConfig, _load_reference, run_pipeline


@dataclass(slots=True)
class GuiState:
    token: str
    temporary: tempfile.TemporaryDirectory[str]
    output_root: Path
    reads: dict[str, Path] = field(default_factory=dict)
    reference: Path | None = None
    output_dir: Path | None = None
    job: dict[str, Any] = field(
        default_factory=lambda: {
            "status": "IDLE",
            "completed": 0,
            "total": 0,
            "sample": "",
            "message": "",
        }
    )
    lock: threading.Lock = field(default_factory=threading.Lock)


def _display_name(filename: str | None) -> str:
    return Path((filename or "read.ab1").replace("\\", "/")).name


def _new_output(root: Path, label: str = "Run") -> Path:
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    candidate = root / f"{label}_{timestamp}"
    if candidate.exists():
        candidate = root / f"{label}_{timestamp}_{uuid4().hex[:6]}"
    return candidate


def _open_directory(path: Path) -> None:
    system = platform.system()
    if system == "Windows":
        os.startfile(path)  # type: ignore[attr-defined]
    elif system == "Darwin":
        subprocess.Popen(["open", str(path)])
    else:
        subprocess.Popen(["xdg-open", str(path)])


def _config_from_payload(payload: dict[str, Any]) -> RunConfig:
    preset = payload.get("preset", "standard")
    settings = payload.get("settings") or {}
    values: dict[str, Any] = {
        "call_mixed_peaks": preset == "mixed",
    }
    integer_fields = {
        "mask_below",
        "min_read_length",
        "quality_delta",
        "min_variant_quality",
        "min_overlap",
        "mixed_peak_min_signal",
    }
    float_fields = {
        "error_cutoff",
        "min_overlap_identity",
        "min_reference_identity",
        "min_reference_coverage",
        "mixed_peak_ratio",
        "orientation_delta",
    }
    if preset == "advanced":
        for name in integer_fields:
            if name in settings:
                values[name] = int(settings[name])
        for name in float_fields:
            if name in settings:
                values[name] = float(settings[name])
    config = RunConfig(**values)
    config.validate()
    return config


def _samples_from_payload(state: GuiState, payload: dict[str, Any]) -> list[SampleInput]:
    grouped: dict[str, dict[str, Path]] = {}
    for row in payload.get("rows") or []:
        file_id = str(row.get("id", ""))
        direction = str(row.get("direction", ""))
        sample = str(row.get("sample", "")).strip()
        if direction == "ignore":
            continue
        if file_id not in state.reads:
            raise ValueError("The uploaded read list changed; preview the files again")
        if not sample:
            raise ValueError("Every included read needs a sample name")
        if direction not in {"forward", "reverse"}:
            raise ValueError(f"Choose forward, reverse, or ignore for {sample}")
        slots = grouped.setdefault(sample, {})
        if direction in slots:
            raise ValueError(f"Sample {sample!r} has more than one {direction} read")
        slots[direction] = state.reads[file_id]
    if not grouped:
        raise ValueError("No reads were selected for analysis")
    return [
        SampleInput(sample, reads.get("forward"), reads.get("reverse"))
        for sample, reads in sorted(grouped.items())
    ]


def create_app(output_root: Path | None = None) -> Flask:
    """Create the local application; exposed separately for automated tests."""
    temporary = tempfile.TemporaryDirectory(prefix="cmhs-sangerflow-gui-")
    state = GuiState(
        token=secrets.token_urlsafe(24),
        temporary=temporary,
        output_root=Path(output_root or Path.home() / "SangerFlow Results").resolve(),
    )
    app = Flask(__name__)
    app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024 * 1024
    app.extensions["sangerflow_state"] = state
    prefix = f"/{state.token}"

    @app.get(prefix + "/")
    def index():
        return GUI_HTML.replace("__TOKEN__", state.token)

    @app.before_request
    def protect_api():
        if (
            request.path.startswith(prefix + "/api/")
            and request.headers.get("X-SangerFlow-Token") != state.token
        ):
            abort(403)

    @app.post(prefix + "/api/preview")
    def preview():
        with state.lock:
            if state.job.get("status") == "RUNNING":
                return jsonify(error="Wait for the current analysis to finish"), 409
        files = request.files.getlist("reads")
        reference = request.files.get("reference")
        if not files:
            return jsonify(error="Select at least one ABI/AB1 chromatogram"), 400
        if reference is None or not reference.filename:
            return jsonify(error="Select a reference FASTA file"), 400

        upload_root = Path(state.temporary.name)
        reads_root = upload_root / "reads"
        if reads_root.exists():
            shutil.rmtree(reads_root)
        reads_root.mkdir()
        state.reads.clear()
        rows = []
        for uploaded in files:
            name = _display_name(uploaded.filename)
            if Path(name).suffix.lower() not in {".ab1", ".abi"}:
                continue
            file_id = uuid4().hex
            destination = reads_root / name
            if destination.exists():
                return jsonify(error=f"More than one selected file is named {name!r}"), 400
            uploaded.save(destination)
            state.reads[file_id] = destination
            try:
                sample, direction = infer_sample_and_direction(Path(name))
            except ValueError:
                sample, direction = Path(name).stem, ""
            rows.append(
                {"id": file_id, "filename": name, "sample": sample, "direction": direction}
            )
        if not rows:
            return jsonify(error="No files ending in .ab1 or .abi were selected"), 400

        reference_name = _display_name(reference.filename)
        reference_path = upload_root / f"reference-{uuid4().hex}.fasta"
        reference.save(reference_path)
        record_name, sequence = _load_reference(reference_path)
        state.reference = reference_path
        with state.lock:
            state.job = {
                "status": "READY",
                "completed": 0,
                "total": len(rows),
                "sample": "",
                "message": "Review the detected pairs, then start the analysis.",
            }
        return jsonify(
            rows=rows,
            reference={"filename": reference_name, "name": record_name, "length": len(sequence)},
        )

    @app.post(prefix + "/api/run")
    def start_run():
        payload = request.get_json(force=True)
        with state.lock:
            if state.job.get("status") == "RUNNING":
                return jsonify(error="An analysis is already running"), 409
        if state.reference is None:
            return jsonify(error="Preview the input files first"), 400
        try:
            samples = _samples_from_payload(state, payload)
            config = _config_from_payload(payload)
            if payload.get("preset", "standard") != "single":
                incomplete = [item.sample for item in samples if not item.forward or not item.reverse]
                if incomplete:
                    raise ValueError(
                        "Choose the Single reads preset or provide both directions for: "
                        + ", ".join(incomplete)
                    )
        except (TypeError, ValueError) as exc:
            return jsonify(error=str(exc)), 400

        output = _new_output(state.output_root)
        state.output_dir = output
        with state.lock:
            state.job = {
                "status": "RUNNING",
                "completed": 0,
                "total": len(samples),
                "sample": "",
                "message": "Preparing analysis…",
                "output": str(output),
            }

        def progress(completed: int, total: int, sample: str, status: str) -> None:
            with state.lock:
                state.job.update(
                    completed=completed,
                    total=total,
                    sample=sample,
                    message=(
                        f"Analyzing {sample}…"
                        if status == "RUNNING"
                        else f"Finished {sample}: {status}"
                    ),
                )

        def worker() -> None:
            try:
                metadata = run_pipeline(
                    samples,
                    state.reference,
                    output,
                    config,
                    progress_callback=progress,
                )
                with state.lock:
                    state.job.update(
                        status="COMPLETE",
                        completed=len(samples),
                        message=(
                            f"Complete: {metadata['passed_samples']} passed, "
                            f"{metadata['failed_samples']} flagged."
                        ),
                        report=f"{prefix}/results/reports/report.html",
                    )
            except Exception as exc:
                with state.lock:
                    state.job.update(status="ERROR", message=str(exc))

        threading.Thread(target=worker, name="sangerflow-analysis", daemon=True).start()
        return jsonify(status="RUNNING", output=str(output))

    @app.get(prefix + "/api/status")
    def status():
        with state.lock:
            return jsonify(dict(state.job))

    @app.post(prefix + "/api/demo")
    def demo():
        with state.lock:
            if state.job.get("status") == "RUNNING":
                return jsonify(error="Wait for the current analysis to finish"), 409
        output = _new_output(state.output_root, "Synthetic_Demo")
        report = create_demo_results(output)
        state.output_dir = output
        with state.lock:
            state.job = {
                "status": "COMPLETE",
                "completed": 1,
                "total": 1,
                "sample": "SYNTHETIC_SAMPLE",
                "message": "Synthetic demonstration created. No patient data were used.",
                "output": str(output),
                "report": f"{prefix}/results/{report.relative_to(output).as_posix()}",
            }
        return jsonify(dict(state.job))

    @app.post(prefix + "/api/open-output")
    def open_output():
        if state.output_dir is None or not state.output_dir.is_dir():
            return jsonify(error="No result directory is available yet"), 404
        _open_directory(state.output_dir)
        return jsonify(status="OPENED")

    @app.get(prefix + "/results/<path:filename>")
    def results(filename: str):
        if state.output_dir is None:
            abort(404)
        return send_from_directory(state.output_dir, filename)

    return app


def launch_gui(*, port: int = 0, open_browser: bool = True) -> int:
    """Start a loopback-only server and optionally open the user's browser."""
    app = create_app()
    state: GuiState = app.extensions["sangerflow_state"]
    server = make_server("127.0.0.1", port, app, threaded=True)
    url = f"http://127.0.0.1:{server.server_port}/{state.token}/"
    print(f"CMHS SangerFlow is running locally at {url}")
    print("Close this window or press Ctrl+C to stop the interface.")
    if open_browser:
        threading.Timer(0.5, webbrowser.open, args=(url,)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.shutdown()
        state.temporary.cleanup()
    return 0


GUI_HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>CMHS SangerFlow Pipeline</title>
<style>
:root{--navy:#17324d;--blue:#235789;--sky:#eaf5ff;--green:#16883e;--red:#b42318;
--line:#d8e1e8;--muted:#617283}*{box-sizing:border-box}body{margin:0;background:#f5f8fa;
color:var(--navy);font:15px system-ui,-apple-system,"Segoe UI",sans-serif}header{background:var(--navy);
color:white;padding:22px max(24px,calc((100% - 1100px)/2))}header h1{margin:0;font-size:25px}
header p{margin:5px 0 0;color:#c9d9e7}.wrap{max-width:1100px;margin:24px auto;padding:0 20px}
.notice{background:#fff7db;border-left:4px solid #d99a00;padding:12px 15px;margin-bottom:18px}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:18px}.card{background:white;border:1px solid var(--line);
border-radius:12px;padding:20px;box-shadow:0 2px 8px #17324d0b}.wide{grid-column:1/-1}h2{font-size:18px;
margin:0 0 14px}.step{display:inline-grid;place-items:center;width:26px;height:26px;border-radius:50%;
background:var(--blue);color:white;margin-right:8px}label{font-weight:650;display:block;margin:12px 0 6px}
input[type=file],select,input[type=text],input[type=number]{width:100%;padding:9px;border:1px solid #b9c7d2;
border-radius:7px;background:white}button,.button{border:0;border-radius:8px;background:var(--blue);color:white;
padding:10px 16px;font-weight:700;cursor:pointer;text-decoration:none;display:inline-block}button.secondary{background:#657789}
button.green{background:var(--green)}button:disabled{opacity:.45;cursor:not-allowed}.actions{display:flex;gap:10px;
flex-wrap:wrap;margin-top:16px}.small{color:var(--muted);font-size:13px}.error{color:var(--red);font-weight:650}
table{width:100%;border-collapse:collapse;margin-top:12px}th,td{text-align:left;padding:8px;border-bottom:1px solid var(--line)}
th{font-size:12px;text-transform:uppercase;color:var(--muted)}td input,td select{margin:0}.filename{word-break:break-all}
.progress{height:15px;background:#e4e9ed;border-radius:9px;overflow:hidden;margin:12px 0}.bar{height:100%;width:0;
background:var(--green);transition:width .25s}.status{padding:12px;background:var(--sky);border-radius:8px;min-height:44px}
details{margin-top:14px}details .advanced{display:grid;grid-template-columns:repeat(3,1fr);gap:10px}
footer{color:var(--muted);padding:24px 0 40px}@media(max-width:760px){.grid{grid-template-columns:1fr}
details .advanced{grid-template-columns:1fr}.wide{grid-column:auto}table{font-size:12px}th,td{padding:5px}}
</style></head><body>
<header><h1>CMHS SangerFlow Pipeline</h1><p>Human Genomics Solutions · local Sanger analysis</p></header>
<main class="wrap"><div class="notice"><strong>Research use only.</strong> Review chromatograms and confirm calls independently.
Files remain on this computer and are not sent over the internet.</div><div class="grid">
<section class="card"><h2><span class="step">1</span>Select files</h2>
<label for="reads">ABI/AB1 chromatograms</label><input id="reads" type="file" accept=".ab1,.abi" multiple>
<label for="reference">Reference FASTA</label><input id="reference" type="file" accept=".fa,.fasta,.fna">
<div class="actions"><button id="preview">Preview pairing</button><button class="secondary" id="demo">Run synthetic demo</button></div>
<p class="small">Select all forward and reverse reads together. The demonstration contains no patient data.</p></section>
<section class="card"><h2><span class="step">2</span>Analysis settings</h2>
<label for="preset">Preset</label><select id="preset"><option value="standard">Standard paired reads</option>
<option value="single">Single reads allowed</option><option value="mixed">Mixed-peak research mode</option>
<option value="advanced">Advanced/custom thresholds</option></select>
<p class="small" id="preset-help">Conservative defaults for forward/reverse sequencing.</p>
<details id="advanced"><summary>Advanced thresholds</summary><div class="advanced">
<label>Error cutoff<input data-setting="error_cutoff" type="number" step="0.01" value="0.05"></label>
<label>Mask below Q<input data-setting="mask_below" type="number" value="20"></label>
<label>Minimum read length<input data-setting="min_read_length" type="number" value="80"></label>
<label>Quality delta<input data-setting="quality_delta" type="number" value="8"></label>
<label>Minimum variant Q<input data-setting="min_variant_quality" type="number" value="20"></label>
<label>Minimum overlap<input data-setting="min_overlap" type="number" value="30"></label>
<label>Overlap identity<input data-setting="min_overlap_identity" type="number" step="0.01" value="0.80"></label>
<label>Reference identity<input data-setting="min_reference_identity" type="number" step="0.01" value="0.80"></label>
<label>Reference coverage<input data-setting="min_reference_coverage" type="number" step="0.01" value="0.80"></label>
</div></details></section>
<section class="card wide" id="pair-card" hidden><h2><span class="step">3</span>Review sample pairing</h2>
<p id="reference-info" class="small"></p><div style="overflow:auto"><table><thead><tr><th>File</th><th>Sample</th>
<th>Direction</th></tr></thead><tbody id="pairs"></tbody></table></div>
<div class="actions"><button class="green" id="run">Start analysis</button></div></section>
<section class="card wide"><h2><span class="step">4</span>Progress and results</h2><div class="status" id="status">Ready.</div>
<div class="progress"><div class="bar" id="bar"></div></div><p class="small" id="output"></p>
<div class="actions"><a class="button" id="report" target="_blank" hidden>Open report</a>
<button class="secondary" id="open-output" hidden>Open results folder</button></div></section></div>
<footer>CMHS SangerFlow runs at 127.0.0.1 only. Closing the launcher stops this interface.</footer></main>
<script>
const TOKEN="__TOKEN__", API=`/${TOKEN}/api`;let previewRows=[];
const el=id=>document.getElementById(id);const headers={"X-SangerFlow-Token":TOKEN};
function message(text,error=false){el("status").textContent=text;el("status").className=error?"status error":"status"}
async function jsonResponse(response){const data=await response.json();if(!response.ok)throw new Error(data.error||"Request failed");return data}
el("preset").onchange=()=>{const p=el("preset").value;el("advanced").open=p==="advanced";
el("preset-help").textContent={standard:"Conservative defaults for forward/reverse sequencing.",single:"Allows samples with only one sequencing direction.",mixed:"Experimental IUPAC calling for strong secondary peaks.",advanced:"Edit validated quality and alignment thresholds below."}[p]};
el("preview").onclick=async()=>{try{const reads=el("reads").files,reference=el("reference").files[0];
if(!reads.length||!reference)throw new Error("Select chromatograms and a reference FASTA first.");message("Uploading files locally and checking the reference…");
const form=new FormData();for(const file of reads)form.append("reads",file);form.append("reference",reference);
const data=await jsonResponse(await fetch(`${API}/preview`,{method:"POST",headers,body:form}));previewRows=data.rows;
el("pairs").replaceChildren(...previewRows.map((row,index)=>{const tr=document.createElement("tr");
const file=document.createElement("td");file.className="filename";file.textContent=row.filename;
const sample=document.createElement("td");const input=document.createElement("input");input.type="text";input.value=row.sample;
input.oninput=()=>previewRows[index].sample=input.value;sample.append(input);const direction=document.createElement("td");
const select=document.createElement("select");for(const [value,label] of [["","Choose…"],["forward","Forward"],["reverse","Reverse"],["ignore","Ignore"]]){
const option=document.createElement("option");option.value=value;option.textContent=label;option.selected=row.direction===value;select.append(option)}
select.onchange=()=>previewRows[index].direction=select.value;direction.append(select);tr.append(file,sample,direction);return tr}));
el("reference-info").textContent=`Reference: ${data.reference.name} · ${data.reference.length} bp · ${data.reference.filename}`;
el("pair-card").hidden=false;message(`Found ${previewRows.length} reads. Review the table before starting.`)}catch(error){message(error.message,true)}};
el("run").onclick=async()=>{try{const settings={};document.querySelectorAll("[data-setting]").forEach(x=>settings[x.dataset.setting]=x.value);
const data=await jsonResponse(await fetch(`${API}/run`,{method:"POST",headers:{...headers,"Content-Type":"application/json"},
body:JSON.stringify({rows:previewRows,preset:el("preset").value,settings})}));message("Analysis started…");el("run").disabled=true;
el("output").textContent=data.output;poll()}catch(error){message(error.message,true)}};
async function poll(){try{const data=await jsonResponse(await fetch(`${API}/status`,{headers}));const pct=data.total?100*data.completed/data.total:0;
el("bar").style.width=`${pct}%`;message(data.message||data.status,data.status==="ERROR");if(data.output)el("output").textContent=`Results: ${data.output}`;
if(data.status==="RUNNING")setTimeout(poll,700);else{el("run").disabled=false;if(data.report){el("report").href=data.report;el("report").hidden=false;
el("open-output").hidden=false}}}catch(error){message(error.message,true)}}
el("demo").onclick=async()=>{try{message("Creating synthetic demonstration…");const data=await jsonResponse(await fetch(`${API}/demo`,{method:"POST",headers}));
message(data.message);el("output").textContent=`Results: ${data.output}`;el("bar").style.width="100%";el("report").href=data.report;
el("report").hidden=false;el("open-output").hidden=false;window.open(data.report,"_blank")}catch(error){message(error.message,true)}};
el("open-output").onclick=async()=>{try{await jsonResponse(await fetch(`${API}/open-output`,{method:"POST",headers}))}catch(error){message(error.message,true)}};
</script></body></html>"""
