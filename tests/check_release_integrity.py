"""Check an existing release without changing data or regenerating figures.

Run from the repository root with Python 3; no third-party packages are needed.
"""
from pathlib import Path
import base64
import csv
import hashlib
import json
import re


ROOT = Path(__file__).resolve().parents[1]
COMPILED = ROOT / "results/compiled_figures"
failures = []


def check(condition, message):
    if not condition:
        failures.append(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


with (ROOT / "docs/manifests/file_manifest.csv").open() as stream:
    manifest = list(csv.DictReader(stream))
for entry in manifest:
    path = ROOT / entry["repository_path"]
    check(path.is_file(), f"Missing manifest input: {entry['repository_path']}")
    if path.is_file():
        check(sha256(path) == entry["sha256"], f"Changed manifest input: {entry['repository_path']}")
        check(path.stat().st_size == int(entry["size_bytes"]), f"Changed input size: {entry['repository_path']}")

figure_names = [f"Figure{n}" for n in range(1, 6)] + [f"FigureS{n}" for n in range(1, 20)]
expected = {f"{name}.{ext}" for name in figure_names for ext in ("png", "pdf")}
expected.update({
    "FigureS11_source.csv", "FigureS11_paired_wilcoxon_statistics.csv",
    "TableS1.xlsx", "TableS2.xlsx", "TableS3.xlsx",
    "compiled_figures.nb.html", "FigureS15_paired_wilcoxon_statistics.csv",
    "FigureS17_source.csv", "FigureS18_source.csv",
    "FigureS1_environment_by_gtdb_class_source.csv",
    "FigureS1_sampling_location_clusters_source.csv",
    "additional_subpanels/05_1_05_2_uptake_rate_paired_tests.csv",
})
for number, scale in ((1, "linear"), (2, "log")):
    expected.update(
        f"additional_subpanels/05_{number}_uptake_rate_by_trait_group_and_added_glucose_{scale}_y.{ext}"
        for ext in ("png", "pdf")
    )
actual = {str(path.relative_to(COMPILED)) for path in COMPILED.rglob("*") if path.is_file()}
check(actual == expected, f"Compiled inventory mismatch; missing={sorted(expected - actual)}, unexpected={sorted(actual - expected)}")
check(all((COMPILED / path).stat().st_size > 0 for path in actual), "Empty compiled output")
extras = {str(path.relative_to(ROOT / "results")) for path in (ROOT / "results").iterdir()} - {"README.md", "compiled_figures"}
check(not extras, f"Unexpected results outside the production output directory: {sorted(extras)}")
temporary = [str(path.relative_to(ROOT)) for path in ROOT.rglob("*")
             if ".git" not in path.parts and path.is_file()
             and (path.name in {".DS_Store", "Rplots.pdf", ".Rhistory", ".RData"}
                  or path.suffix in {".bak", ".tmp"} or path.name.endswith("~"))]
check(not temporary, f"Temporary files in release: {temporary}")

notebook_path = COMPILED / "compiled_figures.nb.html"
if notebook_path.is_file():
    notebook = notebook_path.read_text()
    rmd = (ROOT / "R/compiled_figures/compiled_figures.Rmd").read_text()
    source = re.search(r'<[^>]+id="rmd-source-code"[^>]*>(.*?)</[^>]+>', notebook, re.S)
    check(source is not None, "Missing downloadable Rmd in notebook")
    if source:
        check(base64.b64decode(source[1]).decode() == rmd, "Notebook downloadable Rmd is stale")
    canonical = re.sub(r"\s+", " ", rmd)
    blocks = re.findall(r'<!-- rnb-source-begin (\S+) -->', notebook)
    check(bool(blocks), "Notebook has no displayed R source")
    for encoded in blocks:
        value = json.loads(base64.b64decode(encoded))["data"]
        code = value if isinstance(value, str) else "\n".join(value)
        # Notebook output can split a chunk into several displayed blocks.
        check(re.sub(r"\s+", " ", code).strip() in canonical,
              f"Stale notebook code beginning: {code[:90]}")
    for name in figure_names:
        section_id = "figure-" + name[len("Figure"):].lower()
        if name in {"FigureS12", "FigureS13", "FigureS14"}:
            section_id = "figures-s12-s14-kinetic-model-structure-comparison"
        check(f'id="{section_id}"' in notebook or f'id="{section_id}-' in notebook,
              f"Notebook missing section for {name}")
    # The cached output table records sizes; check each row that it contains.
    cached_sizes = re.findall(r'"1":"([^"\n]+)","2":"TRUE","3":"(\d+)"}', notebook)
    for recorded_path, size in cached_sizes:
        filename = Path(recorded_path).name
        matches = list(COMPILED.rglob(filename))
        check(len(matches) == 1 and matches[0].stat().st_size == int(size),
              f"Stale notebook output size: {filename}")

if failures:
    raise SystemExit("FAIL:\n- " + "\n- ".join(failures))
print(f"PASS: {len(manifest)} input checksums, {len(expected)} compiled files, current notebook source/output sizes, and no temporary/exploratory results.")
