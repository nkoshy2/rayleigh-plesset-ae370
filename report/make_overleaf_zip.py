"""Build report-overleaf.zip (report sources + figure PDFs) for upload to Overleaf.

Run from anywhere:  python report/make_overleaf_zip.py
The zip has main.tex at its top level and the figures in figures/, which
main.tex finds through \\graphicspath.
"""

import pathlib
import zipfile

REPORT = pathlib.Path(__file__).resolve().parent
ROOT = REPORT.parent
OUT = ROOT / "report-overleaf.zip"

with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
    for path in [REPORT / "main.tex", REPORT / "references.bib",
                 *sorted((REPORT / "sections").glob("*.tex"))]:
        z.write(path, path.relative_to(REPORT).as_posix())
    for fig in sorted((ROOT / "figures").glob("*.pdf")):
        z.write(fig, f"figures/{fig.name}")
    names = z.namelist()

print(f"wrote {OUT} with {len(names)} files:")
for n in names:
    print("  " + n)
