"""Konverterer en Quarto-fil (.qmd) til en Jupyter-notebook (.ipynb).

Teksten deles opp i markdown-celler slik at den er lett å redigere i Jupyter:
  - hver overskrift starter en ny celle (med avsnittene under seg)
  - hver `{{< include ... >}}`, hver ```{ojs}-blokk og hver innrammet boks (:::) får sin egen celle
  - YAML-hodet (--- ... ---) blir en rå celle øverst

Bruk:
    python tools/qmd_til_ipynb.py kapitler/01-sentrale-begreper.qmd
    python tools/qmd_til_ipynb.py kapitler/01-sentrale-begreper.qmd --ut kapitler/kap1.ipynb

Krever:  pip install nbformat
"""
import argparse
import re
import sys
from pathlib import Path

import nbformat as nbf


def del_i_blokker(tekst):
    """Del teksten i blokker adskilt av tomme linjer, uten å dele inni kodeblokker eller bokser."""
    blokker, nåværende = [], []
    i_kode = False
    div_dybde = 0
    for linje in tekst.split("\n"):
        strippet = linje.strip()
        if strippet.startswith("```"):
            i_kode = not i_kode
        elif not i_kode:
            if re.match(r"^:{3,}\s*\S", strippet):      # åpner en boks, f.eks. ::: {.callout-tip}
                div_dybde += 1
            elif re.match(r"^:{3,}\s*$", strippet) and div_dybde > 0:
                div_dybde -= 1
        if strippet == "" and not i_kode and div_dybde == 0:
            if nåværende:
                blokker.append("\n".join(nåværende))
                nåværende = []
        else:
            nåværende.append(linje)
    if nåværende:
        blokker.append("\n".join(nåværende))
    return blokker


def type_blokk(blokk):
    første = blokk.lstrip().split("\n", 1)[0]
    if re.match(r"^#{1,6}\s", første):
        return "overskrift"
    if første.startswith("{{< include"):
        return "include"
    if første.startswith("```{ojs}"):
        return "ojs"
    if re.match(r"^:{3,}\s*\S", første):
        return "boks"
    return "tekst"


def lag_celler(brødtekst):
    celler, nåværende = [], []
    for blokk in del_i_blokker(brødtekst):
        t = type_blokk(blokk)
        if t in ("overskrift", "include", "ojs", "boks"):
            if nåværende:
                celler.append("\n\n".join(nåværende))
                nåværende = []
            if t == "overskrift":
                nåværende.append(blokk)
            else:
                celler.append(blokk)
        else:
            nåværende.append(blokk)
    if nåværende:
        celler.append("\n\n".join(nåværende))
    return celler


def konverter(inn: Path, ut: Path):
    tekst = inn.read_text(encoding="utf-8")
    hode, brødtekst = "", tekst
    m = re.match(r"^---\n(.*?)\n---\n", tekst, flags=re.S)
    if m:
        hode, brødtekst = m.group(1), tekst[m.end():]
    if "jupyter:" not in hode:
        hode = (hode + "\n" if hode else "") + "jupyter: python3"
    nb = nbf.v4.new_notebook()
    nb.cells = [nbf.v4.new_raw_cell("---\n" + hode + "\n---")]
    nb.cells += [nbf.v4.new_markdown_cell(c) for c in lag_celler(brødtekst.strip("\n"))]
    nb.metadata = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                   "language_info": {"name": "python"}}
    nbf.write(nb, ut)
    return len(nb.cells)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("qmd", type=Path)
    ap.add_argument("--ut", type=Path, help="filnavn på notebooken (standard: samme navn med .ipynb)")
    args = ap.parse_args()
    ut = args.ut or args.qmd.with_suffix(".ipynb")
    if ut.exists():
        sys.exit(f"{ut} finnes allerede. Slett den først, eller bruk --ut med et nytt navn.")
    n = konverter(args.qmd, ut)
    print(f"Skrev {ut} med {n} celler.")
