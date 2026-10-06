# Biokjemi
En konseptuell tilnærming til biokjemi – skrevet som en [Quarto](https://quarto.org)-bok.

## Skriv og bygg lokalt

```bash
pip install -r requirements.txt
quarto preview        # live forhåndsvisning med auto-oppdatering
quarto render         # bygger hele boka til _book/
```

## Struktur

| Fil | Innhold |
|---|---|
| `_quarto.yml` | Bokoppsett, kapitler og format |
| `custom.scss`, `custom-dark.scss` | Fargepalett og typografi |
| `styles.css` | Egne bokser (`.nokkelpunkt`, `.simulering`, `.laeringsmal`) |
| `kapitler/` | Kapitler (`.ipynb` eller `.qmd`) |
| `_interaktiv/` | Kode til interaktive figurer (Observable JS), hentes inn i kapitlene |
| `_figurer/` | Statiske SVG-illustrasjoner som hentes inn i kapitlene |
| `tools/` | Hjelpeskript (se under) |

Nye kapitler legges til i `chapters:` i `_quarto.yml`.

## Publisering

Hver push til `main` bygger og publiserer boka via GitHub Actions (`.github/workflows/publish.yml`).
Slå på én gang: **Settings → Pages → Source: GitHub Actions**.

## Interaktive figurer

Koden til de interaktive figurene (Observable JS) ligger i `_interaktiv/` og hentes inn i kapitlene med
`{{< include ../_interaktiv/celle.qmd >}}`. Filer og mapper som starter med understrek bygges ikke som egne sider.
Under bygging kommer det advarsler om «OJS block count mismatch». De er ufarlige og påvirker bare linjenumre i feilmeldinger.

## Skrive kapitler som Jupyter-notebook

Et kapittel kan skrives som `.ipynb`. Teksten ligger i markdown-celler (som brytes automatisk i
Jupyter), og interaktive figurer hentes inn med én linje i en egen markdown-celle:

```
{{< include ../_interaktiv/celle.qmd >}}
```

`{ojs}`-blokker fungerer også direkte i markdown-celler. Første celle i notebooken er en rå celle
med YAML-hodet (`title:` og `jupyter: python3`).

## Hjelpeskript i `tools/`

Installer først det som trengs: `pip install -r tools/requirements.txt`.

| Skript | Hva det gjør |
|---|---|
| `qmd_til_ipynb.py` | Gjør en `.qmd`-fil om til en notebook, delt i celler per overskrift, boks og figur. Testet: den bygde siden blir lik. |
| `lag_molekyldata.py` | Regner ut partialladninger og 2D-koordinater til molekylvisningen i kapittel 2 og skriver `_interaktiv/molekyldata.qmd`. Må bare kjøres hvis du endrer listen over molekyler i skriptet. |

