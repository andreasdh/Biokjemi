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

Nye kapitler legges til i `chapters:` i `_quarto.yml`.

## Publisering

Hver push til `main` bygger og publiserer boka via GitHub Actions (`.github/workflows/publish.yml`).
Slå på én gang: **Settings → Pages → Source: GitHub Actions**.
