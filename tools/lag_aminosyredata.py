"""Lager spill/aminosyredata.js: strukturdata for de 20 proteinogene aminosyrene (og noen ekstra).

Dataene brukes både av kapittel 3 (interaktive figurer) og av aminosyrespillet. Koordinatene
kommer fra RDKit (2D, med ryggraden N–Cα–C lagt likt for alle), og R/S-merkingen regnes ut
med CIP-reglene i RDKit, slik at vi kan sjekke at alle er S (unntatt cystein, som er R).

Bruk:
    python tools/lag_aminosyredata.py

Krever:  pip install rdkit
"""
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, rdDepictor
from rdkit.Geometry import Point2D

UT = Path(__file__).resolve().parent.parent / "spill" / "aminosyredata.js"

# id, navn, tre, en, gruppe, essensiell, SMILES (nøytral, L-form), pKa (α-COOH, α-NH3+, sidekjede), hydropati
# (Kyte–Doolittle), kort tekst, tagger. pKa-verdiene er typiske verdier for frie aminosyrer.
A = [
    ("gly", "Glysin", "Gly", "G", "upolar", False, "NCC(=O)O", (2.34, 9.60, None), -0.4,
     "Den minste aminosyren: sidegruppa er bare et hydrogenatom. Glysin er den eneste aminosyren som ikke er kiral, og den gir kjeden stor fleksibilitet.",
     ["ingen sidegruppe", "fleksibel", "ikke kiral"]),
    ("ala", "Alanin", "Ala", "A", "upolar", False, "C[C@H](N)C(=O)O", (2.34, 9.69, None), 1.8,
     "Sidegruppa er en liten metylgruppe (–CH₃). Alanin er upolar og liten, og passer nesten overalt i et protein.",
     ["London-krefter", "liten"]),
    ("val", "Valin", "Val", "V", "upolar", True, "CC(C)[C@H](N)C(=O)O", (2.32, 9.62, None), 4.2,
     "Forgrenet hydrokarbon-sidegruppe. Valin er svært hydrofob og havner gjerne på innsiden av proteiner.",
     ["London-krefter", "hydrofob", "forgrenet"]),
    ("leu", "Leucin", "Leu", "L", "upolar", True, "CC(C)C[C@H](N)C(=O)O", (2.36, 9.60, None), 3.8,
     "Hydrofob og forgrenet. Leucin er en av de vanligste aminosyrene i proteiner og står ofte i den hydrofobe kjernen.",
     ["London-krefter", "hydrofob", "forgrenet"]),
    ("ile", "Isoleucin", "Ile", "I", "upolar", True, "CC[C@H](C)[C@H](N)C(=O)O", (2.36, 9.68, None), 4.5,
     "Isomer av leucin, og den mest hydrofobe av alle. Isoleucin har to kirale sentre (Cα og Cβ).",
     ["London-krefter", "hydrofob", "to kirale sentre"]),
    ("phe", "Fenylalanin", "Phe", "F", "upolar", True, "N[C@@H](Cc1ccccc1)C(=O)O", (1.83, 9.13, None), 2.8,
     "Sidegruppa har en benzenring. Den er stor, flat og hydrofob, og kan stable seg mot andre ringer (π-stabling).",
     ["London-krefter", "aromatisk", "hydrofob"]),
    ("trp", "Tryptofan", "Trp", "W", "upolar", True, "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", (2.38, 9.39, None), -0.9,
     "Den største aminosyren, med en indolring. Mest upolar, men N–H i ringen kan likevel gi en hydrogenbinding. Absorberer UV-lys.",
     ["London-krefter", "aromatisk", "H-binding (donor)", "størst"]),
    ("met", "Metionin", "Met", "M", "upolar", True, "CSCC[C@H](N)C(=O)O", (2.28, 9.21, None), 1.9,
     "Har et svovelatom i en tioeter (–S–CH₃). Metionin er upolar, og AUG som koder for metionin er startkodonet.",
     ["London-krefter", "inneholder svovel", "startkodon"]),
    ("pro", "Prolin", "Pro", "P", "upolar", False, "OC(=O)[C@@H]1CCCN1", (1.99, 10.96, None), -1.6,
     "Sidegruppa er knyttet til aminogruppa og danner en ring. Prolin er stiv og bøyer kjeden, og kan ikke gi en N–H-hydrogenbinding.",
     ["ring", "stiv", "bøyer kjeden"]),
    ("ser", "Serin", "Ser", "S", "polar", False, "N[C@@H](CO)C(=O)O", (2.21, 9.15, None), -0.8,
     "En hydroksylgruppe (–OH) gjør sidegruppa polar. Den kan både gi og ta imot hydrogenbindinger.",
     ["H-binding (donor og akseptor)", "alkohol"]),
    ("thr", "Treonin", "Thr", "T", "polar", True, "C[C@@H](O)[C@H](N)C(=O)O", (2.11, 9.62, None), -0.7,
     "Som serin, men med en ekstra metylgruppe. Treonin har to kirale sentre.",
     ["H-binding (donor og akseptor)", "alkohol", "to kirale sentre"]),
    ("cys", "Cystein", "Cys", "C", "polar", False, "N[C@@H](CS)C(=O)O", (1.96, 10.28, 8.18), 2.5,
     "Tiolgruppen (–SH) kan oksideres, og to cysteiner kan danne en kovalent disulfidbro (–S–S–). Cystein er R, ikke S, fordi S har høyere atomnummer enn O.",
     ["disulfidbro", "inneholder svovel", "R-konfigurasjon"]),
    ("asn", "Asparagin", "Asn", "N", "polar", False, "N[C@@H](CC(N)=O)C(=O)O", (2.02, 8.80, None), -3.5,
     "En amidgruppe (–CONH₂): uladd, men polar, og den danner lett hydrogenbindinger.",
     ["H-binding (donor og akseptor)", "amid"]),
    ("gln", "Glutamin", "Gln", "Q", "polar", False, "N[C@@H](CCC(N)=O)C(=O)O", (2.17, 9.13, None), -3.5,
     "Som asparagin, men med en CH₂-gruppe til. Amidet gjør den uladd og polar.",
     ["H-binding (donor og akseptor)", "amid"]),
    ("tyr", "Tyrosin", "Tyr", "Y", "polar", False, "N[C@@H](Cc1ccc(O)cc1)C(=O)O", (2.20, 9.11, 10.07), -1.3,
     "En benzenring med en OH-gruppe. Tyrosin er både aromatisk og polar, og OH-gruppa kan avgi et proton ved høy pH.",
     ["aromatisk", "H-binding (donor og akseptor)", "fenol"]),
    ("asp", "Asparaginsyre", "Asp", "D", "negativ", False, "N[C@@H](CC(=O)O)C(=O)O", (1.88, 9.60, 3.65), -3.5,
     "Sidegruppa har en karboksylgruppe som er deprotonert (COO⁻) ved fysiologisk pH. I proteiner kalles den ofte aspartat.",
     ["negativ ladning", "ionebinding", "H-binding (akseptor)"]),
    ("glu", "Glutaminsyre", "Glu", "E", "negativ", False, "N[C@@H](CCC(=O)O)C(=O)O", (2.19, 9.67, 4.25), -3.5,
     "Som asparaginsyre, men med en CH₂-gruppe til. I proteiner kalles den ofte glutamat. Natriumsaltet er smaksforsterkeren MSG.",
     ["negativ ladning", "ionebinding", "H-binding (akseptor)"]),
    ("lys", "Lysin", "Lys", "K", "positiv", True, "N[C@@H](CCCCN)C(=O)O", (2.18, 8.95, 10.53), -3.9,
     "Lang sidekjede med en aminogruppe som er positivt ladd (–NH₃⁺) ved fysiologisk pH. Lysin binder gjerne til negativt ladde molekyler som DNA.",
     ["positiv ladning", "ionebinding", "H-binding (donor)"]),
    ("arg", "Arginin", "Arg", "R", "positiv", False, "N[C@@H](CCCNC(N)=N)C(=O)O", (2.17, 9.04, 12.48), -4.5,
     "Guanidiniumgruppa er svært basisk og alltid positivt ladd i kroppen, og kan gi mange hydrogenbindinger samtidig.",
     ["positiv ladning", "ionebinding", "H-binding (donor)"]),
    ("his", "Histidin", "His", "H", "positiv", True, "N[C@@H](Cc1c[nH]cn1)C(=O)O", (1.82, 9.17, 6.00), -3.2,
     "Imidazolringen har pKₐ rundt 6, så histidin kan være både ladd og uladd nær fysiologisk pH. Derfor er den viktig i mange enzymers aktive sete.",
     ["kan være ladd eller uladd", "buffer", "aromatisk"]),
]

# Ikke blant de 20: vises uten pH-avhengig ladning
EKSTRA = [
    ("sec", "Selenocystein", "Sec", "U", "ekstra", False, "N[C@@H]([CH2][SeH])C(=O)O",
     "Den 21. aminosyren: som cystein, men med selen i stedet for svovel. Kodes av et stoppkodon (UGA) i en spesiell sammenheng."),
    ("pyl", "Pyrrolysin", "Pyl", "O", "ekstra", False, "CC1CC=NC1C(=O)NCCCC[C@H](N)C(=O)O",
     "Den 22. aminosyren, funnet i noen arkebakterier. Den er et lysin med en pyrrolinring koblet til."),
    ("orn", "Ornitin", "Orn", "", "ekstra", False, "NCCC[C@H](N)C(=O)O",
     "Ikke kodet av DNA og ikke i proteiner, men et mellomprodukt i urea-syklusen. Som lysin, men med en CH₂-gruppe mindre."),
    ("cit", "Citrullin", "Cit", "", "ekstra", False, "NC(=O)NCCC[C@H](N)C(=O)O",
     "Ikke i proteiner direkte (dannes fra arginin etter at proteinet er laget), og et mellomprodukt i urea-syklusen."),
    ("aspartam", "Aspartam", "", "", "ekstra", False, "COC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](N)CC(=O)O",
     "Søtstoff laget av to aminosyrer: asparaginsyre og fenylalanin (som metylester). Aminosyrene er koblet med en peptidbinding."),
    ("glutamat", "Glutamat (MSG)", "", "", "ekstra", False, "N[C@@H](CCC(=O)[O-])C(=O)O",
     "Natriumglutamat (MSG) er natriumsaltet av glutaminsyre og brukes som smaksforsterker (umami)."),
]

# ryggradens plassering (N, Cα, C) slik at alle aminosyrer tegnes på samme måte
RYGG = {"N": (-1.3, -0.75), "CA": (0.0, 0.0), "C": (1.3, -0.75)}

MOENSTER_RYGG = Chem.MolFromSmarts("[NX3]-[CX4]-[CX3](=O)-[OX2H1,OX1-]")


def finn_ryggrad(mol):
    m = mol.GetSubstructMatches(MOENSTER_RYGG)
    # for peptidlignende ekstra-molekyler kan det være flere treff: bruk det med frie N og C-terminal
    for n, ca, c, o1, o2 in m:
        if mol.GetAtomWithIdx(n).GetTotalNumHs() >= 1 or mol.GetAtomWithIdx(n).GetDegree() <= 2:
            if mol.GetAtomWithIdx(c).GetDegree() == 3:
                return n, ca, c, o1, o2
    return m[0]


def sidekjedeatomer(mol, ca, n, c):
    """Alle atomer som henger på Cα utenom N og C (og det som er koblet til dem), dvs. sidekjeden."""
    if mol.GetBondBetweenAtoms(n, ca) is None:
        return set()
    sc, stack = set(), [a.GetIdx() for a in mol.GetAtomWithIdx(ca).GetNeighbors() if a.GetIdx() not in (n, c)]
    while stack:
        i = stack.pop()
        if i in sc:
            continue
        sc.add(i)
        for nb in mol.GetAtomWithIdx(i).GetNeighbors():
            j = nb.GetIdx()
            if j not in sc and j not in (ca, n, c):
                stack.append(j)
    return sc


def prolin_ring(mol, ca, n):
    """For prolin henger N i samme ring som sidekjeden: regn N som en del av både ryggrad og sidekjede."""
    ri = mol.GetRingInfo()
    return any(ca in r and n in r for r in ri.AtomRings())


def lag_molekyl(id_, smiles, aa=True):
    mol = Chem.MolFromSmiles(smiles)
    n, ca, c, o_dbl, o_oh = finn_ryggrad(mol) if aa else (None,) * 5
    # sørg for at O med dobbeltbinding er O=
    if aa:
        b = mol.GetBondBetweenAtoms(c, o_dbl)
        if b.GetBondTypeAsDouble() != 2:
            o_dbl, o_oh = o_oh, o_dbl
    sc = sidekjedeatomer(mol, ca, n, c) if aa else set()
    cip = ""
    if aa:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        from rdkit.Chem import rdCIPLabeler

        rdCIPLabeler.AssignCIPLabels(mol)
        a = mol.GetAtomWithIdx(ca)
        cip = a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else ""

    # --- protonasjonssteder (fullt protonert form) ---
    site = {}
    q = {}
    if aa:
        mol.GetAtomWithIdx(n).SetFormalCharge(1)
        site[n] = "n"
        site[o_oh] = "a"
        r = id_
        if r == "lys":
            for at in mol.GetAtoms():
                if at.GetIdx() in sc and at.GetSymbol() == "N":
                    at.SetFormalCharge(1)
                    site[at.GetIdx()] = "r"
        elif r == "arg":
            for at in mol.GetAtoms():
                if at.GetIdx() in sc and at.GetSymbol() == "N":
                    b = [bb for bb in at.GetBonds() if bb.GetBondTypeAsDouble() == 2]
                    if b:
                        at.SetFormalCharge(1)
                        site[at.GetIdx()] = "r"
        elif r == "his":
            for at in mol.GetAtoms():
                if at.GetIdx() in sc and at.GetSymbol() == "N" and at.GetTotalNumHs() == 0 and at.GetIsAromatic():
                    at.SetFormalCharge(1)
                    at.SetNumExplicitHs(1)
                    at.SetNoImplicit(True)
                    site[at.GetIdx()] = "r"
        elif r in ("asp", "glu"):
            for at in mol.GetAtoms():
                if at.GetIdx() in sc and at.GetSymbol() == "O" and at.GetTotalNumHs() == 1:
                    site[at.GetIdx()] = "r"
        elif r == "tyr":
            for at in mol.GetAtoms():
                if at.GetIdx() in sc and at.GetSymbol() == "O" and at.GetTotalNumHs() == 1:
                    site[at.GetIdx()] = "r"
        elif r == "cys":
            for at in mol.GetAtoms():
                if at.GetIdx() in sc and at.GetSymbol() == "S":
                    site[at.GetIdx()] = "r"
    Chem.SanitizeMol(mol)

    # --- 2D-koordinater med fast ryggrad ---
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)  # utgangspunkt
    if aa:
        cm = {n: Point2D(*RYGG["N"]), ca: Point2D(*RYGG["CA"]), c: Point2D(*RYGG["C"])}
        rdDepictor.Compute2DCoords(mol, coordMap=cm, canonOrient=False)
    conf = mol.GetConformer()
    Chem.WedgeMolBonds(mol, conf)

    h_tot = {a.GetIdx(): a.GetTotalNumHs() for a in mol.GetAtoms()}
    kek = Chem.Mol(mol)
    Chem.Kekulize(kek, clearAromaticFlags=True)

    atomer = []
    for a in mol.GetAtoms():
        i = a.GetIdx()
        p = conf.GetAtomPosition(i)
        d = {"el": a.GetSymbol(), "x": round(p.x, 3), "y": round(p.y, 3), "h": h_tot[i], "q": a.GetFormalCharge()}
        if i in sc:
            d["sc"] = 1
        if i in site:
            d["site"] = site[i]
        if aa and i == ca:
            d["ca"] = 1
        atomer.append(d)
    bindinger = []
    for b in kek.GetBonds():
        w = 0
        bd = mol.GetBondWithIdx(b.GetIdx()).GetBondDir()
        if bd == Chem.BondDir.BEGINWEDGE:
            w = 1
        elif bd == Chem.BondDir.BEGINDASH:
            w = -1
        bindinger.append([b.GetBeginAtomIdx(), b.GetEndAtomIdx(), int(b.GetBondTypeAsDouble()), w])
    # beregn masse fra nøytral form
    neutral = Chem.MolFromSmiles(smiles)
    mw = round(Descriptors.MolWt(neutral), 2)
    return atomer, bindinger, cip, mw, (n, ca, c)


def main():
    ut = {"aminosyrer": [], "ekstra": []}
    for (id_, navn, tre, en, gruppe, ess, smi, pka, kd, tekst, tagger) in A:
        atomer, bindinger, cip, mw, rygg = lag_molekyl(id_, smi)
        if id_ == "gly":
            cip = ""
        ut["aminosyrer"].append({
            "id": id_, "navn": navn, "tre": tre, "en": en, "gruppe": gruppe, "essensiell": ess,
            "mw": mw, "kd": kd, "pka": list(pka), "cip": cip, "tekst": tekst, "tagger": tagger,
            "atomer": atomer, "bindinger": bindinger,
        })
    for (id_, navn, tre, en, gruppe, ess, smi, tekst) in EKSTRA:
        atomer, bindinger, cip, mw, rygg = lag_molekyl(id_, smi, aa=False)
        ut["ekstra"].append({"id": id_, "navn": navn, "tre": tre, "en": en, "gruppe": gruppe,
                             "mw": mw, "tekst": tekst, "atomer": atomer, "bindinger": bindinger})
    js = ("// Generert av tools/lag_aminosyredata.py. Ikke rediger for hånd.\n"
          "window.BKAA = window.BKAA || {};\n"
          "BKAA.data = " + json.dumps(ut, ensure_ascii=False, separators=(",", ":")) + ";\n")
    UT.parent.mkdir(exist_ok=True)
    UT.write_text(js, encoding="utf-8")
    print(f"Skrev {UT} ({len(js) // 1024} kB).")
    print("CIP:", {a['id']: a['cip'] for a in ut['aminosyrer']})


if __name__ == "__main__":
    main()
