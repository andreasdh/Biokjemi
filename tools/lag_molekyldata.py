"""Lager datagrunnlaget til molekylvisningen i kapittel 2 (_interaktiv/molekyldata.qmd).

Kjør fra prosjektmappen:   python tools/lag_molekyldata.py
Krever:                    pip install rdkit scipy numpy

Partialladningene er Gasteiger-ladninger fra RDKit. De er en forenklet modell som er
nyttig for å vise *hvor* elektrontettheten sitter, ikke for presise beregninger.
"""
import json
import math
import random
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from scipy.optimize import minimize

BOND = 1.5          # standard bindingslengde i RDKit sine 2D-koordinater
H_LENGDE = 0.72     # H-bindinger tegnes litt kortere


# ----------------------------------------------------------------------------------
# Hjelpefunksjoner
# ----------------------------------------------------------------------------------
def rund(x, n=3):
    return round(float(x), n)


def plasser_h(senter, naboer, antall, andre=()):
    """Finn tegnevinkler til `antall` H-atomer rundt et atom, gitt posisjonene til naboene.
    `andre` er posisjoner til alle andre atomer, brukt til å velge den siden med mest plass."""
    vinkler = sorted(math.atan2(n[1] - senter[1], n[0] - senter[0]) for n in naboer)
    if not vinkler:
        return [2 * math.pi * k / antall - math.pi / 2 for k in range(antall)]

    def plass(kandidater):
        # sum av avstand til nærmeste andre atom: større er bedre
        def mål(v):
            p = np.array([senter[0] + BOND * math.cos(v), senter[1] + BOND * math.sin(v)])
            return min([np.linalg.norm(p - np.array(a)) for a in andre] or [9])
        return mål(kandidater) if not isinstance(kandidater, list) else sum(mål(v) for v in kandidater)

    if len(vinkler) == 1:
        v = vinkler[0]
        if antall == 1:
            alt = [v + math.radians(120), v - math.radians(120)]
            return [max(alt, key=plass)]
        if antall == 2:
            return [v + math.radians(120), v - math.radians(120)]
        if antall == 3:
            v += math.pi
            return [v - math.radians(70), v, v + math.radians(70)]
    ut = []
    for _ in range(antall):
        gaps = []
        for i, v in enumerate(vinkler):
            neste = vinkler[(i + 1) % len(vinkler)]
            gaps.append(((neste - v) % (2 * math.pi), v))
        gap, start = max(gaps)
        ny = start + gap / 2
        vinkler.append(ny)
        vinkler.sort()
        ut.append(ny)
    return ut


def bygg(smiles, alle_h=False, manuell=None, symmetriser=None):
    """Bygg en molekylbeskrivelse (atomer, bindinger, ladninger) fra SMILES."""
    mola = Chem.MolFromSmiles(smiles)
    mol0 = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(mol0, clearAromaticFlags=True)
    molh = Chem.AddHs(mol0)
    AllChem.ComputeGasteigerCharges(molh)
    q = [float(a.GetProp("_GasteigerCharge")) for a in molh.GetAtoms()]
    if symmetriser:
        for gruppe in symmetriser:
            snitt = sum(q[i] for i in gruppe) / len(gruppe)
            for i in gruppe:
                q[i] = snitt
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol0)
    conf = mol0.GetConformer()
    xy = {a.GetIdx(): np.array([conf.GetAtomPosition(a.GetIdx()).x,
                                -conf.GetAtomPosition(a.GetIdx()).y]) for a in mol0.GetAtoms()}
    if manuell:
        for i, p in manuell.items():
            xy[i] = np.array(p, dtype=float)

    atomer = []
    for a in mol0.GetAtoms():
        atomer.append({"el": a.GetSymbol(), "x": xy[a.GetIdx()][0], "y": xy[a.GetIdx()][1],
                       "q": q[a.GetIdx()], "fq": a.GetFormalCharge(), "skjul": 0})
    bindinger = [[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), int(b.GetBondTypeAsDouble())]
                 for b in mol0.GetBonds()]

    # hydrogener (indeks fra n_tung og oppover i molh)
    for a in mol0.GetAtoms():
        i = a.GetIdx()
        hs = [n.GetIdx() for n in molh.GetAtomWithIdx(i).GetNeighbors() if n.GetSymbol() == "H"]
        if not hs:
            continue
        vis = alle_h or a.GetSymbol() in ("N", "O", "S", "P")
        naboer = [xy[n.GetIdx()] for n in a.GetNeighbors()]
        andre = [xy[k] for k in xy if k != i and k not in [n.GetIdx() for n in a.GetNeighbors()]]
        vinkler = plasser_h(xy[i], naboer, len(hs), andre) if vis else [0] * len(hs)
        for h_idx, v in zip(hs, vinkler):
            if vis:
                px = xy[i][0] + H_LENGDE * BOND * math.cos(v)
                py = xy[i][1] + H_LENGDE * BOND * math.sin(v)
            else:
                px, py = xy[i]
            atomer.append({"el": "H", "x": px, "y": py, "q": q[h_idx], "fq": 0,
                           "skjul": 0 if vis else 1, "_idx": h_idx})
            bindinger.append([len(atomer) - 1, i, 1])
    # (rekkefølgen i `atomer` er nå: tunge atomer, så H-er)
    # flagg for hydrogenbindingsdonor og -akseptor
    for k, at in enumerate(atomer):
        at["d"] = 0
        at["a"] = 0
    for (i, j, _) in bindinger:
        for h, t in ((i, j), (j, i)):
            if atomer[h]["el"] == "H" and atomer[t]["el"] in ("N", "O"):
                atomer[h]["d"] = 1
    for a in mol0.GetAtoms():
        i = a.GetIdx()
        el = a.GetSymbol()
        aa = mola.GetAtomWithIdx(i)
        if el == "O":
            atomer[i]["a"] = 1
        elif el == "N" and a.GetFormalCharge() == 0:
            if aa.GetIsAromatic():
                # pyridin-lignende N (to bindinger i ringen, ingen H) er akseptor
                if aa.GetTotalDegree() == 2:
                    atomer[i]["a"] = 1
            else:
                amid = any(n.GetSymbol() == "C" and any(
                    b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(n).GetSymbol() == "O" for b in n.GetBonds())
                    for n in a.GetNeighbors())
                anilin = any(mola.GetAtomWithIdx(n.GetIdx()).GetIsAromatic() for n in a.GetNeighbors())
                if not amid and not anilin:
                    atomer[i]["a"] = 1
    # N-H i ringer (T: N3-H) og N-H i amider er donorer (håndtert over), men ikke akseptorer.
    for at in atomer:
        at.pop("_idx", None)
    return {"atomer": atomer, "bindinger": bindinger, "hb": []}


def sentrer(m, skala=1.0):
    xs = [a["x"] for a in m["atomer"]]
    ys = [a["y"] for a in m["atomer"]]
    cx, cy = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2
    for a in m["atomer"]:
        a["x"] = (a["x"] - cx) * skala
        a["y"] = (a["y"] - cy) * skala
    return m


def til_json(m):
    ut = {"atomer": [], "bindinger": m["bindinger"], "hb": m["hb"]}
    for a in m["atomer"]:
        d = {"el": a["el"], "x": rund(a["x"]), "y": rund(a["y"]), "q": rund(a["q"], 2)}
        if a.get("fq"):
            d["fq"] = a["fq"]
        if a.get("skjul"):
            d["skjul"] = 1
        if a.get("d"):
            d["d"] = 1
        if a.get("a"):
            d["a"] = 1
        ut["atomer"].append(d)
    return ut


# ----------------------------------------------------------------------------------
# Enkeltmolekyler
# ----------------------------------------------------------------------------------
def vann_manuell():
    m = bygg("O", alle_h=True, manuell={0: (0, 0)})
    # H-ene: bøyd geometri, 104,5 grader
    a = math.radians(52.25)
    m["atomer"][1]["x"], m["atomer"][1]["y"] = -BOND * 0.8 * math.sin(a), BOND * 0.8 * math.cos(a)
    m["atomer"][2]["x"], m["atomer"][2]["y"] = BOND * 0.8 * math.sin(a), BOND * 0.8 * math.cos(a)
    return m


def metan_manuell():
    m = bygg("C", alle_h=True, manuell={0: (0, 0)})
    pos = [(-1.15, -0.8), (1.15, -0.8), (-0.75, 1.0), (0.75, 1.0)]
    for at, p in zip(m["atomer"][1:], pos):
        at["x"], at["y"] = p
    return m


def hydrert_ion(ion, antall=4):
    """Na+ (vann vender O mot ionet) eller Cl- (vann vender H mot ionet)."""
    q_o, q_h = -0.41, 0.21
    vann = bygg("O", alle_h=True)
    q_o = vann["atomer"][0]["q"]
    q_h = vann["atomer"][1]["q"]
    atomer = [{"el": ion, "x": 0.0, "y": 0.0, "q": 1.0 if ion == "Na" else -1.0, "fq": 0, "d": 0, "a": 0}]
    bindinger, hb = [], []
    avstand = 3.1 if ion == "Na" else 3.4
    vinkler = [math.radians(v) for v in (-135, -45, 45, 135)][:antall]
    for v in vinkler:
        ox, oy = avstand * math.cos(v), avstand * math.sin(v)
        o_idx = len(atomer)
        atomer.append({"el": "O", "x": ox, "y": oy, "q": q_o, "fq": 0, "d": 0, "a": 1})
        # retning bort fra / mot ionet
        if ion == "Na":
            # O mot ionet: H-ene peker bort, ±52,25 grader fra utover-retningen
            for s in (-1, 1):
                h = v + s * math.radians(52.25)
                atomer.append({"el": "H", "x": ox + H_LENGDE * BOND * math.cos(h),
                               "y": oy + H_LENGDE * BOND * math.sin(h), "q": q_h, "fq": 0, "d": 1, "a": 0})
                bindinger.append([len(atomer) - 1, o_idx, 1])
            hb.append([0, o_idx])
        else:
            # én H mot ionet, den andre ±104,5 grader unna
            mot = v + math.pi
            h1 = len(atomer)
            atomer.append({"el": "H", "x": ox + H_LENGDE * BOND * math.cos(mot),
                           "y": oy + H_LENGDE * BOND * math.sin(mot), "q": q_h, "fq": 0, "d": 1, "a": 0})
            bindinger.append([h1, o_idx, 1])
            h2v = mot + math.radians(104.5) * (1 if math.cos(v) > 0 else -1)
            atomer.append({"el": "H", "x": ox + H_LENGDE * BOND * math.cos(h2v),
                           "y": oy + H_LENGDE * BOND * math.sin(h2v), "q": q_h, "fq": 0, "d": 1, "a": 0})
            bindinger.append([len(atomer) - 1, o_idx, 1])
            hb.append([0, h1])
    return {"atomer": atomer, "bindinger": bindinger, "hb": hb}


def vannnettverk():
    """Fem vannmolekyler: ett i midten som danner fire hydrogenbindinger."""
    vann = bygg("O", alle_h=True)
    q_o, q_h = vann["atomer"][0]["q"], vann["atomer"][1]["q"]
    atomer, bindinger, hb = [], [], []

    def lag_vann(ox, oy, h_vinkler):
        o = len(atomer)
        atomer.append({"el": "O", "x": ox, "y": oy, "q": q_o, "fq": 0, "d": 0, "a": 1})
        hs = []
        for v in h_vinkler:
            atomer.append({"el": "H", "x": ox + H_LENGDE * BOND * math.cos(v),
                           "y": oy + H_LENGDE * BOND * math.sin(v), "q": q_h, "fq": 0, "d": 1, "a": 0})
            bindinger.append([len(atomer) - 1, o, 1])
            hs.append(len(atomer) - 1)
        return o, hs

    d = 3.3
    ap = math.radians(104.5)
    # midtmolekylet: H-ene peker nedover mot to naboer
    v1, v2 = math.radians(90 + 52.25), math.radians(90 - 52.25)
    mid, (hm1, hm2) = lag_vann(0, 0, [v1, v2])
    # naboer som mottar H-bindinger fra midten (nede til venstre og nede til høyre)
    for (v, hm) in ((v1, hm1), (v2, hm2)):
        ox, oy = d * math.cos(v), d * math.sin(v)
        ut = v  # retning bort fra midten
        o, hs = lag_vann(ox, oy, [ut - math.radians(35), ut + math.radians(70)])
        hb.append([hm, o])
    # naboer som gir H-binding til midtmolekylets oksygen (oppe til venstre og oppe til høyre)
    for s in (-1, 1):
        v = math.radians(-90 + s * 50)       # retning fra midten ut til naboen
        ox, oy = d * math.cos(v), d * math.sin(v)
        mot = v + math.pi                    # retning fra naboen mot midten
        o, hs = lag_vann(ox, oy, [mot, mot + s * ap])
        hb.append([hs[0], mid])
    return {"atomer": atomer, "bindinger": bindinger, "hb": hb}


# ----------------------------------------------------------------------------------
# Par av molekyler holdt sammen av hydrogenbindinger
# ----------------------------------------------------------------------------------
def h_pa(m, tung_idx):
    """Indekser til alle H-atomer bundet til et tungt atom."""
    ut = []
    for a, b, _ in m["bindinger"]:
        for h, t in ((a, b), (b, a)):
            if t == tung_idx and m["atomer"][h]["el"] == "H":
                ut.append(h)
    return ut


def tilpass_par(m1, m2, lenker, seed=1):
    """Roter/flytt m2 slik at hydrogenbindingene mellom m1 og m2 får realistisk geometri.
    lenker: liste av (side_donor, donoratom, akseptoratom)."""
    rng = random.Random(seed)
    P1 = np.array([[a["x"], a["y"]] for a in m1["atomer"]])
    P2 = np.array([[a["x"], a["y"]] for a in m2["atomer"]])
    c2 = P2.mean(axis=0)
    synlig1 = [i for i, a in enumerate(m1["atomer"]) if not a["skjul"]]
    synlig2 = [i for i, a in enumerate(m2["atomer"]) if not a["skjul"]]

    def transform(params, speil):
        th, tx, ty = params
        P = P2 - c2
        if speil:
            P = P * np.array([1, -1])
        R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
        return P @ R.T + np.array([tx, ty])

    def kostnad(params, speil):
        Q = transform(params, speil)
        k = 0.0
        for (side, donor, akseptor) in lenker:
            if side == 1:
                D, A = P1[donor], Q[akseptor]
                hs = [P1[h] for h in h_pa(m1, donor)]
            else:
                D, A = Q[donor], P1[akseptor]
                hs = [Q[h] for h in h_pa(m2, donor)]
            k += (np.linalg.norm(D - A) - 3.0) ** 2 * 4
            beste_h = min(hs, key=lambda h: np.linalg.norm(h - A))
            k += (np.linalg.norm(beste_h - A) - 2.0) ** 2 * 10
            v1, v2 = D - beste_h, A - beste_h
            cos = float(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9))
            k += (1 + cos) ** 2 * 3
        for a in synlig1:
            for b in synlig2:
                d = np.linalg.norm(P1[a] - Q[b])
                if d < 1.45:
                    k += (1.45 - d) ** 2 * 6
        return k

    beste = (1e9, None, None)
    for speil in (False, True):
        for _ in range(80):
            start = [rng.uniform(0, 2 * math.pi), rng.uniform(-6, 6), rng.uniform(-6, 6)]
            res = minimize(kostnad, start, args=(speil,), method="Nelder-Mead",
                           options={"xatol": 1e-4, "fatol": 1e-7, "maxiter": 3000})
            if res.fun < beste[0]:
                beste = (res.fun, res.x, speil)
    _, params, speil = beste
    Q = transform(params, speil)
    for a, p in zip(m2["atomer"], Q):
        a["x"], a["y"] = p[0], p[1]
    return beste[0]


def lag_par(smiles1, smiles2, lenker):
    """lenker: liste av (side_donor, donoratom, akseptoratom) der side_donor er 1 eller 2.
    Atomindeksene gjelder tunge atomer i hver sin SMILES. Tegner stiplet linje fra
    donor-H (den som peker mot akseptoren) til akseptoratomet."""
    m1 = bygg(smiles1)
    m2 = bygg(smiles2)
    sentrer(m1)
    sentrer(m2)
    kost = tilpass_par(m1, m2, lenker)
    n1 = len(m1["atomer"])
    m = {"atomer": m1["atomer"] + m2["atomer"],
         "bindinger": m1["bindinger"] + [[a + n1, b + n1, o] for a, b, o in m2["bindinger"]],
         "hb": []}
    for (side, donor, akseptor) in lenker:
        if side == 1:
            kandidater, A = h_pa(m1, donor), akseptor + n1
        else:
            kandidater, A = [h + n1 for h in h_pa(m2, donor)], akseptor
        pa = np.array([m["atomer"][A]["x"], m["atomer"][A]["y"]])
        h = min(kandidater, key=lambda k: np.linalg.norm(np.array([m["atomer"][k]["x"], m["atomer"][k]["y"]]) - pa))
        m["hb"].append([h, A])
    m["_kost"] = kost
    print("   tilpasning", round(kost, 3))
    return m


# ----------------------------------------------------------------------------------
# Katalog
# ----------------------------------------------------------------------------------
def katalog():
    ut = []

    def legg(id_, navn, gruppe, tekst, m=None, **ekstra):
        post = {"id": id_, "navn": navn, "gruppe": gruppe, "tekst": tekst}
        if m is not None:
            post.update(til_json(sentrer(m)))
        post.update(ekstra)
        ut.append(post)

    # --- Funksjonelle grupper ---
    legg("metan", "Metan (CH₄)", "gruppe",
         "C–H-bindingen er nesten upolar: elektronene deles omtrent likt, og det blir nesten ingen "
         "delladninger. Slike molekyler kan ikke danne hydrogenbindinger.",
         metan_manuell())
    legg("butan", "Hydrokarbonkjede (butan)", "gruppe",
         "Hydrokarbonkjeder, som i fettsyrer og i sidekjedene til leucin og valin, er upolare. "
         "De trives dårlig i vann og trives godt sammen med hverandre.",
         bygg("CCCC"))
    legg("benzen", "Aromatisk ring (benzen)", "gruppe",
         "Aromatiske ringer, som i fenylalanin, er også upolare. Den store, bevegelige "
         "π-elektronskyen gjør dem gode til London-krefter og stabling.",
         bygg("c1ccccc1"))
    legg("etanol", "Hydroksyl (etanol)", "gruppe",
         "O–H er polar: O er δ⁻ og H er δ⁺. Hydroksylgruppen kan både gi og ta imot hydrogenbindinger "
         "(donor og akseptor). Finnes i serin, treonin og sukker.",
         bygg("CCO"))
    legg("acetaldehyd", "Aldehyd (acetaldehyd)", "gruppe",
         "C=O er sterkt polar: O er δ⁻ og C er δ⁺. Karbonylgruppen kan ta imot hydrogenbindinger, "
         "men har ingen H å gi (akseptor).",
         bygg("CC=O"))
    legg("aceton", "Keton (aceton)", "gruppe",
         "Som aldehydet: en polar C=O som er hydrogenbindings-akseptor.",
         bygg("CC(C)=O"))
    legg("eddiksyre", "Karboksyl (eddiksyre og acetat)", "gruppe",
         "Karboksylgruppen kan gi fra seg et proton. Da flyttes elektrontettheten, og ladningen "
         "deles mellom de to oksygenatomene. Ved fysiologisk pH er den som regel deprotonert og negativ "
         "(aspartat, glutamat).",
         pKa=4.76,
         syre={"navn": "eddiksyre (–COOH)", **til_json(sentrer(bygg("CC(=O)O")))},
         base={"navn": "acetat (–COO⁻)",
               **til_json(sentrer(bygg("CC(=O)[O-]", symmetriser=[[2, 3]])))})
    legg("metylamin", "Amino (metylamin og metylammonium)", "gruppe",
         "Aminogruppen tar lett opp et proton. Da blir den positivt ladd (–NH₃⁺) og kan danne "
         "ionebindinger med negative grupper. Ved fysiologisk pH er amino-gruppen i lysin positiv.",
         pKa=10.6,
         syre={"navn": "metylammonium (–NH₃⁺)", **til_json(sentrer(bygg("C[NH3+]")))},
         base={"navn": "metylamin (–NH₂)", **til_json(sentrer(bygg("CN")))})
    legg("fosfat", "Fosfat (metylfosfat)", "gruppe",
         "Fosfatgrupper (DNA-ryggraden, ATP) bærer typisk to negative ladninger ved fysiologisk pH. "
         "Mange negative ladninger så tett gjør at ATP er et «spent» molekyl.",
         pKa=6.3,
         syre={"navn": "R–O–PO₃H⁻ (monoanion)", **til_json(sentrer(bygg("COP(=O)(O)[O-]", symmetriser=[[3, 5]])))},
         base={"navn": "R–O–PO₃²⁻ (dianion)", **til_json(sentrer(bygg("COP(=O)([O-])[O-]", symmetriser=[[3, 4, 5]])))})
    legg("tiol", "Sulfhydryl (metantiol)", "gruppe",
         "S–H er bare svakt polar (S og H har nesten lik elektronegativitet). Cystein kan likevel "
         "danne disulfidbroer, en kovalent binding som stabiliserer proteiner.",
         bygg("CS"))
    legg("peptid", "Peptidbinding (N-metylacetamid)", "gruppe",
         "Skjelettet i alle proteiner: C=O (δ⁻, akseptor) og N–H (δ⁺, donor) kan danne hydrogenbindinger. "
         "Disse binder sammen α-heliks og β-flak.",
         bygg("CNC(C)=O"))
    legg("vann", "Vann (H₂O)", "gruppe",
         "O er mye mer elektronegativ enn H. O blir δ⁻, H blir δ⁺, og siden molekylet er bøyd, "
         "har hele molekylet en positiv og en negativ side: det er polart.",
         vann_manuell())

    # --- Hydrogenbindinger og ioner ---
    legg("vannnettverk", "Vann: hydrogenbindinger", "binding",
         "Hvert vannmolekyl kan danne opptil fire hydrogenbindinger: to som donor (δ⁺ H) og to som "
         "akseptor (frie elektronpar på O). Nettverket brytes og dannes igjen hele tiden.",
         vannnettverk())
    legg("par_at", "Basepar A–T", "binding",
         "Adenin og tymin kjenner hverandre igjen: to hydrogenbindinger (N–H···O og N–H···N).",
         lag_par("Cn1cnc2c(N)ncnc12", "Cc1cn(C)c(=O)[nH]c1=O",
                 [(1, 6, 9), (2, 7, 7)]))
    legg("par_gc", "Basepar G–C", "binding",
         "Guanin og cytosin danner tre hydrogenbindinger, og et GC-par er derfor litt sterkere enn et AT-par.",
         lag_par("Cn1cnc2c1nc(N)[nH]c2=O", "Cn1ccc(N)nc1=O",
                 [(2, 5, 11), (1, 9, 6), (1, 8, 8)]))
    legg("par_peptid", "To peptidgrupper", "binding",
         "N–H (donor) i én peptidgruppe binder til C=O (akseptor) i en annen. Slike bindinger "
         "gjentas mange ganger i α-heliks og β-flak.",
         lag_par("CNC(C)=O", "CNC(C)=O", [(1, 1, 4)]))
    legg("na_vann", "Na⁺ i vann", "binding",
         "Elektronene er ikke delt i en ionebinding. Når salt løses, vender δ⁻ O i vann seg mot Na⁺ "
         "(ion–dipol-interaksjon) og skjermer ionet.",
         hydrert_ion("Na"))
    legg("cl_vann", "Cl⁻ i vann", "binding",
         "For anioner er det omvendt: δ⁺ H i vann vender seg mot Cl⁻. Slik kan vann bryte "
         "ionegitteret og holde begge ionene i løsning.",
         hydrert_ion("Cl"))
    return ut


def js_streng(tekst, innrykk, bredde=88):
    """Skriv en lang streng som flerlinjes mal-streng (linjeskift blir mellomrom i HTML)."""
    import textwrap
    linjer = textwrap.wrap(tekst, width=bredde)
    if len(linjer) == 1:
        return json.dumps(tekst, ensure_ascii=False)
    pad = " " * (innrykk + 2)
    return "`" + ("\n" + pad).join(linjer) + "`"


def js_liste(liste, innrykk, per_linje=7):
    deler = [json.dumps(x) for x in liste]
    if len(deler) <= per_linje:
        return "[" + ", ".join(deler) + "]"
    pad = " " * (innrykk + 2)
    rader = [", ".join(deler[k:k + per_linje]) for k in range(0, len(deler), per_linje)]
    return "[\n" + ",\n".join(pad + r for r in rader) + "\n" + " " * innrykk + "]"


def main():
    data = katalog()
    linjer = ["```{ojs}", "//| echo: false",
              "// Generert av tools/lag_molekyldata.py. Ikke rediger for hånd.",
              "molekyldata = ["]

    def skriv_geom(g, innrykk):
        pad = " " * innrykk
        ut = [f"{pad}atomer: ["]
        for a in g["atomer"]:
            ut.append(f"{pad}  {json.dumps(a, ensure_ascii=False)},")
        ut.append(f"{pad}],")
        ut.append(f"{pad}bindinger: {js_liste(g['bindinger'], innrykk, 8)},")
        ut.append(f"{pad}hb: {js_liste(g['hb'], innrykk, 8)},")
        return ut

    for post in data:
        linjer.append("  {")
        for nokkel, verdi in post.items():
            if nokkel in ("atomer", "bindinger", "hb", "syre", "base"):
                continue
            if nokkel == "tekst":
                linjer.append(f"    tekst: {js_streng(verdi, 4)},")
            else:
                linjer.append(f"    {nokkel}: {json.dumps(verdi, ensure_ascii=False)},")
        if "atomer" in post:
            linjer += skriv_geom(post, 4)
        for form in ("syre", "base"):
            if form in post:
                linjer.append(f"    {form}: {{")
                linjer.append(f"      navn: {json.dumps(post[form]['navn'], ensure_ascii=False)},")
                linjer += skriv_geom(post[form], 6)
                linjer.append("    },")
        linjer.append("  },")
    linjer += ["]", "```", ""]
    sti = Path(__file__).resolve().parent.parent / "_interaktiv" / "molekyldata.qmd"
    sti.write_text("\n".join(linjer), encoding="utf-8")
    print("skrev", sti, len(linjer), "linjer")


if __name__ == "__main__":
    main()
