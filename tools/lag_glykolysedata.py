"""Lager spill/glykolysedata.js: strukturer, trinn og tall til kapitlet «Stoffskiftet og glykolysen».

Molekylene tegnes som åpne karbonkjeder i sikksakk (lett å se hva som endrer seg), med
fosfatgrupper som en liten oransje sirkel. Koordinatene lages her, uten RDKit.

Bruk:
    python tools/lag_glykolysedata.py

Teksten til trinnene står i listen STEG under, så den er lett å redigere.
"""
import json
import math
from pathlib import Path

UT = Path(__file__).resolve().parent.parent / "spill" / "glykolysedata.js"

SKALA = 1.5  # samme enhet som RDKit-koordinatene i aminosyredata.js (bindingslengde 1,5)


def retn(vinkel, avstand):
    a = math.radians(vinkel)
    return avstand * math.cos(a), avstand * math.sin(a)


def lag(id_, navn, kjede, dobbel=()):
    """kjede: liste med substituenter per karbon: "=O", "OH", "O-", "OP" (ester-O + fosfat).
    dobbel: indekser (1-basert) i på karbon i der bindingen i–(i+1) er en dobbeltbinding."""
    n = len(kjede)
    atomer, bind, tags = [], [], {}

    def ny(el, x, y, h=0, q=0, tag=None):
        atomer.append({"el": el, "x": round(x * SKALA, 3), "y": round(y * SKALA, 3), "h": h, "q": q})
        if tag:
            tags[tag] = len(atomer) - 1
        return len(atomer) - 1

    C = []
    for i in range(1, n + 1):
        C.append(ny("C", (i - 1) * math.cos(math.radians(30)), 0.5 if i % 2 == 1 else 0.0, 0, 0, f"C{i}"))
    for i in range(1, n):
        bind.append([C[i - 1], C[i], 2 if i in dobbel else 1, 0])

    valens = [sum(b[2] for b in bind if c in (b[0], b[1])) for c in C]
    for i in range(1, n + 1):
        subs = kjede[i - 1]
        opp = i % 2 == 1
        if n == 1:
            frie = [90, 270]
        elif i == 1:
            frie = [90, 210]
        elif i == n:
            frie = [90, 330] if opp else [270, 30]
        else:
            frie = [90] if opp else [270]
        for k, s in enumerate(subs):
            v = frie[k]
            c = C[i - 1]
            cx, cy = atomer[c]["x"] / SKALA, atomer[c]["y"] / SKALA
            suff = "" if k == 0 else "b"
            dx, dy = retn(v, 1.0)
            if s == "=O":
                o = ny("O", cx + dx, cy + dy, 0, 0, f"O{i}{suff}")
                bind.append([c, o, 2, 0]); valens[i - 1] += 2
            elif s == "OH":
                o = ny("O", cx + dx, cy + dy, 1, 0, f"O{i}{suff}")
                bind.append([c, o, 1, 0]); valens[i - 1] += 1
            elif s == "O-":
                o = ny("O", cx + dx, cy + dy, 0, -1, f"O{i}{suff}")
                bind.append([c, o, 1, 0]); valens[i - 1] += 1
            elif s == "OP":
                o = ny("O", cx + dx, cy + dy, 0, 0, f"O{i}{suff}")
                dx2, dy2 = retn(v, 2.0)
                p = ny("P", cx + dx2, cy + dy2, 0, 0, f"P{i}{suff}")
                bind.append([c, o, 1, 0]); bind.append([o, p, 1, 0]); valens[i - 1] += 1
            else:
                raise ValueError(s)
    for i, c in enumerate(C):
        atomer[c]["h"] = 4 - valens[i]
    return {"id": id_, "navn": navn, "atomer": atomer, "bindinger": bind, "tags": tags}


M = [
    lag("glc", "Glukose", [["=O"], ["OH"], ["OH"], ["OH"], ["OH"], ["OH"]]),
    lag("g6p", "Glukose-6-fosfat", [["=O"], ["OH"], ["OH"], ["OH"], ["OH"], ["OP"]]),
    lag("f6p", "Fruktose-6-fosfat", [["OH"], ["=O"], ["OH"], ["OH"], ["OH"], ["OP"]]),
    lag("fbp", "Fruktose-1,6-bisfosfat", [["OP"], ["=O"], ["OH"], ["OH"], ["OH"], ["OP"]]),
    lag("dhap", "Dihydroksyacetonfosfat (DHAP)", [["OP"], ["=O"], ["OH"]]),
    lag("g3p", "Glyseraldehyd-3-fosfat (G3P)", [["=O"], ["OH"], ["OP"]]),
    lag("bpg13", "1,3-Bisfosfoglyserat", [["=O", "OP"], ["OH"], ["OP"]]),
    lag("pg3", "3-Fosfoglyserat", [["=O", "O-"], ["OH"], ["OP"]]),
    lag("pg2", "2-Fosfoglyserat", [["=O", "O-"], ["OP"], ["OH"]]),
    lag("pep", "Fosfoenolpyruvat (PEP)", [["=O", "O-"], ["OP"], []], dobbel=(2,)),
    lag("pyr", "Pyruvat", [["=O", "O-"], ["=O"], []]),
    lag("lac", "Laktat", [["=O", "O-"], ["OH"], []]),
]

AKTER = [
    {"nr": 1, "navn": "Investere", "trinn": [1, 3], "farge": "#3d7ea6",
     "tekst": "Cellen bruker 2 ATP for å gjøre glukose klar."},
    {"nr": 2, "navn": "Dele", "trinn": [4, 5], "farge": "#c9962f",
     "tekst": "Sukkeret med seks karbon deles i to like biter med tre karbon."},
    {"nr": 3, "navn": "Høste", "trinn": [6, 10], "farge": "#3f8468",
     "tekst": "Hver bit gir 2 ATP og 1 NADH. Siden vi har to biter, får vi 4 ATP og 2 NADH."},
]

# hl: hvilke atomer (tags) som utheves i substrat og produkt
STEG = [
    {"nr": 1, "akt": 1, "navn": "Heksokinase", "type": "Fosforylering", "sub": ["glc"], "prod": ["g6p"],
     "inn": ["ATP"], "ut": ["ADP"], "dG0": -16.7,
     "rx": {"sub": [["glc", 1], ["atp", 1]], "prod": [["g6p", 1], ["adp", 1]]},
     "hl": {"glc": ["O6"], "g6p": ["O6", "P6"]},
     "hva": "Glukose får en fosfatgruppe på C6. ATP leverer fosfatet og blir til ADP.",
     "hvorfor": [
         "Fosfatgruppen er negativt ladd. Glukose-6-fosfat blir ikke gjenkjent av glukosetransportørene (GLUT), så det kommer seg ikke ut av cellen igjen. Glukosen er «låst inne».",
         "Fosforyleringen holder også konsentrasjonen av fri glukose i cellen lav. Da fortsetter glukose å strømme inn gjennom GLUT, nedover konsentrasjonsgradienten.",
         "Trinnet er svært gunstig og går i praksis bare én vei. Derfor er heksokinase et av tre steder der cellen kan regulere glykolysen."],
     "videre": "Glukose-6-fosfat er et veikryss: det kan også gå til pentosefosfatveien (NADPH og ribose) eller lagres som glykogen.",
     "mer": "Heksokinase lukker seg rundt glukose og ATP («induced fit»). Det holder vann unna, så ATP ikke bare hydrolyseres. I leveren finnes en variant, glukokinase, som ikke hemmes av glukose-6-fosfat.",
     "reg": "Hemmes av produktet glukose-6-fosfat."},
    {"nr": 2, "akt": 1, "navn": "Fosfoglukose-isomerase", "type": "Isomerisering", "sub": ["g6p"], "prod": ["f6p"],
     "inn": [], "ut": [], "dG0": 1.7,
     "rx": {"sub": [["g6p", 1]], "prod": [["f6p", 1]]},
     "hl": {"g6p": ["C1", "O1", "C2", "O2"], "f6p": ["C1", "O1", "C2", "O2"]},
     "hva": "Glukose-6-fosfat (et aldehyd) omdannes til fruktose-6-fosfat (et keton). Fosfatet blir sittende på C6.",
     "hvorfor": [
         "Dette er en forberedelse. Karbonylgruppen flytter fra C1 til C2, og C1 blir en vanlig CH₂OH-gruppe som kan få sitt eget fosfat i neste trinn.",
         "Ketonet på C2 gjør også at kjeden senere kan spaltes midt i molekylet, mellom C3 og C4. Med aldehydet på C1 ville vi fått to biter av ulik størrelse.",
         "ΔG er nær null: trinnet er reversibelt, og retningen avgjøres av konsentrasjonene."],
     "videre": "",
     "mer": "Mekanismen går via en enediol-mellomform, der protoner flyttes mellom C1 og C2.",
     "reg": ""},
    {"nr": 3, "akt": 1, "navn": "Fosfofruktokinase-1 (PFK-1)", "type": "Fosforylering", "sub": ["f6p"], "prod": ["fbp"],
     "inn": ["ATP"], "ut": ["ADP"], "dG0": -14.2,
     "rx": {"sub": [["f6p", 1], ["atp", 1]], "prod": [["fbp", 1], ["adp", 1]]},
     "hl": {"f6p": ["O1"], "fbp": ["O1", "P1"]},
     "hva": "Fruktose-6-fosfat får et nytt fosfat, nå på C1, fra ATP. Det gir fruktose-1,6-bisfosfat.",
     "hvorfor": [
         "Nå har molekylet to fosfatgrupper, og når det deles, får begge halvdelene hver sin. Fosfatene holder dem inne i cellen (ladning) og gir dem evne til å lage ATP.",
         "Dette er det første trinnet som forplikter cellen til glykolysen. Glukose-6-fosfat kan ellers brukes til andre ting, men fruktose-1,6-bisfosfat går bare videre i glykolysen.",
         "Reaksjonen er svært gunstig og går bare én vei. PFK-1 er derfor hovedbryteren for hele glykolysen. Etter dette trinnet har cellen investert 2 ATP."],
     "videre": "",
     "mer": "I leveren er fruktose-2,6-bisfosfat en kraftig aktivator av PFK-1. Det er et eksempel på hormonstyrt regulering (insulin og glukagon).",
     "reg": "Hemmes av ATP og sitrat. Aktiveres av AMP, ADP og fruktose-2,6-bisfosfat."},
    {"nr": 4, "akt": 2, "navn": "Aldolase", "type": "Spalting", "sub": ["fbp"], "prod": ["dhap", "g3p"],
     "inn": [], "ut": [], "dG0": 23.8,
     "rx": {"sub": [["fbp", 1]], "prod": [["dhap", 1], ["g3p", 1]]},
     "hl": {"fbp": ["C3", "C4"], "dhap": ["C3", "O3"], "g3p": ["C1", "O1"]},
     "hva": "Seks-karbonkjeden deles mellom C3 og C4. Resultatet er to tre-karbonmolekyler: dihydroksyacetonfosfat (DHAP) og glyseraldehyd-3-fosfat (G3P).",
     "hvorfor": [
         "Delingen gjør én glukose om til to små biter som kan behandles likt resten av veien. Fra nå av må alt regnes dobbelt.",
         "Under standardbetingelser er trinnet svært ugunstig (ΔG°′ = +24 kJ/mol), og skulle ikke gått. Men i cellen brukes DHAP og G3P opp med en gang, så konsentrasjonene er lave. Da blir Q liten, og ΔG = ΔG°′ + RT ln Q blir nær 0.",
         "Det er et eksempel på at reaksjonene nedstrøms kan «trekke» en ugunstig reaksjon fremover."],
     "videre": "",
     "mer": "I dyr danner aldolase en Schiff-base mellom en lysin-sidekjede (kapittel 3) og ketonet på C2. Den positivt ladde iminiumgruppen trekker til seg elektroner og stabiliserer mellomproduktet.",
     "reg": ""},
    {"nr": 5, "akt": 2, "navn": "Triosefosfat-isomerase (TIM)", "type": "Isomerisering", "sub": ["dhap"], "prod": ["g3p"],
     "inn": [], "ut": [], "dG0": 7.5,
     "rx": {"sub": [["dhap", 1]], "prod": [["g3p", 1]]},
     "hl": {"dhap": ["C2", "O2", "C3", "O3"], "g3p": ["C1", "O1", "C2", "O2"]},
     "hva": "DHAP (et keton) omdannes til G3P (et aldehyd).",
     "hvorfor": [
         "Bare G3P fortsetter i glykolysen. Uten dette trinnet ville halve glukosen blitt stående igjen som DHAP.",
         "Likevekten ligger mot DHAP (ΔG°′ = +7,5 kJ/mol), men G3P brukes opp i neste trinn og trekker reaksjonen fremover. Etter dette trinnet har vi to G3P per glukose."],
     "videre": "DHAP kan også gå en annen vei: som glyserol-3-fosfat er det utgangspunkt for triglyserider og fosfolipider.",
     "mer": "TIM regnes som et «kinetisk perfekt» enzym: hastigheten begrenses bare av hvor fort substratet diffunderer til enzymet.",
     "reg": ""},
    {"nr": 6, "akt": 3, "navn": "Glyseraldehyd-3-fosfat-dehydrogenase (GAPDH)", "type": "Oksidasjon og fosforylering",
     "sub": ["g3p"], "prod": ["bpg13"], "inn": ["NAD⁺", "Pᵢ"], "ut": ["NADH + H⁺"], "dG0": 6.3, "dobbel": True,
     "rx": {"sub": [["g3p", 1], ["pi", 1], ["nadox", 1]], "prod": [["bpg13", 1], ["nadred", 1]]},
     "hl": {"g3p": ["C1", "O1"], "bpg13": ["C1", "O1", "O1b", "P1b"]},
     "hva": "Aldehydet i G3P oksideres, og NAD⁺ tar imot elektronene og blir NADH. Samtidig settes uorganisk fosfat (Pᵢ) på: det dannes 1,3-bisfosfoglyserat.",
     "hvorfor": [
         "Her høstes energi fra glukosen for første gang. Når aldehydkarbonet får en ny binding til oksygen (det oksideres), flyttes elektrontetthet bort fra karbonet. Energien fanges i to former: som NADH (elektroner) og som et acylfosfat som lett gir fra seg fosfatet.",
         "Fosfatet kommer fra fritt Pᵢ, ikke fra ATP. Oksidasjonen «betaler» for fosforyleringen: de to reaksjonene er koblet.",
         "NADH må levere elektronene sine videre, ellers går cellen tom for NAD⁺ og glykolysen stopper (se «Hva skjer videre?»)."],
     "videre": "Elektronene i NADH går til elektrontransportkjeden (hvis oksygen finnes) eller brukes til å lage laktat.",
     "mer": "Enzymet bruker en cystein i det aktive setet (kapittel 3) som angriper aldehydet. Hydridet overføres til NAD⁺, og så angriper Pᵢ den dannede tioesteren.",
     "reg": ""},
    {"nr": 7, "akt": 3, "navn": "Fosfoglyseratkinase", "type": "Fosforylering på substratnivå", "sub": ["bpg13"], "prod": ["pg3"],
     "inn": ["ADP"], "ut": ["ATP"], "dG0": -18.8, "dobbel": True,
     "rx": {"sub": [["bpg13", 1], ["adp", 1]], "prod": [["pg3", 1], ["atp", 1]]},
     "hl": {"bpg13": ["O1b", "P1b"], "pg3": ["C1", "O1", "O1b"]},
     "hva": "Acylfosfatet gir fosfatgruppen sin til ADP. Det blir ATP og 3-fosfoglyserat.",
     "hvorfor": [
         "1,3-bisfosfoglyserat er en bedre fosfatgiver enn ATP (se fosfatstigen). Da er det gunstig at fosfatet går fra 1,3-bisfosfoglyserat til ADP.",
         "Dette kalles fosforylering på substratnivå: ATP lages direkte fra en energirik forbindelse, uten at elektrontransportkjeden er involvert.",
         "Siden vi har to G3P per glukose, får vi 2 ATP her. Da er de 2 ATP som ble investert, betalt tilbake. Fra nå er alt ATP gevinst."],
     "videre": "",
     "mer": "Reaksjonen er nær likevekt i cellen, og enzymet går begge veier (det trengs i glukoneogenesen).",
     "reg": ""},
    {"nr": 8, "akt": 3, "navn": "Fosfoglyseratmutase", "type": "Omplassering av fosfat", "sub": ["pg3"], "prod": ["pg2"],
     "inn": [], "ut": [], "dG0": 4.4, "dobbel": True,
     "rx": {"sub": [["pg3", 1]], "prod": [["pg2", 1]]},
     "hl": {"pg3": ["O3", "P3"], "pg2": ["O2", "P2"]},
     "hva": "Fosfatgruppen flytter fra C3 til C2: 3-fosfoglyserat blir til 2-fosfoglyserat.",
     "hvorfor": [
         "Forberedelse til neste trinn. For å lage en energirik fosfatgiver må fosfatet sitte på C2, ved siden av en OH-gruppe på C3 som kan fjernes som vann.",
         "ΔG er nær null: trinnet er reversibelt."],
     "videre": "3-fosfoglyserat kan også tas ut av glykolysen og brukes til å lage aminosyren serin.",
     "mer": "Fosfatet flyttes via en fosforylert histidin i enzymet (kapittel 3).",
     "reg": ""},
    {"nr": 9, "akt": 3, "navn": "Enolase", "type": "Dehydrering", "sub": ["pg2"], "prod": ["pep"],
     "inn": [], "ut": ["H₂O"], "dG0": 1.8, "dobbel": True,
     "rx": {"sub": [["pg2", 1]], "prod": [["pep", 1]]},
     "hl": {"pg2": ["C2", "C3", "O3"], "pep": ["C2", "C3"]},
     "hva": "Vann fjernes fra 2-fosfoglyserat, og det dannes en dobbeltbinding: fosfoenolpyruvat (PEP).",
     "hvorfor": [
         "Fjerningen av vann «låser» fosfatet på en enolform. PEP blir en av de mest energirike fosfatgiverne i cellen (ΔG°′ for hydrolyse = −62 kJ/mol).",
         "Energien sitter ikke bare «i bindingen». Den ligger i at enolpyruvat (det som blir igjen når fosfatet er borte) straks går over til den mye mer stabile ketoformen, pyruvat. Forskjellen i stabilitet gir drivkraften.",
         "Selve enolasetrinnet er nær likevekt og reversibelt."],
     "videre": "",
     "mer": "Enolase trenger Mg²⁺, som stabiliserer ladningen i mellomproduktet. Fluorid hemmer enolase.",
     "reg": ""},
    {"nr": 10, "akt": 3, "navn": "Pyruvatkinase", "type": "Fosforylering på substratnivå", "sub": ["pep"], "prod": ["pyr"],
     "inn": ["ADP"], "ut": ["ATP"], "dG0": -31.4, "dobbel": True,
     "rx": {"sub": [["pep", 1], ["adp", 1]], "prod": [["pyr", 1], ["atp", 1]]},
     "hl": {"pep": ["O2", "P2", "C3"], "pyr": ["C2", "O2", "C3"]},
     "hva": "PEP gir fosfatgruppen til ADP. Det blir ATP og pyruvat.",
     "hvorfor": [
         "PEP er en mye bedre fosfatgiver enn ATP, så reaksjonen er svært gunstig og går i praksis bare én vei. Pyruvatkinase er derfor den tredje reguleringsbryteren.",
         "Her lages 2 nye ATP (én per G3P). Totalt: 4 ATP laget − 2 ATP brukt = netto 2 ATP, og i tillegg 2 NADH.",
         "Pyruvat har fremdeles mesteparten av karbonenergien. Det er der veien fortsetter."],
     "videre": "Med oksygen går pyruvat til mitokondriene og blir acetyl-CoA. Uten oksygen blir det til laktat.",
     "mer": "I leveren hemmes pyruvatkinase også ved fosforylering (glukagon), slik at leveren ikke bruker glukose når blodsukkeret er lavt.",
     "reg": "Aktiveres av fruktose-1,6-bisfosfat (fremovermating). Hemmes av ATP og alanin."},
]

KONC = {  # typiske konsentrasjoner i røde blodceller (mM); nadox/nadred bare som forhold
    "glc": 5.0, "g6p": 0.083, "f6p": 0.014, "fbp": 0.031, "dhap": 0.14, "g3p": 0.019,
    "bpg13": 0.001, "pg3": 0.12, "pg2": 0.03, "pep": 0.023, "pyr": 0.051,
    "atp": 1.85, "adp": 0.14, "pi": 1.0, "nadox": 700.0, "nadred": 1.0,
}

FOSFATSTIGE = [  # ΔG°′ for hydrolyse av fosfatgruppen (kJ/mol)
    {"id": "pep", "navn": "Fosfoenolpyruvat (PEP)", "g": -61.9},
    {"id": "bpg13", "navn": "1,3-Bisfosfoglyserat", "g": -49.4},
    {"id": "krea", "navn": "Fosfokreatin", "g": -43.1},
    {"id": "atp", "navn": "ATP (→ ADP + Pᵢ)", "g": -30.5},
    {"id": "g1p", "navn": "Glukose-1-fosfat", "g": -20.9},
    {"id": "f6p", "navn": "Fruktose-6-fosfat", "g": -15.9},
    {"id": "g6p", "navn": "Glukose-6-fosfat", "g": -13.8},
    {"id": "gly3p", "navn": "Glyserol-3-fosfat", "g": -9.2},
]

FUNKSJONER = r"""
BKAA.gl.RT = 0.008314 * 310.15; // kJ/mol ved 37 °C
// Q for et trinn: konsentrasjoner i mM (konverteres til M), nadox/nadred bare som forhold
BKAA.gl.Q = function (s, c) {
  const v = (id) => (id === "nadox" || id === "nadred" ? c[id] : c[id] * 1e-3);
  let q = 1;
  s.rx.prod.forEach(([id, n]) => (q *= Math.pow(v(id), n)));
  s.rx.sub.forEach(([id, n]) => (q /= Math.pow(v(id), n)));
  return q;
};
BKAA.gl.dG = function (s, c) {
  const rtlnq = BKAA.gl.RT * Math.log(BKAA.gl.Q(s, c));
  return { dG0: s.dG0, rtlnq, dG: s.dG0 + rtlnq };
};
"""


def main():
    data = {
        "molekyler": {m["id"]: m for m in M},
        "akter": AKTER,
        "steg": STEG,
        "konc": KONC,
        "fosfatstige": FOSFATSTIGE,
    }
    js = ("// Generert av tools/lag_glykolysedata.py. Ikke rediger for hånd.\n"
          "window.BKAA = window.BKAA || {};\n"
          "BKAA.gl = " + json.dumps(data, ensure_ascii=False, separators=(",", ":")) + ";\n" + FUNKSJONER)
    UT.write_text(js, encoding="utf-8")
    print(f"Skrev {UT} ({len(js) // 1024} kB).")


if __name__ == "__main__":
    main()
