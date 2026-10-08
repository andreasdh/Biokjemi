// Delt hjelpekode for kapittel 3 og aminosyrespillet:
//   BKAA.tegn(molekyl, valg)  tegner en skjelettformel som SVG-tekst
//   BKAA.kjemi                ladning, pI, titreringskurve og peptidberegninger
// Krever at aminosyredata.js er lastet først. Ren JavaScript uten avhengigheter.
(function () {
  const BKAA = (window.BKAA = window.BKAA || {});
  const D = BKAA.data;
  BKAA.alle = D.aminosyrer.concat(D.ekstra);
  BKAA.finn = (id) => BKAA.alle.find((a) => a.id === id);

  // ---------- farger ----------
  BKAA.gruppeFarge = { upolar: "#e0b43a", polar: "#4aa57f", positiv: "#3d7ea6", negativ: "#e0584d", ekstra: "#8aa0ad" };
  BKAA.gruppeNavn = { upolar: "upolar", polar: "polar", positiv: "positivt ladd", negativ: "negativt ladd", ekstra: "annen" };
  const ELFARGE = { O: "#e0584d", N: "#4a90c2", S: "#c9962f", Se: "#c9962f", R: "#d9534f" };
  const SUB = "₀₁₂₃₄₅₆₇₈₉";
  const sub = (n) => String(n).split("").map((c) => SUB[+c]).join("");
  const SUP = { 1: "", 2: "²", 3: "³" };

  // ---------- tegning av skjelettformler ----------
  // valg: { bredde, maksBinding, tilstand:{a,n,r} (true = deprotonert), uthev: [{atomer:[...], farge}],
  //         sidegruppe: "farge" (uthev sidekjeden), fs: skriftstørrelse, vis: {cip:true} }
  BKAA.tegn = function (mol, valg = {}) {
    const W = valg.bredde || 320;
    const tilst = valg.tilstand || {};
    const at = mol.atomer.map((a) => {
      const b = { ...a };
      if (a.site && tilst[a.site]) {
        b.h = a.h - 1;
        b.q = a.q - 1;
      }
      return b;
    });
    const naboer = at.map(() => []);
    mol.bindinger.forEach(([a, b, o, w]) => {
      naboer[a].push(b);
      naboer[b].push(a);
    });
    const xs = at.map((a) => a.x),
      ys = at.map((a) => -a.y);
    const minx = Math.min(...xs) - 0.55,
      maxx = Math.max(...xs) + 0.55;
    const miny = Math.min(...ys) - 0.5,
      maxy = Math.max(...ys) + 0.5;
    const pad = 10;
    const maksBinding = valg.maksBinding || 44;
    const s = Math.min((W - 2 * pad) / (maxx - minx), maksBinding / 1.5);
    const H = Math.max(valg.minHoyde || 120, Math.round((maxy - miny) * s + 2 * pad));
    const fs = valg.fs || Math.max(12, Math.min(17, s * 0.5));
    const cw = fs * 0.6;
    const px = (i) => W / 2 + (at[i].x - (minx + maxx) / 2) * s;
    const py = (i) => H / 2 + (-at[i].y - (miny + maxy) / 2) * s;

    // hvilke atomer har etikett, og hvilken side står H-ene på?
    const info = at.map((a, i) => {
      const gradMinst1 = naboer[i].length > 0;
      const metylTerminal = a.el === "C" && naboer[i].length === 1 && a.h > 0;
      const har = a.el !== "C" || a.q !== 0 || naboer[i].length === 0 || metylTerminal;
      if (!har) return { har: false };
      let dx = 0;
      naboer[i].forEach((j) => (dx += px(j) - px(i)));
      return { har: true, venstre: gradMinst1 && dx > 1 }; // naboene til høyre => H-ene til venstre
    });

    const stroke = `stroke="currentColor" stroke-width="${(1.6).toFixed(2)}" stroke-linecap="round"`;
    const korte = (i) => (info[i].har ? fs * (at[i].el.length > 1 ? 0.85 : 0.62) : 0);
    let ut = "";

    // uthevinger (myke flekker bak)
    const utheving = (valg.uthev || []).slice();
    if (valg.sidegruppe) {
      const sc = at.map((a, i) => (a.sc ? i : -1)).filter((i) => i >= 0);
      utheving.push({ atomer: sc, farge: valg.sidegruppe });
    }
    utheving.forEach((u) => {
      const sett = new Set(u.atomer);
      const r = s * 0.42;
      let lag = "";
      mol.bindinger.forEach(([a, b]) => {
        if (sett.has(a) && sett.has(b))
          lag += `<line x1="${px(a)}" y1="${py(a)}" x2="${px(b)}" y2="${py(b)}" stroke="${u.farge}" stroke-width="${r * 2}" stroke-linecap="round"/>`;
      });
      u.atomer.forEach((i) => {
        lag += `<circle cx="${px(i)}" cy="${py(i)}" r="${r}" fill="${u.farge}"/>`;
      });
      ut += `<g opacity="0.28">${lag}</g>`;
    });

    // bindinger
    mol.bindinger.forEach(([a, b, orden, kile]) => {
      const x1 = px(a), y1 = py(a), x2 = px(b), y2 = py(b);
      const dx = x2 - x1, dy = y2 - y1, len = Math.hypot(dx, dy) || 1;
      const ux = dx / len, uy = dy / len, nx = -uy, ny = ux;
      const k1 = korte(a), k2 = korte(b);
      const ax = x1 + ux * k1, ay = y1 + uy * k1, bx = x2 - ux * k2, by = y2 - uy * k2;
      if (kile === 1) {
        const w = s * 0.2;
        ut += `<polygon points="${ax},${ay} ${bx + nx * w},${by + ny * w} ${bx - nx * w},${by - ny * w}" fill="currentColor"/>`;
        return;
      }
      if (kile === -1) {
        let d = "";
        for (let t = 0.12; t <= 1; t += 0.14) {
          const cx = ax + (bx - ax) * t, cy = ay + (by - ay) * t, w = s * 0.2 * t;
          d += `M${cx + nx * w} ${cy + ny * w}L${cx - nx * w} ${cy - ny * w}`;
        }
        ut += `<path d="${d}" ${stroke}/>`;
        return;
      }
      if (orden === 1) {
        ut += `<line x1="${ax}" y1="${ay}" x2="${bx}" y2="${by}" ${stroke}/>`;
      } else if (orden === 2) {
        // ring-dobbeltbinding: bare én side har andre naboer => forskyv mot den siden
        let sum = 0, antall = 0;
        [a, b].forEach((p) =>
          naboer[p].forEach((j) => {
            if (j !== a && j !== b) {
              const side = (px(j) - px(p)) * nx + (py(j) - py(p)) * ny;
              sum += Math.sign(side);
              antall++;
            }
          }),
        );
        const d = 2.8;
        if (antall > 0 && Math.abs(sum) === antall) {
          const sg = Math.sign(sum), inn = Math.min(len * 0.18, 7);
          ut += `<line x1="${ax}" y1="${ay}" x2="${bx}" y2="${by}" ${stroke}/>`;
          ut += `<line x1="${ax + ux * inn + nx * sg * d * 1.7}" y1="${ay + uy * inn + ny * sg * d * 1.7}" x2="${bx - ux * inn + nx * sg * d * 1.7}" y2="${by - uy * inn + ny * sg * d * 1.7}" ${stroke}/>`;
        } else {
          [-1, 1].forEach((sg) => {
            ut += `<line x1="${ax + nx * sg * d}" y1="${ay + ny * sg * d}" x2="${bx + nx * sg * d}" y2="${by + ny * sg * d}" ${stroke}/>`;
          });
        }
      } else {
        [-1, 0, 1].forEach((sg) => {
          ut += `<line x1="${ax + nx * sg * 3}" y1="${ay + ny * sg * 3}" x2="${bx + nx * sg * 3}" y2="${by + ny * sg * 3}" ${stroke}/>`;
        });
      }
    });

    // atometiketter
    at.forEach((a, i) => {
      if (!info[i].har) return;
      const x = px(i), y = py(i);
      const farge = ELFARGE[a.el] || "currentColor";
      const h = a.h > 0 ? "H" + (a.h > 1 ? sub(a.h) : "") : "";
      const venstre = info[i].venstre;
      const tekst = venstre ? h + a.el : a.el + h;
      const bredde = tekst.replace(/[₀-₉]/g, "").length * cw * 0.95 + (tekst.match(/[₀-₉]/g) || []).length * cw * 0.6;
      const ankring = venstre ? "end" : "start";
      const tx = venstre ? x + cw * (a.el.length > 1 ? 0.9 : 0.5) : x - cw * (a.el.length > 1 ? 0.9 : 0.5);
      ut += `<text x="${tx}" y="${y}" dy="0.35em" text-anchor="${ankring}" font-size="${fs}" font-weight="700" fill="${farge}">${tekst}</text>`;
      if (a.q !== 0) {
        const r = fs * 0.34;
        const cx = venstre ? x + cw * 0.35 : tx + bredde + r + 1.5;
        const cy = venstre ? y - fs * 0.95 : y - fs * 0.62;
        const lbl = a.q > 0 ? "+" : "−";
        ut += `<circle cx="${cx}" cy="${cy}" r="${r}" fill="none" stroke="${farge}" stroke-width="1.1"/>`;
        ut += `<path d="M${cx - r * 0.55} ${cy}h${r * 1.1}${a.q > 0 ? `M${cx} ${cy - r * 0.55}v${r * 1.1}` : ""}" stroke="${farge}" stroke-width="1.1" fill="none"/>`;
      }
    });

    (valg.merker || []).forEach((m) => {
      ut += `<text x="${px(m.atom) + (m.dx || 0)}" y="${py(m.atom) + (m.dy || 0)}" font-size="${fs * 0.85}" fill="currentColor" opacity="0.8" text-anchor="${m.anker || "middle"}">${m.tekst}</text>`;
    });
    const cip = valg.cip && mol.cip ? `<text x="${W - 8}" y="${H - 8}" text-anchor="end" font-size="${fs * 0.95}" fill="currentColor" opacity="0.7">(${mol.cip})</text>` : "";
    const svg = `<svg viewBox="0 0 ${W} ${H}" class="bk-aa-svg" role="img" aria-label="${valg.alt || "Strukturformel for " + (mol.navn || "aminosyre")}" style="width:100%;max-width:${W}px;height:auto;display:block;margin:0 auto" font-family="inherit">${ut}${cip}</svg>`;
    return { svg, H, W, px, py, s };
  };

  // ---------- kjemi: ladning, pI og titrering ----------
  const K = (BKAA.kjemi = {});
  // grupper (fullt protonert form): α-COOH (q0=0), α-NH3+ (q0=+1), sidekjede (q0 = 0 eller +1)
  K.grupper = function (aa) {
    const [pa, pn, pr] = aa.pka;
    const g = [
      { s: "a", pK: pa, q0: 0, navn: "α-COOH" },
      { s: "n", pK: pn, q0: 1, navn: "α-NH₃⁺" },
    ];
    if (pr !== null && pr !== undefined) {
      const q0 = ["lys", "arg", "his"].includes(aa.id) ? 1 : 0;
      g.push({ s: "r", pK: pr, q0, navn: "sidegruppe" });
    }
    return g;
  };
  const andelDeprot = (pK, pH) => 1 / (1 + Math.pow(10, pK - pH));
  K.andelDeprot = andelDeprot;
  K.ladning = (aa, pH) => K.grupper(aa).reduce((q, g) => q + g.q0 - andelDeprot(g.pK, pH), 0);
  K.tilstand = (aa, pH) => {
    const t = {};
    K.grupper(aa).forEach((g) => (t[g.s] = pH > g.pK));
    return t;
  };
  K.heltallsladning = (aa, pH) => K.grupper(aa).reduce((q, g) => q + g.q0 - (pH > g.pK ? 1 : 0), 0);
  K.bisect = (f, lo = -2, hi = 16) => {
    // finner pH der f(pH) = 0 for en avtakende funksjon f
    for (let i = 0; i < 60; i++) {
      const m = (lo + hi) / 2;
      if (f(m) > 0) lo = m;
      else hi = m;
    }
    return (lo + hi) / 2;
  };
  // pI = gjennomsnittet av de to pKₐ-verdiene som ligger på hver side av den nøytrale formen
  K.pI = (aa) => {
    const g = K.grupper(aa).sort((a, b) => a.pK - b.pK);
    const q0 = g.reduce((s, x) => s + x.q0, 0);
    return (g[q0 - 1].pK + g[q0].pK) / 2;
  };
  // pH som funksjon av antall ekvivalenter OH⁻ (0 = fullt protonert)
  K.titrering = (aa) => {
    const g = K.grupper(aa);
    const pts = [];
    for (let pH = 0; pH <= 14.001; pH += 0.05) {
      const x = g.reduce((s, gr) => s + andelDeprot(gr.pK, pH), 0);
      pts.push([x, pH]);
    }
    return pts;
  };

  // ---------- peptider: typiske pKa-verdier i proteiner ----------
  K.PKA_PEPTID = { nTerm: 8.0, cTerm: 3.1, D: 4.1, E: 4.1, H: 6.0, C: 8.3, Y: 10.9, K: 10.8, R: 12.5 };
  K.peptidGrupper = (seq) => {
    const P = K.PKA_PEPTID;
    const g = [{ pK: P.nTerm, q0: 1 }, { pK: P.cTerm, q0: 0 }];
    seq.split("").forEach((c) => {
      if (P[c] !== undefined && c.length === 1) g.push({ pK: P[c], q0: "KRH".includes(c) ? 1 : 0 });
    });
    return g;
  };
  K.peptidLadning = (seq, pH) => K.peptidGrupper(seq).reduce((q, g) => q + g.q0 - andelDeprot(g.pK, pH), 0);
  K.peptidPI = (seq) => K.bisect((pH) => K.peptidLadning(seq, pH));
  K.peptidMasse = (seq) => {
    const vann = 18.015;
    return seq.split("").reduce((m, c) => {
      const aa = D.aminosyrer.find((a) => a.en === c);
      return m + (aa ? aa.mw - vann : 0);
    }, vann);
  };
  K.fraEn = (c) => D.aminosyrer.find((a) => a.en === c);

  // ---------- tilfeldig tall med frø (for spill-utfordringer) ----------
  BKAA.rng = function (frø) {
    let t = [...String(frø)].reduce((h, c) => (Math.imul(h ^ c.charCodeAt(0), 2654435761) >>> 0), 2166136261) || 1;
    return function () {
      t += 0x6d2b79f5;
      let r = Math.imul(t ^ (t >>> 15), 1 | t);
      r ^= r + Math.imul(r ^ (r >>> 7), 61 | r);
      return ((r ^ (r >>> 14)) >>> 0) / 4294967296;
    };
  };
})();
