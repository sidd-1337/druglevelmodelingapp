"use client";

import React, { useEffect, useMemo, useState } from "react";
import { useSearchParams } from "next/navigation";
import {
  Perpetrator,
  PresetSummary,
  RegimenDrugIn,
  RegimenResult,
  SubstanceIn,
  multiDrugApi,
  perpetratorsApi,
  presetsApi,
} from "@/lib/api";
import { num } from "@/lib/format";
import { cz } from "@/lib/typo";
import OverlayChart from "@/components/OverlayChart";
import EnzymeChart from "@/components/EnzymeChart";
import SectionHead from "@/components/ui/SectionHead";
import Slider from "@/components/ui/Slider";
import Status from "@/components/ui/Status";

const NBSP = " ";
const MAX_DRUGS = 6; // one colour each in the overlay
const MAX_ENTRIES = 8; // the API's limit for drugs and substances together
const SUBSTANCE = "substance:";

const EFFECT: Record<string, string> = {
  serotonergic: "Serotonergní účinek",
  cns_depressant: "Útlum CNS",
  qt_prolonging: "Prodloužení QT",
  sympathomimetic: "Sympatomimetický účinek",
};

type Phase = "loading" | "ready" | "offline";

/** Preset ids are safe to join with commas; substance names are not ("ACE inhibitor (lisinopril, enalapril)"). */
function idList(param: string | null): string[] {
  return param ? param.split(",").map((s) => s.trim()).filter(Boolean) : [];
}

export default function CombinationTool() {
  const params = useSearchParams();
  // Read the link once: replaceState below also updates useSearchParams, and re-reading
  // it would reset every dose the user has set.
  const [initial] = useState(() => ({ lek: params.get("lek"), latka: params.getAll("latka") }));
  const [presets, setPresets] = useState<PresetSummary[] | null>(null);
  const [perpetrators, setPerpetrators] = useState<Perpetrator[]>([]);
  const [phase, setPhase] = useState<Phase>("loading");
  const [drugs, setDrugs] = useState<RegimenDrugIn[]>([]);
  const [substances, setSubstances] = useState<SubstanceIn[]>([]);
  const [days, setDays] = useState(10);
  const [weightKg, setWeightKg] = useState(70);
  const [result, setResult] = useState<RegimenResult | null>(null);
  const [selected, setSelected] = useState<number | null>(null);

  // Catalogues, then the regimen in the link (?lek=lithium,mdma&latka=Paroxetine 20 mg/day&latka=...).
  useEffect(() => {
    let active = true;
    Promise.all([presetsApi(), perpetratorsApi()]).then(([ps, perps]) => {
      if (!active) return;
      if (!ps) {
        setPhase("offline");
        return;
      }
      setPresets(ps);
      setPerpetrators(perps ?? []);
      const byId = new Map(ps.map((p) => [p.id, p]));
      setDrugs(
        idList(initial.lek)
          .filter((id) => byId.has(id))
          .slice(0, MAX_DRUGS)
          .map((id) => ({ drugId: id, doseMg: byId.get(id)!.defaultDoseMg, intervalH: byId.get(id)!.defaultIntervalH }))
      );
      const known = new Set((perps ?? []).map((p) => p.name));
      setSubstances(
        initial.latka
          .filter((name) => known.has(name))
          .map((name) => ({ name, startH: null, stopH: null }))
      );
      setPhase("ready");
    });
    return () => {
      active = false;
    };
  }, [initial]);

  // Keep the link shareable without navigating.
  useEffect(() => {
    if (phase !== "ready") return;
    const q = new URLSearchParams();
    if (drugs.length) q.set("lek", drugs.map((d) => d.drugId).join(","));
    substances.forEach((s) => q.append("latka", s.name));
    const query = q.toString();
    window.history.replaceState(null, "", query ? `?${query}` : window.location.pathname);
  }, [phase, drugs, substances]);

  // Simulate, debounced so dragging a slider does not fire a request per step.
  useEffect(() => {
    if (phase === "offline" || drugs.length === 0) {
      setResult(null);
      return;
    }
    let active = true;
    const timer = setTimeout(() => {
      // A stop on or before the start is shown as a warning and treated as "throughout".
      const valid = substances.map((x) =>
        x.stopH !== null && x.startH !== null && x.stopH <= x.startH ? { ...x, stopH: null } : x
      );
      multiDrugApi({ drugs, substances: valid, durationH: days * 24, bodyWeightKg: weightKg }).then((r) => {
        if (!active) return;
        if (r) setResult(r);
        else setPhase("offline");
      });
    }, 250);
    return () => {
      active = false;
      clearTimeout(timer);
    };
  }, [phase, drugs, substances, days, weightKg]);

  const presetById = useMemo(() => new Map((presets ?? []).map((p) => [p.id, p])), [presets]);
  const perpByName = useMemo(() => new Map(perpetrators.map((p) => [p.name, p])), [perpetrators]);
  const shown = (name: string) => perpByName.get(name)?.nameCs ?? name;
  const participant = (key: string) =>
    key.startsWith(SUBSTANCE) ? shown(key.slice(SUBSTANCE.length)) : presetById.get(key)?.name ?? key;

  const entries = drugs.length + substances.length;
  const canAddDrug = drugs.length < MAX_DRUGS && entries < MAX_ENTRIES;
  const canAddSubstance = entries < MAX_ENTRIES;
  const addablePresets = (presets ?? []).filter((p) => p.id !== "custom" && !drugs.some((d) => d.drugId === p.id));
  const addablePerps = perpetrators
    .filter((p) => !substances.some((s) => s.name === p.name))
    .sort((a, b) => a.nameCs.localeCompare(b.nameCs, "cs"));

  const addDrug = (id: string) => {
    const p = presetById.get(id);
    if (p) setDrugs((ds) => [...ds, { drugId: id, doseMg: p.defaultDoseMg, intervalH: p.defaultIntervalH }]);
  };
  const updateDrug = (id: string, change: Partial<RegimenDrugIn>) =>
    setDrugs((ds) => ds.map((d) => (d.drugId === id ? { ...d, ...change } : d)));
  const updateSubstance = (name: string, change: Partial<SubstanceIn>) =>
    setSubstances((ss) => ss.map((s) => (s.name === name ? { ...s, ...change } : s)));

  if (phase === "offline") {
    return (
      <p className="max-w-prose text-body">
        <Status tone="warn">
          {cz(
            "Výpočet kombinací teď není dostupný. Sdílené enzymy a časový průběh interakcí počítá jen simulační služba, která neodpovídá. Zkuste to prosím později."
          )}
        </Status>
      </p>
    );
  }

  return (
    <div className="space-y-20">
      {/* ------------------------------------------------------------ regimen */}
      <section aria-labelledby="rezim" className="grid gap-10 lg:grid-cols-[minmax(0,1fr)_19rem] lg:gap-14">
        <div className="space-y-10">
          <div>
            <h2 id="rezim" className="text-sub">Léky</h2>
            <ul className="mt-4 divide-y divide-line-soft border-y border-line-soft">
              {drugs.map((d, i) => {
                const p = presetById.get(d.drugId);
                if (!p) return null;
                const step = p.defaultDoseMg > 50 ? 5 : 0.5;
                return (
                  <li key={d.drugId} className="grid gap-4 py-5 sm:grid-cols-[12rem_minmax(0,1fr)_minmax(0,1fr)_auto] sm:items-end sm:gap-6">
                    <span className="flex items-center gap-3 text-ui font-medium">
                      <span aria-hidden className="h-[3px] w-5 rounded-full"
                            style={{ background: `rgb(var(--series-${(i % 6) + 1}))` }} />
                      {p.name}
                    </span>
                    <Slider label="Dávka" display={`${num(d.doseMg, Number.isInteger(d.doseMg) ? 0 : 1)}${NBSP}mg`}
                            value={d.doseMg} min={Math.max(step, Math.round(p.defaultDoseMg * 0.2))}
                            max={Math.round(p.defaultDoseMg * 3)} step={step}
                            onChange={(v) => updateDrug(d.drugId, { doseMg: v })} />
                    <Slider label="Interval" display={`každých${NBSP}${d.intervalH}${NBSP}h`} value={d.intervalH}
                            min={4} max={48} step={2} onChange={(v) => updateDrug(d.drugId, { intervalH: v })} />
                    <button type="button" className="link text-ui sm:pb-1"
                            onClick={() => setDrugs((ds) => ds.filter((x) => x.drugId !== d.drugId))}>
                      Odebrat<span className="sr-only"> {p.name}</span>
                    </button>
                  </li>
                );
              })}
            </ul>
            {canAddDrug ? (
              <label className="mt-4 flex flex-wrap items-center gap-3 text-ui">
                <span className="text-ink-2">Přidat lék</span>
                <select className="rounded-ctl border border-line bg-bg px-3 py-2 text-ui" value=""
                        onChange={(e) => e.target.value && addDrug(e.target.value)}>
                  <option value="">Vyberte…</option>
                  {addablePresets.map((p) => (
                    <option key={p.id} value={p.id}>{p.name}</option>
                  ))}
                </select>
              </label>
            ) : (
              <p className="mt-4 text-small text-ink-2">{cz(`Najednou lze porovnat nejvýše ${MAX_DRUGS} léků.`)}</p>
            )}
          </div>

          <div>
            <h2 className="text-sub">Další látky</h2>
            <p className="mt-2 max-w-prose text-ui text-ink-2">
              {cz("Látky, které ovlivňují hladinu jiných léků (například inhibitory nebo induktory enzymů). Samy se v grafu nekreslí.")}
            </p>
            <ul className="mt-4 divide-y divide-line-soft border-y border-line-soft empty:hidden">
              {substances.map((s) => {
                const p = perpByName.get(s.name);
                const already = s.startH === null;
                const stopBeforeStart = s.stopH !== null && s.startH !== null && s.stopH <= s.startH;
                return (
                  <li key={s.name} className="space-y-4 py-5">
                    <div className="flex flex-wrap items-baseline justify-between gap-x-6 gap-y-1">
                      <span className="text-ui font-medium">
                        {shown(s.name)}
                        {p?.enzyme && (
                          <span className="ml-2 font-normal text-ink-2">
                            {p.induces ? "indukuje" : "tlumí"} {p.enzyme}
                          </span>
                        )}
                      </span>
                      <button type="button" className="link text-ui"
                              onClick={() => setSubstances((ss) => ss.filter((x) => x.name !== s.name))}>
                        Odebrat<span className="sr-only"> {shown(s.name)}</span>
                      </button>
                    </div>
                    {p?.timed ? (
                      <div className="grid gap-4 sm:grid-cols-[auto_minmax(0,1fr)_minmax(0,1fr)] sm:items-end sm:gap-6">
                        <label className="flex items-center gap-2 text-ui">
                          <input type="checkbox" checked={already}
                                 onChange={(e) => updateSubstance(s.name, { startH: e.target.checked ? null : 0 })} />
                          Už užívá
                        </label>
                        {!already ? (
                          <Slider label="Začíná" display={`den${NBSP}${num((s.startH ?? 0) / 24)}`}
                                  value={(s.startH ?? 0) / 24} min={0} max={days} step={1}
                                  onChange={(v) => updateSubstance(s.name, { startH: v * 24 })} />
                        ) : (
                          <span className="text-small text-ink-2">{cz("Působí plně od začátku simulace.")}</span>
                        )}
                        <Slider label="Končí" display={s.stopH === null ? "neukončeno" : `den${NBSP}${num(s.stopH / 24)}`}
                                value={s.stopH === null ? 0 : s.stopH / 24} min={0} max={days} step={1}
                                onChange={(v) => updateSubstance(s.name, { stopH: v === 0 ? null : v * 24 })} />
                      </div>
                    ) : (
                      <p className="text-small text-ink-2">
                        {cz("Bez vlastních farmakokinetických údajů: působí stejně po celou dobu simulace.")}
                      </p>
                    )}
                    {stopBeforeStart && (
                      <p className="text-small">
                        <Status tone="warn">{cz("Ukončení je dřív než začátek, proto se počítá s užíváním po celou dobu.")}</Status>
                      </p>
                    )}
                  </li>
                );
              })}
            </ul>
            {canAddSubstance && (
              <label className="mt-4 flex flex-wrap items-center gap-3 text-ui">
                <span className="text-ink-2">Přidat látku</span>
                <select className="rounded-ctl border border-line bg-bg px-3 py-2 text-ui" value=""
                        onChange={(e) =>
                          e.target.value &&
                          setSubstances((ss) => [...ss, { name: e.target.value, startH: null, stopH: null }])
                        }>
                  <option value="">Vyberte…</option>
                  {addablePerps.map((p) => (
                    <option key={p.name} value={p.name}>{p.nameCs}</option>
                  ))}
                </select>
              </label>
            )}
          </div>
        </div>

        <div className="panel space-y-5 p-6 lg:sticky lg:top-24 lg:self-start">
          <Slider label="Délka simulace" display={`${days}${NBSP}dní`} value={days} min={2} max={28} step={1}
                  onChange={setDays} />
          <Slider label="Tělesná hmotnost" display={`${weightKg}${NBSP}kg`} value={weightKg} min={40} max={140}
                  step={1} onChange={setWeightKg} />
        </div>
      </section>

      {drugs.length === 0 && phase === "ready" && (
        <p className="max-w-prose text-body text-ink-2">{cz("Přidejte alespoň jeden lék a uvidíte jeho hladinu spolu s ostatními.")}</p>
      )}

      {result && result.series.length > 0 && (
        <>
          {/* ------------------------------------------------------------ overlay */}
          <section aria-labelledby="hladiny">
            <SectionHead id="hladiny" title="Všechny léky na jedné škále">
              Každá křivka ukazuje, kde se lék nachází ve svém vlastním terapeutickém okně: 0 je MEC, 1 je MTC. Jednotky i okna se mezi léky liší, proto je to jediná společná škála.
            </SectionHead>
            <div className="mt-8">
              <OverlayChart
                timeGrid={result.timeGrid}
                series={result.series.map((s) => ({
                  key: s.drugId, name: s.name, unit: s.unit, mec: s.mec, mtc: s.mtc, values: s.concentrations,
                }))}
                selected={selected}
                onSelect={setSelected}
              />
            </div>
          </section>

          {/* ------------------------------------------------------------ clashes */}
          <section aria-labelledby="stretry">
            <SectionHead id="stretry" title="Střety účinků">
              Některé kombinace jsou nebezpečné, i když se žádná hladina nezmění: účinky se sčítají.
            </SectionHead>
            <ul className="mt-6 max-w-prose space-y-3 text-body">
              {result.clashes.map((c) => {
                const who = c.participants.map(participant).join(", ");
                return (
                  <li key={c.effect}>
                    {c.overlapH > 0 ? (
                      <Status tone="bad">
                        <strong className="font-medium">{EFFECT[c.effect] ?? c.effect}:</strong>{" "}
                        {cz(`${who} působí současně ${num(c.overlapH, 0)} h (poprvé v ${num(c.firstH ?? 0, 1)} h). Jejich účinky se sčítají.`)}
                      </Status>
                    ) : (
                      <Status tone="warn">
                        <strong className="font-medium">{EFFECT[c.effect] ?? c.effect}:</strong>{" "}
                        {cz(`${who} sdílejí tento účinek, ale v této simulaci nepůsobí současně. Jiné rozvržení dávek je může svést dohromady.`)}
                      </Status>
                    )}
                  </li>
                );
              })}
              {result.clashes.length === 0 && (
                <li className="text-ink-2">{cz("Posouzené léky a látky nesdílejí žádný ze sledovaných účinků.")}</li>
              )}
            </ul>
            {result.pdUnassessed.length > 0 && (
              <p className="mt-4 max-w-prose text-small text-ink-2">
                {cz(`Neposouzeno: ${result.pdUnassessed.map(participant).join(", ")}. To neznamená, že je kombinace bezpečná.`)}
              </p>
            )}
          </section>

          {/* ------------------------------------------------------- interactions */}
          <section aria-labelledby="interakce">
            <SectionHead id="interakce" title="Vliv na hladiny">
              Jak se expozice každého léku mění oproti samostatnému užití, v posledním dávkovacím intervalu.
            </SectionHead>
            <dl className="mt-6 border-b border-line-soft">
              {result.series.map((s) => {
                const acting = s.applied.filter((a) => !a.skipped || a.skipped === "iv");
                const pctChange = Math.round(100 * (s.exposureFactor - 1));
                return (
                  <div key={s.drugId} className="grid gap-2 border-t border-line-soft py-5 sm:grid-cols-[15rem_minmax(0,1fr)] sm:gap-12">
                    <dt className="text-ui font-medium">{s.name}</dt>
                    <dd className="text-ui">
                      {acting.length === 0 && s.drivers.length === 0 ? (
                        <span className="text-ink-2">{cz("Žádná doložená interakce. To neznamená, že žádná neexistuje.")}</span>
                      ) : (
                        <>
                          <p>
                            <span className="figures font-semibold">
                              ×{num(s.exposureFactor, 2)}
                            </span>{" "}
                            <span className="text-ink-2">
                              ({pctChange > 0 ? "+" : ""}{pctChange}{NBSP}%)
                            </span>
                            {s.drivers.length > 0 && (
                              <span className="text-ink-2">
                                {" "}· {cz(`nastupuje a odeznívá s: ${s.drivers.map(shown).join(", ")}`)}
                              </span>
                            )}
                          </p>
                          <ul className="mt-1.5 space-y-0.5 text-small text-ink-2">
                            {acting.map((a) => (
                              <li key={`${a.perpetrator}-${a.site}`}>
                                {shown(a.perpetrator)}:{" "}
                                {a.skipped === "iv"
                                  ? "u nitrožilní dávky bez vlivu"
                                  : a.viaEnzyme
                                  ? "přes sdílený enzym (simulováno)"
                                  : `×${num(a.factor ?? 1, 2)} ${a.site === "clearance" ? "přes clearance" : "přes vstřebání"}`}
                                {a.pooled && ` · sloučeno v${NBSP}${a.pooled}`}
                              </li>
                            ))}
                          </ul>
                          {s.pooled.map((p) => (
                            <p key={p.enzyme} className="mt-1.5 max-w-prose text-small text-ink-2">
                              {cz(`${p.enzyme}: ${acting.filter((a) => a.pooled === p.enzyme).map((a) => shown(a.perpetrator)).join(", ")} působí na stejný enzym, proto se skládají v něm: ×${num(p.combined, 2)}, ne ×${num(p.multiplied, 2)} vynásobením faktorů.`)}
                            </p>
                          ))}
                        </>
                      )}
                    </dd>
                  </div>
                );
              })}
            </dl>
            {result.unmatched.length > 0 && (
              <p className="mt-4 max-w-prose text-small text-ink-2">
                {cz(`${result.unmatched.map(shown).join(", ")}: s vybranými léky nemá doloženou ani odvozenou interakci, takže jejich hladiny nemění.`)}
              </p>
            )}
          </section>

          {/* ------------------------------------------------------------ enzymes */}
          {Object.keys(result.enzymeActivity).length > 0 && (
            <section aria-labelledby="enzymy">
              <SectionHead id="enzymy" title="Aktivita enzymů">
                Zničený nebo indukovaný enzym se mění jen tvorbou či odbouráváním, tedy během dnů. Účinek reverzibilního inhibitoru sleduje jeho hladiny.
              </SectionHead>
              <div className="mt-8">
                <EnzymeChart
                  timeGrid={result.timeGrid}
                  activity={result.enzymeActivity}
                  drivers={Object.fromEntries(
                    Object.entries(result.enzymeDrivers).map(([e, names]) => [e, names.map(shown)])
                  )}
                  selected={selected}
                  onSelect={setSelected}
                />
              </div>
            </section>
          )}

          <p className="max-w-prose text-small text-ink-2">
            {cz("Modelový výpočet pro vzdělávání, ne doporučení k dávkování. Střety účinků vycházejí z SPC a jsou hlášeny, dokud obě léčiva působí (u léku nad MEC); skutečné riziko na hladině přímo nezávisí.")}
          </p>
        </>
      )}
    </div>
  );
}
