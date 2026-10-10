import React from "react";
import { DrugPreset, SuklPackageDetail } from "@/lib/types";
import { FileText, ExternalLink, ShieldCheck } from "lucide-react";

interface AboutDrugTabProps {
  drug: DrugPreset | null;
  suklDetail: SuklPackageDetail;
}

export default function AboutDrugTab({ drug, suklDetail }: AboutDrugTabProps) {
  return (
    <div className="space-y-6">
      {/* Top Summary Card */}
      <div className="apple-card p-6 space-y-6">
        <div className="flex flex-wrap items-start justify-between gap-4 pb-5 border-b border-zinc-200/80 dark:border-zinc-800">
          <div>
            <div className="flex items-center gap-2.5">
              <h3 className="text-xl font-semibold tracking-tight text-zinc-950 dark:text-zinc-50">
                {suklDetail.nazev}
              </h3>
              <span className="text-xs px-2.5 py-1 rounded-full bg-zinc-100 dark:bg-zinc-800 text-zinc-700 dark:text-zinc-300 border border-zinc-200/80 dark:border-zinc-700 font-medium">
                SÚKL {suklDetail.kodSukl}
              </span>
            </div>
            <p className="text-sm text-zinc-500 dark:text-zinc-400 mt-1 font-normal">
              {suklDetail.doplnek || "Registrováno v databázi SÚKL"}
            </p>
          </div>

          <div className="flex items-center gap-2">
            <span
              className={`text-xs px-3 py-1 rounded-full font-medium border ${
                suklDetail.vydejKod === "blue_stripe"
                  ? "bg-blue-50 text-blue-700 border-blue-200/80 dark:bg-blue-950/40 dark:text-blue-300 dark:border-blue-800/60"
                  : suklDetail.vydejKod === "otc"
                  ? "bg-emerald-50 text-emerald-700 border-emerald-200/80 dark:bg-emerald-950/40 dark:text-emerald-300 dark:border-emerald-800/60"
                  : "bg-zinc-100 text-zinc-700 border-zinc-200 dark:bg-zinc-800 dark:text-zinc-300 dark:border-zinc-700"
              }`}
            >
              {suklDetail.vydejText}
            </span>
          </div>
        </div>

        {/* Technical Data Grid */}
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3.5 text-sm">
          <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60">
            <span className="text-xs text-zinc-500 dark:text-zinc-400 block font-normal">Účinná látka</span>
            <span className="font-semibold text-zinc-900 dark:text-zinc-100 mt-0.5 block">
              {suklDetail.lecivaLatka || drug?.substance || suklDetail.atcNazev}
            </span>
            <span className="text-xs text-blue-600 dark:text-blue-400 font-medium">
              {suklDetail.atcNazev && suklDetail.atcNazev !== (suklDetail.lecivaLatka || "") ? `${suklDetail.atcNazev} ` : ""}
              {suklDetail.atcKod ? `(${suklDetail.atcKod})` : ""}
            </span>
          </div>

          <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60">
            <span className="text-xs text-zinc-500 dark:text-zinc-400 block font-normal">Léková forma & Síla</span>
            <span className="font-semibold text-zinc-900 dark:text-zinc-100 mt-0.5 block">{suklDetail.sila}</span>
            <span className="text-xs text-zinc-500 dark:text-zinc-400">{suklDetail.forma}</span>
          </div>

          <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60">
            <span className="text-xs text-zinc-500 dark:text-zinc-400 block font-normal">Cesta podání</span>
            <span className="font-semibold text-zinc-900 dark:text-zinc-100 mt-0.5 block">{suklDetail.cesta}</span>
          </div>

          <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60">
            <span className="text-xs text-zinc-500 dark:text-zinc-400 block font-normal">Registrační číslo</span>
            <span className="text-zinc-900 dark:text-zinc-100 mt-0.5 block font-medium">{suklDetail.registracniCislo}</span>
          </div>

          <div className="p-4 rounded-xl bg-zinc-50 dark:bg-zinc-800/60 border border-zinc-200/80 dark:border-zinc-700/60 md:col-span-2">
            <span className="text-xs text-zinc-500 dark:text-zinc-400 block font-normal">Držitel rozhodnutí o registraci</span>
            <span className="font-semibold text-zinc-900 dark:text-zinc-100 mt-0.5 block">
              {suklDetail.drzitel || "Neuvedeno v databázi"}
            </span>
          </div>
        </div>
      </div>

      {/* Official SÚKL PDF Documents Section */}
      <div className="apple-card p-6 space-y-4">
        <div>
          <h4 className="font-semibold text-base text-zinc-950 dark:text-zinc-50 flex items-center gap-2">
            <FileText className="w-4 h-4 text-blue-600 dark:text-blue-400" />
            <span>Oficiální dokumentace SÚKL</span>
          </h4>
          <p className="text-xs text-zinc-500 dark:text-zinc-400 mt-0.5">
            Schválené dokumenty Státním ústavem pro kontrolu léčiv ve formátu PDF.
          </p>
        </div>

        <div className="grid grid-cols-1 sm:grid-cols-3 gap-3.5">
          {suklDetail.spcUrl && (
            <a
              href={suklDetail.spcUrl}
              target="_blank"
              rel="noopener noreferrer"
              className="p-4 rounded-xl bg-zinc-50 hover:bg-zinc-100/90 dark:bg-zinc-800/60 dark:hover:bg-zinc-800 border border-zinc-200/80 dark:border-zinc-700/60 transition-all flex items-center justify-between group"
            >
              <div>
                <div className="font-semibold text-sm text-zinc-900 dark:text-zinc-100 group-hover:text-blue-600 dark:group-hover:text-blue-400 transition-colors">
                  Stáhnout SPC (PDF)
                </div>
                <div className="text-xs text-zinc-500 dark:text-zinc-400 mt-0.5 font-normal">Souhrn údajů pro lékaře</div>
              </div>
              <ExternalLink className="w-4 h-4 text-zinc-400 group-hover:text-blue-600 dark:group-hover:text-blue-400 group-hover:translate-x-0.5 transition-all" />
            </a>
          )}

          {suklDetail.pilUrl && (
            <a
              href={suklDetail.pilUrl}
              target="_blank"
              rel="noopener noreferrer"
              className="p-4 rounded-xl bg-zinc-50 hover:bg-zinc-100/90 dark:bg-zinc-800/60 dark:hover:bg-zinc-800 border border-zinc-200/80 dark:border-zinc-700/60 transition-all flex items-center justify-between group"
            >
              <div>
                <div className="font-semibold text-sm text-zinc-900 dark:text-zinc-100 group-hover:text-blue-600 dark:group-hover:text-blue-400 transition-colors">
                  Stáhnout PIL (PDF)
                </div>
                <div className="text-xs text-zinc-500 dark:text-zinc-400 mt-0.5 font-normal">Příbalová informace pro pacienta</div>
              </div>
              <ExternalLink className="w-4 h-4 text-zinc-400 group-hover:text-blue-600 dark:group-hover:text-blue-400 group-hover:translate-x-0.5 transition-all" />
            </a>
          )}

          <a
            href={suklDetail.suklUrl}
            target="_blank"
            rel="noopener noreferrer"
            className="p-4 rounded-xl bg-zinc-50 hover:bg-zinc-100/90 dark:bg-zinc-800/60 dark:hover:bg-zinc-800 border border-zinc-200/80 dark:border-zinc-700/60 transition-all flex items-center justify-between group"
          >
            <div>
              <div className="font-semibold text-sm text-zinc-900 dark:text-zinc-100 group-hover:text-blue-600 dark:group-hover:text-blue-400 transition-colors">
                Karta v databázi SÚKL
              </div>
              <div className="text-xs text-zinc-500 dark:text-zinc-400 mt-0.5 font-normal">Oficiální web prehledy.sukl.cz</div>
            </div>
            <ExternalLink className="w-4 h-4 text-zinc-400 group-hover:text-blue-600 dark:group-hover:text-blue-400 group-hover:translate-x-0.5 transition-all" />
          </a>
        </div>
      </div>

      {/* Clinical Notes Card */}
      {drug && drug.notes && (
        <div className="apple-card p-6">
          <div className="flex items-center gap-2 mb-2">
            <ShieldCheck className="w-4 h-4 text-emerald-600 dark:text-emerald-400" />
            <h4 className="font-semibold text-sm text-zinc-950 dark:text-zinc-50">
              Klinický farmakologický souhrn
            </h4>
          </div>
          <p className="text-sm text-zinc-600 dark:text-zinc-300 leading-relaxed font-normal">
            {drug.notes}
          </p>
        </div>
      )}
    </div>
  );
}
