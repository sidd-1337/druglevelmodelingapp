"use client";

import React from "react";
import { RouteType } from "@/lib/types";

interface RouteSelectorProps {
  value: RouteType;
  onChange: (route: RouteType) => void;
}

const ROUTES: Array<{ id: RouteType; label: string; icon: string; badge: string }> = [
  { id: "oral", label: "Oral", icon: "💊", badge: "Perorální" },
  { id: "iv_bolus", label: "IV Bolus", icon: "💉", badge: "Nitrožilní" },
  { id: "insufflated", label: "Insufflated", icon: "👃", badge: "Nazální (nosní)" },
  { id: "smoked", label: "Smoked", icon: "💨", badge: "Inhalační (kouření)" },
];

export default function RouteSelector({ value, onChange }: RouteSelectorProps) {
  return (
    <div className="inline-flex p-1 bg-zinc-100 dark:bg-zinc-900 rounded-2xl w-full border border-zinc-200/80 dark:border-zinc-800">
      {ROUTES.map((r) => {
        const isActive = value === r.id;
        return (
          <button
            key={r.id}
            type="button"
            onClick={() => onChange(r.id)}
            className={`flex-1 flex flex-col items-center justify-center py-2.5 px-2 rounded-xl text-center transition-all ${
              isActive
                ? "bg-white dark:bg-zinc-800 text-zinc-950 dark:text-zinc-50 shadow-sm border border-black/5 dark:border-white/10 font-semibold"
                : "text-zinc-600 dark:text-zinc-400 hover:text-zinc-900 dark:hover:text-zinc-100 font-medium"
            }`}
          >
            <div className="flex items-center gap-1.5 text-xs sm:text-sm">
              <span>{r.icon}</span>
              <span>{r.label}</span>
            </div>
            <span className="text-[10px] text-zinc-400 dark:text-zinc-500 mt-0.5">{r.badge}</span>
          </button>
        );
      })}
    </div>
  );
}
