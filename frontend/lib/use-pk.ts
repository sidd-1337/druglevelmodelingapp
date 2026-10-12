"use client";

import { useEffect, useState } from "react";

import {
  PkHarmReduction,
  PkSimulation,
  harmReductionApi,
  simulatePkApi,
} from "./api";
import { RouteType } from "./types";

export interface PkRequest {
  drugId: string | null;
  doseMg: number;
  intervalH: number;
  durationH: number;
  bodyWeightKg: number;
  route: RouteType;
  elapsedH: number;
}

export type PkSource = "api" | "local";

export interface PkResult {
  simulation: PkSimulation | null;
  harmReduction: PkHarmReduction | null;
  /** "api" once the backend has answered; "local" while loading or if it is unreachable. */
  source: PkSource;
}

/**
 * Pharmacokinetics from the backend, which is the single source of truth.
 *
 * The Python core models active metabolites, enzyme inactivation, elimination
 * pathways and genotype - none of which the local TypeScript engine can represent.
 * Callers fall back to that engine only while this is loading or when the backend
 * cannot be reached, and should say so when they do.
 */
export function usePk(req: PkRequest): PkResult {
  const [result, setResult] = useState<PkResult>({
    simulation: null,
    harmReduction: null,
    source: "local",
  });

  const { drugId, doseMg, intervalH, durationH, bodyWeightKg, route, elapsedH } = req;

  useEffect(() => {
    if (!drugId) {
      setResult({ simulation: null, harmReduction: null, source: "local" });
      return;
    }
    let active = true;
    Promise.all([
      simulatePkApi({ drugId, doseMg, intervalH, durationH, bodyWeightKg, route }),
      harmReductionApi({ drugId, doseMg, bodyWeightKg, route, elapsedH }),
    ])
      .then(([simulation, harmReduction]) => {
        if (!active) return;
        setResult({
          simulation,
          harmReduction,
          source: simulation && harmReduction ? "api" : "local",
        });
      })
      .catch(() => {
        if (active) setResult({ simulation: null, harmReduction: null, source: "local" });
      });
    return () => {
      active = false; // a later request must not be overwritten by an earlier reply
    };
  }, [drugId, doseMg, intervalH, durationH, bodyWeightKg, route, elapsedH]);

  return result;
}
