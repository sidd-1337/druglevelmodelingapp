"use client";

import React from "react";
import { RouteType } from "@/lib/types";
import Segmented, { SegmentOption } from "@/components/ui/Segmented";

const ROUTES: SegmentOption<RouteType>[] = [
  { value: "oral", label: "Perorální" },
  { value: "iv_bolus", label: "Nitrožilní" },
  { value: "insufflated", label: "Nazální" },
  { value: "smoked", label: "Inhalační" },
];

interface RouteSelectorProps {
  value: RouteType;
  onChange: (route: RouteType) => void;
}

export default function RouteSelector({ value, onChange }: RouteSelectorProps) {
  return <Segmented label="Cesta podání" value={value} options={ROUTES} onChange={onChange} />;
}
