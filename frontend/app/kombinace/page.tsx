import React, { Suspense } from "react";
import type { Metadata } from "next";
import CombinationTool from "@/components/CombinationTool";

export const metadata: Metadata = {
  title: "Kombinace léků – hladiny, interakce a střety účinků",
  description:
    "Několik léků najednou: hladiny v jednom grafu, interakce v čase, aktivita enzymů a sčítání účinků podle SPC.",
};

export default function KombinacePage() {
  return (
    <div className="page pb-16 pt-14 sm:pt-20">
      <h1 className="text-title">Kombinace léků</h1>
      <p className="mt-4 max-w-prose text-lead text-ink-2">
        Hladiny několika léků najednou, jak se navzájem ovlivňují v čase a kdy se jejich účinky sčítají.
      </p>
      <div className="mt-14 sm:mt-16">
        {/* useSearchParams reads the regimen from the link, which needs a Suspense boundary. */}
        <Suspense fallback={null}>
          <CombinationTool />
        </Suspense>
      </div>
    </div>
  );
}
