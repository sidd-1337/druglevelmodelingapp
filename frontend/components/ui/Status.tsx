import React from "react";

export type Tone = "ok" | "warn" | "bad" | "neutral";

const DOT: Record<Tone, string> = {
  ok: "bg-ok-fill",
  warn: "bg-warn-fill",
  bad: "bg-bad-fill",
  neutral: "bg-ink-2",
};

/** A status is a dot and a sentence. The words carry the meaning; the dot only repeats it. */
export default function Status({
  tone = "neutral",
  children,
}: {
  tone?: Tone;
  children: React.ReactNode;
}) {
  return (
    <span className="inline-flex items-start gap-2">
      <span aria-hidden className={`mt-[0.45em] h-2 w-2 shrink-0 rounded-full ${DOT[tone]}`} />
      <span>{children}</span>
    </span>
  );
}
