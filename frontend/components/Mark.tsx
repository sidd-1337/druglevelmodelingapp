import React from "react";

/** The concentration curve: quick rise, slow decay. Used as the logo glyph. */
export default function Mark({ className = "h-4 w-6" }: { className?: string }) {
  return (
    <svg viewBox="0 0 24 16" fill="none" aria-hidden className={className}>
      <path
        d="M1.5 14 C5 14 6.5 2.5 10 2.5 C14 2.5 16 12.5 22.5 14"
        stroke="currentColor"
        strokeWidth="2"
        strokeLinecap="round"
      />
    </svg>
  );
}
