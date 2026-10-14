"use client";

import React from "react";
import Link from "next/link";
import { usePathname } from "next/navigation";

const LINKS = [
  { href: "/", label: "Hledání" },
  { href: "/slovnik", label: "Slovník" },
];

export default function SiteNav() {
  const pathname = usePathname();

  return (
    <nav aria-label="Hlavní navigace" className="flex items-center">
      {LINKS.map((link) => {
        // The drug pages belong to "Hledání".
        const current =
          link.href === "/" ? pathname === "/" || pathname.startsWith("/leky") : pathname.startsWith(link.href);
        return (
          <Link
            key={link.href}
            href={link.href}
            aria-current={current ? "page" : undefined}
            className={`rounded-ctl px-3 py-1.5 text-ui transition-colors ${
              current ? "font-medium text-ink" : "text-ink-2 hover:text-ink"
            }`}
          >
            {link.label}
          </Link>
        );
      })}
    </nav>
  );
}
