"use client";

import React, { useEffect, useState } from "react";
import { Sun, Moon } from "lucide-react";

export default function ThemeToggle() {
  const [isDark, setIsDark] = useState<boolean>(false);
  const [mounted, setMounted] = useState<boolean>(false);

  useEffect(() => {
    setMounted(true);
    const stored = localStorage.getItem("theme");
    const systemPrefersDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
    const activeIsDark = stored ? stored === "dark" : systemPrefersDark;

    setIsDark(activeIsDark);
    if (activeIsDark) {
      document.documentElement.classList.add("dark");
    } else {
      document.documentElement.classList.remove("dark");
    }
  }, []);

  const toggleTheme = () => {
    const nextIsDark = !isDark;
    setIsDark(nextIsDark);

    if (nextIsDark) {
      document.documentElement.classList.add("dark");
      localStorage.setItem("theme", "dark");
    } else {
      document.documentElement.classList.remove("dark");
      localStorage.setItem("theme", "light");
    }
  };

  // Reserve the space so the header does not shift when the icon appears.
  if (!mounted) return <div className="h-8 w-8" aria-hidden />;

  return (
    <button
      type="button"
      onClick={toggleTheme}
      aria-label={isDark ? "Přepnout na světlý režim" : "Přepnout na tmavý režim"}
      title={isDark ? "Světlý režim" : "Tmavý režim"}
      className="inline-flex h-8 w-8 items-center justify-center rounded-ctl text-ink-2 transition-colors hover:bg-ink/[0.07] hover:text-ink"
    >
      {isDark ? <Sun className="h-[18px] w-[18px]" strokeWidth={1.75} /> : <Moon className="h-[18px] w-[18px]" strokeWidth={1.75} />}
    </button>
  );
}
