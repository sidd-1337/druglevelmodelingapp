"use client";

import React, { KeyboardEvent, useEffect, useId, useRef, useState, useTransition } from "react";
import { useRouter } from "next/navigation";
import { Search, X, Loader2 } from "lucide-react";
import { SearchResultItem, searchDrugsApi } from "@/lib/api";
import { packageText, sentenceCase } from "@/lib/format";
import { searchSuklCatalog } from "@/lib/sukl-data";

export default function SuklSearch() {
  const [query, setQuery] = useState("");
  const [results, setResults] = useState<SearchResultItem[]>([]);
  const [isOpen, setIsOpen] = useState(false);
  const [isSearching, setIsSearching] = useState(false);
  const [active, setActive] = useState(-1);
  const router = useRouter();
  const [, startTransition] = useTransition();
  const searchTimeoutRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  const requestRef = useRef(0);
  const rootRef = useRef<HTMLDivElement>(null);
  const listId = useId();

  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const val = e.target.value;
    setQuery(val);
    setActive(-1);
    const request = ++requestRef.current;

    if (val.trim().length >= 2) {
      // 1. Instant local hits
      setResults(searchSuklCatalog(val));
      setIsOpen(true);

      // 2. Debounced query to the FastAPI backend; a newer keystroke makes this one stale.
      if (searchTimeoutRef.current) {
        clearTimeout(searchTimeoutRef.current);
      }
      setIsSearching(true);
      searchTimeoutRef.current = setTimeout(async () => {
        try {
          const apiHits = await searchDrugsApi(val);
          if (request === requestRef.current && apiHits.length > 0) {
            setResults(apiHits);
          }
        } finally {
          if (request === requestRef.current) setIsSearching(false);
        }
      }, 250);
    } else {
      setResults([]);
      setIsOpen(false);
      setIsSearching(false);
    }
  };

  useEffect(() => {
    return () => {
      if (searchTimeoutRef.current) {
        clearTimeout(searchTimeoutRef.current);
      }
    };
  }, []);

  // Close the list when the pointer goes down anywhere outside the search.
  useEffect(() => {
    const onPointerDown = (e: PointerEvent) => {
      if (rootRef.current && !rootRef.current.contains(e.target as Node)) setIsOpen(false);
    };
    document.addEventListener("pointerdown", onPointerDown);
    return () => document.removeEventListener("pointerdown", onPointerDown);
  }, []);

  const handleSelect = (kodSukl: string) => {
    setIsOpen(false);
    startTransition(() => {
      router.push(`/leky/${kodSukl}`);
    });
  };

  const clear = () => {
    requestRef.current++;
    setQuery("");
    setResults([]);
    setIsOpen(false);
    setIsSearching(false);
  };

  const onKeyDown = (e: KeyboardEvent<HTMLInputElement>) => {
    if (e.key === "Escape") {
      setIsOpen(false);
      return;
    }
    if (!results.length) return;
    if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      e.preventDefault();
      setIsOpen(true);
      const step = e.key === "ArrowDown" ? 1 : -1;
      setActive((i) => (i + step + results.length) % results.length);
    } else if (e.key === "Enter" && isOpen) {
      e.preventDefault();
      handleSelect(results[active >= 0 ? active : 0].kodSukl);
    }
  };

  const showEmpty = isOpen && !isSearching && query.trim().length >= 2 && results.length === 0;
  const showList = isOpen && results.length > 0;

  return (
    <div ref={rootRef} className="relative mx-auto w-full max-w-xl text-left">
      <div className="relative flex items-center">
        <Search
          className="pointer-events-none absolute left-4 h-5 w-5 text-ink-2"
          strokeWidth={1.75}
          aria-hidden
        />
        <input
          type="text"
          role="combobox"
          aria-expanded={showList}
          aria-controls={listId}
          aria-autocomplete="list"
          aria-activedescendant={active >= 0 ? `${listId}-${active}` : undefined}
          aria-label="Hledat lék"
          autoComplete="off"
          spellCheck={false}
          value={query}
          onChange={handleInputChange}
          onKeyDown={onKeyDown}
          onFocus={() => query.trim().length >= 2 && setIsOpen(true)}
          placeholder="Název léku, látka nebo kód SÚKL"
          className="h-[52px] w-full rounded-ctl bg-ink/[0.06] pl-12 pr-12 text-body text-ink outline-none transition-shadow placeholder:text-ink-2 focus:bg-bg focus:ring-2 focus:ring-accent"
        />
        {isSearching ? (
          <Loader2 className="absolute right-4 h-4 w-4 animate-spin text-ink-2" aria-hidden />
        ) : (
          query && (
            <button
              type="button"
              onClick={clear}
              aria-label="Vymazat hledání"
              className="absolute right-2.5 flex h-8 w-8 items-center justify-center rounded-full text-ink-2 hover:text-ink"
            >
              <X className="h-4 w-4" aria-hidden />
            </button>
          )
        )}
      </div>

      <div role="status" className="sr-only">
        {showList ? `Nalezeno ${results.length} léčiv` : showEmpty ? "Nic nenalezeno" : ""}
      </div>

      {showList && (
        <ul
          id={listId}
          role="listbox"
          aria-label="Nalezená léčiva"
          className="absolute inset-x-0 top-full z-30 mt-2 max-h-80 overflow-y-auto rounded-panel bg-bg py-2 shadow-pop ring-1 ring-line-soft"
        >
          {results.map((item, i) => (
            <li
              key={item.kodSukl}
              id={`${listId}-${i}`}
              role="option"
              aria-selected={i === active}
              onClick={() => handleSelect(item.kodSukl)}
              onMouseMove={() => setActive(i)}
              className={`flex cursor-pointer items-baseline justify-between gap-4 px-5 py-2.5 ${
                i === active ? "bg-ink/[0.06]" : ""
              }`}
            >
              <span className="min-w-0">
                <span className="block text-ui font-semibold text-ink">{sentenceCase(item.nazev)}</span>
                <span className="line-clamp-1 block text-small text-ink-2">
                  {packageText(item.doplnek) || "Léčivý přípravek"}
                </span>
              </span>
              <span className="figures shrink-0 text-small text-ink-2">{item.kodSukl}</span>
            </li>
          ))}
        </ul>
      )}

      {showEmpty && (
        <div className="absolute inset-x-0 top-full z-30 mt-2 rounded-panel bg-bg px-5 py-4 text-ui text-ink-2 shadow-pop ring-1 ring-line-soft">
          Nic jsme nenašli. Zkuste název účinné látky nebo kód SÚKL.
        </div>
      )}
    </div>
  );
}
