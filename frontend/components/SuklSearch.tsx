"use client";

import React, { useEffect, useRef, useState, useTransition } from "react";
import { useRouter } from "next/navigation";
import { Search, X, Loader2 } from "lucide-react";
import { SearchResultItem, searchDrugsApi } from "@/lib/api";
import { searchSuklCatalog } from "@/lib/sukl-data";

export default function SuklSearch() {
  const [query, setQuery] = useState("");
  const [results, setResults] = useState<SearchResultItem[]>([]);
  const [isOpen, setIsOpen] = useState(false);
  const [isSearching, setIsSearching] = useState(false);
  const router = useRouter();
  const [, startTransition] = useTransition();
  const searchTimeoutRef = useRef<NodeJS.Timeout | null>(null);

  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const val = e.target.value;
    setQuery(val);

    if (val.trim().length >= 2) {
      // 1. Instant local hits
      const localHits = searchSuklCatalog(val);
      setResults(localHits);
      setIsOpen(true);

      // 2. Debounced query to FastAPI backend
      if (searchTimeoutRef.current) {
        clearTimeout(searchTimeoutRef.current);
      }
      setIsSearching(true);
      searchTimeoutRef.current = setTimeout(async () => {
        try {
          const apiHits = await searchDrugsApi(val);
          if (apiHits.length > 0) {
            setResults(apiHits);
          }
        } finally {
          setIsSearching(false);
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

  const handleSelect = (kodSukl: string) => {
    setIsOpen(false);
    startTransition(() => {
      router.push(`/leky/${kodSukl}`);
    });
  };

  return (
    <div className="relative w-full max-w-2xl mx-auto">
      <div className="relative flex items-center">
        <div className="absolute left-4 text-zinc-400 dark:text-zinc-500 pointer-events-none">
          {isSearching ? (
            <Loader2 className="w-5 h-5 animate-spin text-blue-600 dark:text-blue-400" />
          ) : (
            <Search className="w-5 h-5" />
          )}
        </div>
        <input
          type="text"
          value={query}
          onChange={handleInputChange}
          onFocus={() => query.length >= 2 && setIsOpen(true)}
          placeholder="Vyhledat lék podle názvu, účinné látky nebo kódu SÚKL..."
          className="w-full pl-12 pr-10 py-3.5 bg-white dark:bg-zinc-900 border border-zinc-200/90 dark:border-zinc-800 rounded-2xl text-zinc-900 dark:text-zinc-100 placeholder-zinc-400 dark:placeholder-zinc-500 shadow-[0_2px_12px_rgba(0,0,0,0.04)] dark:shadow-none focus:outline-none focus:border-blue-500 dark:focus:border-blue-400 focus:ring-4 focus:ring-blue-500/10 transition-all text-sm sm:text-base font-normal"
        />
        {query && (
          <button
            onClick={() => {
              setQuery("");
              setResults([]);
              setIsOpen(false);
            }}
            className="absolute right-4 text-zinc-400 hover:text-zinc-700 dark:hover:text-zinc-200 transition-colors"
          >
            <X className="w-4 h-4" />
          </button>
        )}
      </div>

      {isOpen && results.length > 0 && (
        <div className="absolute left-0 right-0 mt-2 bg-white/95 dark:bg-zinc-900/95 backdrop-blur-xl border border-zinc-200/90 dark:border-zinc-800 rounded-2xl shadow-xl overflow-hidden z-50 divide-y divide-zinc-100 dark:divide-zinc-800/80">
          <div className="px-4 py-2.5 text-xs font-medium text-zinc-500 dark:text-zinc-400 bg-zinc-50/70 dark:bg-zinc-950/40 flex items-center justify-between">
            <span>Nalezená léčiva ({results.length})</span>
            {isSearching && (
              <span className="text-[11px] text-blue-600 dark:text-blue-400 font-medium">Hledám v SÚKL...</span>
            )}
          </div>
          <div className="max-h-72 overflow-y-auto">
            {results.map((item) => (
              <button
                key={item.kodSukl}
                onClick={() => handleSelect(item.kodSukl)}
                className="w-full text-left px-4 py-3 hover:bg-zinc-50 dark:hover:bg-zinc-800/60 transition-colors flex items-center justify-between gap-3 group"
              >
                <div>
                  <div className="font-semibold text-sm text-zinc-900 dark:text-zinc-100 group-hover:text-blue-600 dark:group-hover:text-blue-400 transition-colors">
                    {item.nazev}
                  </div>
                  <div className="text-xs text-zinc-500 dark:text-zinc-400 line-clamp-1">
                    {item.doplnek || "Léčivý přípravek"}
                  </div>
                </div>
                <div className="text-right shrink-0">
                  <span className="inline-block font-medium text-[11px] px-2 py-0.5 rounded-md bg-zinc-100 dark:bg-zinc-800 text-zinc-600 dark:text-zinc-300 border border-zinc-200/60 dark:border-zinc-700/60">
                    {item.kodSukl}
                  </span>
                </div>
              </button>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
