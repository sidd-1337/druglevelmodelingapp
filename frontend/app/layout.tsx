import type { Metadata } from "next";
import Link from "next/link";
import ThemeToggle from "@/components/ThemeToggle";
import "./globals.css";

export const metadata: Metadata = {
  title: "Lékový Portál ČR – Informace, Úhrady, Výpadky & Farmakokinetika",
  description:
    "Jednoduchý a moderní lékový portál s otevřenými daty SÚKL ČR, cenami a úhradami pojišťoven, výpadky na trhu a klinickou farmakokinetikou.",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="cs" suppressHydrationWarning>
      <head>
        {/* FOUC prevention script */}
        <script
          dangerouslySetInnerHTML={{
            __html: `
              try {
                const stored = localStorage.getItem('theme');
                const isDark = stored ? stored === 'dark' : window.matchMedia('(prefers-color-scheme: dark)').matches;
                if (isDark) {
                  document.documentElement.classList.add('dark');
                } else {
                  document.documentElement.classList.remove('dark');
                }
              } catch (e) {}
            `,
          }}
        />
      </head>
      <body className="flex flex-col min-h-screen text-zinc-900 dark:text-zinc-100 antialiased selection:bg-blue-500/20 selection:text-blue-700 dark:selection:text-blue-300">
        {/* Apple Style Frosted Header */}
        <header className="sticky top-0 z-50 backdrop-blur-xl bg-white/80 dark:bg-black/80 border-b border-zinc-200/80 dark:border-zinc-800 px-4 sm:px-8 py-3 transition-colors">
          <div className="max-w-6xl mx-auto flex items-center justify-between">
            <Link href="/" className="flex items-center gap-2.5 group">
              <div className="w-8 h-8 rounded-full bg-zinc-900 dark:bg-zinc-100 text-white dark:text-zinc-950 flex items-center justify-center text-sm font-semibold shadow-sm transition-transform group-hover:scale-105">
                💊
              </div>
              <div className="flex items-baseline gap-2">
                <span className="font-semibold text-base tracking-tight text-zinc-950 dark:text-zinc-50">
                  Lékový portál ČR
                </span>
                <span className="hidden sm:inline-block text-[11px] font-medium text-zinc-500 dark:text-zinc-400">
                  SÚKL & PK
                </span>
              </div>
            </Link>

            <div className="flex items-center gap-3 sm:gap-6">
              <nav className="flex items-center gap-1 sm:gap-4 text-sm font-medium">
                <Link
                  href="/"
                  className="px-2.5 py-1.5 rounded-lg text-zinc-600 hover:text-zinc-950 dark:text-zinc-400 dark:hover:text-zinc-100 transition-colors"
                >
                  Hledání
                </Link>
                <Link
                  href="/slovnik"
                  className="px-2.5 py-1.5 rounded-lg text-zinc-600 hover:text-zinc-950 dark:text-zinc-400 dark:hover:text-zinc-100 transition-colors"
                >
                  Slovník
                </Link>
              </nav>

              <div className="pl-2 sm:pl-3 border-l border-zinc-200 dark:border-zinc-800">
                <ThemeToggle />
              </div>
            </div>
          </div>
        </header>

        {/* Main View */}
        <main className="flex-1 max-w-6xl w-full mx-auto px-4 sm:px-8 py-8 sm:py-12">
          {children}
        </main>

        {/* Apple Style Restrained Footer */}
        <footer className="border-t border-zinc-200/80 dark:border-zinc-800 bg-white/50 dark:bg-black/40 px-4 sm:px-8 py-8 text-xs text-zinc-500 dark:text-zinc-400 transition-colors">
          <div className="max-w-6xl mx-auto flex flex-col sm:flex-row items-center justify-between gap-4">
            <p>
              Zdravotnické upozornění: Data vycházejí z otevřených registrů SÚKL ČR. Aplikace slouží pro výukové, modelovací a harm reduction účely.
            </p>
            <div className="flex items-center gap-3 text-zinc-400 dark:text-zinc-500">
              <span>Next.js 15</span>
              <span>•</span>
              <span>FastAPI</span>
              <span>•</span>
              <span>Data SÚKL</span>
            </div>
          </div>
        </footer>
      </body>
    </html>
  );
}
