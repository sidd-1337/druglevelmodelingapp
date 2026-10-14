import type { Metadata, Viewport } from "next";
import Link from "next/link";
import "@fontsource-variable/inter/opsz.css";
import Mark from "@/components/Mark";
import SiteNav from "@/components/SiteNav";
import ThemeToggle from "@/components/ThemeToggle";
import "./globals.css";

export const metadata: Metadata = {
  title: "Lékový portál ČR – cena, úhrada, dostupnost a hladina léku v těle",
  description:
    "Vyhledejte lék v registru SÚKL a zjistěte jeho cenu a úhradu, dostupnost v lékárnách a modelovanou hladinu v krvi.",
};

export const viewport: Viewport = {
  themeColor: [
    { media: "(prefers-color-scheme: light)", color: "#fbfbfd" },
    { media: "(prefers-color-scheme: dark)", color: "#000000" },
  ],
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="cs" suppressHydrationWarning>
      <head>
        {/* Applies the saved or system theme before first paint to avoid a flash. */}
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
      <body className="flex min-h-screen flex-col font-sans">
        <a
          href="#obsah"
          className="sr-only focus:not-sr-only focus:fixed focus:left-4 focus:top-3 focus:z-[60] focus:rounded-ctl focus:bg-bg focus:px-4 focus:py-2 focus:text-ui focus:shadow-pop"
        >
          Přeskočit na obsah
        </a>

        <header className="sticky top-0 z-50 border-b border-line-soft bg-bg/80 backdrop-blur-xl backdrop-saturate-150">
          <div className="page flex h-12 items-center justify-between">
            <Link href="/" className="flex items-center gap-2.5 text-ui font-semibold">
              <Mark />
              <span>Lékový portál ČR</span>
            </Link>
            <div className="flex items-center gap-1">
              <SiteNav />
              <ThemeToggle />
            </div>
          </div>
        </header>

        <main id="obsah" className="flex-1">
          {children}
        </main>

        <footer className="bg-band">
          <div className="page py-12 text-small text-ink-2">
            <p className="max-w-prose">
              Zdravotnické upozornění. Data vycházejí z otevřených registrů SÚKL ČR. Portál slouží k výuce,
              modelování a snižování rizik při užívání látek. Simulace jsou orientační, nejsou zdravotnickým
              prostředkem a nenahrazují radu lékaře ani lékárníka.
            </p>
            <div className="mt-8 flex flex-wrap items-center justify-between gap-x-8 gap-y-2 border-t border-line pt-4">
              <span>Zdroj dat: SÚKL, otevřené registry</span>
              <Link href="/slovnik" className="hover:text-ink">
                Slovník pojmů
              </Link>
            </div>
          </div>
        </footer>
      </body>
    </html>
  );
}
