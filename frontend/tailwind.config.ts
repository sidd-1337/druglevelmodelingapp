import type { Config } from "tailwindcss";

// Every colour is a CSS variable (space-separated RGB channels) so light and dark
// share one set of class names and alpha modifiers (bg-ink/10) keep working.
const token = (name: string) => `rgb(var(--${name}) / <alpha-value>)`;

const config: Config = {
  darkMode: "class",
  content: [
    "./pages/**/*.{js,ts,jsx,tsx,mdx}",
    "./components/**/*.{js,ts,jsx,tsx,mdx}",
    "./app/**/*.{js,ts,jsx,tsx,mdx}",
  ],
  theme: {
    // Replaces Tailwind's palette on purpose: the site has one neutral ramp,
    // one blue and three status hues, nothing else.
    colors: {
      transparent: "transparent",
      current: "currentColor",
      bg: token("bg"),
      band: token("band"),
      ink: token("ink"),
      "ink-2": token("ink-2"),
      line: token("line"),
      "line-soft": token("line-soft"),
      link: token("link"),
      accent: token("accent"),
      ok: token("ok"),
      warn: token("warn"),
      bad: token("bad"),
      "ok-fill": token("ok-fill"),
      "warn-fill": token("warn-fill"),
      "bad-fill": token("bad-fill"),
      thumb: token("thumb"),
    },
    borderColor: ({ theme }) => ({
      ...theme("colors"),
      DEFAULT: "rgb(var(--line-soft))",
    }),
    fontSize: {
      small: ["0.8125rem", { lineHeight: "1.4", letterSpacing: "-0.003em" }],
      ui: ["0.9375rem", { lineHeight: "1.4", letterSpacing: "-0.008em" }],
      body: ["1.0625rem", { lineHeight: "1.47", letterSpacing: "-0.011em" }],
      h3: ["1.3125rem", { lineHeight: "1.24", letterSpacing: "-0.012em", fontWeight: "600" }],
      sub: ["1.5rem", { lineHeight: "1.2", letterSpacing: "-0.018em", fontWeight: "600" }],
      lead: [
        "clamp(1.1875rem, 2vw, 1.5rem)",
        { lineHeight: "1.38", letterSpacing: "-0.014em", fontWeight: "400" },
      ],
      figure: [
        "clamp(2rem, 4.2vw, 3rem)",
        { lineHeight: "1.05", letterSpacing: "-0.03em", fontWeight: "600" },
      ],
      headline: [
        "clamp(1.75rem, 3.2vw, 2.5rem)",
        { lineHeight: "1.1", letterSpacing: "-0.025em", fontWeight: "600" },
      ],
      title: [
        "clamp(2.25rem, 5vw, 3.5rem)",
        { lineHeight: "1.06", letterSpacing: "-0.032em", fontWeight: "600" },
      ],
      display: [
        "clamp(2.75rem, 7.2vw, 5.25rem)",
        { lineHeight: "1.02", letterSpacing: "-0.04em", fontWeight: "600" },
      ],
    },
    borderRadius: {
      none: "0",
      // Two radii only: controls, and the surfaces that hold a tool.
      ctl: "12px",
      panel: "20px",
      full: "9999px",
    },
    extend: {
      boxShadow: {
        thumb: "0 0 0 0.5px rgb(0 0 0 / 0.08), 0 1px 4px rgb(0 0 0 / 0.18)",
        pop: "0 12px 40px -8px rgb(0 0 0 / 0.22), 0 0 0 0.5px rgb(0 0 0 / 0.1)",
      },
      fontFamily: {
        sans: [
          "-apple-system",
          "BlinkMacSystemFont",
          '"SF Pro Display"',
          '"SF Pro Text"',
          '"Inter Variable"',
          '"Segoe UI"',
          "system-ui",
          "sans-serif",
        ],
      },
      maxWidth: {
        page: "1080px",
        prose: "40rem",
      },
    },
  },
  plugins: [],
};

export default config;
