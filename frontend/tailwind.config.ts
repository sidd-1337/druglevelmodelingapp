import type { Config } from "tailwindcss";

const config: Config = {
  darkMode: "class",
  content: [
    "./pages/**/*.{js,ts,jsx,tsx,mdx}",
    "./components/**/*.{js,ts,jsx,tsx,mdx}",
    "./app/**/*.{js,ts,jsx,tsx,mdx}",
  ],
  theme: {
    extend: {
      colors: {
        zinc: {
          150: "#ececee",
          750: "#323236",
          850: "#1f1f23",
        },
        apple: {
          canvas: "#f5f5f7",
          darkcanvas: "#000000",
          card: "#ffffff",
          darkcard: "#18181b",
          border: "#e4e4e7",
          darkborder: "#27272a",
          blue: "#0071e3",
          accent: "#2997ff",
          subtext: "#86868b",
          darksubtext: "#a1a1a6",
        },
      },
      boxShadow: {
        "apple-sm": "0 1px 2px 0 rgba(0, 0, 0, 0.05)",
        "apple-card": "0 2px 8px -2px rgba(0, 0, 0, 0.05), 0 0 0 1px rgba(0, 0, 0, 0.06)",
        "apple-card-dark": "0 4px 20px -2px rgba(0, 0, 0, 0.5), 0 0 0 1px rgba(255, 255, 255, 0.08)",
        "apple-float": "0 20px 40px -15px rgba(0, 0, 0, 0.1), 0 0 0 1px rgba(0, 0, 0, 0.06)",
      },
      fontFamily: {
        sans: [
          "-apple-system",
          "BlinkMacSystemFont",
          '"SF Pro Text"',
          '"SF Pro Display"',
          '"Segoe UI"',
          "Inter",
          "system-ui",
          "sans-serif",
        ],
        mono: ['"SF Mono"', "Menlo", "Monaco", "Consolas", "monospace"],
      },
    },
  },
  plugins: [],
};

export default config;
