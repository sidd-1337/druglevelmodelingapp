import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  reactStrictMode: true,
  // Support standalone output for clean container deployment
  output: "standalone",
  // Single-origin: prohlížeč volá /api/v1/* na vlastní doméně a Next to proxuje na backend.
  async rewrites() {
    const api = process.env.INTERNAL_API_URL || "http://localhost:8000";
    return [{ source: "/api/v1/:path*", destination: `${api}/api/v1/:path*` }];
  },
  // Allow SVGs and images from official portals if needed
  images: {
    remotePatterns: [
      {
        protocol: "https",
        hostname: "www.sukl.cz",
      },
    ],
  },
};

export default nextConfig;
