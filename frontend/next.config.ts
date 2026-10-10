import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  reactStrictMode: true,
  // Support standalone output for clean container deployment
  output: "standalone",
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
