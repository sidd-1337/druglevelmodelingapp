import React from "react";
import { cz } from "@/lib/typo";

interface SectionHeadProps {
  id?: string;
  title: string;
  size?: "headline" | "sub";
  children?: React.ReactNode;
}

export default function SectionHead({ id, title, size = "sub", children }: SectionHeadProps) {
  return (
    <div className="max-w-prose">
      <h2 id={id} className={size === "headline" ? "text-headline" : "text-sub"}>{title}</h2>
      {children && <p className="mt-2 text-body text-ink-2">{typeof children === "string" ? cz(children) : children}</p>}
    </div>
  );
}
