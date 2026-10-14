const NBSP = " ";

// One-letter preposition or conjunction (k, s, v, z, o, u, a, i) followed by a space.
const ONE_LETTER = /(^|[\s(])([KkSsVvZzOoUuAaIi])\s+(?=\S)/g;

/**
 * Czech typographic rule: a one-letter preposition or conjunction never ends a line,
 * and a number stays with its unit. Apply to running text.
 */
export function cz(text: string): string {
  // Two passes: the first consumes the space before the next one-letter word in "a v lékárně".
  return text
    .replace(ONE_LETTER, `$1$2${NBSP}`)
    .replace(ONE_LETTER, `$1$2${NBSP}`)
    .replace(/(\d)\s+(mg|µg|g|kg|ml|mL|l|h|%|Kč|let)(?![\p{L}])/gu, `$1${NBSP}$2`)
    .replace(/§\s+(\d)/g, `§${NBSP}$1`)
    .replace(/č\.\s+(\d)/g, `č.${NBSP}$1`);
}
