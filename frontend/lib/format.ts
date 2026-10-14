// Czech number and text conventions: decimal comma, space as thousands separator,
// and a non-breaking space between a number and its unit.

const NBSP = " ";

const fixed: Record<number, Intl.NumberFormat> = {};

/** Number with a fixed count of decimals: 3.5 -> "3,50". */
export function num(value: number, digits = 0): string {
  fixed[digits] ??= new Intl.NumberFormat("cs-CZ", {
    minimumFractionDigits: digits,
    maximumFractionDigits: digits,
  });
  return fixed[digits].format(value);
}

const significant = new Intl.NumberFormat("cs-CZ", { maximumSignificantDigits: 3 });

/** Concentration or other value spanning orders of magnitude: 0.0183 -> "0,0183", 12.34 -> "12,3". */
export function sig(value: number): string {
  return significant.format(value);
}

/** Koruna amount: 5000 -> "5 000 Kč". */
export function kc(value: number, digits = 0): string {
  return `${num(value, digits)}${NBSP}Kč`;
}

/** Value with unit: ("12,3", "h") -> "12,3 h". */
export function withUnit(value: string, unit: string): string {
  return `${value}${NBSP}${unit}`;
}

/** SÚKL stores names in capitals ("RITALIN"). Show them as words: "Ritalin". Mixed-case text is left alone. */
export function sentenceCase(text: string): string {
  if (!text || text !== text.toLocaleUpperCase("cs")) return text;
  return text
    .toLocaleLowerCase("cs")
    .replace(/(^|[\s\-(/])(\p{L})/gu, (_, lead: string, ch: string) => lead + ch.toLocaleUpperCase("cs"));
}

/** Package text in capitals ("10MG TBL NOB 30") becomes "10 mg tbl nob 30". Mixed-case text is left alone. */
export function packageText(text: string): string {
  if (!text || text !== text.toLocaleUpperCase("cs")) return text;
  return text.toLocaleLowerCase("cs").replace(/(\d)(mg|mcg|µg|ml|iu|g)(?![\p{L}])/gu, `$1${NBSP}$2`);
}

/** Texts built by the simulation use decimal points; Czech uses a comma. */
export function czDecimals(text: string): string {
  return text.replace(/(\d)\.(\d)/g, "$1,$2");
}

/** Removes pictographs (✅ ⚠️) the simulation puts in headlines; status colour and wording carry the meaning. */
export function stripEmoji(text: string): string {
  return text.replace(/[\p{Extended_Pictographic}️‍]/gu, "").replace(/\s{2,}/g, " ").trim();
}

/** Capital letter first, the rest as written: "na lékařský předpis" -> "Na lékařský předpis". */
export function upperFirst(text: string): string {
  return text ? text.charAt(0).toLocaleUpperCase("cs") + text.slice(1) : text;
}

/** Substance names arrive in capitals ("METHYLFENIDÁT-HYDROCHLORID"); show them as one sentence-case name. */
export function substanceCase(text: string): string {
  if (!text || text !== text.toLocaleUpperCase("cs")) return text;
  return upperFirst(text.toLocaleLowerCase("cs"));
}
