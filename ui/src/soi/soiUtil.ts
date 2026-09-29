import type { Alphabet, TargetFeature } from "@platforma-open/milaboratories.lineage-trees.model";
import type { ListOption } from "@platforma-sdk/ui-vue";

export const alphabetOptions: ListOption<Alphabet>[] = [
  { value: "nucleotide", label: "Nucleotide" },
  { value: "amino-acid", label: "Amino acid" },
];

export const targetFeatureOptions: ListOption<TargetFeature>[] = [
  { value: "CDR3", label: "CDR3" },
  { value: "VDJRegion", label: "Full sequence" },
];

/** The first name from the constructor that is not taken. */
export function inferNewName(existingNames: string[], nameConstructor: (i: number) => string) {
  const names = new Set(existingNames);
  let i = 1;
  let name = nameConstructor(i);
  while (names.has(name)) {
    i++;
    name = nameConstructor(i);
  }
  return name;
}

/**
 * Sequences whose length misses the target feature; end-to-end search never hits them.
 * Loose bounds: a CDR3 is under 120 nt or 40 aa, a full VDJ region over 100 nt or 35 aa.
 */
export function lengthMismatches(
  sequences: { name: string; sequence: string }[],
  alphabet: Alphabet,
  targetFeature: TargetFeature,
): string[] {
  const nt = alphabet === "nucleotide";
  const out: string[] = [];
  for (const { name, sequence } of sequences) {
    const length = sequence.replace(/[\s.-]/g, "").length;
    if (targetFeature === "CDR3" && length > (nt ? 120 : 40)) {
      out.push(
        `${name || "(unnamed)"}: ${length} ${nt ? "nt" : "aa"} is longer than a CDR3; it belongs in a full-sequence list`,
      );
    } else if (targetFeature === "VDJRegion" && length < (nt ? 100 : 35)) {
      out.push(
        `${name || "(unnamed)"}: ${length} ${nt ? "nt" : "aa"} is too short for a full sequence; a CDR3 belongs in a CDR3 list`,
      );
    }
  }
  return out;
}
