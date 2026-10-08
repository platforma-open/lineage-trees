import type { AlignmentSource, BlockArgs, DatasetRun } from "./index";

/** "a, b and c". */
export function joinLabels(labels: string[]): string {
  if (labels.length <= 1) return labels.join("");
  return `${labels.slice(0, -1).join(", ")} and ${labels[labels.length - 1]}`;
}

/** The results page's summary of what ran. `methods` and `routes` are per donor, as the workflow reports. */
export function describeRun(
  args: BlockArgs,
  runs: DatasetRun[],
  methods: string[],
  routes: (string | undefined)[],
): string {
  const builder =
    args.igPhyMLScope === "all"
      ? "IgPhyML HLP19"
      : args.igPhyMLScope === "known"
        ? "IgPhyML HLP19 for lineages holding a known antibody and FastTree+RAxML GTR elsewhere"
        : "FastTree+RAxML GTR";
  const knownSets = runs.filter((run) => run.isKnown).length;
  const knownPart =
    knownSets === 0
      ? ""
      : ` ${knownSets === 1 ? "1 dataset" : `${knownSets} datasets`} as known antibody datasets.`;
  const single = runs.some((run) => run.modality === "paired-sc");
  const bulk = runs.some((run) => run.modality === "bulk-heavy");
  const data =
    single && bulk ? "Single cell and bulk" : single ? "Paired single cell" : "Bulk heavy chain";
  const cdr3Only = methods.filter((m) => m === "cdr3").length;
  const clustering =
    args.clusteringMode === "adaptive"
      ? "Lineages by HILARy's adaptive clustering with the shared-mutation test" +
        (cdr3Only > 0
          ? `, CDR3 only for ${cdr3Only === 1 ? "1 donor" : `${cdr3Only} donors`} with clonotypes lacking alignments`
          : "")
      : `Lineages from heavy chain CDR3 by single linkage at ${args.clusteringThreshold} of the CDR3 length`;
  // Light chains are used wherever a dataset carries them, known antibody datasets included.
  const lineages = single
    ? `${clustering}, split by light chain V, J and junction length;` +
      ` trees over heavy and light with ${builder} and a separate model per chain.`
    : `${clustering}; trees with ${builder}.`;
  const count = (source: AlignmentSource) =>
    runs.filter((run) => run.alignmentSource === source).length;
  const datasets = (n: number) => (n === 1 ? "1 dataset" : `${n} datasets`);
  const alignments: string[] = [];
  if (count("mixcr") > 0) alignments.push(`${datasets(count("mixcr"))} from the MiXCR clns files`);
  if (count("upstream") > 0) {
    alignments.push(`${datasets(count("upstream"))} on its own annotation and alignments`);
  }
  if (count("none") > 0) alignments.push(`${datasets(count("none"))} without any, so no trees`);
  // Only the allele stage knows if a donor supported inference, so use its report.
  const byRoute = (route: string) => routes.filter((r) => r === route).length;
  const inferredCount = routes.length - byRoute("reference");
  const alleles =
    routes.length === 0
      ? "V alleles are inferred for every donor with TIgGER; which donors it could infer for is reported once the run reaches that stage."
      : inferredCount === 0
        ? "Every donor fell back to the reference alleles; the allele log says why for each."
        : byRoute("reference") === 0
          ? `V alleles inferred with TIgGER for every donor, pooled over its datasets.`
          : `V alleles inferred with TIgGER for ${inferredCount} of ${routes.length} donors; reference alleles for the rest.`;
  return `${data}.${knownPart} ${lineages} Alignments: ${joinLabels(alignments)}. ${alleles}`;
}
