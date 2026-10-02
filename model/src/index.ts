import type {
  AxisSpec,
  InferOutputsType,
  PColumnSpec,
  PlDataTableModel,
  PlDataTableStateV2,
  PObjectSpec,
  PlRef,
  ResultPool,
  TreeNodeAccessor,
} from "@platforma-sdk/model";
import {
  BlockModelV3,
  DataColumn,
  DataModelBuilder,
  canonicalizeAxisId,
  createPFrameForGraphs,
  createPlDataTableStateV2,
  createPlDataTableV3,
  getAxisId,
  getUniquePartitionKeys,
  isPColumnSpec,
  parseResourceMap,
} from "@platforma-sdk/model";
import type { GraphMakerState } from "@milaboratories/graph-maker";
import { kind } from "@platforma-open/milaboratories.lineage-trees.kind";
import type { SOIList } from "./soi";
import { describeRun } from "./runStatement";
import { byDepth, lineageFilter, nodeTable, nodeTableParts, nodesFilter } from "./nodeTables";

export * from "./soi";

// Default from HILARy
export const DEFAULT_CLUSTERING_THRESHOLD = 0.2;

/** What the tools prefix a progress line with; the last such line is what the UI shows. */
export const PROGRESS_PREFIX = "[==PROGRESS==]";

/** Member clonotypes per lineage, and the lineage table's default sort key. */
const LINEAGE_SIZE_COLUMN = "pl7.app/clustering/clusterSize";

/** Content id of the run's trees. Lineage and node ids are valid only within it. Undefined until a run settles. */
function runIdOf(outputs: TreeNodeAccessor | undefined): string | undefined {
  return outputs
    ?.resolve({ field: "runId", allowPermanentAbsence: true, stableIfNotFound: true })
    ?.getDataAsJson<{ runId: string }>()?.runId;
}

/** Whether a saved entry belongs to the result on show. Entries without a run key never match. */
export const isCurrent = (entry: { runKey?: string }, runKey: string | undefined) =>
  runKey !== undefined && entry.runKey === runKey;

/** Views opened on the result on show. */
function currentViews<V extends { runKey?: string }>(
  views: V[] | undefined,
  runKey: string | undefined,
): V[] {
  return (views ?? []).filter((view) => isCurrent(view, runKey));
}

// "Anchor" means a characterised antibody dataset. SDK anchors appear only as the bundle's "ds0", "ds1".

export type BlockData = {
  /** Upstream clonotyping datasets in picker order; `args` sorts them. */
  datasets: PlRef[];
  /** Anchor sets. Their members survive tree filters and set each member's tree distance. `args` drops refs not in `datasets`. */
  anchorDatasets: PlRef[];
  /** Sample metadata naming each sample's donor. Clustering never crosses donors. */
  donorColumn?: PlRef;
  /** `fixed`: single linkage per V, J and CDR3 length class. `adaptive`: HILARy's full method, human only. */
  clusteringMode: ClusteringMode;
  /** Fixed mode: a fraction of the CDR3 length. */
  clusteringThreshold: number;
  /** Adaptive mode: desired precision, HILARy's default 0.99. */
  precision: number;
  /** Adaptive mode: desired sensitivity, HILARy's default 0.9. */
  sensitivity: number;
  /** Larger lineages are subsampled to this, anchors kept. Left-out clonotypes lose only tree columns. */
  maxTipsPerTree?: number;
  /** Smaller lineages get no tree. Unset means the floor of two. */
  minTipsPerTree?: number;
  /** Lineages built with IgPhyML instead of FastTree+RAxML. More accurate but slow. */
  igPhyMLScope: IgPhyMLScope;
  /** Lists of sequences to find in the trees; only non-empty lists reach the workflow. */
  sequencesOfInterest: SOIList[];
  /** UI-only, never projected. */
  treesTableState: PlDataTableStateV2;
  /** One per opened tree section. UI-only, never projected. `tab` and `tableState` are absent on older projects. */
  treeViews: {
    id: string;
    lineageId: string;
    /** Run it was opened on. A rerun changes lineage ids and the lineage axis. */
    runKey?: string;
    state: GraphMakerState;
    tab?: "graph" | "table";
    tableState?: PlDataTableStateV2;
  }[];
  /** One per opened path section. Written only on a user gesture. UI-only, never projected. */
  pathViews: {
    id: string;
    lineageId: string;
    lineageLabel: string;
    /** Run it was opened on, as for tree views. */
    runKey?: string;
    nodeId?: string;
    nodeLabel?: string;
    nodeIds: string[];
    tableState: PlDataTableStateV2;
  }[];
  /** UI-only, never projected. */
  expansionGraphState: GraphMakerState;
  /** Named sets of collected nodes, one section each. UI-only, never projected. Absent on older projects. */
  baskets: NodeBasket[];
};

/** Ids point into the `runKey` run; the rest is copied so the entry survives a rerun. */
export type BasketNode = {
  lineageId: string;
  nodeId: string;
  runKey: string;
  lineageLabel: string;
  nodeLabel: string;
  heavySequence?: string;
  lightSequence?: string;
};

export type NodeBasket = {
  id: string;
  name: string;
  nodes: BasketNode[];
  tableState: PlDataTableStateV2;
};

export type ClusteringMode = "fixed" | "adaptive";

export const CLUSTERING_MODE_OPTIONS = [
  { value: "fixed", label: "Fixed threshold" },
  { value: "adaptive", label: "Adaptive, HILARy (human only)" },
] as const satisfies readonly { value: ClusteringMode; label: string }[];

/** HILARy's own defaults for its adaptive methods. */
export const DEFAULT_PRECISION = 0.99;
export const DEFAULT_SENSITIVITY = 0.9;

/** Data shape before several datasets were accepted. */
type BlockDataV1 = Omit<
  BlockData,
  | "datasets"
  | "anchorDatasets"
  | "igPhyMLScope"
  | "clusteringMode"
  | "precision"
  | "sensitivity"
  | "maxTipsPerTree"
  | "minTipsPerTree"
  | "sequencesOfInterest"
> & {
  inputAnchor?: PlRef;
  useIgPhyML: boolean;
  useLightChains: boolean;
  overviewTableState: PlDataTableStateV2;
};

export type IgPhyMLScope = "none" | "anchored" | "all";

export const IGPHYML_SCOPE_OPTIONS = [
  { value: "none", label: "None (FastTree and RAxML for every lineage)" },
  { value: "anchored", label: "Lineages holding an anchor" },
  { value: "all", label: "Every lineage" },
] as const satisfies readonly { value: IgPhyMLScope; label: string }[];

export type BlockArgs = {
  /** Sorted and deduplicated, so reordering the picker changes nothing. */
  datasets: PlRef[];
  /** Sorted, and only refs that are also in `datasets`. */
  anchorDatasets: PlRef[];
  donorColumn?: PlRef;
  clusteringMode: ClusteringMode;
  clusteringThreshold: number;
  precision: number;
  sensitivity: number;
  maxTipsPerTree?: number;
  minTipsPerTree?: number;
  /** Lists with at least one sequence, sorted by id. */
  sequencesOfInterest: SOIList[];
  igPhyMLScope: IgPhyMLScope;
};

/** Which shape of clonotyping output a dataset is. */
export type Modality = "bulk-heavy" | "paired-sc";

/** Where the tree step gets alignments: `mixcr` from clns, `upstream` from imported columns, `none` gives no trees. */
export type AlignmentSource = "mixcr" | "upstream" | "none";

/** Per donor group, what the overview shows beside its progress. */
export type DonorStats = { donor: string; clonotype_count: number; lineage_count: number };

/** Picked datasets with samples lacking a donor value, and whether no sample has one. */
export type SamplesWithoutDonor = {
  datasets: { dataset: string; missing: number; total: number }[];
  noneNamed: boolean;
};

/** What the workflow reports about each dataset of the last run. */
export type DatasetRun = {
  runId: string;
  modality: Modality;
  alignmentSource: AlignmentSource;
  isAnchor: boolean;
};

const RUN_ID_DOMAIN = "pl7.app/vdj/clonotypingRunId";

/** A dataset's alignment route from pool specs, for the UI before the run. Must mirror the workflow. */
export function alignmentRouteFor(anchor: PColumnSpec, specs: PObjectSpec[]): AlignmentSource {
  const sampleAxis = anchor.axesSpec[0];
  const clonotypeAxis = anchor.axesSpec[1];
  const runId = clonotypeAxis?.domain?.[RUN_ID_DOMAIN];
  const columns = specs.filter(isPColumnSpec);

  const perSample = (spec: PColumnSpec) =>
    spec.axesSpec.length === 1 &&
    spec.axesSpec[0]?.name === sampleAxis?.name &&
    spec.domain?.[RUN_ID_DOMAIN] === runId;
  const perClonotype = (spec: PColumnSpec) =>
    spec.axesSpec.length === 1 &&
    spec.axesSpec[0]?.name === clonotypeAxis?.name &&
    spec.axesSpec[0]?.domain?.[RUN_ID_DOMAIN] === runId;
  const nucleotide = (spec: PColumnSpec) => spec.domain?.["pl7.app/alphabet"] === "nucleotide";

  if (columns.some((spec) => perSample(spec) && spec.name === "mixcr.com/clns")) return "mixcr";
  const carries = (name: string, feature?: string) =>
    columns.some(
      (spec) =>
        perClonotype(spec) &&
        nucleotide(spec) &&
        spec.name === name &&
        (feature === undefined || spec.domain?.["pl7.app/vdj/feature"] === feature),
    );
  // Both are needed: a sequence without its germline is not an alignment.
  if (carries("pl7.app/vdj/sequenceAlignment") && carries("pl7.app/vdj/germlineAlignment")) {
    return "upstream";
  }
  return "none";
}

/** A dataset's modality, or undefined if not accepted. Only IG: T cells have no hypermutation. */
export function datasetModality(spec: PObjectSpec): Modality | undefined {
  if (!isPColumnSpec(spec)) return undefined;
  if (spec.annotations?.["pl7.app/isAnchor"] !== "true") return undefined;
  if (spec.axesSpec.length < 2) return undefined;
  if (spec.axesSpec[0]?.name !== "pl7.app/sampleId") return undefined;

  const keyAxis = spec.axesSpec[1];
  if (keyAxis === undefined) return undefined;
  if (keyAxis.name === "pl7.app/vdj/clonotypeKey") {
    return keyAxis.domain?.["pl7.app/vdj/chain"] === "IGHeavy" ? "bulk-heavy" : undefined;
  }
  if (keyAxis.name === "pl7.app/vdj/scClonotypeKey") {
    return keyAxis.domain?.["pl7.app/vdj/receptor"] === "IG" ? "paired-sc" : undefined;
  }
  return undefined;
}

/** Axis identity as a string: name plus sorted domain. Tells which datasets a donor column reaches. */
function axisKey(axis: AxisSpec): string {
  const id = getAxisId(axis);
  const domain = Object.entries(id.domain ?? {}).sort(([a], [b]) => a.localeCompare(b));
  return JSON.stringify([id.name, domain]);
}

/** One string per ref, for sets and maps; the workflow keys anchors the same way. */
export const refKey = (ref: PlRef) => `${ref.blockId}/${ref.name}`;

/** Picked anchor sets, deduplicated and sorted; none with one dataset. Shared by args and settings. */
export function effectiveAnchors(data: Pick<BlockData, "datasets" | "anchorDatasets">): PlRef[] {
  const datasets = canonicalRefs(data.datasets ?? []);
  if (datasets.length < 2) return [];
  const picked = new Set(datasets.map(refKey));
  return canonicalRefs(data.anchorDatasets ?? []).filter((ref) => picked.has(refKey(ref)));
}

function canonicalRefs(refs: PlRef[]): PlRef[] {
  const byKey = new Map<string, PlRef>();
  for (const ref of refs) byKey.set(refKey(ref), ref);
  return [...byKey.entries()].sort(([a], [b]) => a.localeCompare(b)).map(([, ref]) => ref);
}

/** What the UI needs to know about one picked dataset, from the pool alone. */
export type DatasetInfo = {
  ref: PlRef;
  label: string;
  runId: string | undefined;
  modality: Modality;
  alignmentRoute: AlignmentSource;
  sampleAxisKey: string;
  isAnchor: boolean;
};

/** A per-donor resource map among the workflow's outputs. */
const donorResourceMap = (outputs: TreeNodeAccessor | undefined, field: string) =>
  outputs?.resolve({ field, assertFieldType: "Input", allowPermanentAbsence: true });

/** One parsed value per donor. */
function donorMap<T>(
  outputs: TreeNodeAccessor | undefined,
  field: string,
  read: (acc: TreeNodeAccessor) => T | undefined,
) {
  return parseResourceMap(donorResourceMap(outputs, field), read, false);
}

/** A stage's per-donor live stdout, kept for donors whose log has not started. */
function stageStream<T>(
  outputs: TreeNodeAccessor | undefined,
  field: string,
  read: (acc: TreeNodeAccessor) => T | undefined,
) {
  return parseResourceMap(donorResourceMap(outputs, field), read, true);
}

/** Log handles for the log view. */
const stageLogs = (outputs: TreeNodeAccessor | undefined, field: string) =>
  stageStream(outputs, field, (acc) => acc.getLogHandle());

/** The last progress line per donor, for the bar. */
const stageProgress = (outputs: TreeNodeAccessor | undefined, field: string) =>
  stageStream(outputs, field, (acc) => acc.getProgressLog(PROGRESS_PREFIX));

const collectStream = (outputs: TreeNodeAccessor | undefined) =>
  outputs?.resolve({ field: "collectLog", assertFieldType: "Input", allowPermanentAbsence: true });

const readAlleleRoute = (acc: TreeNodeAccessor) =>
  acc.getDataAsJson<{ route: string; reason: string }>();
const readClusteringMethod = (acc: TreeNodeAccessor) =>
  acc.getDataAsJson<{ method: string; reason: string }>();

/** Clonotyping datasets: bulk and single-cell anchor columns keyed by sample. */
const DATASET_QUERY = [
  {
    axes: [{ name: "pl7.app/sampleId" }, { name: "pl7.app/vdj/clonotypeKey" }],
    annotations: { "pl7.app/isAnchor": "true" },
  },
  {
    axes: [{ name: "pl7.app/sampleId" }, { name: "pl7.app/vdj/scClonotypeKey" }],
    annotations: { "pl7.app/isAnchor": "true" },
  },
];

/** The picker's label for each dataset, keyed by `refKey`, so every message names it alike. */
function datasetLabels(resultPool: ResultPool): Map<string, string> {
  return new Map(resultPool.getOptions(DATASET_QUERY).map((o) => [refKey(o.ref), o.label]));
}

/** A new block's data; tests start from it too. A template can seed only `clusteringThreshold`. */
export function defaultBlockData(
  clusteringThreshold: number = DEFAULT_CLUSTERING_THRESHOLD,
): BlockData {
  return {
    datasets: [],
    anchorDatasets: [],
    igPhyMLScope: "none",
    clusteringMode: "fixed",
    precision: DEFAULT_PRECISION,
    sensitivity: DEFAULT_SENSITIVITY,
    sequencesOfInterest: [],
    clusteringThreshold,
    treesTableState: createPlDataTableStateV2(),
    treeViews: [],
    pathViews: [],
    baskets: [],
    expansionGraphState: {
      title: "Lineage expansion",
      template: "curve",
      // Log axes: rank against abundance is a power law.
      axesSettings: {
        axisX: { scale: "log", gridlines: false },
        axisY: { scale: "log", gridlines: true },
      },
    },
  };
}

const dataModel = new DataModelBuilder({ kind })
  .from<BlockDataV1>("v1")
  // v2: dataset becomes a list, IgPhyML switch a scope; two old fields dropped.
  .migrate<BlockData>(
    "v2",
    ({ inputAnchor, useIgPhyML, useLightChains: _l, overviewTableState: _, ...rest }) => ({
      ...rest,
      datasets: inputAnchor === undefined ? [] : [inputAnchor],
      anchorDatasets: [],
      igPhyMLScope: useIgPhyML ? "all" : "none",
      clusteringMode: "fixed",
      precision: DEFAULT_PRECISION,
      sensitivity: DEFAULT_SENSITIVITY,
      sequencesOfInterest: [],
    }),
  )
  .init(({ params }) => defaultBlockData(params?.clusteringThreshold));

export const platforma = BlockModelV3.create({ dataModel, kind })
  .args<BlockArgs>((data) => {
    const datasets = canonicalRefs(data.datasets ?? []);
    if (datasets.length === 0) throw new Error("Pick at least one input dataset");
    // Anchor refs no longer in `datasets` are dropped rather than failing the run.
    const anchorDatasets = effectiveAnchors(data);
    const clusteringMode = data.clusteringMode ?? "fixed";
    if (
      clusteringMode === "fixed" &&
      !(data.clusteringThreshold > 0 && data.clusteringThreshold < 1)
    ) {
      throw new Error("Clustering threshold must be a fraction between 0 and 1");
    }
    if (clusteringMode === "adaptive") {
      const fraction = (v: number) => v > 0 && v <= 1;
      if (
        !fraction(data.precision ?? DEFAULT_PRECISION) ||
        !fraction(data.sensitivity ?? DEFAULT_SENSITIVITY)
      ) {
        throw new Error("Precision and sensitivity must be fractions between 0 and 1");
      }
    }
    // Unset means no cap and a floor of two.
    const maxTipsPerTree = data.maxTipsPerTree;
    // The tree builders need three tips; below that no cap makes sense.
    if (
      maxTipsPerTree !== undefined &&
      !(Number.isInteger(maxTipsPerTree) && maxTipsPerTree >= 3)
    ) {
      throw new Error("Maximum tips per tree must be a whole number of at least 3");
    }
    const minTipsPerTree = data.minTipsPerTree;
    if (
      minTipsPerTree !== undefined &&
      !(
        Number.isInteger(minTipsPerTree) &&
        minTipsPerTree >= 2 &&
        (maxTipsPerTree === undefined || minTipsPerTree <= maxTipsPerTree)
      )
    ) {
      throw new Error("Minimum tips per tree must be a whole number between 2 and the maximum");
    }
    // Drop empty lists and sort, so adding or reordering lists does not stale the block.
    const sequencesOfInterest = (data.sequencesOfInterest ?? [])
      .filter((list) => list.sequences.length > 0)
      .map((list) => ({
        parameters: list.parameters,
        sequences: [...list.sequences].sort((a, b) => a.id.localeCompare(b.id)),
      }))
      .sort((a, b) => a.parameters.id.localeCompare(b.parameters.id));
    return {
      maxTipsPerTree,
      minTipsPerTree,
      sequencesOfInterest,
      // "anchored" without anchors builds no IgPhyML tree but labels distances as mixed.
      igPhyMLScope:
        data.igPhyMLScope === "anchored" && anchorDatasets.length === 0
          ? "none"
          : (data.igPhyMLScope ?? "none"),
      datasets,
      anchorDatasets,
      // Absent means one donor; the workflow treats it as a single group.
      donorColumn: data.donorColumn,
      clusteringMode,
      // Only the active mode's parameters pass, so editing the other mode does not stale.
      clusteringThreshold:
        clusteringMode === "fixed" ? data.clusteringThreshold : DEFAULT_CLUSTERING_THRESHOLD,
      precision:
        clusteringMode === "adaptive" ? (data.precision ?? DEFAULT_PRECISION) : DEFAULT_PRECISION,
      sensitivity:
        clusteringMode === "adaptive"
          ? (data.sensitivity ?? DEFAULT_SENSITIVITY)
          : DEFAULT_SENSITIVITY,
    };
  })
  // Inverse of `init`, so an exported template round-trips.
  .templateParams((data) => ({ clusteringThreshold: data.clusteringThreshold }))

  // The receptor is in the axis domain, which the query cannot express, so IG is checked here.
  .output("inputOptions", (ctx) =>
    ctx.resultPool.getOptions(DATASET_QUERY).filter((option) => {
      const spec = ctx.resultPool.getPColumnSpecByRef(option.ref);
      return spec !== undefined && datasetModality(spec) !== undefined;
    }),
  )

  /** Sample columns that can name a donor: `pl7.app/metadata` or `pl7.app/label`. */
  .output("donorOptions", (ctx) =>
    ctx.resultPool.getOptions((spec) => {
      if (!isPColumnSpec(spec)) return false;
      if (spec.axesSpec.length !== 1) return false;
      if (spec.axesSpec[0]?.name !== "pl7.app/sampleId") return false;
      return spec.name === "pl7.app/metadata" || spec.name === "pl7.app/label";
    }),
  )

  /** Picked datasets in `args` order; picks whose spec left the pool are dropped. */
  .output("datasets", (ctx): DatasetInfo[] => {
    const specs = ctx.resultPool.getSpecs().entries.map((entry) => entry.obj);
    const labels = datasetLabels(ctx.resultPool);
    const anchors = new Set((ctx.data.anchorDatasets ?? []).map(refKey));
    const out: DatasetInfo[] = [];
    for (const ref of canonicalRefs(ctx.data.datasets ?? [])) {
      const spec = ctx.resultPool.getPColumnSpecByRef(ref);
      if (spec === undefined) continue;
      const modality = datasetModality(spec);
      if (modality === undefined) continue;
      out.push({
        ref,
        label: labels.get(refKey(ref)) ?? ref.name,
        runId: spec.axesSpec[1]?.domain?.[RUN_ID_DOMAIN],
        modality,
        alignmentRoute: alignmentRouteFor(spec, specs),
        sampleAxisKey: axisKey(spec.axesSpec[0]),
        isAnchor: anchors.has(refKey(ref)),
      });
    }
    return out;
  })

  /** Datasets on another sample axis than the donor column. They get no donor and are not clustered. */
  .output("datasetsOutsideDonorColumn", (ctx): string[] => {
    if (ctx.data.donorColumn === undefined) return [];
    const donorSpec = ctx.resultPool.getPColumnSpecByRef(ctx.data.donorColumn);
    const donorAxis = donorSpec?.axesSpec[0];
    if (donorAxis === undefined) return [];
    const donorKey = axisKey(donorAxis);
    const labels = datasetLabels(ctx.resultPool);
    const out: string[] = [];
    for (const ref of canonicalRefs(ctx.data.datasets ?? [])) {
      const spec = ctx.resultPool.getPColumnSpecByRef(ref);
      if (spec === undefined || axisKey(spec.axesSpec[0]) === donorKey) continue;
      out.push(labels.get(refKey(ref)) ?? ref.name);
    }
    return out;
  })

  /**
   * Per picked dataset on the donor column's samples, how many of its samples have no donor
   * value (absent, null or blank). Those samples are left out of clustering; when no sample
   * has one, nothing can be clustered and the workflow stops with the same message.
   * Undefined while the donor column or a dataset's samples are not readable yet.
   */
  .output("samplesWithoutDonor", (ctx): SamplesWithoutDonor | undefined => {
    if (ctx.data.donorColumn === undefined) return { datasets: [], noneNamed: false };
    const donor = ctx.resultPool.getPColumnByRef(ctx.data.donorColumn);
    const donorAxis = donor?.spec.axesSpec[0];
    const values = donor?.data?.getDataAsJsonOrUndefined<{ data?: Record<string, unknown> }>()
      ?.data;
    if (donorAxis === undefined || values === undefined) return undefined;
    const named = (sample: string) => {
      const v = values[JSON.stringify([sample])];
      return v !== undefined && v !== null && String(v).trim() !== "";
    };
    const donorKey = axisKey(donorAxis);
    const labels = datasetLabels(ctx.resultPool);
    const out: SamplesWithoutDonor["datasets"] = [];
    let covered = 0;
    let namedSamples = 0;
    for (const ref of canonicalRefs(ctx.data.datasets ?? [])) {
      const column = ctx.resultPool.getPColumnByRef(ref);
      // Datasets on another sample axis are reported by datasetsOutsideDonorColumn.
      if (column === undefined || axisKey(column.spec.axesSpec[0]) !== donorKey) continue;
      const samples = getUniquePartitionKeys(column.data)?.[0];
      if (samples === undefined) return undefined;
      const missing = samples.filter((s) => !named(String(s))).length;
      covered++;
      namedSamples += samples.length - missing;
      if (missing > 0) {
        out.push({ dataset: labels.get(refKey(ref)) ?? ref.name, missing, total: samples.length });
      }
    }
    return { datasets: out, noneNamed: covered > 0 && namedSamples === 0 };
  })

  /** What ran. Reads `activeArgs` so it matches the results on screen, not unrun edits. */
  .output("modeStatement", (ctx) => {
    const args = ctx.activeArgs;
    if (args === undefined) return undefined;
    // Absent on projects run before this output existed.
    const runs = ctx.outputs
      ?.resolve({ field: "datasets", allowPermanentAbsence: true, stableIfNotFound: true })
      ?.getDataAsJson<DatasetRun[]>();
    if (runs === undefined) return undefined;
    const methods = donorMap(ctx.outputs, "clusteringMethods", readClusteringMethod).data.map(
      (entry) => entry.value.method,
    );
    const routes = donorMap(ctx.outputs, "alleleRoutes", readAlleleRoute).data.map(
      (entry) => entry.value.route,
    );
    return describeRun(args, runs, methods, routes);
  })

  /** Donors where every clonotype is its own lineage: no V, J and CDR3-length group held two. */
  .output("singletonDonors", (ctx): string[] =>
    donorMap(ctx.outputs, "clusteringMethods", readClusteringMethod)
      .data.filter((entry) => entry.value.method === "singletons")
      .map((entry) => String(entry.key[0]))
      .sort(),
  )

  /** Donors with no clonotypes in the picked datasets, once clustering has run. */
  .output("emptyDonors", (ctx): string[] =>
    (
      ctx.outputs
        ?.resolve({ field: "donorStats", allowPermanentAbsence: true, stableIfNotFound: true })
        ?.getDataAsJson<DonorStats[]>() ?? []
    )
      .filter((stats) => stats.clonotype_count === 0)
      .map((stats) => stats.donor)
      .sort(),
  )

  /** Per-donor clonotype and lineage counts, once clustering has run. */
  .output("donorStats", (ctx) =>
    ctx.outputs
      ?.resolve({ field: "donorStats", allowPermanentAbsence: true, stableIfNotFound: true })
      ?.getDataAsJson<DonorStats[]>(),
  )

  /** One row per lineage. */
  .outputWithStatus("treesTable", (ctx) => {
    const own = ctx.outputs?.resolve("treesTable")?.getPColumns();
    if (own === undefined || own.length === 0) return undefined;
    // Best hit per lineage for each sequence list.
    const soi = (
      ctx.outputs
        ?.resolve({ field: "soiTreesResults", allowPermanentAbsence: true, stableIfNotFound: true })
        ?.mapFields((_, v) => v?.getPColumns() ?? []) ?? []
    ).flat();
    const columns = [...own, ...soi];
    const recipes = columns.map((column) => DataColumn.fromColumn(column));
    // Largest lineages first by default.
    const size = recipes.find((recipe) => recipe.getSpec().name === LINEAGE_SIZE_COLUMN);
    // Size is primary: every lineage has one. The label column must not be primary.
    const primary = size ?? recipes[0];
    return createPlDataTableV3(ctx, {
      primaryColumns: [primary],
      columns: recipes.filter((recipe) => recipe !== primary),
      tableState: ctx.data.treesTableState,
      sorting:
        size === undefined
          ? undefined
          : [
              {
                column: { type: "column", id: size.id },
                ascending: false,
                naAndAbsentAreLeastValues: true,
              },
            ],
    });
  })

  // createPFrameForGraphs adds the related pool columns GraphMaker needs.
  .outputWithStatus("expansionPlot", (ctx) => {
    const columns = ctx.outputs?.resolve("expansionPlot")?.getPColumns();
    if (columns === undefined || columns.length === 0) return undefined;
    return createPFrameForGraphs(ctx, columns);
  })

  // Raw stdout, not JSON. getDataAsString gives undefined while empty.
  .output("clusteringLog", (ctx) => ctx.outputs?.resolve("clusteringLog")?.getDataAsString())

  // How many clonotypes and abundance rows each dataset brought to the merge.
  .output("mergeLog", (ctx) =>
    ctx.outputs
      ?.resolve({ field: "mergeLog", allowPermanentAbsence: true, stableIfNotFound: true })
      ?.getDataAsString(),
  )

  // Join counts, light chain split and skipped lineages.
  .output("treesLog", (ctx) => ctx.outputs?.resolve("treesLog")?.getDataAsString())

  /** Why no trees were built, read from the tree step's log. */
  .output("noTreesReason", (ctx) => {
    const log = ctx.outputs?.resolve("treesLog")?.getDataAsString();
    if (log === undefined || log === "") return undefined;
    const built = [...log.matchAll(/built (\d+) trees/g)].reduce(
      (total, match) => total + Number(match[1]),
      0,
    );
    if (built > 0) return undefined;
    const reasons = [...log.matchAll(/^(.+): no trees for this group$/gm)].map((m) => m[1]);
    return reasons.length > 0
      ? `No trees were built: ${[...new Set(reasons)].join("; ")}.`
      : "No trees were built.";
  })

  .outputWithStatus("treeNodesPf", (ctx) => {
    const columns = ctx.outputs?.resolve("treeNodes")?.getPColumns();
    if (columns === undefined || columns.length === 0) return undefined;
    // Sequence list hits per node.
    const soi = (
      ctx.outputs
        ?.resolve({ field: "soiNodesResults", allowPermanentAbsence: true, stableIfNotFound: true })
        ?.mapFields((_, v) => v?.getPColumns() ?? []) ?? []
    ).flat();
    return createPFrameForGraphs(ctx, [...columns, ...soi]);
  })

  /** Per sequence list, how many nodes and lineages the last run's search hit. */
  .output("soiHits", (ctx): Record<string, { nodes: number; lineages: number }> | undefined => {
    const acc = ctx.outputs?.resolve({
      field: "soiHitCounts",
      allowPermanentAbsence: true,
      stableIfNotFound: true,
    });
    if (acc === undefined) return undefined;
    const entries = acc.mapFields((id, counts) => {
      const value = counts?.getDataAsJson<{ nodes: number; lineages: number }>();
      return value === undefined ? undefined : ([id, value] as const);
    });
    return Object.fromEntries((entries ?? []).filter((e) => e !== undefined));
  })

  /** Whether every sequence search of the last run has finished. */
  .output("soiReady", (ctx) =>
    ctx.outputs
      ?.resolve({ field: "soiNodesResults", allowPermanentAbsence: true, stableIfNotFound: true })
      ?.getIsReadyOrError(),
  )

  // Lineage axis with domain, to match a clicked row. By name: the linker has the clonotype axis first.
  .output("lineageAxisSpec", (ctx) => {
    const columns = ctx.outputs?.resolve("treeNodes")?.getPColumns();
    return columns
      ?.flatMap((column) => column.spec.axesSpec)
      .find((axis) => axis.name === "pl7.app/clusterId");
  })

  /** Columns keyed on lineage and node, without the three-axis node-to-clonotype linker. */
  .output("treeNodeColumns", (ctx) => {
    const columns = ctx.outputs?.resolve("treeNodes")?.getPColumns();
    if (columns === undefined) return undefined;
    const nodeScoped = columns.filter((column) => column.spec.axesSpec.length === 2);
    if (nodeScoped.length === 0) return undefined;
    const topology = nodeScoped.find(
      (column) => column.spec.name === "pl7.app/dendrogram/topology",
    );
    if (topology === undefined) return undefined;
    const label = nodeScoped.find((column) => column.spec.name === "pl7.app/label");
    const sequenceOf = (chain: string) =>
      nodeScoped.find(
        (column) =>
          column.spec.name === "pl7.app/vdj/sequenceAlignment" &&
          column.spec.domain?.["pl7.app/vdj/chain"] === chain,
      );
    const heavySequence = sequenceOf("IGHeavy");
    const lightSequence = sequenceOf("IGLight");
    return {
      // For turning a node into a path: each node's parent, and its label.
      topologyId: topology.id,
      labelId: label?.id,
      lineageAxis: getAxisId(topology.spec.axesSpec[0]),
      nodeAxis: getAxisId(topology.spec.axesSpec[1]),
      // The tree page colours tips by anchor when set.
      hasAnchorProperty: nodeScoped.some(
        (column) => column.spec.name === "pl7.app/dendrogram/isAnchor",
      ),
      // Emitted only on runs over several datasets; the tree page colours tips by it.
      hasDatasetProperty: nodeScoped.some(
        (column) => column.spec.name === "pl7.app/dendrogram/dataset",
      ),
      // Only runs with light chains emit the light reconstructed sequence.
      hasLightSequence: lightSequence !== undefined,
      // What a basket copies from a node when it is added.
      heavySequenceId: heavySequence?.id,
      lightSequenceId: lightSequence?.id,
    };
  })

  /** See `runIdOf`. */
  .output("runKey", (ctx) => runIdOf(ctx.outputs))

  /** Node table per opened tree, keyed by view id. */
  .outputWithStatus("treeNodeTables", (ctx) => {
    const views = currentViews(ctx.data.treeViews, runIdOf(ctx.outputs));
    if (views.length === 0) return undefined;
    const parts = nodeTableParts(ctx.outputs?.resolve("treeNodes")?.getPColumns());
    if (parts === undefined) return undefined;

    const tables: Record<string, PlDataTableModel> = {};
    for (const view of views) {
      const table = nodeTable(ctx, parts, view.tableState ?? createPlDataTableStateV2(), {
        type: "and",
        filters: [lineageFilter(parts, view.lineageId)],
      });
      if (table !== undefined) tables[view.id] = table;
    }
    return tables;
  })

  /** Table per opened path, keyed by view id. A record: the model cannot see which section is shown. */
  .outputWithStatus("mutationalPaths", (ctx) => {
    const views = currentViews(ctx.data.pathViews, runIdOf(ctx.outputs)).filter(
      (view) => view.nodeIds.length > 0,
    );
    if (views.length === 0) return undefined;
    const parts = nodeTableParts(ctx.outputs?.resolve("treeNodes")?.getPColumns());
    if (parts === undefined) return undefined;

    const tables: Record<string, PlDataTableModel> = {};
    for (const view of views) {
      const table = nodeTable(ctx, parts, view.tableState, {
        type: "and",
        filters: [lineageFilter(parts, view.lineageId), nodesFilter(parts, view.nodeIds)],
      });
      if (table === undefined) continue;
      // Hide the lineage axis, the same on every row. Visibility rules match columns only.
      const axesMeta = table.columnsMeta?.axes;
      if (axesMeta !== undefined)
        axesMeta[canonicalizeAxisId(parts.lineageAxis)] = { hidden: true };
      tables[view.id] = table;
    }
    return tables;
  })

  /** Node table per basket, from its nodes in the run on show. None if it has no such node. */
  .outputWithStatus("basketTables", (ctx) => {
    const baskets = ctx.data.baskets ?? [];
    const runKey = runIdOf(ctx.outputs);
    if (baskets.length === 0 || runKey === undefined) return undefined;
    const parts = nodeTableParts(ctx.outputs?.resolve("treeNodes")?.getPColumns());
    if (parts === undefined) return undefined;

    const tables: Record<string, PlDataTableModel> = {};
    for (const basket of baskets) {
      // One branch per lineage: the node ids are only unique inside one.
      const byLineage = new Map<string, string[]>();
      for (const node of basket.nodes) {
        if (!isCurrent(node, runKey)) continue;
        const ids = byLineage.get(node.lineageId) ?? [];
        ids.push(node.nodeId);
        byLineage.set(node.lineageId, ids);
      }
      if (byLineage.size === 0) continue;
      const table = nodeTable(
        ctx,
        parts,
        basket.tableState ?? createPlDataTableStateV2(),
        {
          type: "or",
          filters: [...byLineage.entries()]
            .sort(([a], [b]) => a.localeCompare(b))
            .map(([lineageId, nodeIds]) => ({
              type: "and" as const,
              filters: [lineageFilter(parts, lineageId), nodesFilter(parts, [...nodeIds].sort())],
            })),
        },
        // Lineage by lineage, germline first inside each.
        [
          {
            column: { type: "axis", id: parts.lineageAxis },
            ascending: true,
            naAndAbsentAreLeastValues: true,
          },
          byDepth(parts),
        ],
      );
      if (table !== undefined) tables[basket.id] = table;
    }
    return tables;
  })

  /** Per-donor stage logs and progress. Absent on projects run before they existed. */
  .output("allelesLogs", (ctx) => stageLogs(ctx.outputs, "allelesLogs"))
  .output("allelesProgress", (ctx) => stageProgress(ctx.outputs, "allelesLogs"))
  /** Per donor, the route the allele stage took, "tigger" or "reference", and why. */
  .output("alleleRoutes", (ctx) => donorMap(ctx.outputs, "alleleRoutes", readAlleleRoute))
  // Alignment shows as part of allele inference; its log is read only for liveness.
  .output("alignmentsLogs", (ctx) => stageLogs(ctx.outputs, "alignmentsLogs"))
  .output("clusteringLogs", (ctx) => stageLogs(ctx.outputs, "clusteringLogs"))
  .output("clusteringProgress", (ctx) => stageProgress(ctx.outputs, "clusteringLogs"))
  .output("treesLogs", (ctx) => stageLogs(ctx.outputs, "treesLogs"))
  .output("treesProgress", (ctx) => stageProgress(ctx.outputs, "treesLogs"))
  /** The run-wide collect step after the last donor's trees; absent on older projects. */
  .output("collectLog", (ctx) => collectStream(ctx.outputs)?.getLogHandle())
  .output("collectProgress", (ctx) => collectStream(ctx.outputs)?.getProgressLog(PROGRESS_PREFIX))

  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .sections((ctx) => {
    const trees = currentViews(ctx.data.treeViews, runIdOf(ctx.outputs)).map((v) => ({
      type: "link" as const,
      href: `/tree?id=${v.id}` as const,
      label: v.state.title,
    }));
    // The lineage it belongs to, and the node it ends at once one is chosen.
    const paths = currentViews(ctx.data.pathViews, runIdOf(ctx.outputs)).map((v) => ({
      type: "link" as const,
      href: `/path?id=${v.id}` as const,
      label: v.nodeLabel ? `Path / ${v.lineageLabel} / ${v.nodeLabel}` : `Path / ${v.lineageLabel}`,
    }));
    const baskets = (ctx.data.baskets ?? []).map((b) => ({
      type: "link" as const,
      href: `/basket?id=${b.id}` as const,
      label: b.name,
    }));
    return [
      { type: "link" as const, href: "/" as const, label: "Overview" },
      { type: "link" as const, href: "/expansion" as const, label: "Lineage expansion plot" },
      { type: "link" as const, href: "/trees" as const, label: "Lineage table" },
      { type: "link" as const, href: "/soi" as const, label: "Sequence search" },
      ...(trees.length ? [{ type: "delimiter" as const }, ...trees] : []),
      ...(paths.length ? [{ type: "delimiter" as const }, ...paths] : []),
      ...(baskets.length ? [{ type: "delimiter" as const }, ...baskets] : []),
    ];
  })
  .done();

export type BlockOutputs = InferOutputsType<typeof platforma>;
