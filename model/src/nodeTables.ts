// Node tables for a lineage, a path or a basket, with shared leading columns and order.
import type { PColumn, TreeNodeAccessor } from "@platforma-sdk/model";
import {
  DataColumn,
  createPlDataTableV3,
  getAxisId,
  upgradePlDataTableStateV2,
} from "@platforma-sdk/model";

/** Distance from the germline in nodes; orders path tables. */
export const NODE_DEPTH_COLUMN = "pl7.app/dendrogram/nodeDepth";

type Leading = { name: string; label?: string; last?: boolean }[];

/** Default columns in display order per table, matched by name and label (heavy and light share a name). `last` puts wide ones at the end. */
const LEADING_COLUMNS = {
  // The tree's own table lists observed sequences only, so it leads with what was observed.
  tree: [
    { name: "pl7.app/vdj/geneHitWithAllele", label: "V allele" },
    { name: "pl7.app/vdj/geneHitWithAllele", label: "D allele" },
    { name: "pl7.app/vdj/geneHitWithAllele", label: "J allele" },
    { name: "pl7.app/vdj/sequence", label: "CDR3 aa" },
    { name: "pl7.app/vdj/sequence", label: "VDJRegion nt", last: true },
  ],
  // Paths and baskets hold inferred nodes too: steps and changes since the parent and the MRCA.
  path: [
    { name: NODE_DEPTH_COLUMN },
    { name: "pl7.app/dendrogram/isObserved" },
    { name: "pl7.app/dendrogram/mutationCount", label: "Heavy #AA mutations" },
    { name: "pl7.app/dendrogram/mutationsAcquired", label: "Heavy AA mutations" },
    { name: "pl7.app/dendrogram/mutationCount", label: "Heavy #NT mutations" },
    { name: "pl7.app/dendrogram/mutationsAcquired", label: "Heavy NT mutations" },
    { name: "pl7.app/dendrogram/mutationCount", label: "Heavy #AA mutations from MRCA" },
    { name: "pl7.app/dendrogram/mutationsAcquired", label: "Heavy AA mutations from MRCA" },
    { name: "pl7.app/dendrogram/mutationCount", label: "Heavy #NT mutations from MRCA" },
    { name: "pl7.app/dendrogram/mutationsAcquired", label: "Heavy NT mutations from MRCA" },
    { name: "pl7.app/vdj/sequence", label: "CDR3 aa" },
    { name: "pl7.app/vdj/sequenceAlignment", label: "Heavy reconstructed sequence", last: true },
    { name: "pl7.app/vdj/sequenceAlignment", label: "Light reconstructed sequence", last: true },
  ],
} satisfies Record<string, Leading>;
export type NodeTableKind = keyof typeof LEADING_COLUMNS;

// Selectors are regexes by default; exact match keeps "V gene" from matching "Light V gene".
const selectorsOf = (leading: Leading) =>
  leading.map((wanted) => ({
    name: { type: "exact" as const, value: wanted.name },
    ...(wanted.label === undefined
      ? {}
      : { annotations: { "pl7.app/label": { type: "exact" as const, value: wanted.label } } }),
  }));

// Leading columns rank above the workflow's orderPriority, `last` ones below it; the rest start hidden.
const displayOptionsOf = (leading: Leading) => {
  const selectors = selectorsOf(leading);
  return {
    ordering: selectors.map((match, i) => ({
      match,
      priority: (leading[i].last ? -1_000_000 : 1_000_000) - i,
    })),
    visibility: [
      ...selectors.map((match) => ({ match, visibility: "default" as const })),
      { match: { name: { type: "regex" as const, value: ".*" } }, visibility: "optional" as const },
    ],
  };
};

/** Columns split into leading and the rest. Undefined until a tree exists: no depth, no order. */
export function nodeTableParts(columns: PColumn<TreeNodeAccessor>[] | undefined) {
  if (columns === undefined || columns.length === 0) return undefined;
  const nodeScoped = columns.filter((column) => column.spec.axesSpec.length === 2);
  if (nodeScoped.length === 0) return undefined;

  const recipes = nodeScoped.map((column) => DataColumn.fromColumn(column));
  const depth = recipes.find((recipe) => recipe.getSpec().name === NODE_DEPTH_COLUMN);
  if (depth === undefined) return undefined;
  const observed = recipes.find(
    (recipe) => recipe.getSpec().name === "pl7.app/dendrogram/isObserved",
  );
  const axes = nodeScoped[0].spec.axesSpec;
  return {
    recipes,
    depth,
    observed,
    lineageAxis: getAxisId(axes[0]),
    nodeAxis: getAxisId(axes[1]),
  };
}

/** A table's leading columns, those the run has, and the rest. */
function split(parts: NodeTableParts, kind: NodeTableKind) {
  const leading = LEADING_COLUMNS[kind]
    .map((wanted: Leading[number]) =>
      parts.recipes.find((recipe) => {
        const spec = recipe.getSpec();
        return (
          spec.name === wanted.name &&
          (wanted.label === undefined || spec.annotations?.["pl7.app/label"] === wanted.label)
        );
      }),
    )
    .filter((recipe) => recipe !== undefined);
  // A run without any of the kind's columns still needs a primary column.
  if (leading.length === 0) leading.push(parts.depth);
  return { leading, rest: parts.recipes.filter((recipe) => !leading.includes(recipe)) };
}

export type NodeTableParts = NonNullable<ReturnType<typeof nodeTableParts>>;
type TableOptions = Parameters<typeof createPlDataTableV3>[1];

/** The nodes of one lineage. */
export const lineageFilter = (parts: NodeTableParts, lineageId: string) => ({
  type: "patternEquals" as const,
  column: { type: "axis" as const, id: parts.lineageAxis },
  value: lineageId,
});

/** Nodes by id. Ids are unique only within a lineage. */
export const nodesFilter = (parts: NodeTableParts, nodeIds: string[]) => ({
  type: "inSet" as const,
  column: { type: "axis" as const, id: parts.nodeAxis },
  value: nodeIds,
});

/** Observed nodes only, for the tree's own table. No filter if the run has no such column. */
export const observedFilter = (parts: NodeTableParts) =>
  parts.observed === undefined
    ? []
    : [
        {
          type: "patternEquals" as const,
          column: { type: "column" as const, id: parts.observed.id },
          value: "true",
        },
      ];

/** Germline first. */
export const byDepth = (parts: NodeTableParts) => ({
  column: { type: "column" as const, id: parts.depth.id },
  ascending: true,
  naAndAbsentAreLeastValues: true,
});

/**
 * A node table with its kind's default columns and ordering. `filters` is the view itself (its
 * lineage, its nodes), not a user filter: the SDK would show it in the filter panel as a default
 * filter, where clearing it shows every lineage. So a saved view's default filters are ignored,
 * and the panel gets none.
 */
export function nodeTable(
  ctx: Parameters<typeof createPlDataTableV3>[0],
  parts: NodeTableParts,
  kind: NodeTableKind,
  tableState: TableOptions["tableState"],
  filters: TableOptions["filters"],
  sorting: TableOptions["sorting"] = [byDepth(parts)],
) {
  const state = upgradePlDataTableStateV2(tableState);
  const { leading, rest } = split(parts, kind);
  const table = createPlDataTableV3(ctx, {
    primaryColumns: leading,
    columns: rest,
    tableState:
      state.pTableParams.sourceId === null
        ? state
        : { ...state, pTableParams: { ...state.pTableParams, defaultFilters: null } },
    displayOptions: displayOptionsOf(LEADING_COLUMNS[kind]),
    filters,
    sorting,
  });
  return table === undefined ? undefined : { ...table, defaultFilters: undefined };
}
