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

/** Default columns in display order, matched by name and label (heavy and light share a name). `last` puts wide ones at the end. */
const PATH_LEADING_COLUMNS: { name: string; label?: string; last?: boolean }[] = [
  { name: NODE_DEPTH_COLUMN },
  { name: "pl7.app/dendrogram/isObserved" },
  { name: "pl7.app/dendrogram/mutationCount", label: "Heavy #AA mutations" },
  { name: "pl7.app/dendrogram/mutationsAcquired", label: "Heavy AA mutations" },
  { name: "pl7.app/dendrogram/mutationCount", label: "Heavy #NT mutations" },
  { name: "pl7.app/dendrogram/mutationsAcquired", label: "Heavy NT mutations" },
  { name: "pl7.app/vdj/sequence", label: "CDR3 aa" },
  { name: "pl7.app/vdj/sequenceAlignment", label: "Heavy reconstructed sequence", last: true },
  { name: "pl7.app/vdj/sequenceAlignment", label: "Light reconstructed sequence", last: true },
];

// Selectors are regexes by default; exact match keeps "V gene" from matching "Light V gene".
const PATH_LEADING_SELECTORS = PATH_LEADING_COLUMNS.map((wanted) => ({
  name: { type: "exact" as const, value: wanted.name },
  ...(wanted.label === undefined
    ? {}
    : { annotations: { "pl7.app/label": { type: "exact" as const, value: wanted.label } } }),
}));

// Leading columns rank above the workflow's orderPriority, `last` ones below it.
const PATH_DISPLAY_OPTIONS = {
  ordering: PATH_LEADING_SELECTORS.map((match, i) => ({
    match,
    priority: (PATH_LEADING_COLUMNS[i].last ? -1_000_000 : 1_000_000) - i,
  })),
  visibility: [
    ...PATH_LEADING_SELECTORS.map((match) => ({ match, visibility: "default" as const })),
    { match: { name: { type: "regex" as const, value: ".*" } }, visibility: "optional" as const },
  ],
};

/** Columns split into leading and the rest. Undefined until a tree exists: no depth, no order. */
export function nodeTableParts(columns: PColumn<TreeNodeAccessor>[] | undefined) {
  if (columns === undefined || columns.length === 0) return undefined;
  const nodeScoped = columns.filter((column) => column.spec.axesSpec.length === 2);
  if (nodeScoped.length === 0) return undefined;

  const recipes = nodeScoped.map((column) => DataColumn.fromColumn(column));
  const depth = recipes.find((recipe) => recipe.getSpec().name === NODE_DEPTH_COLUMN);
  if (depth === undefined) return undefined;
  // Columns the run lacks are skipped.
  const leading = PATH_LEADING_COLUMNS.map((wanted) =>
    recipes.find((recipe) => {
      const spec = recipe.getSpec();
      return (
        spec.name === wanted.name &&
        (wanted.label === undefined || spec.annotations?.["pl7.app/label"] === wanted.label)
      );
    }),
  ).filter((recipe) => recipe !== undefined);
  const rest = recipes.filter((recipe) => !leading.includes(recipe));
  const axes = nodeScoped[0].spec.axesSpec;
  return { leading, rest, depth, lineageAxis: getAxisId(axes[0]), nodeAxis: getAxisId(axes[1]) };
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

/** Germline first. */
export const byDepth = (parts: NodeTableParts) => ({
  column: { type: "column" as const, id: parts.depth.id },
  ascending: true,
  naAndAbsentAreLeastValues: true,
});

/**
 * A node table with the path table's columns and ordering. `filters` is the view itself (its
 * lineage, its nodes), not a user filter: the SDK would show it in the filter panel as a default
 * filter, where clearing it shows every lineage. So a saved view's default filters are ignored,
 * and the panel gets none.
 */
export function nodeTable(
  ctx: Parameters<typeof createPlDataTableV3>[0],
  parts: NodeTableParts,
  tableState: TableOptions["tableState"],
  filters: TableOptions["filters"],
  sorting: TableOptions["sorting"] = [byDepth(parts)],
) {
  const state = upgradePlDataTableStateV2(tableState);
  const table = createPlDataTableV3(ctx, {
    primaryColumns: parts.leading,
    columns: parts.rest,
    tableState:
      state.pTableParams.sourceId === null
        ? state
        : { ...state, pTableParams: { ...state.pTableParams, defaultFilters: null } },
    displayOptions: PATH_DISPLAY_OPTIONS,
    filters,
    sorting,
  });
  return table === undefined ? undefined : { ...table, defaultFilters: undefined };
}
