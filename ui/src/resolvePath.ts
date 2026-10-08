import type { BlockOutputs } from "@platforma-open/milaboratories.lineage-trees.model";
import type { PFrameHandle } from "@platforma-sdk/model";
import { readLineageColumns } from "./lineageColumns";

/** What an output holds once it has settled, which is what the UI side sees. */
type Settled<T> = T extends { ok: true; value: infer V } ? V : never;

/** Taken from the model rather than restated, so the axis ids keep their types. */
export type TreeNodeColumns = NonNullable<Settled<BlockOutputs["treeNodeColumns"]>>;

export type ResolvedPath = { nodeIds: string[]; nodeLabel: string };

/** Node ids from the germline down to the node, and its label. Run on click; the page only reads. */
export async function resolvePath(
  frame: PFrameHandle,
  columns: TreeNodeColumns,
  lineageId: string,
  nodeId: string,
): Promise<ResolvedPath> {
  const wanted = columns.labelId ? [columns.topologyId, columns.labelId] : [columns.topologyId];
  const read = await readLineageColumns(frame, columns, lineageId, wanted);
  // The germline has no parent, so it gets no entry and a walk stops there.
  const parentOf = read.get(columns.topologyId);
  if (!parentOf) throw new Error("this tree's topology is not available");
  const labelOf = (columns.labelId && read.get(columns.labelId)) || new Map<string, string>();
  // The root has no parent and no label, but is some node's parent.
  const isRoot = !parentOf.has(nodeId) && [...parentOf.values()].includes(nodeId);
  if (!parentOf.has(nodeId) && !labelOf.has(nodeId) && !isRoot) {
    throw new Error("that node is not part of this lineage's tree");
  }

  const nodeIds: string[] = [];
  const seen = new Set<string>();
  let current = nodeId;
  // Guards against a corrupted tree with a cycle, which would hang the UI.
  while (!seen.has(current)) {
    nodeIds.push(current);
    seen.add(current);
    const parent = parentOf.get(current);
    if (parent === undefined) break;
    current = parent;
  }
  nodeIds.reverse();

  // Observed nodes are named by clonotype; inferred ones by their step on the path.
  const nodeLabel = labelOf.get(nodeId) || (isRoot ? "Root" : `Step ${nodeIds.length - 1}`);
  return { nodeIds, nodeLabel };
}
