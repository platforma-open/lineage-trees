import type { BasketNode } from "@platforma-open/milaboratories.lineage-trees.model";
import type { PFrameHandle, PlSelectionModel, PObjectId } from "@platforma-sdk/model";
import { readLineageColumns } from "./lineageColumns";
import type { TreeNodeColumns } from "./resolvePath";

/** One entry per node and run: the same node added twice is kept once. */
export const basketNodeKey = (node: Pick<BasketNode, "runKey" | "lineageId" | "nodeId">) =>
  `${node.runKey}|${node.lineageId}|${node.nodeId}`;

export function mergeBasketNodes(into: BasketNode[], added: BasketNode[]): BasketNode[] {
  const seen = new Set(into.map(basketNodeKey));
  const merged = [...into];
  for (const node of added) {
    const key = basketNodeKey(node);
    if (seen.has(key)) continue;
    seen.add(key);
    merged.push(node);
  }
  return merged;
}

/** Lineage and node ids of the selected rows. Key parts follow `axesSpec` order. */
export function selectedNodes(
  selection: PlSelectionModel,
  columns: TreeNodeColumns,
): { lineageId: string; nodeId: string }[] {
  const lineageAt = selection.axesSpec.findIndex((a) => a.name === columns.lineageAxis.name);
  const nodeAt = selection.axesSpec.findIndex((a) => a.name === columns.nodeAxis.name);
  if (lineageAt < 0 || nodeAt < 0) return [];
  return selection.selectedKeys.flatMap((key) => {
    const lineage = key[lineageAt];
    const node = key[nodeAt];
    if (lineage == null || node == null) return [];
    return [{ lineageId: String(lineage), nodeId: String(node) }];
  });
}

/** Basket entries with a copy of each node's label and sequences, so they survive a rerun. */
export async function snapshotNodes(
  frame: PFrameHandle,
  columns: TreeNodeColumns,
  runKey: string,
  lineageId: string,
  lineageLabel: string,
  nodeIds: string[],
): Promise<BasketNode[]> {
  const wanted = [columns.labelId, columns.heavySequenceId, columns.lightSequenceId].filter(
    (id): id is PObjectId => id !== undefined,
  );
  const read = await readLineageColumns(frame, columns, lineageId, wanted);
  const valueOf = (id: PObjectId | undefined, nodeId: string) =>
    id === undefined ? undefined : read.get(id)?.get(nodeId);
  return nodeIds.map((nodeId) => {
    return {
      lineageId,
      nodeId,
      runKey,
      lineageLabel,
      // An observed node is named by its clonotype; an inferred one has no name.
      nodeLabel: valueOf(columns.labelId, nodeId) ?? `Inferred node ${nodeId}`,
      heavySequence: valueOf(columns.heavySequenceId, nodeId),
      lightSequence: valueOf(columns.lightSequenceId, nodeId),
    };
  });
}
