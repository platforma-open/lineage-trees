import type { BasketNode } from "@platforma-open/milaboratories.lineage-trees.model";
import { ref } from "vue";
import { useApp } from "./app";
import { snapshotNodes } from "./baskets";

/**
 * Backs every "Add to basket" button: copy the nodes, then open the dialog. The dialog's Add
 * is then a single write, and a copy failure shows beside the button.
 */
export function useAddToBasket() {
  const app = useApp();
  const pending = ref<BasketNode[]>([]);
  const error = ref<string | undefined>();

  const start = async (lineageId: string, lineageLabel: string, nodeIds: string[]) => {
    error.value = undefined;
    if (nodeIds.length === 0) return;
    const columns = app.model.outputs.treeNodeColumns;
    const runKey = app.model.outputs.runKey;
    const frameStatus = app.model.outputs.treeNodesPf;
    const frame = frameStatus?.ok === true ? frameStatus.value : undefined;
    if (!columns || !frame || runKey === undefined) {
      error.value = "The trees are still loading. Try again in a moment.";
      return;
    }
    try {
      pending.value = await snapshotNodes(frame, columns, runKey, lineageId, lineageLabel, nodeIds);
    } catch (caught) {
      error.value = caught instanceof Error ? caught.message : String(caught);
    }
  };

  return { pending, error, start, close: () => (pending.value = []) };
}
