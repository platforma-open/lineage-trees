import { isCurrent } from "@platforma-open/milaboratories.lineage-trees.model";
import { computed, watch, type Ref } from "vue";
import { useApp } from "./app";

/** Drops views opened on another run. Call on a user gesture; writes only if something drops. */
export function useDropStaleViews() {
  const app = useApp();
  return (runKey: string | undefined) => {
    if (runKey === undefined) return;
    const current = (view: { runKey?: string }) => isCurrent(view, runKey);
    const trees = app.model.data.treeViews ?? [];
    if (!trees.every(current)) app.model.data.treeViews = trees.filter(current);
    const paths = app.model.data.pathViews ?? [];
    if (!paths.every(current)) app.model.data.pathViews = paths.filter(current);
  };
}

/**
 * Leaves for the lineage table once the page's view is gone for good: closed, or from another
 * run. Until a run settles there is no run key, so a view only looks gone and the page waits.
 * Moves only this client's page; no data is written.
 */
export function useLeaveWhenGone(index: Ref<number>) {
  const app = useApp();
  const gone = computed(
    () =>
      index.value < 0 &&
      (app.model.outputs.runKey !== undefined || app.model.outputs.runFinished === true),
  );
  watch(
    gone,
    (isGone) => {
      if (isGone) void app.navigateTo("/trees");
    },
    { immediate: true },
  );
}
