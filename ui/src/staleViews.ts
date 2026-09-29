import { isCurrent } from "@platforma-open/milaboratories.lineage-trees.model";
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
