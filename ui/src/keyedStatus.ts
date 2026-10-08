/** One id's entry from a status-carrying record output. Keeps the status so grids show loading. */
export function keyedStatus<T, E extends { ok: false }>(
  status: { ok: true; value: Record<string, T> | undefined; stable: boolean } | E | undefined,
  id: string | undefined,
): { ok: true; value: T | undefined; stable: boolean } | E {
  if (status === undefined) return { ok: true, value: undefined, stable: false };
  if (status.ok === false) return status;
  return {
    ok: true,
    value: id === undefined ? undefined : status.value?.[id],
    stable: status.stable,
  };
}

import { isCurrent } from "@platforma-open/milaboratories.lineage-trees.model";

/** Index of the view named by the query id, -1 once closed. */
export const indexById = (items: { id: string }[] | undefined, id: string | undefined) =>
  (items ?? []).findIndex((item) => item.id === id);

/** Like `indexById`, but a view from another run reads as closed: its ids may mean other nodes. */
export const indexOfCurrent = (
  items: { id: string; runKey?: string }[] | undefined,
  id: string | undefined,
  runKey: string | undefined,
) => (items ?? []).findIndex((item) => item.id === id && isCurrent(item, runKey));
