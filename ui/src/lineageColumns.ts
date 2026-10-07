import type { PFrameHandle, PObjectId } from "@platforma-sdk/model";
import { getRawPlatformaInstance, isValueNA } from "@platforma-sdk/model";
import type { TreeNodeColumns } from "./resolvePath";

/** One lineage's rows of the wanted columns. */
async function queryLineage(
  frame: PFrameHandle,
  columns: TreeNodeColumns,
  lineageId: string,
  wanted: PObjectId[],
) {
  const request = {
    src: {
      type: "full" as const,
      entries: wanted.map((column) => ({ type: "column" as const, column })),
    },
    filters: [
      {
        type: "bySingleColumnV2" as const,
        column: { type: "axis" as const, id: columns.lineageAxis },
        predicate: { operator: "Equal" as const, reference: lineageId },
      },
    ],
    sorting: [],
  };
  // The request is structure-cloned, and the axis id is a reactive proxy that cannot be.
  // A JSON round trip makes it plain data.
  return await getRawPlatformaInstance().pFrameDriver.calculateTableData(
    frame,
    JSON.parse(JSON.stringify(request)) as typeof request,
  );
}

/** One lineage's values per column, by node id. Missing columns and NA values are left out. */
export async function readLineageColumns(
  frame: PFrameHandle,
  columns: TreeNodeColumns,
  lineageId: string,
  wanted: PObjectId[],
): Promise<Map<PObjectId, Map<string, string>>> {
  const out = new Map<PObjectId, Map<string, string>>();
  if (wanted.length === 0) return out;
  const data = await queryLineage(frame, columns, lineageId, wanted);
  const nodeAxis = data.find(
    (entry) => entry.spec.type === "axis" && entry.spec.spec.name === columns.nodeAxis.name,
  );
  if (!nodeAxis) return out;
  for (const id of wanted) {
    const entry = data.find((e) => e.spec.type === "column" && e.spec.id === id);
    if (!entry) continue;
    const values = new Map<string, string>();
    for (let row = 0; row < nodeAxis.data.data.length; row++) {
      if (!isValueNA(entry.data, row)) {
        values.set(String(nodeAxis.data.data[row]), String(entry.data.data[row]));
      }
    }
    out.set(id, values);
  }
  return out;
}

/** One lineage's clonotype keys by node id, from a node-to-clonotype link column. */
export async function readLineageLinks(
  frame: PFrameHandle,
  columns: TreeNodeColumns,
  lineageId: string,
  link: TreeNodeColumns["clonotypeLinks"][number],
): Promise<Map<string, string[]>> {
  const data = await queryLineage(frame, columns, lineageId, [link.id]);
  const axis = (name: string) =>
    data.find((entry) => entry.spec.type === "axis" && entry.spec.spec.name === name);
  const nodeAxis = axis(columns.nodeAxis.name);
  const clonotypeAxis = axis(link.clonotypeAxis);
  const out = new Map<string, string[]>();
  if (!nodeAxis || !clonotypeAxis) return out;
  for (let row = 0; row < nodeAxis.data.data.length; row++) {
    const node = String(nodeAxis.data.data[row]);
    out.set(node, [...(out.get(node) ?? []), String(clonotypeAxis.data.data[row])]);
  }
  return out;
}
