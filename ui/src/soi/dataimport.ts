// TSV or CSV reader, with a header row, for sequence import. Text only, to avoid an xlsx dep.

export type ImportDataColumn = { readonly header: string };
export type ImportDataRow = (string | undefined)[];
export type ImportData = { readonly columns: ImportDataColumn[]; readonly rows: ImportDataRow[] };

// A field may be quoted, and a quoted field may hold the delimiter or a doubled quote.
function splitLine(line: string, delimiter: string): string[] {
  const out: string[] = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (quoted) {
      if (c === '"') {
        if (line[i + 1] === '"') {
          field += '"';
          i++;
        } else {
          quoted = false;
        }
      } else {
        field += c;
      }
    } else if (c === '"') {
      quoted = true;
    } else if (c === delimiter) {
      out.push(field);
      field = "";
    } else {
      field += c;
    }
  }
  out.push(field);
  return out;
}

export function readFileForImport(data: Uint8Array, fileName: string): ImportData {
  const text = new TextDecoder().decode(data).replace(/^﻿/, "");
  const lines = text.split(/\r?\n/);
  const delimiter = fileName.toLowerCase().endsWith(".csv") ? "," : "\t";
  const header = lines[0] === undefined ? [] : splitLine(lines[0], delimiter).map((h) => h.trim());
  const columns: ImportDataColumn[] = header.map((h, i) => ({ header: h || `Column ${i + 1}` }));
  const rows: ImportDataRow[] = [];
  for (const line of lines.slice(1)) {
    if (line.trim() === "") continue;
    const cells = splitLine(line, delimiter).map((c) => c.trim());
    rows.push(
      columns.map((_, i) => (cells[i] === "" || cells[i] === undefined ? undefined : cells[i])),
    );
  }
  return { columns, rows };
}
