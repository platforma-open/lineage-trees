<script setup lang="ts">
import { AgGridVue } from "ag-grid-vue3";
import type { ColDef, GridOptions } from "ag-grid-enterprise";
import { ClientSideRowModelModule, ModuleRegistry } from "ag-grid-enterprise";
import type { PlAgHeaderComponentParams } from "@platforma-sdk/ui-vue";
import {
  AgGridTheme,
  PlAgOverlayLoading,
  PlAgOverlayNoRows,
  PlAgTextAndButtonCell,
  createAgGridColDef,
  makeRowNumberColDef,
} from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "./app";
import type { DonorRow } from "./donorRows";
import { STAGES, useDonorRows } from "./donorRows";

// One row per donor group. Double-clicking a row opens its logs.
const app = useApp();

const emit = defineEmits<{ open: [donor: string] }>();
const donorRows = useDonorRows();

ModuleRegistry.registerModules([ClientSideRowModelModule]);

const defaultColDef: ColDef = {
  suppressHeaderMenuButton: true,
  lockPinned: true,
  sortable: false,
};

const columnDefs: ColDef<DonorRow>[] = [
  makeRowNumberColDef(),
  createAgGridColDef<DonorRow, string>({
    colId: "donor",
    field: "donor",
    headerName: "Donor",
    headerComponentParams: { type: "Text" } satisfies PlAgHeaderComponentParams,
    pinned: "left",
    lockPinned: true,
    sortable: true,
    cellRenderer: PlAgTextAndButtonCell,
    cellRendererParams: { invokeRowsOnDoubleClick: true },
  }),
  ...STAGES.map(({ key, label }) =>
    createAgGridColDef<DonorRow, string>({
      colId: `${key}Progress`,
      headerName: label,
      headerComponentParams: { type: "Progress" } satisfies PlAgHeaderComponentParams,
      flex: 1,
      minWidth: 160,
      valueGetter: (p) => p.data?.progress[key].text,
      progress(_value, cellData) {
        const stage = cellData.data?.progress[key];
        if (stage === undefined) return { status: "not_started", text: "Queued" };
        return {
          status: stage.status,
          percent: stage.percent,
          text: stage.text,
          suffix: stage.percent === undefined ? "" : `${stage.percent.toFixed(0)}%`,
        };
      },
    }),
  ),
  createAgGridColDef<DonorRow, number | undefined>({
    colId: "clonotypes",
    headerName: "Clonotypes",
    headerComponentParams: { type: "Number" } satisfies PlAgHeaderComponentParams,
    width: 130,
    valueGetter: (p) => p.data?.clonotypes,
    valueFormatter: (p) => (p.value == null ? "" : p.value.toLocaleString()),
  }),
  createAgGridColDef<DonorRow, number | undefined>({
    colId: "lineages",
    headerName: "Lineages",
    headerComponentParams: { type: "Number" } satisfies PlAgHeaderComponentParams,
    width: 130,
    valueGetter: (p) => p.data?.lineages,
    valueFormatter: (p) => (p.value == null ? "" : p.value.toLocaleString()),
  }),
];

const gridOptions: GridOptions<DonorRow> = {
  getRowId: (row) => row.data.donor,
  onRowDoubleClicked: (e) => {
    if (e.data) emit("open", e.data.donor);
  },
  components: { PlAgTextAndButtonCell },
};

const loadingOverlayParams = computed(() =>
  app.model.outputs.isRunning
    ? { variant: "running" as const, runningText: "Starting" }
    : { variant: "not-ready" as const },
);
</script>

<template>
  <AgGridVue
    :theme="AgGridTheme"
    :style="{ height: '100%' }"
    :rowData="donorRows"
    :defaultColDef="defaultColDef"
    :columnDefs="columnDefs"
    :grid-options="gridOptions"
    :loadingOverlayComponentParams="loadingOverlayParams"
    :loadingOverlayComponent="PlAgOverlayLoading"
    :noRowsOverlayComponent="PlAgOverlayNoRows"
  />
</template>
