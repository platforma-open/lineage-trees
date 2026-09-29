<script setup lang="ts">
import type { ColDef, GridApi, GridOptions } from "ag-grid-enterprise";
import { AgGridVue } from "ag-grid-vue3";
import type { SequenceOfInterest } from "@platforma-open/milaboratories.lineage-trees.model";
import { AgGridTheme } from "@platforma-sdk/ui-vue";
import { ref, watch } from "vue";

// The sequences of one list, editable in place.
const model = defineModel<SequenceOfInterest[]>({ required: true });

const columnDefs: ColDef<SequenceOfInterest>[] = [
  { colId: "name", field: "name", headerName: "Name", editable: true },
  { colId: "sequence", field: "sequence", headerName: "Sequence", editable: true, flex: 1 },
];

const gridApi = ref<GridApi<SequenceOfInterest>>();

const gridOptions: GridOptions<SequenceOfInterest> = {
  getRowId: (row) => String(row.data.id),
  onGridReady: (p) => (gridApi.value = p.api),
};

watch(
  () => [model.value, gridApi.value] as const,
  ([, api]) => api?.autoSizeColumns(["name"]),
  { deep: true },
);
</script>

<template>
  <AgGridVue
    :theme="AgGridTheme"
    :style="{ height: '100%' }"
    :rowData="model"
    :columnDefs="columnDefs"
    :grid-options="gridOptions"
  />
</template>
