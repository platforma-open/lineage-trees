<script setup lang="ts">
import type { BasketNode } from "@platforma-open/milaboratories.lineage-trees.model";
import { isCurrent } from "@platforma-open/milaboratories.lineage-trees.model";
import type { PlSelectionModel } from "@platforma-sdk/model";
import { createPlSelectionModel } from "@platforma-sdk/model";
import {
  AgGridTheme,
  PlAgDataTableV2,
  PlAlert,
  PlBlockPage,
  PlBtnGhost,
  PlBtnPrimary,
  PlBtnSecondary,
  PlDialogModal,
  PlEditableTitle,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import type { ColDef } from "ag-grid-enterprise";
import { ClientSideRowModelModule, ModuleRegistry } from "ag-grid-enterprise";
import { AgGridVue } from "ag-grid-vue3";
import { computed, ref } from "vue";
import { useApp } from "./app";
import { indexById, keyedStatus } from "./keyedStatus";
import { selectedNodes } from "./baskets";

// One page per basket, picked from `baskets` by the query id. Writes only on user clicks.
const app = useApp<`/basket?id=${string}` | "/trees">();

const index = computed(() => indexById(app.model.data.baskets, app.queryParams.id));
const basket = computed(() => (index.value < 0 ? undefined : app.model.data.baskets[index.value]));

// Ids from another run may name other nodes now, so those nodes show their saved copy.
const runKey = computed(() => app.model.outputs.runKey);
const fromEarlierRun = computed(() =>
  runKey.value === undefined
    ? []
    : (basket.value?.nodes ?? []).filter((node) => !isCurrent(node, runKey.value)),
);
const inThisRun = computed(() =>
  runKey.value === undefined
    ? 0
    : (basket.value?.nodes ?? []).filter((node) => isCurrent(node, runKey.value)).length,
);

// Only clonotypes are exported: inferred nodes have none, and entries from before the export stored none.
const inferred = computed(
  () =>
    (basket.value?.nodes ?? []).filter(
      (node) => node.clonotypes !== undefined && Object.keys(node.clonotypes).length === 0,
    ).length,
);
const withoutKeys = computed(
  () => (basket.value?.nodes ?? []).filter((node) => node.clonotypes === undefined).length,
);

const tableStatus = computed(() => keyedStatus(app.model.outputs.basketTables, basket.value?.id));

const settings = usePlDataTableSettingsV2({
  model: () => tableStatus.value,
  sourceId: () => basket.value?.id,
});

const selection = ref<PlSelectionModel>(createPlSelectionModel());
const selected = computed(() => {
  const columns = app.model.outputs.treeNodeColumns;
  return columns ? selectedNodes(selection.value, columns) : [];
});

const removeSelected = () => {
  const b = basket.value;
  if (!b || runKey.value === undefined) return;
  const gone = new Set(selected.value.map((n) => `${n.lineageId}|${n.nodeId}`));
  b.nodes = b.nodes.filter(
    (node) => node.runKey !== runKey.value || !gone.has(`${node.lineageId}|${node.nodeId}`),
  );
  selection.value = createPlSelectionModel();
};

const removeEarlier = () => {
  const b = basket.value;
  if (!b || runKey.value === undefined) return;
  b.nodes = b.nodes.filter((node) => node.runKey === runKey.value);
};

const confirmDelete = ref(false);
const deleteBasket = async () => {
  confirmDelete.value = false;
  const at = index.value;
  if (at < 0) return;
  app.model.data.baskets.splice(at, 1);
  await app.navigateTo("/trees");
};

ModuleRegistry.registerModules([ClientSideRowModelModule]);

// The light column only when some copied node carries a light sequence.
const earlierColumns = computed<ColDef<BasketNode>[]>(() => [
  { field: "lineageLabel", headerName: "Lineage" },
  { field: "nodeLabel", headerName: "Node" },
  { field: "heavySequence", headerName: "Heavy reconstructed sequence", flex: 1 },
  ...(fromEarlierRun.value.some((n) => n.lightSequence)
    ? [{ field: "lightSequence" as const, headerName: "Light reconstructed sequence", flex: 1 }]
    : []),
]);
</script>

<template>
  <PlBlockPage>
    <template v-if="basket" #title>
      <PlEditableTitle v-model="basket.name" max-width="600px" :max-length="60" />
    </template>
    <template v-if="basket" #append>
      <PlBtnGhost v-if="selected.length > 0" icon="close" @click.stop="removeSelected">
        Remove selected
      </PlBtnGhost>
      <PlBtnGhost icon="delete-bin" @click.stop="confirmDelete = true">Delete basket</PlBtnGhost>
    </template>

    <div v-if="!basket">This basket no longer exists.</div>
    <template v-else>
      <div v-if="basket.nodes.length === 0">
        This basket is empty. Add nodes from a tree's table or from a mutational path.
      </div>
      <PlAlert v-else type="info">
        Once the block runs, downstream blocks such as lead selection can filter each dataset to
        this basket's clonotypes.
        <template v-if="inferred > 0">
          {{ inferred }} inferred {{ inferred === 1 ? "node has" : "nodes have" }} no clonotype and
          {{ inferred === 1 ? "is" : "are" }} not exported.
        </template>
        <template v-if="withoutKeys > 0">
          {{ withoutKeys }} {{ withoutKeys === 1 ? "node was" : "nodes were" }} added before baskets
          were exported and {{ withoutKeys === 1 ? "is" : "are" }} not exported until added again.
        </template>
      </PlAlert>
      <PlAgDataTableV2
        v-if="inThisRun > 0"
        v-model="basket.tableState"
        v-model:selection="selection"
        :settings="settings"
        show-export-button
        show-columns-panel
      />

      <template v-if="fromEarlierRun.length > 0">
        <PlAlert type="info">
          {{ fromEarlierRun.length }}
          {{ fromEarlierRun.length === 1 ? "node was" : "nodes were" }} added from an earlier run.
          The trees have been rebuilt since, so they are listed as they were when added.
        </PlAlert>
        <AgGridVue
          :theme="AgGridTheme"
          :style="{ height: '300px' }"
          :row-data="fromEarlierRun"
          :column-defs="earlierColumns"
          :default-col-def="{ sortable: true, resizable: true }"
        />
        <div>
          <PlBtnSecondary @click="removeEarlier">Remove nodes from earlier runs</PlBtnSecondary>
        </div>
      </template>
    </template>

    <PlDialogModal v-model="confirmDelete">
      <template #title>Delete "{{ basket?.name }}"?</template>
      <template #actions>
        <PlBtnPrimary @click="deleteBasket">Delete</PlBtnPrimary>
        <PlBtnSecondary @click="confirmDelete = false">Cancel</PlBtnSecondary>
      </template>
    </PlDialogModal>
  </PlBlockPage>
</template>
