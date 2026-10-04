<script setup lang="ts">
import type { PlSelectionModel } from "@platforma-sdk/model";
import { createPlSelectionModel } from "@platforma-sdk/model";
import {
  PlAgDataTableV2,
  PlAlert,
  PlBlockPage,
  PlBtnGhost,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { computed, ref, watch } from "vue";
import AddToBasketModal from "./AddToBasketModal.vue";
import { useApp } from "./app";
import { indexOfCurrent, keyedStatus } from "./keyedStatus";
import { selectedNodes } from "./baskets";
import { useAddToBasket } from "./useAddToBasket";
import { useLeaveWhenGone } from "./staleViews";

// One page per opened path, picked from `pathViews` by query id. Read only: the path was
// resolved on click.
const app = useApp<`/path?id=${string}` | "/trees">();

const index = computed(() =>
  indexOfCurrent(app.model.data.pathViews, app.queryParams.id, app.model.outputs.runKey),
);
const view = computed(() => (index.value < 0 ? undefined : app.model.data.pathViews[index.value]));
useLeaveWhenGone(index);

// Keeps the output's status so the table shows loading and errors instead of going blank.
const tableStatus = computed(() => keyedStatus(app.model.outputs.mutationalPaths, view.value?.id));
const table = computed(() => (tableStatus.value.ok ? tableStatus.value.value : undefined));

const settings = usePlDataTableSettingsV2({
  model: () => tableStatus.value,
  // Node ids repeat across lineages, so the section id, which names both, is the key.
  sourceId: () => view.value?.id,
});

// A path goes into a basket whole, or only the rows picked from it.
const selection = ref<PlSelectionModel>(createPlSelectionModel());
watch(
  () => view.value?.id,
  () => (selection.value = createPlSelectionModel()),
);
const selected = computed(() => {
  const columns = app.model.outputs.treeNodeColumns;
  const v = view.value;
  if (!columns || !v) return [];
  return selectedNodes(selection.value, columns)
    .filter((n) => n.lineageId === v.lineageId)
    .map((n) => n.nodeId);
});
const basket = useAddToBasket();
const addToBasket = () => {
  const v = view.value;
  if (!v) return;
  void basket.start(
    v.lineageId,
    v.lineageLabel,
    selected.value.length > 0 ? selected.value : v.nodeIds,
  );
};

// Drops only this section; the tree it came from stays open.
const close = async () => {
  const at = index.value;
  if (at < 0) return;
  app.model.data.pathViews.splice(at, 1);
  await app.navigateTo("/trees");
};
</script>

<template>
  <PlBlockPage>
    <template #title>
      {{ view ? `Path / ${view.lineageLabel} / ${view.nodeLabel}` : "Path" }}
    </template>
    <template v-if="view" #append>
      <PlBtnGhost icon="add" @click.stop="addToBasket">
        {{ selected.length > 0 ? "Add selected to basket" : "Add path to basket" }}
      </PlBtnGhost>
      <PlBtnGhost icon="close" @click.stop="close">Close</PlBtnGhost>
    </template>

    <PlAlert v-if="basket.error.value" type="error">{{ basket.error.value }}</PlAlert>
    <div v-if="!view">The path is loading.</div>
    <template v-else>
      <PlAgDataTableV2
        v-if="table"
        v-model="view.tableState"
        v-model:selection="selection"
        :settings="settings"
        show-export-button
        show-columns-panel
      />
    </template>
    <AddToBasketModal
      v-if="basket.pending.value.length > 0"
      :nodes="basket.pending.value"
      @close="basket.close"
    />
  </PlBlockPage>
</template>
