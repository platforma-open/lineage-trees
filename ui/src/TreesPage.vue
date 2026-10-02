<script setup lang="ts">
import type { AxisId, PTableKey } from "@platforma-sdk/model";
import {
  PlAgDataTableV2,
  PlBlockPage,
  PlBtnGhost,
  PlMaskIcon24,
  PlSlideModal,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { computed, reactive } from "vue";
import SettingsPanel from "./SettingsPanel.vue";
import { useApp } from "./app";
import { useDropStaleViews } from "./staleViews";

const app = useApp();
const dropStaleViews = useDropStaleViews();

const view = reactive({ settingsOpen: false });

const treesSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.treesTable,
});

// Matched with whole-object equality, so the domain has to come from the emitted spec.
const lineageAxis = computed<AxisId>(() => {
  const spec = app.model.outputs.lineageAxisSpec;
  return { type: "String", name: "pl7.app/clusterId", domain: spec?.domain ?? {} };
});

// Id is `<donor>/<clone>` with a donor column, else the clone. Split at the last "/":
// a donor name may hold one, a clone id never does.
const lineageTitle = (lineageId: string) => {
  const cut = lineageId.lastIndexOf("/");
  return cut < 0 ? `Tree ${lineageId}` : `${lineageId.slice(0, cut)} / ${lineageId.slice(cut + 1)}`;
};

const onLineageClicked = async (key?: PTableKey) => {
  if (!key) return;
  // Absent on blocks created before this field.
  if (!app.model.data.treeViews) app.model.data.treeViews = [];
  const lineageId = String(key[0]);
  const id = `tree-${lineageId}`;
  // Views from an earlier run name other lineages now; dropped on this click.
  const runKey = app.model.outputs.runKey;
  dropStaleViews(runKey);
  if (!app.model.data.treeViews.some((v) => v.id === id)) {
    app.model.data.treeViews.push({
      id,
      lineageId,
      runKey,
      state: { title: lineageTitle(lineageId), template: "dendro" },
    });
  }
  await app.navigateTo(`/tree?id=${encodeURIComponent(id)}`);
};
</script>

<template>
  <PlBlockPage>
    <template #title>Lineage table</template>
    <template #append>
      <PlBtnGhost @click.stop="() => (view.settingsOpen = true)">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>

    <PlAgDataTableV2
      v-model="app.model.data.treesTableState"
      :settings="treesSettings"
      :show-cell-button-for-axis-id="lineageAxis"
      show-columns-panel
      @cell-button-clicked="onLineageClicked"
    />
  </PlBlockPage>

  <PlSlideModal v-model="view.settingsOpen" :shadow="true">
    <template #title>Settings</template>
    <SettingsPanel />
  </PlSlideModal>
</template>
