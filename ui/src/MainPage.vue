<script setup lang="ts">
import {
  PlAlert,
  PlBlockPage,
  PlBtnGhost,
  PlMaskIcon24,
  PlSlideModal,
} from "@platforma-sdk/ui-vue";
import { reactive, watch } from "vue";
import DonorReportPanel from "./DonorReportPanel.vue";
import DonorTable from "./DonorTable.vue";
import SettingsPanel from "./SettingsPanel.vue";
import { useApp } from "./app";

const app = useApp();

// Local, not `app.model.data`: a watcher writing drawer state into block data is a hairpin.
const view = reactive<{ settingsOpen: boolean; reportOpen: boolean; donor?: string }>({
  settingsOpen: (app.model.data.datasets ?? []).length === 0,
  reportOpen: false,
});

const openReport = (donor: string) => {
  view.donor = donor;
  view.reportOpen = true;
};

// Close the drawer once clustering starts.
watch(
  () => app.model.outputs.isRunning,
  (isRunning, wasRunning) => {
    if (wasRunning === false && isRunning === true) view.settingsOpen = false;
  },
);
</script>

<template>
  <PlBlockPage>
    <template #title>Lineage Trees</template>

    <template #append>
      <PlBtnGhost @click.stop="() => (view.settingsOpen = true)">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>

    <PlAlert v-if="(app.model.data.datasets ?? []).length === 0" type="info">
      Pick one or more IG clonotyping datasets in Settings to infer lineages from. Bulk heavy chain
      and paired single cell are both accepted, imported or from MiXCR; TCR is not.
    </PlAlert>

    <!-- A lineage with no tree draws the same empty plot as a broken one, so say which. -->
    <PlAlert v-if="app.model.outputs.noTreesReason" type="warn">
      {{ app.model.outputs.noTreesReason }} Double-click a donor to open its logs.
    </PlAlert>

    <div :style="{ flex: 1 }">
      <DonorTable @open="openReport" />
    </div>
  </PlBlockPage>

  <PlSlideModal v-model="view.settingsOpen" :shadow="true">
    <template #title>Settings</template>
    <SettingsPanel />
  </PlSlideModal>
  <PlSlideModal v-model="view.reportOpen" width="80%">
    <template #title>Logs for {{ view.donor ?? "..." }}</template>
    <DonorReportPanel v-model="view.donor" />
  </PlSlideModal>
</template>
