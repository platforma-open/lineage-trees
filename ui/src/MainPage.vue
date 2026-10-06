<script setup lang="ts">
import {
  PlAlert,
  PlBlockPage,
  PlBtnGhost,
  PlMaskIcon24,
  PlSlideModal,
} from "@platforma-sdk/ui-vue";
import { computed, reactive, watch } from "vue";
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

const hasDonorColumn = computed(() => app.model.data.donorColumn !== undefined);
// Every clonotype is its own lineage: no V, J and CDR3-length group held two.
// Five genes are enough to see the pattern; the share is of that source's clonotypes.
const geneMismatch = computed(() => app.model.outputs.geneMismatch ?? []);
const geneMismatchText = computed(() =>
  geneMismatch.value
    .slice(0, 5)
    .map((g) => `${g.gene} (${g.source}, ${(100 * g.share).toFixed(0)}%)`)
    .join(", "),
);
const singletonText = computed(() => {
  const donors = app.model.outputs.singletonDonors ?? [];
  if (donors.length === 0) return undefined;
  return hasDonorColumn.value
    ? `No lineage in ${donors.length === 1 ? "donor" : "donors"} ${donors.join(", ")} has more than one clonotype.`
    : "No lineage has more than one clonotype.";
});
const emptyText = computed(() => {
  const donors = app.model.outputs.emptyDonors ?? [];
  if (donors.length === 0 || !hasDonorColumn.value) return undefined;
  return donors.length === 1
    ? `Donor ${donors[0]} has no clonotypes in the picked datasets.`
    : `Donors ${donors.join(", ")} have no clonotypes in the picked datasets.`;
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

    <!-- Run is disabled with no reason given; say which setting holds it back. -->
    <PlAlert v-if="app.model.outputs.settingsProblem" type="warn">
      {{ app.model.outputs.settingsProblem }} Run stays disabled until it is fixed in Settings.
    </PlAlert>

    <PlAlert v-if="app.model.outputs.samplesWithoutDonor?.noneNamed" type="error">
      No sample in the picked datasets has a value in the donor column, so nothing can be clustered.
      Fill in the donor column or clear it in Settings.
    </PlAlert>

    <PlAlert v-if="singletonText" type="warn">{{ singletonText }}</PlAlert>
    <PlAlert v-if="emptyText" type="info">{{ emptyText }}</PlAlert>

    <!-- A lineage with no tree draws the same empty plot as a broken one, so say which. -->
    <!-- With single-clonotype lineages that line already says why no tree was built. -->
    <PlAlert v-if="app.model.outputs.noTreesReason && !singletonText" type="warn">
      {{ app.model.outputs.noTreesReason }} Double-click a donor to open its logs.
    </PlAlert>

    <!-- Common genes one source never uses: the references likely name them differently. -->
    <PlAlert v-if="geneMismatch.length > 0" type="warn">
      Some V/J genes are common in the {{ geneMismatch[0].source }} datasets but never appear in the
      {{ geneMismatch[0].source === "MiXCR" ? "imported" : "MiXCR" }} ones, which usually means the
      two were annotated against different references. Clonotypes only join a lineage when their V
      gene, J gene and CDR3 length match, so these did not cluster across the two:
      {{ geneMismatchText }}.
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
