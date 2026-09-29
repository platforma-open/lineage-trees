<script setup lang="ts">
import type { SimpleOption } from "@platforma-sdk/ui-vue";
import { PlBtnGroup, PlLogView } from "@platforma-sdk/ui-vue";
import { computed, ref, watch } from "vue";
import type { Stage } from "./donorRows";
import { STAGES, useDonorRows } from "./donorRows";
import { useApp } from "./app";

// A donor's logs, one stage at a time.
const donor = defineModel<string | undefined>();
const donorRows = useDonorRows();
const app = useApp();

const row = computed(() => donorRows.value.find((r) => r.donor === donor.value));

type Tab = "preparation" | Stage;
const tabs: SimpleOption<Tab>[] = [
  { value: "preparation", text: "Preparation" },
  ...STAGES.map(({ key, label }) => ({ value: key, text: label })),
];

// Opens on the stage the donor is at; local, since it is view state.
const tab = ref<Tab>("preparation");
watch(donor, () => (tab.value = row.value?.currentStage ?? "preparation"), { immediate: true });

// Run-wide steps, the same on every donor, headed by the mode statement.
const preparationLog = computed(() =>
  [app.model.outputs.modeStatement, app.model.outputs.mergeLog, app.model.outputs.clusteringLog]
    .filter((l) => l)
    .join("\n"),
);

const handle = computed(() =>
  tab.value === "preparation" ? undefined : row.value?.logs[tab.value],
);
</script>

<template>
  <div v-if="row !== undefined" class="panel">
    <PlBtnGroup v-model="tab" :options="tabs" />
    <div class="log">
      <PlLogView v-if="tab === 'preparation'" :value="preparationLog || 'Not started.'" />
      <PlLogView v-else-if="handle" :log-handle="handle" />
      <PlLogView v-else value="Not started." />
    </div>
  </div>
  <div v-else>No donor selected</div>
</template>

<style lang="css" scoped>
.panel {
  display: flex;
  flex-direction: column;
  gap: 12px;
  height: 100%;
}
/* The log fills what the tabs leave. */
.log {
  flex: 1;
  min-height: 0;
}
</style>
