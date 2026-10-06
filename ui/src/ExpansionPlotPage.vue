<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import { PlAlert, PlBlockPage } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "./app";

const app = useApp();

const SIZE_RANK = "pl7.app/clustering/sizeRank";
const FREQUENCY = "pl7.app/clustering/abundancePercent";
const SAMPLE_AXIS = "pl7.app/sampleId";

// With few lineages, each sample's frequencies rest on a handful of points and read oddly.
const FEW_LINEAGES = 100;
const fewLineages = computed(() => {
  const stats = app.model.outputs.donorStats;
  if (stats === undefined) return false;
  return stats.reduce((n, s) => n + s.lineage_count, 0) < FEW_LINEAGES;
});

// Defaults only; GraphMaker's own controls stay in charge. Empty axesSpec matches on name alone.
const defaultOptions = computed((): PredefinedGraphOption<"scatterplot">[] => [
  {
    inputName: "x",
    selectedSource: { kind: "PColumn", valueType: "Int", name: SIZE_RANK, axesSpec: [] },
  },
  {
    inputName: "y",
    selectedSource: { kind: "PColumn", valueType: "Double", name: FREQUENCY, axesSpec: [] },
  },
  {
    inputName: "grouping",
    selectedSource: { type: "String", name: SAMPLE_AXIS },
  },
]);
</script>

<!-- GraphMaker renders its own title and tabs; only a warning may sit above it, as on the tree page. -->
<template>
  <PlBlockPage no-body-gutters>
    <PlAlert v-if="fewLineages" type="warn">
      The data clustered into fewer than 100 lineages, which may mean that the lineage expansion
      plot is not reliable.
    </PlAlert>
    <GraphMaker
      v-model="app.model.data.expansionGraphState"
      chart-type="scatterplot"
      :p-frame="app.model.outputs.expansionPlot"
      :default-options="defaultOptions"
      :status-text="{
        noPframe: { title: 'Pick an input dataset on the Overview page and run the block' },
        notReady: { title: 'Choose columns for the axes under Data Mapping' },
        empty: { title: 'No lineages to plot' },
      }"
    />
  </PlBlockPage>
</template>
