<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import { PlBlockPage } from "@platforma-sdk/ui-vue";
import { computed } from "vue";
import { useApp } from "./app";

const app = useApp();

const SIZE_RANK = "pl7.app/clustering/sizeRank";
const FREQUENCY = "pl7.app/clustering/abundancePercent";
const SAMPLE_AXIS = "pl7.app/sampleId";

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

<!-- GraphMaker must be the only child: it renders its own title and tabs, and needs the full height. -->
<template>
  <PlBlockPage no-body-gutters>
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
