<script setup lang="ts">
import type {
  SearchParameters,
  SOIListParameters,
} from "@platforma-open/milaboratories.lineage-trees.model";
import type { ListOption } from "@platforma-sdk/ui-vue";
import { PlDropdown, PlTextField } from "@platforma-sdk/ui-vue";
import { alphabetOptions, targetFeatureOptions } from "./soiUtil";

// Alphabet and target feature are fixed at creation: the sequences were checked against them.
const model = defineModel<SOIListParameters>({ required: true });

const searchSettings: ListOption<SearchParameters>[] = [
  {
    value: { type: "preset_alignment_search_top", dissimilarityPercent: 1 },
    label: "Max dissimilarity 1%",
  },
  {
    value: { type: "preset_alignment_search_top", dissimilarityPercent: 2 },
    label: "Max dissimilarity 2%",
  },
  {
    value: { type: "preset_alignment_search_top", dissimilarityPercent: 5 },
    label: "Max dissimilarity 5%",
  },
  {
    value: { type: "preset_alignment_search_top", dissimilarityPercent: 10 },
    label: "Max dissimilarity 10%",
  },
  {
    value: { type: "preset_alignment_search_top", dissimilarityPercent: 15 },
    label: "Max dissimilarity 15%",
  },
  {
    value: { type: "tree_search_top", parameters: "oneMismatch" },
    label: "Fuzzy exact: 1 mismatch",
  },
  {
    value: { type: "tree_search_top", parameters: "oneMismatchOrIndel" },
    label: "Fuzzy exact: 1 mismatch or indel",
  },
  {
    value: { type: "tree_search_top", parameters: "twoMismatches" },
    label: "Fuzzy exact: 2 mismatches",
  },
  {
    value: { type: "tree_search_top", parameters: "twoMismatchesOrIndels" },
    label: "Fuzzy exact: 2 mismatches or indels",
  },
  {
    value: { type: "tree_search_top", parameters: "threeMismatches" },
    label: "Fuzzy exact: 3 mismatches",
  },
  {
    value: { type: "tree_search_top", parameters: "threeMismatchesOrIndels" },
    label: "Fuzzy exact: 3 mismatches or indels",
  },
  {
    value: { type: "tree_search_top", parameters: "fourMismatches" },
    label: "Fuzzy exact: 4 mismatches",
  },
  {
    value: { type: "tree_search_top", parameters: "fourMismatchesOrIndels" },
    label: "Fuzzy exact: 4 mismatches or indels",
  },
];
</script>

<template>
  <PlDropdown v-model="model.type" :options="alphabetOptions" label="Alphabet" :disabled="true" />
  <PlDropdown
    v-model="model.targetFeature"
    :options="targetFeatureOptions"
    label="Target feature"
    :disabled="true"
  />
  <PlTextField v-model="model.name" label="List name" />
  <PlDropdown v-model="model.searchParameters" :options="searchSettings" label="Search parameters">
    <template #tooltip>
      Max dissimilarity aligns each sequence against every node and keeps hits within the given
      share of mismatches. Fuzzy exact keeps hits within a fixed number of mismatches, or mismatches
      and indels.
    </template>
  </PlDropdown>
</template>
