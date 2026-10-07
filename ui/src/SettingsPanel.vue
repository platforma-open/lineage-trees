<script setup lang="ts">
import {
  CLUSTERING_MODE_OPTIONS,
  IGPHYML_SCOPE_OPTIONS,
  effectiveKnown,
  refKey,
} from "@platforma-open/milaboratories.lineage-trees.model";
import {
  PlAccordion,
  PlAccordionSection,
  PlAlert,
  PlDropdown,
  PlDropdownMultiRef,
  PlDropdownRef,
  PlNumberField,
  PlRow,
} from "@platforma-sdk/ui-vue";
import { computed, ref } from "vue";
import { useApp } from "./app";

const clusteringModeOptions = [...CLUSTERING_MODE_OPTIONS];

const app = useApp();

// Read from data, not outputs, so the option does not flicker while outputs load.
const hasKnown = computed(() => effectiveKnown(app.model.data).length > 0);
const igPhyMLOptions = computed(() =>
  IGPHYML_SCOPE_OPTIONS.filter((o) => o.value !== "known" || hasKnown.value),
);
// A pick that leaves no known antibody resets "known", so the stored scope matches the screen.
const settleScope = () => {
  if (!hasKnown.value && app.model.data.igPhyMLScope === "known") {
    app.model.data.igPhyMLScope = "none";
  }
};
const pickedDatasets = computed({
  get: () => app.model.data.datasets,
  set: (refs) => {
    app.model.data.datasets = refs;
    settleScope();
  },
});
const pickedKnown = computed({
  get: () => app.model.data.knownDatasets,
  set: (refs) => {
    app.model.data.knownDatasets = refs;
    settleScope();
  },
});
const advancedOpen = ref(false);

// A known antibody dataset is one of the picked datasets; read from inputOptions so the list comes fast.
const pickedKeys = computed(() => new Set(app.model.data.datasets.map(refKey)));
const knownOff = computed(() => pickedKeys.value.size < 2);
const knownOptions = computed(() =>
  app.model.outputs.inputOptions
    ?.filter((o) => pickedKeys.value.has(refKey(o.ref)))
    .sort((a, b) => refKey(a.ref).localeCompare(refKey(b.ref))),
);
// Picked datasets with neither MiXCR alignments nor alignment columns: no germline, so no trees.
const unaligned = computed(() => {
  const routes = app.model.outputs.alignmentRoutes ?? {};
  return (app.model.outputs.inputOptions ?? [])
    .filter((o) => pickedKeys.value.has(refKey(o.ref)) && routes[refKey(o.ref)] === "none")
    .map((o) => o.label);
});
const outsideDonor = computed(() => app.model.outputs.datasetsOutsideDonorColumn ?? []);
const withoutDonor = computed(() => app.model.outputs.samplesWithoutDonor);
// View state, so closing the note is not a block edit; it shows again when the panel reopens.
const adaptiveNoteOpen = ref(true);
</script>

<template>
  <PlAlert v-if="app.model.outputs.settingsProblem" type="warn">
    {{ app.model.outputs.settingsProblem }}
  </PlAlert>
  <PlDropdownMultiRef
    v-model="pickedDatasets"
    :options="app.model.outputs.inputOptions"
    label="Datasets"
    placeholder="Pick one or more IG datasets"
    :required="true"
    :error="pickedDatasets.length === 0 ? 'Input dataset is required' : undefined"
  >
  </PlDropdownMultiRef>
  <PlDropdownMultiRef
    v-model="pickedKnown"
    :options="knownOptions"
    :disabled="knownOff"
    label="Known antibody datasets"
    :placeholder="knownOff ? 'Pick two or more datasets' : 'None'"
  >
    <template #tooltip>
      Datasets of antibodies whose relatives you are looking for in the others, such as
      characterised leads or hits from a screen. Their clonotypes are clustered with everything
      else: every lineage says how many it holds, and every other member of such a lineage gets its
      distance along the tree to the nearest known antibody, which is exported for ranking.
    </template>
  </PlDropdownMultiRef>
  <PlAlert v-if="unaligned.length > 0" type="warn">
    {{ unaligned.join(", ") }} {{ unaligned.length === 1 ? "has" : "have" }} no germline alignment,
    so {{ unaligned.length === 1 ? "its" : "their" }} clonotypes join lineages but get no trees.
  </PlAlert>
  <PlDropdownRef
    v-model="app.model.data.donorColumn"
    :options="app.model.outputs.donorOptions"
    label="Donor column"
    placeholder="Single donor"
    clearable
  >
    <template #tooltip>
      Sample metadata naming the subject each sample came from. Clustering runs inside a donor.
    </template>
  </PlDropdownRef>
  <PlAlert v-if="outsideDonor.length > 0" type="warn">
    The donor column has no samples from {{ outsideDonor.join(", ") }}, so that data is left out of
    clustering. Clear the donor column to cluster everything as one donor.
  </PlAlert>
  <PlAlert v-if="withoutDonor?.noneNamed" type="error">
    No sample in the picked datasets has a value in the donor column, so nothing can be clustered.
    Fill in the donor column or clear it.
  </PlAlert>
  <PlAlert v-else-if="(withoutDonor?.datasets.length ?? 0) > 0" type="warn">
    <div v-for="d in withoutDonor?.datasets" :key="d.dataset">
      {{ d.dataset }}: {{ d.missing }} of {{ d.total }} samples have no donor, so they are left out
      of clustering.
    </div>
  </PlAlert>
  <PlAccordion :multiple="true">
    <PlAccordionSection v-model="advancedOpen" label="Advanced">
      <PlDropdown
        v-model="app.model.data.clusteringMode"
        :options="clusteringModeOptions"
        label="Clustering"
      >
        <template #tooltip>
          How clonotypes are grouped into lineages. Only clonotypes with the same V gene, J gene and
          CDR3 length are compared.
          <br /><br />
          <b>Fixed threshold</b>: two clonotypes join one lineage when their heavy-chain CDR3s are
          closer than the Clustering threshold below. <br /><br />
          <b>Adaptive</b>: HILARy picks a threshold for each group from the data itself, aiming at
          the Precision and Sensitivity below. When alignments are available it also uses mutations
          the clonotypes share outside the CDR3.
        </template>
      </PlDropdown>
      <PlAlert
        v-if="app.model.data.clusteringMode === 'adaptive'"
        v-model="adaptiveNoteOpen"
        type="warn"
        closeable
      >
        HILARy's adaptive mode is calibrated on human repertoires; on other species prefer the fixed
        threshold.
      </PlAlert>
      <PlNumberField
        v-if="app.model.data.clusteringMode !== 'adaptive'"
        v-model="app.model.data.clusteringThreshold"
        label="Clustering threshold"
        :min-value="0.01"
        :max-value="0.99"
        :step="0.01"
      >
        <template #tooltip>
          How different two heavy-chain CDR3s may be and still belong to one lineage, as a share of
          the CDR3's length: at 0.2 they join when no more than about 20% of their positions differ.
          Lower values give smaller, stricter lineages; higher values join more distant relatives
          but risk joining unrelated clonotypes. 0.2 is HILARy's default.
        </template>
      </PlNumberField>
      <PlRow v-if="app.model.data.clusteringMode === 'adaptive'">
        <PlNumberField
          v-model="app.model.data.precision"
          label="Precision"
          :min-value="0.5"
          :max-value="1"
          :step="0.01"
        >
          <template #tooltip>
            How cautious the grouping is: of all the clonotype pairs put in the same lineage, the
            share that should really be related. Raise it to keep lineages pure, at the cost of
            splitting some real ones. The default is 0.99.
          </template>
        </PlNumberField>
        <PlNumberField
          v-model="app.model.data.sensitivity"
          label="Sensitivity"
          :min-value="0.5"
          :max-value="1"
          :step="0.01"
        >
          <template #tooltip>
            How complete the grouping is: of all the truly related clonotype pairs, the share that
            should end up in the same lineage. Raise it to keep real lineages together, at the cost
            of joining in some unrelated clonotypes. The default is 0.9.
          </template>
        </PlNumberField>
      </PlRow>
      <PlDropdown v-model="app.model.data.igPhyMLScope" :options="igPhyMLOptions" label="IgPhyML">
        <template #tooltip>
          Which lineages IgPhyML builds instead of FastTree and RAxML. Its HLP19 codon model
          accounts for the context-dependent hotspot biases of somatic hypermutation and resolves
          topology better, but it is far slower.
        </template>
      </PlDropdown>
      <PlNumberField
        v-model="app.model.data.maxTipsPerTree"
        label="Maximum tips per tree"
        :min-value="3"
        :max-value="100000"
        :step="50"
        clearable
      >
        <template #tooltip>
          Lineage with more than this many sequences will be subsampled to that many at random.
        </template>
      </PlNumberField>
      <PlNumberField
        v-model="app.model.data.minTipsPerTree"
        label="Minimum tips per tree"
        :min-value="2"
        :max-value="100000"
        :step="1"
        clearable
      >
        <template #tooltip>
          Lineages with fewer than this many sequences will not be built a tree.
        </template>
      </PlNumberField>
    </PlAccordionSection>
  </PlAccordion>
</template>
