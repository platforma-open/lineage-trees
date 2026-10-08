<script setup lang="ts">
import type {
  SequenceOfInterest,
  SOIList,
  SOIListParameters,
} from "@platforma-open/milaboratories.lineage-trees.model";
import type { LocalImportFileHandle, PlId } from "@platforma-sdk/model";
import { getRawPlatformaInstance, uniquePlId } from "@platforma-sdk/model";
import type { ListOption } from "@platforma-sdk/ui-vue";
import {
  PlAlert,
  PlBlockPage,
  PlBtnGhost,
  PlBtnPrimary,
  PlBtnSecondary,
  PlDialogModal,
  PlDropdown,
  PlDropdownLine,
  PlSlideModal,
} from "@platforma-sdk/ui-vue";
import { computed, reactive, watch } from "vue";
import { useApp } from "../app";
import SOIImportModal from "./SOIImportModal.vue";
import SOISettingsPanel from "./SOISettingsPanel.vue";
import SOITable from "./SOITable.vue";
import { alphabetOptions, inferNewName, lengthMismatches, targetFeatureOptions } from "./soiUtil";

// Sequence lists to find in the trees. Top hits per node and lineage show in the tree view
// and lineage table. From mixcr-shm-trees.
const app = useApp();

const lists = computed(() => app.model.data.sequencesOfInterest ?? []);

const listOptions = computed(
  () =>
    [
      ...lists.value.map((l) => ({ value: l.parameters.id as string, label: l.parameters.name })),
      { value: "", label: "+ Add new list" },
    ] as ListOption<string>[],
);

type NewListSettings = Pick<SOIListParameters, "type" | "targetFeature">;

// View state, not block state.
const view = reactive<{
  currentListId?: PlId;
  importFile?: LocalImportFileHandle;
  settingsOpen: boolean;
  newList?: NewListSettings;
  listToDelete?: PlId;
}>({ settingsOpen: false });

watch(
  () => lists.value.map((v) => v.parameters.id),
  (ids) => {
    if (view.currentListId === undefined && ids.length > 0) view.currentListId = ids[0];
    else if (!ids.includes(view.currentListId as PlId)) view.currentListId = ids[0];
  },
  { immediate: true },
);

function startAddingNewList() {
  view.newList = { type: "nucleotide", targetFeature: "CDR3" };
}

function addNewList() {
  if (!view.newList) return;
  const id = uniquePlId();
  const prefix = `${view.newList.targetFeature === "CDR3" ? "CDR3" : "Full sequence"} ${
    view.newList.type === "nucleotide" ? "nt" : "aa"
  }`;
  const name = inferNewName(
    lists.value.map((l) => l.parameters.name),
    (i) => `${prefix} list (${i})`,
  );
  if (app.model.data.sequencesOfInterest === undefined) app.model.data.sequencesOfInterest = [];
  app.model.data.sequencesOfInterest.push({
    sequences: [],
    parameters: {
      id,
      name,
      type: view.newList.type,
      targetFeature: view.newList.targetFeature,
      searchParameters: {
        type: "preset_alignment_search_top",
        dissimilarityPercent: view.newList.type === "nucleotide" ? 10 : 2,
      },
    },
  });
  view.currentListId = id;
  view.newList = undefined;
  view.settingsOpen = true;
}

const currentIndex = computed(() =>
  lists.value.findIndex((l) => l.parameters.id === view.currentListId),
);

function deleteList() {
  const idx = currentIndex.value;
  if (idx < 0) return;
  app.model.data.sequencesOfInterest.splice(idx, 1);
  view.listToDelete = undefined;
}

async function importFile() {
  const file = await getRawPlatformaInstance().lsDriver.showOpenSingleFileDialog({
    title: "Select a sequence list file",
    buttonLabel: "Import",
    filters: [{ name: "FASTA or table", extensions: ["tsv", "csv", "txt", "fa", "fasta"] }],
  });
  if (!file.file) return;
  view.importFile = file.file;
}

function onImport(records: SequenceOfInterest[]) {
  const idx = currentIndex.value;
  if (idx < 0) return;
  const list = app.model.data.sequencesOfInterest[idx];
  list.sequences = [...list.sequences, ...records];
  view.importFile = undefined;
}

const currentListIdForDropdown = computed<string>({
  get: () => view.currentListId ?? "",
  set: (value) => {
    if (value === "") startAddingNewList();
    else view.currentListId = value as PlId;
  },
});

const currentList = computed<SOIList | undefined>(() =>
  currentIndex.value < 0 ? undefined : lists.value[currentIndex.value],
);

// Sequences in the open list that cannot hit its target feature.
const mismatches = computed(() =>
  currentList.value === undefined
    ? []
    : lengthMismatches(
        currentList.value.sequences,
        currentList.value.parameters.type,
        currentList.value.parameters.targetFeature,
      ),
);

// What the last run found for the open list; undefined until a run searched it.
const hits = computed(() => {
  const id = currentList.value?.parameters.id;
  return id === undefined ? undefined : app.model.outputs.soiHits?.[id];
});
</script>

<template>
  <PlBlockPage>
    <template #title>
      <PlDropdownLine
        v-if="view.currentListId"
        v-model="currentListIdForDropdown"
        :options="listOptions"
      />
      <PlBtnGhost v-else icon="add" @click.stop="startAddingNewList">Create new list</PlBtnGhost>
    </template>
    <template #append>
      <template v-if="view.currentListId">
        <PlBtnGhost icon="delete-bin" @click.stop="() => (view.listToDelete = view.currentListId)"
          >Delete current list</PlBtnGhost
        >
        <PlBtnGhost icon="dna-import" @click.stop="importFile">Import sequences</PlBtnGhost>
        <PlBtnGhost icon="settings" @click.stop="() => (view.settingsOpen = true)"
          >Settings</PlBtnGhost
        >
      </template>
    </template>

    <PlAlert v-if="currentList && currentList.sequences.length === 0" type="info">
      Import sequences into this list. Only lists with sequences are searched, so adding a list does
      not require a new run until it has some.
    </PlAlert>
    <PlAlert v-if="mismatches.length > 0" type="warn">
      {{ mismatches.length }} of {{ currentList?.sequences.length }} sequences do not fit this
      list's target feature and will not be found. {{ mismatches[0]
      }}{{ mismatches.length > 1 ? " (and more)" : "" }}
    </PlAlert>
    <PlAlert
      v-if="currentList && currentList.sequences.length > 0 && hits !== undefined"
      :type="hits.lineages === 0 ? 'warn' : 'success'"
    >
      <template v-if="hits.lineages === 0">
        The last run found no lineage matching this list. Check the alphabet and target feature,
        that CDR3s include the conserved Cys and Trp or Phe, and the dissimilarity setting.
      </template>
      <template v-else>
        The last run matched {{ hits.nodes }} nodes in {{ hits.lineages }} lineages. The hits are
        the "Top Hit {{ currentList.parameters.name }}" columns in the lineage table and the tree
        view.
      </template>
    </PlAlert>

    <div v-if="currentList" :style="{ flex: 1 }">
      <SOITable v-model="currentList.sequences" />
    </div>

    <SOIImportModal
      v-if="view.importFile && currentList"
      :file="view.importFile"
      :alphabet="currentList.parameters.type"
      :target-feature="currentList.parameters.targetFeature"
      @on-close="view.importFile = undefined"
      @on-import="onImport"
    />

    <PlDialogModal
      v-if="view.newList !== undefined"
      :model-value="true"
      @update:model-value="
        (v) => {
          if (!v) view.newList = undefined;
        }
      "
    >
      <template #title>New list</template>
      <PlDropdown v-model="view.newList.type" :options="alphabetOptions" label="Alphabet" />
      <PlDropdown
        v-model="view.newList.targetFeature"
        :options="targetFeatureOptions"
        label="Target feature"
      />
      <p>
        The list is searched against the observed clonotypes' CDR3 or full sequence in this
        alphabet. Nucleotide sequences are translated when added to an amino acid list.
      </p>
      <template #actions>
        <PlBtnPrimary @click="addNewList">Create</PlBtnPrimary>
        <PlBtnSecondary @click="() => (view.newList = undefined)">Cancel</PlBtnSecondary>
      </template>
    </PlDialogModal>

    <PlDialogModal
      v-if="view.listToDelete !== undefined"
      :model-value="true"
      @update:model-value="
        (v) => {
          if (!v) view.listToDelete = undefined;
        }
      "
    >
      <template #title>Delete this list?</template>
      <template #actions>
        <PlBtnPrimary @click="deleteList">Delete</PlBtnPrimary>
        <PlBtnSecondary @click="() => (view.listToDelete = undefined)">Cancel</PlBtnSecondary>
      </template>
    </PlDialogModal>

    <PlSlideModal v-if="currentList" v-model="view.settingsOpen">
      <template #title>List settings</template>
      <SOISettingsPanel v-model="currentList.parameters" />
    </PlSlideModal>
  </PlBlockPage>
</template>
