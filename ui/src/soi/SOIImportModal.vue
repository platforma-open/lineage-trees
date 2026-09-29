<script setup lang="ts">
import type {
  Alphabet,
  SequenceOfInterest,
  TargetFeature,
} from "@platforma-open/milaboratories.lineage-trees.model";
import type { LocalImportFileHandle } from "@platforma-sdk/model";
import { getFileNameFromHandle, getRawPlatformaInstance, uniquePlId } from "@platforma-sdk/model";
import type { ListOption } from "@platforma-sdk/ui-vue";
import {
  PlBtnGhost,
  PlBtnPrimary,
  PlDialogModal,
  PlDropdown,
  PlLogView,
} from "@platforma-sdk/ui-vue";
import { computedAsync } from "@vueuse/core";
import { computed, reactive } from "vue";
import { detectAlphabet, translate } from "./alphabets";
import { readFileForImport } from "./dataimport";
import { lengthMismatches } from "./soiUtil";

// Reads a FASTA or table, checks its alphabet against the list's, and translates nucleotides
// for an amino acid list. From mixcr-shm-trees.
const props = defineProps<{
  file: LocalImportFileHandle;
  alphabet: Alphabet;
  targetFeature: TargetFeature;
}>();

const emit = defineEmits<{ onClose: []; onImport: [data: SequenceOfInterest[]] }>();

type ErrorMessage = { title: string; message?: string };
const data = reactive<{
  errorMessage?: ErrorMessage;
  importing: boolean;
  sequenceColumn?: number;
  nameColumn?: number;
}>({ importing: false });

const fileName = computed(() => getFileNameFromHandle(props.file));
const fileType = computed<"table" | "fasta">(() =>
  /\.(tsv|csv|txt)$/i.test(fileName.value) ? "table" : "fasta",
);

const fileContent = computedAsync(async () => {
  const pl = getRawPlatformaInstance();
  if ((await pl.lsDriver.getLocalFileSize(props.file)) > 5_000_000) {
    data.errorMessage = { title: "File is too big" };
    return undefined;
  }
  try {
    return await pl.lsDriver.getLocalFileContent(props.file);
  } catch (e: unknown) {
    console.error(e);
    data.errorMessage = { title: "Error reading file", message: String(e) };
    return undefined;
  }
});

const tableData = computed(() => {
  if (fileType.value !== "table" || fileContent.value === undefined) return undefined;
  try {
    return readFileForImport(fileContent.value, fileName.value);
  } catch (e: unknown) {
    console.error(e);
    data.errorMessage = { title: "Error reading table", message: String(e) };
    return undefined;
  }
});

const columnOptions = computed<ListOption<number>[]>(
  () => tableData.value?.data.columns.map((c, idx) => ({ value: idx, label: c.header })) ?? [],
);

type FastaRecord = { readonly description: string; readonly sequence: string };

function parseFasta(content: string): FastaRecord[] {
  const records: FastaRecord[] = [];
  let description = "";
  let sequence = "";
  for (let line of content.split("\n")) {
    line = line.trim();
    if (line.startsWith(">")) {
      if (description || sequence) records.push({ description, sequence });
      description = line.replace(/^>\s*/, "");
      sequence = "";
    } else if (line.length > 0) {
      sequence += line;
    }
  }
  if (description || sequence) records.push({ description, sequence });
  return records;
}

const fastaData = computed(() => {
  if (fileType.value !== "fasta" || fileContent.value === undefined) return undefined;
  try {
    return parseFasta(new TextDecoder().decode(fileContent.value));
  } catch (e: unknown) {
    console.error(e);
    data.errorMessage = { title: "Error reading FASTA file", message: String(e) };
    return undefined;
  }
});

const recordsToImport = computed<FastaRecord[] | undefined>(() => {
  if (fileType.value === "table") {
    const sc = data.sequenceColumn;
    const nc = data.nameColumn;
    if (!tableData.value || sc === undefined || nc === undefined) return undefined;
    return tableData.value.data.rows
      .filter((r) => r[nc] !== undefined && r[sc] !== undefined)
      .map((r) => ({ description: String(r[nc]), sequence: String(r[sc]) }));
  }
  return fastaData.value;
});

const sequencesToImport = computed(() => {
  if (fileType.value === "table") {
    const sc = data.sequenceColumn;
    if (!tableData.value || sc === undefined) return undefined;
    return tableData.value.data.rows.filter((r) => r[sc] !== undefined).map((r) => String(r[sc]));
  }
  return fastaData.value?.map((r) => r.sequence);
});

const detectedAlphabet = computed(() => {
  if (sequencesToImport.value === undefined || sequencesToImport.value.length === 0)
    return undefined;
  try {
    return detectAlphabet(sequencesToImport.value);
  } catch (e: unknown) {
    console.log(e);
    return undefined;
  }
});

type TextAndRecords = { text: string; records?: FastaRecord[] };

const importData = computed<TextAndRecords>(() => {
  if (data.errorMessage) return { text: `Error:\n${data.errorMessage.title}` };
  if (sequencesToImport.value === undefined) return { text: "Please select the sequence column" };
  const da = detectedAlphabet.value;
  if (!da) return { text: "Can't detect alphabet" };
  if (da.type === "amino-acid" && props.alphabet === "nucleotide") {
    return {
      text: "Amino acid sequences cannot go into a nucleotide list; create an amino acid list.",
    };
  }
  let records = recordsToImport.value;
  if (records === undefined) return { text: "Please select the name column" };
  let text = "";
  if (da.type === "nucleotide" && props.alphabet === "amino-acid") {
    try {
      records = records.map(({ description, sequence }) => ({
        description,
        sequence: translate(sequence),
      }));
    } catch (e: unknown) {
      return { text: `Can't translate sequences: ${e instanceof Error ? e.message : String(e)}` };
    }
    text = "Sequences are translated to fit the list's alphabet.\n";
  }
  text += `Records to import: ${records.length}`;
  // Warn now: a CDR3 never hits a full sequence, nor the reverse, and the search is silent.
  const mismatches = lengthMismatches(
    records.map((r) => ({ name: r.description, sequence: r.sequence })),
    props.alphabet,
    props.targetFeature,
  );
  if (mismatches.length > 0) {
    text += `\n\nWarning, ${mismatches.length} of ${records.length} sequences do not fit a ${
      props.targetFeature === "CDR3" ? "CDR3" : "full sequence"
    } list and will not be found:\n${mismatches.slice(0, 10).join("\n")}${
      mismatches.length > 10 ? "\n..." : ""
    }`;
  }
  return { text, records };
});

const canImport = computed(() => (importData.value.records?.length ?? 0) > 0);

function runImport() {
  const records = importData.value.records;
  if (!records) return;
  data.importing = true;
  emit(
    "onImport",
    records.map((r) => ({ id: uniquePlId(), name: r.description, sequence: r.sequence })),
  );
}
</script>

<template>
  <PlDialogModal
    width="800px"
    :model-value="true"
    @update:model-value="
      (v) => {
        if (!v) emit('onClose');
      }
    "
  >
    <template #title>Import sequences</template>
    <template v-if="fileType === 'table'">
      <PlDropdown
        v-model="data.nameColumn"
        :options="columnOptions"
        clearable
        label="Name column"
      />
      <PlDropdown
        v-model="data.sequenceColumn"
        :options="columnOptions"
        clearable
        label="Sequence column"
      />
    </template>
    <PlLogView :value="importData.text" label="Import information" />
    <template #actions>
      <PlBtnPrimary :loading="data.importing" :disabled="!canImport" @click="runImport">
        Import
      </PlBtnPrimary>
      <PlBtnGhost @click="emit('onClose')">Cancel</PlBtnGhost>
    </template>
  </PlDialogModal>
</template>
