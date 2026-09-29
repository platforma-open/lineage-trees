<script setup lang="ts">
import type { BasketNode } from "@platforma-open/milaboratories.lineage-trees.model";
import { createPlDataTableStateV2, uniquePlId } from "@platforma-sdk/model";
import {
  PlBtnGhost,
  PlBtnPrimary,
  PlDialogModal,
  PlDropdown,
  PlTextField,
} from "@platforma-sdk/ui-vue";
import { computed, ref } from "vue";
import { useApp } from "./app";
import { mergeBasketNodes } from "./baskets";

// The nodes arrive already copied, so adding is a single write on the click.
const props = defineProps<{ nodes: BasketNode[] }>();
const emit = defineEmits<{ close: [] }>();

const app = useApp();

const baskets = computed(() => app.model.data.baskets ?? []);

// Read once, when the dialog opens: a name typed by the user is theirs to keep.
const firstFreeName = () => {
  const taken = new Set(baskets.value.map((b) => b.name));
  let i = baskets.value.length + 1;
  while (taken.has(`Basket ${i}`)) i++;
  return `Basket ${i}`;
};

// An empty target means a new basket.
const target = ref<string>("");
const newName = ref(firstFreeName());

const options = computed(() => [
  ...baskets.value.map((b) => ({ value: b.id, label: b.name })),
  { value: "", label: "New basket..." },
]);

const canAdd = computed(() => target.value !== "" || newName.value.trim() !== "");

const add = () => {
  // Absent on blocks created before baskets.
  if (!app.model.data.baskets) app.model.data.baskets = [];
  const existing = app.model.data.baskets.find((b) => b.id === target.value);
  if (existing) {
    existing.nodes = mergeBasketNodes(existing.nodes, props.nodes);
  } else {
    app.model.data.baskets.push({
      id: uniquePlId(),
      name: newName.value.trim(),
      nodes: mergeBasketNodes([], props.nodes),
      tableState: createPlDataTableStateV2(),
    });
  }
  emit("close");
};
</script>

<template>
  <PlDialogModal :model-value="true" @update:model-value="(open) => !open && emit('close')">
    <template #title>
      Add {{ nodes.length }} {{ nodes.length === 1 ? "node" : "nodes" }} to a basket
    </template>
    <PlDropdown v-model="target" :options="options" label="Basket" />
    <PlTextField v-if="target === ''" v-model="newName" label="New basket name" />
    <template #actions>
      <PlBtnPrimary :disabled="!canAdd" @click="add">Add</PlBtnPrimary>
      <PlBtnGhost @click="emit('close')">Cancel</PlBtnGhost>
    </template>
  </PlDialogModal>
</template>
