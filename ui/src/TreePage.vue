<script setup lang="ts">
import type { PredefinedGraphOption } from "@milaboratories/graph-maker";
import { GraphMaker } from "@milaboratories/graph-maker";
import type { PlSelectionModel } from "@platforma-sdk/model";
import { createPlDataTableStateV2, createPlSelectionModel } from "@platforma-sdk/model";
import {
  PlAgDataTableV2,
  PlAlert,
  PlBlockPage,
  PlBtnGhost,
  usePlDataTableSettingsV2,
} from "@platforma-sdk/ui-vue";
import { computed, ref, watch } from "vue";
import AddToBasketModal from "./AddToBasketModal.vue";
import { selectedNodes } from "./baskets";
import { resolvePath } from "./resolvePath";
import { useApp } from "./app";
import { indexOfCurrent, keyedStatus } from "./keyedStatus";
import { useAddToBasket } from "./useAddToBasket";

const app = useApp<`/tree?id=${string}` | `/path?id=${string}` | "/trees">();

const index = computed(() =>
  indexOfCurrent(app.model.data.treeViews, app.queryParams.id, app.model.outputs.runKey),
);
const missing = computed(() => index.value < 0);
const view = computed({
  get: () => app.model.data.treeViews[index.value],
  set: (v) => (app.model.data.treeViews[index.value] = v),
});

// Tooltip columns. Clonotype-axis columns stall the dendrogram, so the workflow copies them
// onto tree axes under this marker.
const NODE_PROPERTY = { "pl7.app/dendrogram/nodeProperty": "true" };

const nodeProperty = (name: string, label: string) => ({
  inputName: "tableContent",
  selectedSource: {
    kind: "PColumn",
    name,
    valueType: "String",
    annotations: { ...NODE_PROPERTY, "pl7.app/label": label },
    axesSpec: [],
  },
});

const hasDataset = computed(() => app.model.outputs.treeNodeColumns?.hasDatasetProperty === true);

const defaultOptions = computed(
  () =>
    [
      // Several datasets: tips colored by dataset. Anchors then take the shape, else the color.
      ...(hasDataset.value
        ? [
            {
              inputName: "nodeColor",
              selectedSource: {
                kind: "PColumn",
                name: "pl7.app/dendrogram/dataset",
                valueType: "String",
                annotations: { ...NODE_PROPERTY, "pl7.app/label": "Dataset" },
                axesSpec: [],
              },
            },
          ]
        : []),
      ...(app.model.outputs.treeNodeColumns?.hasAnchorProperty
        ? [
            {
              inputName: hasDataset.value ? "nodeShape" : "nodeColor",
              selectedSource: {
                kind: "PColumn",
                name: "pl7.app/dendrogram/isAnchor",
                valueType: "String",
                annotations: { ...NODE_PROPERTY, "pl7.app/label": "Anchor" },
                axesSpec: [],
              },
            },
          ]
        : []),
      {
        inputName: "value",
        selectedSource: {
          kind: "PColumn",
          name: "pl7.app/dendrogram/topology",
          valueType: "Long",
          axesSpec: [],
        },
      },
      {
        inputName: "height",
        selectedSource: {
          kind: "PColumn",
          name: "pl7.app/dendrogram/distance",
          valueType: "Double",
          annotations: { "pl7.app/dendrogram/distance/from": "parent" },
          axesSpec: [],
        },
      },
      // Of the sequences, only the heavy CDR3 aa; the rest are a column pick away.
      nodeProperty("pl7.app/vdj/sequence", "CDR3 aa"),
      {
        inputName: "tableContent",
        selectedSource: {
          kind: "PColumn",
          name: "pl7.app/dendrogram/isObserved",
          valueType: "String",
          axesSpec: [],
        },
      },
      nodeProperty("pl7.app/vdj/geneHit", "V gene"),
      nodeProperty("pl7.app/vdj/geneHit", "J gene"),
      // The light chain, only on runs that used light chains.
      ...(app.model.outputs.treeNodeColumns?.hasLightSequence
        ? [
            nodeProperty("pl7.app/vdj/geneHit", "Light V gene"),
            nodeProperty("pl7.app/vdj/geneHit", "Light J gene"),
          ]
        : []),
    ] as PredefinedGraphOption<"dendro">[],
);

// Paths are resolved on click so the path page only reads. GraphMaker's dendrogram passes no
// node id yet (`info[0].id` off `rawIndexes` is always undefined), so this errors until fixed.
const pathError = ref<string | undefined>();
// The page is reused across trees, so an error belongs to the tree it was raised on.
watch(
  () => view.value?.id,
  () => (pathError.value = undefined),
);

const openPath = async (clicked: unknown) => {
  pathError.value = undefined;
  const opened = view.value;
  if (!opened) return;
  const nodeId =
    typeof clicked === "number" || typeof clicked === "string" ? String(clicked) : undefined;
  if (nodeId === undefined) {
    pathError.value = "This node carries no id, so its path cannot be traced.";
    return;
  }
  const columns = app.model.outputs.treeNodeColumns;
  // A status-carrying output, so the handle is only there once it has settled.
  const frameStatus = app.model.outputs.treeNodesPf;
  const frame = frameStatus?.ok === true ? frameStatus.value : undefined;
  if (!columns || !frame) {
    pathError.value = "The tree is still loading. Try again in a moment.";
    return;
  }
  // Absent on blocks created before this field.
  if (!app.model.data.pathViews) app.model.data.pathViews = [];
  // One section per lineage and node; reopening an open path goes back to it.
  const id = `path-${opened.lineageId}-${nodeId}`;
  if (!app.model.data.pathViews.some((v) => v.id === id)) {
    try {
      const { nodeIds, nodeLabel } = await resolvePath(frame, columns, opened.lineageId, nodeId);
      app.model.data.pathViews.push({
        id,
        lineageId: opened.lineageId,
        runKey: app.model.outputs.runKey,
        // The same name the tree's own section carries, so the two read as a pair.
        lineageLabel: opened.state.title,
        nodeId,
        nodeLabel,
        nodeIds,
        tableState: createPlDataTableStateV2(),
      });
    } catch (caught) {
      pathError.value = caught instanceof Error ? caught.message : String(caught);
      return;
    }
  }
  await app.navigateTo(`/path?id=${id}`);
};

const TOOLTIP_BUTTONS = [
  { id: "path", label: "Mutational path" },
  { id: "basket", label: "Add to basket" },
];

const onTooltipButton = (clicked: unknown, buttonId: string) => {
  if (buttonId !== "basket") return openPath(clicked);
  const v = view.value;
  if (!v || (typeof clicked !== "number" && typeof clicked !== "string")) return;
  return basket.start(v.lineageId, v.state.title, [String(clicked)]);
};

// Drops only the view section; the tree result is untouched.
const close = async () => {
  const at = index.value;
  if (at < 0) return;
  app.model.data.treeViews.splice(at, 1);
  await app.navigateTo("/trees");
};

// Stored per tree. A tree saved before the table existed opens as a graph.
const tab = computed({
  get: () => view.value?.tab ?? "graph",
  set: (value) => {
    if (view.value) view.value = { ...view.value, tab: value };
  },
});

const tableState = computed({
  get: () => view.value?.tableState ?? createPlDataTableStateV2(),
  set: (value) => {
    if (view.value) view.value = { ...view.value, tableState: value };
  },
});

// Keeps the output's status so the grid shows loading and errors instead of going blank.
const tableStatus = computed(() => keyedStatus(app.model.outputs.treeNodeTables, view.value?.id));

const tableSettings = usePlDataTableSettingsV2({
  model: () => tableStatus.value,
  sourceId: () => view.value?.id,
});

// Local on purpose: a selection belongs to this view, not to the block's data.
const selection = ref<PlSelectionModel>(createPlSelectionModel());
const selected = computed(() => {
  const columns = app.model.outputs.treeNodeColumns;
  return columns ? selectedNodes(selection.value, columns) : [];
});
const basket = useAddToBasket();
const addSelected = () => {
  const v = view.value;
  if (!v) return;
  // Only this tree's rows: the page is reused, so a selection may hold another lineage's ids.
  const ids = selected.value.filter((n) => n.lineageId === v.lineageId).map((n) => n.nodeId);
  void basket.start(v.lineageId, v.state.title, ids);
};
watch(
  () => view.value?.id,
  () => (selection.value = createPlSelectionModel()),
);

// Must be plural `selectedFilterValues`: the singular is ignored and the dendro draws every lineage.
const fixedOptions = computed(() => {
  const axis = app.model.outputs.lineageAxisSpec;
  const v = view.value;
  if (!axis || !v) return undefined;
  return [
    { inputName: "filters", selectedSource: axis, selectedFilterValues: [v.lineageId] },
  ] as PredefinedGraphOption<"dendro">[];
});
</script>

<template>
  <!-- The page header needs a title and GraphMaker draws its own, so the graph puts its
  buttons in its title line and the table in the page's. -->
  <PlBlockPage :no-body-gutters="tab === 'graph'">
    <template v-if="view && tab === 'table'" #title>{{ view.state.title }}</template>
    <template v-if="view && tab === 'table'" #append>
      <PlBtnGhost v-if="selected.length > 0" icon="add" @click.stop="addSelected">
        Add to basket
      </PlBtnGhost>
      <PlBtnGhost icon="graph" @click.stop="tab = 'graph'">Go to Graph</PlBtnGhost>
      <PlBtnGhost icon="close" @click.stop="close">Close</PlBtnGhost>
    </template>
    <div v-if="missing">This tree is no longer open.</div>
    <PlAlert v-if="basket.error.value" type="error">{{ basket.error.value }}</PlAlert>
    <PlAlert v-if="pathError" type="error">{{ pathError }}</PlAlert>
    <template v-if="view">
      <GraphMaker
        v-if="tab === 'graph'"
        v-model="view.state"
        chart-type="dendro"
        :p-frame="app.model.outputs.treeNodesPf"
        :default-options="defaultOptions"
        :fixed-options="fixedOptions"
        :tooltip-buttons="TOOLTIP_BUTTONS"
        @tooltip-btn-click="onTooltipButton"
      >
        <template #titleLineSlot>
          <PlBtnGhost :style="{ marginLeft: '12px' }" icon="table" @click.stop="tab = 'table'">
            Go to Table
          </PlBtnGhost>
          <PlBtnGhost icon="close" @click.stop="close">Close</PlBtnGhost>
        </template>
      </GraphMaker>
      <PlAgDataTableV2
        v-else
        v-model="tableState"
        v-model:selection="selection"
        :settings="tableSettings"
        show-export-button
        show-columns-panel
      />
    </template>
    <AddToBasketModal
      v-if="basket.pending.value.length > 0"
      :nodes="basket.pending.value"
      @close="basket.close"
    />
  </PlBlockPage>
</template>
