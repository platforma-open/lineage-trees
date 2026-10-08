import { platforma } from "@platforma-open/milaboratories.lineage-trees.model";
import { defineAppV3 } from "@platforma-sdk/ui-vue";
import BasketPage from "./BasketPage.vue";
import MainPage from "./MainPage.vue";
import MutationalPathPage from "./MutationalPathPage.vue";
import SOIPage from "./soi/SOIPage.vue";
import TreePage from "./TreePage.vue";
import TreesPage from "./TreesPage.vue";

export const sdkPlugin = defineAppV3(platforma, (app) => ({
  // The green bar while the block runs.
  progress: () => app.model.outputs.isRunning,
  routes: {
    "/": () => MainPage,
    "/trees": () => TreesPage,
    "/tree": () => TreePage,
    "/path": () => MutationalPathPage,
    "/soi": () => SOIPage,
    "/basket": () => BasketPage,
  },
}));

export const useApp = sdkPlugin.useApp;
