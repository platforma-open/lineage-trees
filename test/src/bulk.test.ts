/* Bulk inputs: each shape must run without error and give sane outputs. No single cell. */

import type { DonorStats, SOIList } from "@platforma-open/milaboratories.lineage-trees.model";
import { uniquePlId } from "@platforma-sdk/model";
import { blockTest } from "@platforma-sdk/test";
import type { Outputs, TestCtx } from "./helpers";
import {
  ALPACA,
  ALPACA_PRESET,
  HUMAN_TCR_BULK,
  addSamples,
  alpacaClonotyped,
  clonotype,
  importVdj,
  runLineageTrees,
  value,
} from "./helpers";

const TIMEOUT = 2_400_000;

type TreeNodeColumns = { hasAnchorProperty: boolean; topologyId: string };

function expectTrees({ expect }: TestCtx, outputs: Outputs, label?: string) {
  expect(value(outputs, "noTreesReason"), label).toBeUndefined();
  expect(value<TreeNodeColumns>(outputs, "treeNodeColumns"), label).toBeDefined();
}

function expectClustered({ expect }: TestCtx, outputs: Outputs, label?: string) {
  const stats = value<DonorStats[]>(outputs, "donorStats") ?? [];
  expect(stats.length, label).toBeGreaterThan(0);
  expect(
    stats.reduce((n, s) => n + s.clonotype_count, 0),
    label,
  ).toBeGreaterThan(0);
  expect(
    stats.reduce((n, s) => n + s.lineage_count, 0),
    label,
  ).toBeGreaterThan(0);
}

// IGH AIRR tables with duplicate_count, which Import V(D)J requires.
const AIRR_BULK = [
  { label: "airr-mixcr", sample: "airr-mixcr", asset: "airr-mixcr.tsv" },
  { label: "airr-igblast", sample: "airr-igblast", asset: "airr-igblast.tsv" },
];

blockTest("no datasets picked", { timeout: 120_000 }, async ({ rawPrj, helpers, expect, ml }) => {
  const ctx = { rawPrj, helpers, expect, ml };
  const outputs = await runLineageTrees(ctx, { from: [] });
  ctx.expect(value(outputs, "inputOptions")).toEqual([]);
  ctx.expect(value(outputs, "datasets")).toEqual([]);
});

blockTest(
  "alpaca twice, one an anchor set",
  { timeout: TIMEOUT },
  async ({ rawPrj, helpers, expect, ml }) => {
    const ctx = { rawPrj, helpers, expect, ml };
    const { mixcrBlockIds } = await alpacaClonotyped(ctx, 2);
    const outputs = await runLineageTrees(ctx, {
      from: mixcrBlockIds,
      anchors: [mixcrBlockIds[0]],
    });
    expectTrees(ctx, outputs);
    ctx.expect(value<TreeNodeColumns>(outputs, "treeNodeColumns")?.hasAnchorProperty).toBe(true);
    ctx.expect(value<string>(outputs, "modeStatement")).toContain("1 dataset as anchor sets");
  },
);

blockTest("imported AIRR", { timeout: TIMEOUT }, async ({ rawPrj, helpers, expect, ml }) => {
  const ctx = { rawPrj, helpers, expect, ml };
  const { blockId: sndBlockId } = await addSamples(ctx, { tsv: AIRR_BULK });
  const imports = [];
  for (const ds of AIRR_BULK) {
    imports.push(await importVdj(ctx, { sndBlockId, datasetLabel: ds.label, format: "airr" }));
  }
  const outputs = await runLineageTrees(ctx, { from: imports });
  // Import V(D)J drops alignment columns, so there are no trees.
  ctx
    .expect(value<string>(outputs, "modeStatement"))
    .toContain("2 datasets without any, so no trees");
  expectClustered(ctx, outputs);
});

// A bare antibody set keys on pl7.app/variantKey, not a clonotype key, so it is not offered.
blockTest(
  "imported bare heavy set is not offered",
  { timeout: TIMEOUT },
  async ({ rawPrj, helpers, expect, ml }) => {
    const ctx = { rawPrj, helpers, expect, ml };
    const { blockId: sndBlockId } = await addSamples(ctx, {
      tsv: [{ label: "bare-heavy", sample: "bare-heavy", asset: "bare-heavy-only.tsv" }],
    });
    const importId = await importVdj(ctx, {
      sndBlockId,
      datasetLabel: "bare-heavy",
      format: "custom",
      bareSet: {
        identity: "mAb ID",
        chainSelection: "IGHeavy",
        sequences: { IGHeavy: "VH" },
        scheme: "imgt",
      },
    });
    const outputs = await runLineageTrees(ctx, { from: [] });
    const options = value<{ ref: { blockId: string } }[]>(outputs, "inputOptions") ?? [];
    ctx.expect(options.filter((o) => o.ref.blockId === importId)).toEqual([]);
  },
);

blockTest(
  "alpaca MiXCR plus imported AIRR",
  { timeout: TIMEOUT },
  async ({ rawPrj, helpers, expect, ml }) => {
    const ctx = { rawPrj, helpers, expect, ml };
    const { blockId: sndBlockId } = await addSamples(ctx, {
      fastq: [{ label: "Alpaca", samples: [ALPACA] }],
      tsv: [AIRR_BULK[0]],
    });
    const mixcrId = await clonotype(ctx, {
      sndBlockId,
      datasetLabel: "Alpaca",
      preset: ALPACA_PRESET,
      species: "alpaca",
      chains: ["IGHeavy"],
      label: "Alpaca",
    });
    const importId = await importVdj(ctx, {
      sndBlockId,
      datasetLabel: AIRR_BULK[0].label,
      format: "airr",
    });
    const outputs = await runLineageTrees(ctx, { from: [mixcrId, importId] });
    const statement = value<string>(outputs, "modeStatement");
    ctx.expect(statement).toContain("from the MiXCR clns files");
    expectClustered(ctx, outputs);
  },
);

blockTest("near-empty IG input", { timeout: TIMEOUT }, async ({ rawPrj, helpers, expect, ml }) => {
  const ctx = { rawPrj, helpers, expect, ml };
  const { blockId: sndBlockId } = await addSamples(ctx, {
    fastq: [{ label: "Human bulk", samples: [HUMAN_TCR_BULK] }],
  });
  const mixcrId = await clonotype(ctx, {
    sndBlockId,
    datasetLabel: "Human bulk",
    preset: { type: "name", name: "neb-human-rna-xcr-umi-nebnext" },
    chains: ["IGHeavy"],
    label: "Human bulk",
  });
  const outputs = await runLineageTrees(ctx, { from: [mixcrId] });
  ctx
    .expect(value(outputs, "donorStats"))
    .toMatchObject([{ clonotype_count: 1, lineage_count: 1 }]);
});

blockTest(
  "donor column from sample metadata",
  { timeout: TIMEOUT },
  async ({ rawPrj, helpers, expect, ml }) => {
    const ctx = { rawPrj, helpers, expect, ml };
    const { sndBlockId, mixcrBlockIds } = await alpacaClonotyped(ctx, 1, {
      label: "Donor",
      values: { [ALPACA.label]: "Ty1" },
    });
    const outputs = await runLineageTrees(ctx, {
      from: mixcrBlockIds,
      donor: { blockId: sndBlockId, label: "Donor" },
    });
    expectTrees(ctx, outputs);
    const stats = value<DonorStats[]>(outputs, "donorStats") ?? [];
    ctx.expect(stats.map((s) => s.donor).join(",")).toContain("Ty1");
  },
);

// A CDR3 shared by five alpaca clonotypes, searched in both modes.
const SHARED_CDR3 = "TGCGCGGCAGATGGCATCCCCCTGGGTAACTGTCTGGATTACAAGGACATGGACTATTGG";

function soiList(name: string, searchParameters: SOIList["parameters"]["searchParameters"]) {
  return {
    parameters: {
      id: uniquePlId(),
      name,
      type: "nucleotide" as const,
      targetFeature: "CDR3" as const,
      searchParameters,
    },
    sequences: [{ id: uniquePlId(), name: "shared", sequence: SHARED_CDR3 }],
  };
}

// One clonotyping, then one Lineage Trees block per setting.
blockTest(
  "alpaca, one clonotyping under several settings",
  { timeout: 2 * TIMEOUT },
  async ({ rawPrj, helpers, expect, ml }) => {
    const ctx = { rawPrj, helpers, expect, ml };
    const { mixcrBlockIds } = await alpacaClonotyped(ctx);
    const from = mixcrBlockIds;

    let label = "default settings";
    let outputs = await runLineageTrees(ctx, { from, label });
    expectTrees(ctx, outputs, label);
    expectClustered(ctx, outputs, label);
    ctx
      .expect(value<string>(outputs, "modeStatement"), label)
      .toContain("from the MiXCR clns files");

    label = "min tips above every lineage";
    outputs = await runLineageTrees(ctx, { from, label, data: { minTipsPerTree: 100 } });
    ctx.expect(value<string>(outputs, "noTreesReason"), label).toBeDefined();
    expectClustered(ctx, outputs, label);

    label = "max tips 3";
    outputs = await runLineageTrees(ctx, { from, label, data: { maxTipsPerTree: 3 } });
    expectTrees(ctx, outputs, label);

    label = "adaptive clustering";
    outputs = await runLineageTrees(ctx, { from, label, data: { clusteringMode: "adaptive" } });
    ctx.expect(value<string>(outputs, "modeStatement"), label).toContain("adaptive clustering");
    expectClustered(ctx, outputs, label);

    label = "sequence search in both modes";
    const alignment = soiList("alignment", {
      type: "preset_alignment_search_top",
      dissimilarityPercent: 2,
    });
    const tree = soiList("tree", { type: "tree_search_top", parameters: "oneMismatch" });
    outputs = await runLineageTrees(ctx, {
      from,
      label,
      data: { sequencesOfInterest: [alignment, tree] },
    });
    ctx.expect(value<boolean>(outputs, "soiReady"), label).toBe(true);
    const hits = value<Record<string, { nodes: number; lineages: number }>>(outputs, "soiHits");
    ctx.expect(hits?.[alignment.parameters.id]?.lineages, label).toBeGreaterThan(0);
    ctx.expect(hits?.[tree.parameters.id]?.lineages, label).toBeGreaterThan(0);
  },
);
