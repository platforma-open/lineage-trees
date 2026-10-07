import type { BlockData } from "@platforma-open/milaboratories.lineage-trees.model";
import { defaultBlockData } from "@platforma-open/milaboratories.lineage-trees.model";
import { ImportVdjBlockPointer } from "@platforma-open/milaboratories.import-vdj";
import { MixcrClonotyping2BlockPointer } from "@platforma-open/milaboratories.mixcr-clonotyping-2";
import { SamplesAndDataBlockPointer } from "@platforma-open/milaboratories.samples-and-data";
import type { PlRef } from "@platforma-sdk/model";
import { createPlDataTableStateV2, uniquePlId } from "@platforma-sdk/model";
import { awaitStableState, blockTest } from "@platforma-sdk/test";
import { LineageTreesBlockPointer } from "this-block";

/** The pieces of a blockTest's fixture the helpers need. */
export type TestCtx = Pick<
  Parameters<NonNullable<Parameters<typeof blockTest>[2]>>[0],
  "rawPrj" | "helpers" | "expect" | "ml"
>;

type Output = { ok: boolean; value?: unknown; errors?: unknown[] };
export type Outputs = Record<string, Output | undefined>;

export const UPSTREAM_TIMEOUT = 600_000;
export const LINEAGE_TIMEOUT = 1_500_000;

/** Unwraps a model output; undefined when absent or errored. */
export function value<T>(outputs: Outputs, name: string): T | undefined {
  const out = outputs[name];
  return out?.ok ? (out.value as T | undefined) : undefined;
}

async function stableOutputs(ctx: TestCtx, blockId: string, timeout = 100_000): Promise<Outputs> {
  const state = (await awaitStableState(ctx.rawPrj.getBlockState(blockId), timeout)) as {
    outputs?: Outputs;
  };
  return state.outputs ?? {};
}

export type FastqSample = { label: string; r1: string; r2: string };

/** A Samples & Data block with FASTQ and TSV datasets and optional metadata. */
export async function addSamples(
  ctx: TestCtx,
  opts: {
    fastq?: { label: string; samples: FastqSample[] }[];
    tsv?: { label: string; sample: string; asset: string }[];
    metadata?: { label: string; values: Record<string, string> };
  },
): Promise<{ blockId: string; sampleIds: Record<string, string> }> {
  const { rawPrj: project, helpers } = ctx;
  const blockId = await project.addBlock("Samples & Data", SamplesAndDataBlockPointer);
  const sampleIds: Record<string, string> = {};
  const idOf = (label: string) => (sampleIds[label] ??= uniquePlId());
  const datasets = [];
  for (const ds of opts.fastq ?? []) {
    const data: Record<string, unknown> = {};
    for (const s of ds.samples) {
      data[idOf(s.label)] = {
        R1: await helpers.getLocalFileHandle(`./assets/${s.r1}`),
        R2: await helpers.getLocalFileHandle(`./assets/${s.r2}`),
      };
    }
    datasets.push({
      id: uniquePlId(),
      label: ds.label,
      content: { type: "Fastq", readIndices: ["R1", "R2"], gzipped: true, data },
    });
  }
  for (const ds of opts.tsv ?? []) {
    datasets.push({
      id: uniquePlId(),
      label: ds.label,
      content: {
        type: "Xsv",
        xsvType: "tsv",
        data: { [idOf(ds.sample)]: await helpers.getLocalFileHandle(`./assets/${ds.asset}`) },
      },
    });
  }
  const metadata =
    opts.metadata === undefined
      ? []
      : [
          {
            id: uniquePlId(),
            label: opts.metadata.label,
            global: false,
            valueType: "String",
            data: Object.fromEntries(
              Object.entries(opts.metadata.values).map(([sample, v]) => [idOf(sample), v]),
            ),
          },
        ];
  // Facade PlId brands do not unify with ours, so the literal is untyped.
  await project.mutateBlockStorage(blockId, {
    operation: "update-block-data",
    value: {
      metadata,
      sampleIds: Object.values(sampleIds),
      sampleLabelColumnLabel: "Sample Name",
      sampleLabels: Object.fromEntries(Object.entries(sampleIds).map(([l, id]) => [id, l])),
      datasets,
      h5adFilesToPreprocess: [],
      seuratFilesToPreprocess: [],
      suggestedImport: false,
    },
  });
  await project.runBlock(blockId);
  await helpers.awaitBlockDoneAndGetStableBlockState(blockId, UPSTREAM_TIMEOUT);
  return { blockId, sampleIds };
}

/** The option `blockId` produced, optionally matched by label. */
function pickOption(
  ctx: TestCtx,
  options: { ref: PlRef; label: string }[] | undefined,
  blockId: string,
  label?: string,
): PlRef {
  const found = (options ?? []).find(
    (o) => o.ref.blockId === blockId && (label === undefined || o.label.includes(label)),
  );
  ctx.expect(found, `option from ${blockId} among ${JSON.stringify(options)}`).toBeDefined();
  return found!.ref;
}

export type MixcrPreset = { type: "name"; name: string } | { type: "file"; asset: string };

/** A MiXCR Clonotyping block over one Samples & Data dataset, run to completion. */
export async function clonotype(
  ctx: TestCtx,
  opts: {
    sndBlockId: string;
    datasetLabel: string;
    preset: MixcrPreset;
    species?: string;
    /** Germline library uploaded with the block, from `./assets`. */
    libraryAsset?: string;
    chains: string[];
    label: string;
  },
): Promise<string> {
  const { rawPrj: project, helpers } = ctx;
  const blockId = await project.addBlock("MiXCR Clonotyping", MixcrClonotyping2BlockPointer);
  const before = await stableOutputs(ctx, blockId);
  const input = pickOption(
    ctx,
    value<{ ref: PlRef; label: string }[]>(before, "inputOptions"),
    opts.sndBlockId,
    opts.datasetLabel,
  );
  const preset =
    opts.preset.type === "name"
      ? opts.preset
      : { type: "file", file: await helpers.getLocalFileHandle(`./assets/${opts.preset.asset}`) };
  await project.mutateBlockStorage(blockId, {
    operation: "update-block-data",
    value: {
      defaultBlockLabel: opts.label,
      customBlockLabel: opts.label,
      input,
      preset,
      species: opts.species,
      libraryFile:
        opts.libraryAsset === undefined
          ? undefined
          : await helpers.getLocalFileHandle(`./assets/${opts.libraryAsset}`),
      // An uploaded library takes its species from here, not from `species`.
      customSpecies: opts.libraryAsset === undefined ? undefined : opts.species,
      isLibraryFileGzipped: opts.libraryAsset?.endsWith(".gz"),
      chains: opts.chains,
      cloneClusteringMode: "default",
      runMode: "full",
      tableState: createPlDataTableStateV2(),
    },
  });
  await project.runBlock(blockId);
  await helpers.awaitBlockDoneAndGetStableBlockState(blockId, UPSTREAM_TIMEOUT);
  return blockId;
}

/** An Import V(D)J Data block over one Samples & Data dataset, run to completion. */
export async function importVdj(
  ctx: TestCtx,
  opts: {
    sndBlockId: string;
    datasetLabel: string;
    format: "airr" | "custom";
    bareSet?: Record<string, unknown> & { identity: string };
  },
): Promise<string> {
  const { rawPrj: project, helpers } = ctx;
  const blockId = await project.addBlock("Import V(D)J Data", ImportVdjBlockPointer);
  const before = await stableOutputs(ctx, blockId);
  const datasetRef = pickOption(
    ctx,
    value<{ ref: PlRef; label: string }[]>(before, "datasetOptions"),
    opts.sndBlockId,
    opts.datasetLabel,
  );
  // Stands in for the UI watcher that writes prerun's verdict.
  const prerunCheck =
    opts.bareSet === undefined
      ? {
          check: "dataset",
          subject: [datasetRef.blockId, datasetRef.name, opts.format].join("\0"),
          columnsPresent: true,
        }
      : { check: "columns", subject: opts.bareSet.identity, identityCollides: false };
  await project.mutateBlockStorage(blockId, {
    operation: "update-block-data",
    value: {
      defaultBlockLabel: opts.datasetLabel,
      customBlockLabel: "",
      datasetRef,
      format: opts.format,
      chains: ["IGHeavy"],
      primaryCountType: "read",
      bareSet: opts.bareSet,
      tableState: createPlDataTableStateV2(),
      settingsOpen: true,
      prerunCheck,
    },
  });
  await project.runBlock(blockId);
  await helpers.awaitBlockDoneAndGetStableBlockState(blockId, UPSTREAM_TIMEOUT);
  return blockId;
}

/** Full block data, since `update-block-data` replaces it whole. */
export function lineageData(fields: Partial<BlockData>): BlockData {
  return { ...defaultBlockData(), ...fields };
}

/** Names of outputs that carry an error. */
export function erroredOutputs(outputs: Outputs): string[] {
  return Object.entries(outputs)
    .filter(([, out]) => out !== undefined && !out.ok)
    .map(([name, out]) => `${name}: ${JSON.stringify(out!.errors).slice(0, 2000)}`);
}

/** Last lines of every per-donor stage log, for the failure report. */
async function stageLogs(ctx: TestCtx, outputs: Outputs): Promise<string> {
  const parts: string[] = [];
  for (const name of ["allelesLogs", "alignmentsLogs", "clusteringLogs", "treesLogs"]) {
    const map = value<{ data: { key: unknown[]; value?: string }[] }>(outputs, name);
    for (const entry of map?.data ?? []) {
      if (entry.value === undefined) continue;
      try {
        const res = (await ctx.ml.driverKit.logDriver.lastLines(
          entry.value as Parameters<typeof ctx.ml.driverKit.logDriver.lastLines>[0],
          40,
        )) as { logs?: string };
        parts.push(`--- ${name} ${JSON.stringify(entry.key)}\n${res.logs ?? ""}`);
      } catch (e) {
        parts.push(`--- ${name} ${JSON.stringify(entry.key)}: unreadable (${String(e)})`);
      }
    }
  }
  for (const name of ["mergeLog", "clusteringLog", "treesLog"]) {
    const text = value<string>(outputs, name);
    if (text) parts.push(`--- ${name}\n${text.slice(-4000)}`);
  }
  return parts.join("\n");
}

/** Adds and runs Lineage Trees on the given datasets. On error, prints stage logs and fails. */
export async function runLineageTrees(
  ctx: TestCtx,
  opts: {
    /** Upstream block ids whose datasets to pick, in order. */
    from: string[];
    /** Upstream block ids whose datasets are known antibody datasets. */
    known?: string[];
    /** Upstream block id and label of the donor metadata column. */
    donor?: { blockId: string; label: string };
    data?: Partial<BlockData>;
    /** Tells apart several blocks in one project. */
    label?: string;
  },
): Promise<Outputs> {
  const { rawPrj: project, helpers, expect } = ctx;
  const label = opts.label ?? "Lineage Trees";
  const blockId = await project.addBlock(label, LineageTreesBlockPointer);
  const before = await stableOutputs(ctx, blockId);
  expect(erroredOutputs(before), label).toEqual([]);
  const options = value<{ ref: PlRef; label: string }[]>(before, "inputOptions");
  const datasets = opts.from.map((id) => pickOption(ctx, options, id));
  const knownDatasets = (opts.known ?? []).map((id) => pickOption(ctx, options, id));
  const donorColumn =
    opts.donor === undefined
      ? undefined
      : pickOption(
          ctx,
          value<{ ref: PlRef; label: string }[]>(before, "donorOptions"),
          opts.donor.blockId,
          opts.donor.label,
        );

  await project.mutateBlockStorage(blockId, {
    operation: "update-block-data",
    value: lineageData({ ...opts.data, datasets, knownDatasets, donorColumn }),
  });
  if (datasets.length === 0) return await stableOutputs(ctx, blockId);

  await project.runBlock(blockId);
  let runError: unknown;
  try {
    await helpers.awaitBlockDone(blockId, LINEAGE_TIMEOUT);
  } catch (e) {
    runError = e;
  }
  const outputs = await stableOutputs(ctx, blockId, 300_000);
  const errors = erroredOutputs(outputs);
  if (runError !== undefined || errors.length > 0) {
    console.error(`${label} run failed: ${String(runError)}\n${errors.join("\n")}`);
    console.error(await stageLogs(ctx, outputs));
  }
  expect(runError, label).toBeUndefined();
  expect(errors, label).toEqual([]);
  return outputs;
}

// Fixtures: alpaca VHH bulk (SRA PRJNA638614) and a mostly TCR human bulk.
export const ALPACA: FastqSample = {
  label: "SRR11974622",
  r1: "SRR11974622_head1000_1.fastq.gz",
  r2: "SRR11974622_head1000_2.fastq.gz",
};
export const HUMAN_TCR_BULK: FastqSample = {
  label: "SRR11233652",
  r1: "SRR11233652_sampledBulk_R1.fastq.gz",
  r2: "SRR11233652_sampledBulk_R2.fastq.gz",
};
export const ALPACA_PRESET: MixcrPreset = { type: "file", asset: "mixcr-preset-alpaca.yaml" };

/** The alpaca sample, clonotyped `times` times. */
export async function alpacaClonotyped(
  ctx: TestCtx,
  times = 1,
  metadata?: { label: string; values: Record<string, string> },
): Promise<{ sndBlockId: string; mixcrBlockIds: string[] }> {
  const { blockId: sndBlockId } = await addSamples(ctx, {
    fastq: [{ label: "Alpaca", samples: [ALPACA] }],
    metadata,
  });
  const mixcrBlockIds: string[] = [];
  for (let i = 0; i < times; i++) {
    mixcrBlockIds.push(
      await clonotype(ctx, {
        sndBlockId,
        datasetLabel: "Alpaca",
        preset: ALPACA_PRESET,
        species: "alpaca",
        chains: ["IGHeavy"],
        label: `Alpaca ${i + 1}`,
      }),
    );
  }
  return { sndBlockId, mixcrBlockIds };
}
