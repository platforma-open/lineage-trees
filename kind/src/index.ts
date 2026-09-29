import { assertParamsObject, defineBlockKind } from "@platforma-sdk/block-kind";
import { name, version } from "../package.json" with { type: "json" };

/**
 * This block's init-params contract.
 *
 * Only the clustering threshold is templatable. The input dataset is not: it is
 * a `PlRef` into one project's result pool, so it carries no meaning in a
 * template meant to seed a fresh project.
 */
export type BlockParams = {
  clusteringThreshold?: number;
};

function parseInitializationParams(value: unknown): BlockParams {
  assertParamsObject(value);

  const { clusteringThreshold } = value;
  if (clusteringThreshold === undefined) return {};
  if (typeof clusteringThreshold !== "number" || !Number.isFinite(clusteringThreshold)) {
    throw new Error("'clusteringThreshold' must be a finite number when present.");
  }

  return { clusteringThreshold };
}

// Identity (`name`/`version`) comes from this package's own `package.json`, so
// the on-wire `{name}@{version}` reference can never drift from what npm
// publishes; the bundler inlines the JSON import.
export const kind = defineBlockKind<BlockParams>({
  name,
  version,
  parseInitializationParams,
});
