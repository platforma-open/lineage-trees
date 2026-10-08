import type { PlId } from "@platforma-sdk/model";

// Sequences of interest, searched with mitool like mixcr-shm-trees: top hit per node and lineage.

export type SequenceOfInterest = {
  id: PlId;
  name: string;
  sequence: string;
};

export type KnownTreeSearchParameters =
  | "oneMismatch"
  | "oneMismatchOrIndel"
  | "twoMismatches"
  | "twoMismatchesOrIndels"
  | "threeMismatches"
  | "threeMismatchesOrIndels"
  | "fourMismatches"
  | "fourMismatchesOrIndels";

export type SearchParametersTreeSearchTop = {
  type: "tree_search_top";
  parameters: KnownTreeSearchParameters;
};

export type SearchPresetTop = {
  type: "preset_alignment_search_top";
  dissimilarityPercent: number;
};

export type SearchParameters = SearchParametersTreeSearchTop | SearchPresetTop;

export type Alphabet = "nucleotide" | "amino-acid";
export type TargetFeature = "CDR3" | "VDJRegion";

export type SOIListParameters = {
  id: PlId;
  name: string;
  type: Alphabet;
  targetFeature: TargetFeature;
  searchParameters: SearchParameters;
};

export type SOIList = {
  parameters: SOIListParameters;
  sequences: SequenceOfInterest[];
};
