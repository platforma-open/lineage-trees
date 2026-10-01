import type { AnyLogHandle } from "@platforma-sdk/model";
import { isLiveLog } from "@platforma-sdk/model";
import { PROGRESS_PREFIX } from "@platforma-open/milaboratories.lineage-trees.model";
import { computed } from "vue";
import { useApp } from "./app";

// The overview's rows, one per donor group. Shared by the grid and the log drawer.

/** The stages a donor shows, each a bar of its own and a log of its own. */
export type Stage = "alleles" | "clustering" | "trees" | "results";
export const STAGES: { key: Stage; label: string }[] = [
  { key: "alleles", label: "Allele inference" },
  { key: "clustering", label: "Lineage clustering" },
  { key: "trees", label: "Tree reconstruction" },
  // Run-wide: merged once every donor has its trees, then saved as the block's outputs.
  { key: "results", label: "Results" },
];

/** The stages each donor runs on its own; `results` is shared by all of them. */
type DonorStage = Exclude<Stage, "results">;

export type StageProgress = {
  status: "not_started" | "running" | "done";
  /** What the bar says: the tool's current step while running. */
  text: string;
  percent?: number;
};

export type DonorRow = {
  donor: string;
  progress: Record<Stage, StageProgress>;
  /** The stage whose log is the one to open first. */
  currentStage: Stage;
  logs: Partial<Record<Stage, AnyLogHandle>>;
  clonotypes?: number;
  lineages?: number;
};

// "[0:01:23] Trees: 42.5%": the elapsed-time stamp, the step, then its percentage if any.
const PROGRESS_LINE = /^(?:\[[0-9:]+\]\s*)?(?<step>.*?)(?::\s*(?<progress>[0-9.]+)%)?$/;

const byDonor = <T>(data: { key: unknown[]; value: T }[] | undefined) =>
  new Map((data ?? []).map((entry) => [String(entry.key[0]), entry.value]));

function parseLine(line: string | undefined, fallback: string) {
  const text = (line ?? "").replace(PROGRESS_PREFIX, "").trim();
  const match = text.match(PROGRESS_LINE)?.groups;
  const step = match?.step?.trim() || fallback;
  const percent = match?.progress === undefined ? undefined : Number(match.progress);
  return { step, percent };
}

export function useDonorRows() {
  const app = useApp();
  return computed<DonorRow[]>(() => {
    const out = app.model.outputs;
    const logs = {
      alleles: byDonor(out.allelesLogs?.data),
      alignments: byDonor(out.alignmentsLogs?.data),
      clustering: byDonor(out.clusteringLogs?.data),
      trees: byDonor(out.treesLogs?.data),
    };
    const lines = {
      alleles: byDonor(out.allelesProgress?.data),
      clustering: byDonor(out.clusteringProgress?.data),
      trees: byDonor(out.treesProgress?.data),
    };
    const collectLog = out.collectLog;
    const isRunning = out.isRunning === true;
    const stats = new Map((out.donorStats ?? []).map((s) => [s.donor, s]));
    const donors = [
      ...new Set([...Object.values(logs).flatMap((m) => [...m.keys()]), ...stats.keys()]),
    ].sort();
    const allTreesDone = donors.every((d) => {
      const h = logs.trees.get(d);
      return h !== undefined && !isLiveLog(h);
    });

    return donors.map((donor) => {
      const handle = (stage: keyof typeof logs) => logs[stage].get(donor);
      const live = (stage: keyof typeof logs) => {
        const h = handle(stage);
        return h !== undefined && isLiveLog(h);
      };

      // `waiting`: the prior stage is done but this log has not started, so the backend is
      // scheduling it; show it as running, not queued.
      const stage = (key: DonorStage, waiting: boolean, idleText: string): StageProgress => {
        const h = handle(key);
        if (h !== undefined && !isLiveLog(h)) return { status: "done", text: "Done", percent: 100 };
        if (h === undefined && !waiting) return { status: "not_started", text: "Queued" };
        if (h === undefined) return { status: "running", text: idleText };
        const { step, percent } = parseLine(lines[key].get(donor), idleText);
        return { status: "running", text: step, percent };
      };

      const done = (key: keyof typeof logs) => handle(key) !== undefined && !live(key);
      const results = (): StageProgress => {
        if (collectLog !== undefined && !isLiveLog(collectLog)) {
          // The backend then imports the tables as columns, which reports no progress.
          return isRunning
            ? { status: "running", text: "Saving results" }
            : { status: "done", text: "Done", percent: 100 };
        }
        if (collectLog !== undefined) {
          const { step, percent } = parseLine(out.collectProgress, "Collecting results");
          return { status: "running", text: step, percent };
        }
        if (!done("trees")) return { status: "not_started", text: "Queued" };
        // Projects run before the collect log existed have none to show.
        if (!isRunning) return { status: "done", text: "Done", percent: 100 };
        return { status: "running", text: allTreesDone ? "Starting" : "Waiting for other donors" };
      };

      const progress: Record<Stage, StageProgress> = {
        // Alignments are joined first and allele inference reads them, so joining shows here.
        alleles: stage("alleles", handle("alignments") !== undefined, "Joining alignments"),
        clustering: stage("clustering", done("alleles"), "Lineage clustering"),
        trees: stage("trees", done("clustering"), "Starting"),
        results: results(),
      };
      // The tree tool labels its build "Trees"; say what it is doing.
      if (progress.trees.status === "running") {
        progress.trees.text = progress.trees.text.replace(/^Trees\b/, "Building trees");
      }

      const handles: Partial<Record<Stage, AnyLogHandle>> = {};
      for (const { key } of STAGES) {
        const h = key === "results" ? collectLog : handle(key);
        if (h !== undefined) handles[key] = h;
      }
      const running = STAGES.find(({ key }) => progress[key].status === "running");

      return {
        donor,
        progress,
        currentStage: running?.key ?? (done("trees") ? "trees" : "alleles"),
        logs: handles,
        clonotypes: stats.get(donor)?.clonotype_count,
        lineages: stats.get(donor)?.lineage_count,
      } satisfies DonorRow;
    });
  });
}
