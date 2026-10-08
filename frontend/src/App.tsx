import { useEffect, useMemo, useRef, useState } from "react";
import { api, download } from "./client";
import { latestSettings, object } from "./projectSettings";
import { BenchmarkResult } from "./BenchmarkResult";
import type {
  Analysis,
  BenchmarkReport,
  Catalog,
  Feature,
  Help,
  Job,
  JobPayload,
  Json,
  Policy,
  Project,
} from "./types";

type Page =
  | "Input"
  | "Sequence"
  | "Design"
  | "RBS"
  | "Assembly"
  | "Reports"
  | "Reliability";
const pages: Page[] = [
  "Input",
  "Sequence",
  "Design",
  "RBS",
  "Assembly",
  "Reports",
  "Reliability",
];
const fmt = (value: number | null | undefined, digits = 3) =>
  value == null
    ? "Unassessed"
    : value.toLocaleString(undefined, { maximumFractionDigits: digits });
const emptyFeature = {
  id: "",
  name: "",
  feature_type: "protected region",
  start_1based: 1,
  end_1based: 1,
  strand: 0,
  mutability: "LOCKED",
  translation_table: 11,
  cds_mode: "complete",
  frame_offset: 0,
  rbs_mask: "",
};

function HelpBox({ entry }: { entry: Help | undefined }) {
  return entry ? (
    <details className="help">
      <summary>{entry.title} · help</summary>
      <p>{entry.definition}</p>
      <p>{entry.effect}</p>
      <p>{entry.limitations}</p>
      {entry.reference.startsWith("https:") ? (
        <a href={entry.reference} target="_blank" rel="noreferrer">
          Scientific reference ↗
        </a>
      ) : (
        <span>Reference: {entry.reference}</span>
      )}
    </details>
  ) : null;
}
function Field({
  label,
  children,
}: {
  label: string;
  children: React.ReactNode;
}) {
  return (
    <label className="field">
      <span>{label}</span>
      {children}
    </label>
  );
}
function ScientificResult({ payload }: { payload: JobPayload | null }) {
  if (!payload)
    return (
      <p className="muted">
        Run a design to see model results and residual target distance.
      </p>
    );
  const result = payload.result as Record<string, Json>;
  const rows = Array.isArray(payload.result)
    ? payload.result
    : [payload.result];
  return (
    <div>
      {rows.map((item, i) => {
        const row = item as Record<string, Json>;
        const predictions =
          (row.predictions as Record<string, Json>[] | undefined) ??
          (row.relative_tir ? [row] : []);
        return (
          <div className="model-result" key={i}>
            {row.target_relative_tir && (
              <h3>Target {fmt(Number(row.target_relative_tir))}</h3>
            )}
            {predictions.map((prediction, j) => (
              <div key={j}>
                <div className="metric-strip">
                  <div>
                    <small>Relative model TIR</small>
                    <strong>{fmt(Number(prediction.relative_tir))}</strong>
                  </div>
                  <div>
                    <small>Spacing</small>
                    <strong>{String(prediction.spacing_bp)} bp</strong>
                  </div>
                  <div>
                    <small>Log-fold target distance</small>
                    <strong>
                      {fmt((row.log_fold_errors as number[] | undefined)?.[j])}
                    </strong>
                  </div>
                </div>
                <table>
                  <tbody>
                    {Object.entries(
                      (prediction.energy_components as Record<
                        string,
                        number
                      >) ?? {},
                    ).map(([key, value]) => (
                      <tr key={key}>
                        <td>{key}</td>
                        <td>{fmt(value)} kcal/mol</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
                <details>
                  <summary>Unbound RNA structure and model provenance</summary>
                  <pre>{String(prediction.unbound_mrna_structure ?? "")}</pre>
                  <pre>
                    {JSON.stringify(prediction.model_versions, null, 2)}
                  </pre>
                </details>
                {((prediction.warnings as string[]) ?? []).map((note) => (
                  <p className="note" key={note}>
                    {note}
                  </p>
                ))}
              </div>
            ))}
          </div>
        );
      })}
      {result.status && (
        <p className="note">Search status: {String(result.status)}</p>
      )}
      <details>
        <summary>Complete result and settings</summary>
        <pre>{JSON.stringify(payload.result, null, 2)}</pre>
      </details>
    </div>
  );
}

export default function App() {
  const [page, setPage] = useState<Page>("Input");
  const [catalog, setCatalog] = useState<Catalog | null>(null);
  const [project, setProject] = useState<Project | null>(null);
  const [analysis, setAnalysis] = useState<Analysis | null>(null);
  const [host, setHost] = useState("ecoli-mg1655");
  const [error, setError] = useState("");
  const [notice, setNotice] = useState("");
  const [loading, setLoading] = useState(false);
  const [job, setJob] = useState<Job | null>(null);
  const [payload, setPayload] = useState<JobPayload | null>(null);
  const [advanced, setAdvanced] = useState(false);
  const [text, setText] = useState("");
  const [format, setFormat] = useState("raw");
  const [name, setName] = useState("My construct");
  const [topology, setTopology] = useState("linear");
  const [wholeCDS, setWholeCDS] = useState(false);
  const [seed, setSeed] = useState(7);
  const [iterations, setIterations] = useState(500);
  const [alternatives, setAlternatives] = useState(1);
  const [nTerminal, setNTerminal] = useState(30);
  const [algorithm, setAlgorithm] = useState("simulated_annealing");
  const [gcTarget, setGCTarget] = useState("");
  const [rareThreshold, setRareThreshold] = useState("");
  const [codonWeight, setCodonWeight] = useState(1);
  const [synthesisWeight, setSynthesisWeight] = useState(1);
  const [customProfile, setCustomProfile] = useState<Record<
    string,
    Json
  > | null>(null);
  const [gcWeight, setGCWeight] = useState(0);
  const [rareWeight, setRareWeight] = useState(0);
  const [preset, setPreset] = useState("Balanced");
  const [policyID, setPolicyID] = useState("generic");
  const [loadedPolicy, setLoadedPolicy] = useState<Policy | null>(null);
  const [homopolymer, setHomopolymer] = useState("");
  const [repeatLength, setRepeatLength] = useState("");
  const [invertedLength, setInvertedLength] = useState("");
  const [gcMin, setGCMin] = useState("");
  const [gcMax, setGCMax] = useState("");
  const [gcWindow, setGCWindow] = useState(50);
  const [entropyWindow, setEntropyWindow] = useState("");
  const [entropyMin, setEntropyMin] = useState("");
  const [selection, setSelection] = useState<[number, number] | null>(null);
  const [viewerStart, setViewerStart] = useState(0);
  const anchor = useRef(0);
  const dragging = useRef(false);
  const [feature, setFeature] = useState(emptyFeature);
  const [motif, setMotif] = useState("");
  const [enzyme, setEnzyme] = useState("");
  const [motifMode, setMotifMode] = useState("AVOID");
  const [motifQuery, setMotifQuery] = useState("");
  const [selectedCDS, setSelectedCDS] = useState("");
  const [selectedRBS, setSelectedRBS] = useState("");
  const [target, setTarget] = useState(1000);
  const [targets, setTargets] = useState("250,500,1000");
  const [joint, setJoint] = useState(0);
  const [rbsSeed, setRBSSeed] = useState(7);
  const [rbsIterations, setRBSIterations] = useState(500);
  const [rbsAlternatives, setRBSAlternatives] = useState(1);
  const [protectDownstream, setProtectDownstream] = useState(30);
  const [upstreamNT, setUpstreamNT] = useState(35);
  const [downstreamNT, setDownstreamNT] = useState(35);
  const [antiSD, setAntiSD] = useState("");
  const [antiSource, setAntiSource] = useState("");
  const [uncalibrated, setUncalibrated] = useState(false);
  const [maximum, setMaximum] = useState(2500);
  const [overlap, setOverlap] = useState(30);
  const [minimum, setMinimum] = useState(60);
  const [preferredTm, setPreferredTm] = useState("");
  const [sodium, setSodium] = useState(50);
  const [magnesium, setMagnesium] = useState(0);
  const [strandConcentration, setStrandConcentration] = useState(25);
  const [strand2Concentration, setStrand2Concentration] = useState(25);
  const [potassium, setPotassium] = useState(0);
  const [tris, setTris] = useState(0);
  const [dntps, setDNTPs] = useState(0);
  const [saltCorrection, setSaltCorrection] = useState(5);
  const [avoidProtected, setAvoidProtected] = useState(false);
  const [overlapGCMin, setOverlapGCMin] = useState("");
  const [overlapGCMax, setOverlapGCMax] = useState("");
  const [boundaryWindow, setBoundaryWindow] = useState(100);
  const [benchmarkText, setBenchmarkText] = useState("");
  const [benchmarkManifest, setBenchmarkManifest] = useState("");
  const [benchmarkResult, setBenchmarkResult] = useState<{
    csv: string;
    manifest: string;
    report: BenchmarkReport;
  } | null>(null);
  // A returned report belongs to the exact inputs evaluated, including while
  // an earlier evaluation is still in flight when the inputs change.
  const currentBenchmark =
    benchmarkResult?.csv === benchmarkText &&
    benchmarkResult.manifest === benchmarkManifest
      ? benchmarkResult.report
      : null;
  const busy = loading || job?.status === "running" || job?.status === "queued";
  const help = (key: string) => <HelpBox entry={catalog?.help[key]} />;
  const guard = async (work: () => Promise<void>) => {
    setError("");
    setLoading(true);
    try {
      await work();
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    api<Catalog>("/catalog")
      .then(setCatalog)
      .catch((e) => setError(e.message));
  }, []);
  const policy = useMemo(() => {
    const base =
      loadedPolicy ??
      catalog?.synthesis_profiles.find((row) => row.id === policyID) ??
      catalog?.synthesis_profiles[0];
    if (!base) return null;
    const updates: Record<string, Json> = {};
    if (homopolymer) updates.homopolymer_limit = Number(homopolymer);
    if (repeatLength) updates.repeat_min_length = Number(repeatLength);
    if (invertedLength)
      updates.inverted_repeat_min_length = Number(invertedLength);
    if (gcMin && gcMax)
      updates.global_gc_range = {
        minimum: Number(gcMin) / 100,
        maximum: Number(gcMax) / 100,
      };
    if (entropyWindow && entropyMin) {
      updates.complexity_window = Number(entropyWindow);
      updates.minimum_entropy = Number(entropyMin);
    }
    return Object.keys(updates).length
      ? ({
          ...base,
          ...updates,
          id: "custom",
          version: "user-policy-v1",
          source: "Explicit local policy adapted from " + String(base.source),
        } as Policy)
      : base;
  }, [
    catalog,
    loadedPolicy,
    policyID,
    homopolymer,
    repeatLength,
    invertedLength,
    gcMin,
    gcMax,
    entropyWindow,
    entropyMin,
  ]);
  useEffect(() => {
    if (!project) return;
    setAnalysis(null);
    let cancelled = false;
    const timer = setTimeout(
      () =>
        api<Analysis>("/analyze", {
          project,
          host,
          codon_profile: customProfile,
          rare_threshold: rareThreshold ? Number(rareThreshold) : null,
          gc_window: gcWindow,
          synthesis_profile: policy,
          motif_search: motifQuery || null,
        })
          .then((row) => {
            if (!cancelled) setAnalysis(row);
          })
          .catch((e) => {
            if (!cancelled) setError(e.message);
          }),
      120,
    );
    return () => {
      cancelled = true;
      clearTimeout(timer);
    };
  }, [
    project,
    host,
    customProfile,
    rareThreshold,
    gcWindow,
    policy,
    motifQuery,
  ]);
  useEffect(() => {
    if (selection)
      setFeature((row) => ({
        ...row,
        start_1based: selection[0],
        end_1based: selection[1],
      }));
  }, [selection]);
  useEffect(() => {
    const end = () => {
      dragging.current = false;
    };
    window.addEventListener("pointerup", end);
    return () => window.removeEventListener("pointerup", end);
  }, []);
  useEffect(() => {
    if (!job || !["running", "queued"].includes(job.status)) return;
    let cancelled = false;
    const timer = setTimeout(
      () =>
        api<Job>("/jobs/" + job.id)
          .then((next) => {
            if (cancelled) return;
            setJob(next);
            if (next.status === "completed" && next.result) {
              setPayload(next.result);
              setProject(next.result.project);
              setNotice(
                "Local run completed. Final validation is shown in Reports.",
              );
            }
            if (next.status === "failed")
              setError(next.error ?? "The local run failed.");
            if (next.status === "cancelled")
              setNotice(
                "Run cancelled. Your previous construct was preserved.",
              );
          })
          .catch((e) => {
            if (!cancelled) {
              setError(e.message);
              setJob(null);
            }
          }),
      300,
    );
    return () => {
      cancelled = true;
      clearTimeout(timer);
    };
  }, [job]);

  const adopt = (next: Project) => {
    setProject(next);
    setAnalysis(null);
    setPayload(null);
    setJob(null);
    setSelection(null);
    setViewerStart(0);
    setHomopolymer("");
    setRepeatLength("");
    setInvertedLength("");
    setGCMin("");
    setGCMax("");
    setEntropyWindow("");
    setEntropyMin("");
    setGCWindow(50);
    setFeature(emptyFeature);
    setMotif("");
    setEnzyme("");
    setMotifQuery("");
    setSelectedCDS(
      next.original.features.find((f) => f.feature_type === "CDS")?.id ?? "",
    );
    setSelectedRBS(
      next.original.features.find((f) => f.feature_type === "RBS")?.id ?? "",
    );
    const saved = latestSettings(next, "search");
    const rbsSaved = latestSettings(next, "rbs_context");
    const context = object(rbsSaved.rbs_context);
    const rbsHost = object(context.host);
    setHost(
      typeof saved.host === "string"
        ? saved.host
        : typeof rbsHost.id === "string"
          ? rbsHost.id
          : (next.original.host_profile?.id ?? "ecoli-mg1655"),
    );
    setCustomProfile(
      saved.codon_profile
        ? (saved.codon_profile as Record<string, Json>)
        : null,
    );
    setLoadedPolicy(next.synthesis_policy);
    if (next.synthesis_policy) {
      setPolicyID(next.synthesis_policy.id);
    } else setPolicyID("generic");
    const number = (value: Json | undefined, fallback: number) =>
      typeof value === "number" ? value : fallback;
    const config = object(saved.search);
    const weights = object(config.weights);
    setSeed(number(config.seed, 7));
    setIterations(number(config.iterations, 500));
    setAlternatives(number(config.alternatives, 1));
    setNTerminal(number(config.n_terminal_lock_nt, 30));
    setAlgorithm(
      typeof config.algorithm === "string"
        ? config.algorithm
        : "simulated_annealing",
    );
    setCodonWeight(number(weights.codon_adaptation, 1));
    setGCWeight(number(weights.global_gc, 0));
    setRareWeight(number(weights.rare_codons, 0));
    setPreset(Object.keys(weights).length ? "Custom" : "Balanced");
    setGCTarget(
      typeof config.gc_target === "number"
        ? String(config.gc_target * 100)
        : "",
    );
    setRareThreshold(
      typeof config.rare_threshold === "number"
        ? String(config.rare_threshold)
        : "",
    );
    setSynthesisWeight(number(saved.synthesis_weight, 1));
    const assembly = object(
      next.gibson?.settings ??
        latestSettings(next, "gibson_settings").gibson_settings,
    );
    setMaximum(number(assembly.maximum_ordered_length, 2500));
    setOverlap(number(assembly.overlap_length, 30));
    setMinimum(number(assembly.minimum_ordered_length, 60));
    setBoundaryWindow(number(assembly.search_window, 100));
    setAvoidProtected(Boolean(assembly.avoid_protected_overlaps));
    setPreferredTm(
      typeof assembly.preferred_tm_c === "number"
        ? String(assembly.preferred_tm_c)
        : "",
    );
    const range = object(assembly.overlap_gc_range);
    setOverlapGCMin(
      typeof range.minimum === "number" ? String(range.minimum * 100) : "",
    );
    setOverlapGCMax(
      typeof range.maximum === "number" ? String(range.maximum * 100) : "",
    );
    const solution = object(assembly.conditions);
    setStrandConcentration(number(solution.strand1_nm, 25));
    setStrand2Concentration(number(solution.strand2_nm, 25));
    setSodium(number(solution.sodium_mm, 50));
    setMagnesium(number(solution.magnesium_mm, 0));
    setPotassium(number(solution.potassium_mm, 0));
    setTris(number(solution.tris_mm, 0));
    setDNTPs(number(solution.dntps_mm, 0));
    setSaltCorrection(number(solution.salt_correction, 5));
    if (typeof context.cds_id === "string") setSelectedCDS(context.cds_id);
    if (typeof context.rbs_id === "string") setSelectedRBS(context.rbs_id);
    const upstream = object(context.upstream_region);
    // Display the length of the backend's stored interval; no sequence extraction here.
    setUpstreamNT(
      typeof upstream.end === "number" && typeof upstream.start === "number"
        ? upstream.end - upstream.start
        : 35,
    );
    setDownstreamNT(number(context.downstream_nt, 35));
    setAntiSD(typeof rbsHost.anti_sd === "string" ? rbsHost.anti_sd : "");
    setAntiSource(typeof rbsHost.source === "string" ? rbsHost.source : "");
    setUncalibrated(Boolean(context.allow_uncalibrated));
    const rbsConfig = object(latestSettings(next, "rbs_settings").rbs_settings);
    setRBSSeed(number(rbsConfig.seed, 7));
    setRBSIterations(number(rbsConfig.iterations, 500));
    setRBSAlternatives(number(rbsConfig.alternatives, 1));
    setJoint(number(rbsConfig.joint_cds_nt, 0));
    setProtectDownstream(number(rbsConfig.protect_downstream_nt, 30));
    setTarget(number(rbsSaved.target_relative_tir, 1000));
    setTargets(
      Array.isArray(rbsSaved.library_targets) && rbsSaved.library_targets.length
        ? rbsSaved.library_targets.join(",")
        : "250,500,1000",
    );
  };
  const importText = () =>
    guard(async () => {
      adopt(
        await api<Project>("/import", {
          format,
          text,
          name,
          topology,
          annotate_cds: wholeCDS,
          seed,
          host,
        }),
      );
      setPage("Sequence");
      setNotice("Reference loaded. Review or edit annotations before design.");
    });
  const fileInput = async (file: File | undefined) => {
    if (file) setText(await file.text());
  };
  const saveAnnotation = (remove = false) =>
    guard(async () => {
      if (!project) return;
      adopt(
        await api<Project>("/annotations", {
          project,
          ...feature,
          id: feature.id || "feature-" + crypto.randomUUID(),
          delete: remove,
          rbs_mask:
            feature.mutability === "DESIGNABLE_RBS" ? feature.rbs_mask : null,
        }),
      );
      setFeature(emptyFeature);
      setNotice(
        "Annotations updated. Previous designs were cleared; the original DNA remains the reference.",
      );
    });
  const editFeature = (row: Feature) => {
    setFeature({
      id: row.id,
      name: row.name,
      feature_type: row.feature_type,
      start_1based: row.start + 1,
      end_1based: row.end,
      strand: row.strand,
      mutability: row.mutability,
      translation_table: row.translation_table,
      cds_mode: row.cds_mode,
      frame_offset: row.frame_offset,
      rbs_mask: row.rbs_mask ?? "",
    });
    setSelection([row.start + 1, row.end]);
  };
  const setPresetChoice = (value: string) => {
    setPreset(value);
    if (value === "Expression priority") {
      setCodonWeight(2);
      setSynthesisWeight(0.25);
      setGCWeight(0);
      setRareWeight(0);
    } else if (value === "Synthesis priority") {
      setCodonWeight(0.25);
      setSynthesisWeight(2);
      setGCWeight(0);
      setRareWeight(0);
    } else if (value === "Balanced") {
      setCodonWeight(1);
      setSynthesisWeight(1);
      setGCWeight(0);
      setRareWeight(0);
    }
  };
  const search = () => ({
    ...object(project ? latestSettings(project, "search").search : undefined),
    seed,
    iterations,
    alternatives,
    algorithm,
    n_terminal_lock_nt: nTerminal,
    weights: {
      codon_adaptation: codonWeight,
      global_gc: gcWeight,
      rare_codons: rareWeight,
      custom: 0,
    },
    gc_target: gcTarget ? Number(gcTarget) / 100 : null,
    rare_threshold: rareThreshold ? Number(rareThreshold) : null,
  });
  const changeHost = (value: string) => {
    setHost(value);
    setCustomProfile(null);
    setAntiSD("");
    setAntiSource("");
    setUncalibrated(false);
    setNotice(
      "Host changed. Codon scoring uses the selected built-in profile; previous anti-SD settings were cleared.",
    );
  };
  const run = (operation: string) =>
    guard(async () => {
      if (!project) return;
      if (
        Boolean(gcMin) !== Boolean(gcMax) ||
        Boolean(entropyWindow) !== Boolean(entropyMin)
      )
        throw new Error(
          "Complete both bounds/window and threshold for the selected synthesis criterion.",
        );
      let request: Record<string, unknown> = { operation, project, host };
      if (operation === "codon" || operation === "synthesis")
        request = {
          ...request,
          search: search(),
          synthesis_profile: policy,
          codon_profile: customProfile,
          synthesis_weight: synthesisWeight,
        };
      if (operation.startsWith("rbs")) {
        const context = await api<Json>("/rbs-context", {
          project,
          cds_id: selectedCDS,
          rbs_id: selectedRBS,
          upstream_nt: upstreamNT,
          downstream_nt: downstreamNT,
          host,
          anti_sd: antiSD || null,
          source: antiSource || null,
          allow_uncalibrated: uncalibrated,
        });
        request = {
          ...request,
          rbs_context: context,
          target_relative_tir: target,
          library_targets: targets.split(",").map(Number),
          rbs_settings: {
            ...object(latestSettings(project, "rbs_settings").rbs_settings),
            seed: rbsSeed,
            iterations: rbsIterations,
            alternatives: rbsAlternatives,
            joint_cds_nt: joint,
            protect_downstream_nt: protectDownstream,
          },
        };
      }
      if (operation === "gibson")
        if (Boolean(overlapGCMin) !== Boolean(overlapGCMax))
          throw new Error("Provide both overlap GC bounds.");
      if (operation === "gibson")
        request = {
          ...request,
          gibson_settings: {
            ...object(
              project.gibson?.settings ??
                latestSettings(project, "gibson_settings").gibson_settings,
            ),
            maximum_ordered_length: maximum,
            overlap_length: overlap,
            minimum_ordered_length: minimum,
            search_window: boundaryWindow,
            overlap_gc_range:
              overlapGCMin && overlapGCMax
                ? {
                    minimum: Number(overlapGCMin) / 100,
                    maximum: Number(overlapGCMax) / 100,
                  }
                : null,
            preferred_tm_c: preferredTm ? Number(preferredTm) : null,
            avoid_protected_overlaps: avoidProtected,
            synthesis_profile: policy,
            conditions: {
              strand1_nm: strandConcentration,
              strand2_nm: strand2Concentration,
              sodium_mm: sodium,
              magnesium_mm: magnesium,
              potassium_mm: potassium,
              tris_mm: tris,
              dntps_mm: dntps,
              salt_correction: saltCorrection,
            },
          },
        };
      setJob(await api<Job>("/jobs", request));
      setNotice("");
    });
  const exportDesign = (kind: string) =>
    guard(async () => {
      if (project)
        await download(
          "/export",
          { project, format: kind },
          "dna-design." +
            { genbank: "gb", fasta: "fasta", csv: "csv", json: "json" }[kind],
        );
    });
  const changed = new Set(
    analysis?.changes.map((row) => row.position_1based) ?? [],
  );
  const matchBases = new Set(
    analysis?.motif_matches.flatMap((row) => row.positions.map((p) => p + 1)) ??
      [],
  );
  const permissionAt = (position: number) =>
    analysis?.permissions.find(
      (row) => row.start <= position && position < row.end,
    )?.state ?? "LOCKED";

  return (
    <div className="app">
      <aside className="sidebar">
        <div className="brand">
          <span className="brand-mark">
            AT
            <br />
            GC
          </span>
          <div>
            DNA Design<span>Local laboratory workbench</span>
          </div>
        </div>
        <nav aria-label="Workspace">
          {pages.map((item, index) => (
            <button
              key={item}
              className={page === item ? "active" : ""}
              onClick={() => setPage(item)}
            >
              <span>{String(index + 1).padStart(2, "0")}</span>
              {item}
            </button>
          ))}
        </nav>
        <div className="privacy">
          <span className="dot" /> Local sequence processing
          <p>
            No sequence uploads to prediction services. Save a project file to
            keep your work.
          </p>
        </div>
      </aside>
      <main>
        <header className="topbar">
          <div>
            <p className="eyebrow">CONSTRUCT WORKSPACE</p>
            <h1>{project?.original.name ?? "Start a DNA design"}</h1>
          </div>
          <div className="header-actions">
            {project && (
              <button
                disabled={busy}
                onClick={() =>
                  guard(() =>
                    download("/save-draft", project, "dna-project.json"),
                  )
                }
              >
                Save project
              </button>
            )}
            <span
              className={"status " + (analysis?.valid ? "pass" : "neutral")}
            >
              {project
                ? analysis?.valid
                  ? "Validated"
                  : "Review pending"
                : "No construct"}
            </span>
          </div>
        </header>
        {error && (
          <div role="alert" className="alert error">
            <strong>Action needs attention</strong>
            <p>{error}</p>
            <button onClick={() => setError("")}>Dismiss</button>
          </div>
        )}
        {notice && (
          <div role="status" className="alert notice">
            {notice}
          </div>
        )}
        {job && ["queued", "running"].includes(job.status) && (
          <div className="job-progress" role="status">
            <div>
              <strong>{job.stage}</strong>
              <span>{fmt(job.progress * 100, 0)}%</span>
            </div>
            <progress value={job.progress} max={1} />
            <button
              onClick={() =>
                guard(async () => {
                  await api("/jobs/" + job.id, undefined, "DELETE");
                })
              }
            >
              Cancel run
            </button>
          </div>
        )}
        <div className="page-heading">
          <div>
            <p className="eyebrow">
              {page === "Input" ? "IMPORT & REFERENCE" : "DNA DESIGN PLATFORM"}
            </p>
            <h2>
              {
                {
                  Input: "Bring your construct",
                  Sequence: "Annotate and protect",
                  Design: "Codon and synthesis design",
                  RBS: "Translation initiation",
                  Assembly: "Gibson assembly",
                  Reports: "Validation and exports",
                  Reliability: "Reliability workspace",
                }[page]
              }
            </h2>
          </div>
          {project && (
            <span className="muted">
              {project.original.sequence.length.toLocaleString()} bp ·{" "}
              {project.original.topology} · {project.original.features.length}{" "}
              features
            </span>
          )}
        </div>

        {page === "Input" && (
          <div className="grid two">
            <section className="card">
              <h3>Sequence input</h3>
              <p className="muted">
                Import DNA, annotations, a protein, or a saved project. Only
                synthetic/public DNA is used by the demo.
              </p>
              <div className="form-grid">
                <Field label="Input format">
                  <select
                    value={format}
                    onChange={(e) => setFormat(e.target.value)}
                  >
                    {["raw", "fasta", "genbank", "protein", "project"].map(
                      (row) => (
                        <option key={row}>{row}</option>
                      ),
                    )}
                  </select>
                </Field>
                <Field label="Construct name">
                  <input
                    value={name}
                    onChange={(e) => setName(e.target.value)}
                  />
                </Field>
                <Field label="Topology">
                  <select
                    value={topology}
                    onChange={(e) => setTopology(e.target.value)}
                  >
                    <option>linear</option>
                    <option>circular</option>
                  </select>
                </Field>
                <Field label="Local file">
                  <input
                    type="file"
                    accept=".txt,.fasta,.fa,.gb,.gbk,.json"
                    onChange={(e) =>
                      guard(() => fileInput(e.target.files?.[0]))
                    }
                  />
                </Field>
              </div>
              <Field label="Sequence or file contents">
                <textarea
                  className="sequence-input"
                  value={text}
                  onChange={(e) => setText(e.target.value)}
                  placeholder="Paste DNA, FASTA, GenBank, protein or project JSON"
                  spellCheck={false}
                />
              </Field>
              {["raw", "fasta"].includes(format) && (
                <label className="checkbox">
                  <input
                    type="checkbox"
                    checked={wholeCDS}
                    onChange={(e) => setWholeCDS(e.target.checked)}
                  />
                  Explicitly annotate the entire sequence as a forward, complete
                  table-11 CDS
                </label>
              )}
              <div className="actions">
                <button
                  className="primary"
                  disabled={busy || !text}
                  onClick={importText}
                >
                  Load reference
                </button>
                <button
                  disabled={busy}
                  onClick={() =>
                    guard(async () => {
                      adopt(await api<Project>("/demo"));
                      setPage("Sequence");
                    })
                  }
                >
                  Load synthetic demo
                </button>
              </div>
            </section>
            <section className="card">
              <h3>Start with explicit biology</h3>
              <p>
                Unannotated DNA and functional non-coding features begin
                protected. CDS changes must preserve every original protein,
                including overlaps and reverse genes.
              </p>
              <p>
                Coordinates shown here are 1-based inclusive. DNA and
                annotations are processed by the local application.
              </p>
              {help("locked")}
              {help("validation")}
              <h3>Continuing a project</h3>
              <p>
                Choose <strong>project</strong> to restore the original
                reference, final DNA, seeds, history and designed initiation
                locks. Draft files can contain unresolved checks; final exports
                require fresh validation.
              </p>
            </section>
          </div>
        )}

        {!project && page !== "Input" && (
          <section className="card empty">
            <h3>Load a construct first</h3>
            <p>Use a local sequence file or explore the synthetic demo.</p>
            <button className="primary" onClick={() => setPage("Input")}>
              Go to input
            </button>
          </section>
        )}
        {project && page === "Sequence" && (
          <>
            <section className="card">
              <div className="section-heading">
                <h3>Feature map</h3>
                <span className="muted">
                  Click a feature to inspect it; drag bases or shift-click to
                  select an interval.
                </span>
              </div>
              <svg
                className="feature-map"
                viewBox={`0 0 1000 ${Math.max(70, project.original.features.length * 30 + 30)}`}
                role="img"
                aria-label="Annotated sequence map"
              >
                <line x1={20} y1={16} x2={980} y2={16} stroke="#cdd9d4" />
                {project.original.features.map((row, i) => (
                  <g
                    key={row.id}
                    onClick={() => editFeature(row)}
                    className="feature-bar"
                    role="button"
                    tabIndex={0}
                    aria-label={row.name}
                    onKeyDown={(e) => {
                      if (e.key === "Enter") editFeature(row);
                    }}
                  >
                    <rect
                      x={
                        20 +
                        (row.start / project.original.sequence.length) * 960
                      }
                      y={32 + i * 30}
                      width={Math.max(
                        4,
                        ((row.end - row.start) /
                          project.original.sequence.length) *
                          960,
                      )}
                      height={18}
                      rx={3}
                      fill={
                        row.feature_type === "CDS"
                          ? "#2f7569"
                          : row.mutability === "LOCKED"
                            ? "#8b9c94"
                            : "#cc9744"
                      }
                    />
                    <text
                      x={
                        24 +
                        (row.start / project.original.sequence.length) * 960
                      }
                      y={45 + i * 30}
                      fill="white"
                      fontSize={11}
                    >
                      {row.name}{" "}
                      {row.strand === -1 ? "←" : row.strand === 1 ? "→" : ""}
                    </text>
                  </g>
                ))}
              </svg>
              <div className="legend">
                <span>
                  <i className="cds" />
                  CDS
                </span>
                <span>
                  <i className="locked" />
                  Protected
                </span>
                <span>
                  <i className="mutable" />
                  Designable
                </span>
                <span>
                  <i className="mutation" />
                  Changed base
                </span>
              </div>
            </section>
            <section className="card">
              <div className="section-heading">
                <h3>Sequence selection</h3>
                <div className="actions">
                  <button
                    disabled={viewerStart === 0}
                    onClick={() =>
                      setViewerStart(Math.max(0, viewerStart - 700))
                    }
                  >
                    Previous
                  </button>
                  <span>
                    {viewerStart + 1}–
                    {Math.min(viewerStart + 700, project.sequence.length)}
                  </span>
                  <button
                    disabled={viewerStart + 700 >= project.sequence.length}
                    onClick={() => setViewerStart(viewerStart + 700)}
                  >
                    Next
                  </button>
                </div>
              </div>
              <div className="bases" aria-label="DNA sequence">
                {Array.from(
                  {
                    length: Math.ceil(
                      project.sequence.slice(viewerStart, viewerStart + 700)
                        .length / 70,
                    ),
                  },
                  (_, line) => {
                    const start = viewerStart + line * 70;
                    return (
                      <div className="base-line" key={start}>
                        <span className="line-number">{start + 1}</span>
                        <div>
                          {project.sequence
                            .slice(start, start + 70)
                            .split("")
                            .map((base, j) => {
                              const p = start + j;
                              const selected =
                                selection &&
                                selection[0] <= p + 1 &&
                                p + 1 <= selection[1];
                              return (
                                <button
                                  key={p}
                                  title={`${p + 1}: ${permissionAt(p)}`}
                                  aria-label={`Base ${p + 1} ${base}`}
                                  className={[
                                    selected ? "selected" : "",
                                    changed.has(p + 1) ? "changed" : "",
                                    matchBases.has(p + 1) ? "matched" : "",
                                    permissionAt(p) === "LOCKED"
                                      ? "protected"
                                      : "",
                                  ].join(" ")}
                                  onPointerDown={(e) => {
                                    e.preventDefault();
                                    anchor.current =
                                      e.shiftKey && selection
                                        ? selection[0]
                                        : p + 1;
                                    dragging.current = true;
                                    setSelection([
                                      Math.min(anchor.current, p + 1),
                                      Math.max(anchor.current, p + 1),
                                    ]);
                                  }}
                                  onPointerEnter={() => {
                                    if (dragging.current)
                                      setSelection([
                                        Math.min(anchor.current, p + 1),
                                        Math.max(anchor.current, p + 1),
                                      ]);
                                  }}
                                  onKeyDown={(e) => {
                                    if (e.key === "Enter")
                                      setSelection([p + 1, p + 1]);
                                  }}
                                >
                                  {base}
                                </button>
                              );
                            })}
                        </div>
                      </div>
                    );
                  },
                )}
              </div>
              <p className="muted">
                {selection
                  ? `Selected ${selection[0]}–${selection[1]} (inclusive).`
                  : "Select an interval to prefill annotation coordinates."}
              </p>
              <details>
                <summary>Full DNA for copying</summary>
                <textarea
                  readOnly
                  value={project.sequence}
                  className="sequence-input"
                  aria-label="Full designed DNA"
                />
              </details>
              <div className="form-grid">
                <Field label="IUPAC motif search">
                  <input
                    value={motifQuery}
                    onChange={(e) =>
                      setMotifQuery(e.target.value.toUpperCase())
                    }
                    placeholder="GAATTC or GCNNGC"
                  />
                </Field>
              </div>
              {analysis?.motif_matches.map((row, i) => (
                <button
                  className="match-button"
                  key={i}
                  onClick={() => {
                    setSelection([row.start + 1, row.start + row.length]);
                    setViewerStart(Math.floor(row.start / 700) * 700);
                  }}
                >
                  Match {row.start + 1} · strand {row.strand}
                </button>
              ))}
            </section>
            <div className="grid two">
              <section className="card">
                <h3>{feature.id ? "Edit annotation" : "Create annotation"}</h3>
                <p className="muted">
                  Editing reference annotations clears earlier designs.
                  Overlapping permissions compose.
                </p>
                <div className="form-grid">
                  <Field label="Feature name">
                    <input
                      value={feature.name}
                      onChange={(e) =>
                        setFeature({ ...feature, name: e.target.value })
                      }
                    />
                  </Field>
                  <Field label="Feature type">
                    <select
                      value={feature.feature_type}
                      onChange={(e) => {
                        const type = e.target.value;
                        setFeature({
                          ...feature,
                          feature_type: type,
                          strand: type === "CDS" || type === "RBS" ? 1 : 0,
                          mutability:
                            type === "CDS"
                              ? "SYNONYMOUS_ONLY"
                              : type === "RBS"
                                ? "DESIGNABLE_RBS"
                                : "LOCKED",
                        });
                      }}
                    >
                      {[
                        "CDS",
                        "RBS",
                        "promoter",
                        "terminator",
                        "protected region",
                        "mutable non-coding region",
                        "custom feature",
                      ].map((row) => (
                        <option key={row}>{row}</option>
                      ))}
                    </select>
                  </Field>
                  <Field label="Start (1-based)">
                    <input
                      type="number"
                      min={1}
                      value={feature.start_1based}
                      onChange={(e) =>
                        setFeature({
                          ...feature,
                          start_1based: Number(e.target.value),
                        })
                      }
                    />
                  </Field>
                  <Field label="End (inclusive)">
                    <input
                      type="number"
                      min={1}
                      value={feature.end_1based}
                      onChange={(e) =>
                        setFeature({
                          ...feature,
                          end_1based: Number(e.target.value),
                        })
                      }
                    />
                  </Field>
                  <Field label="Strand">
                    <select
                      value={feature.strand}
                      onChange={(e) =>
                        setFeature({
                          ...feature,
                          strand: Number(e.target.value),
                        })
                      }
                    >
                      <option value={0}>Unspecified</option>
                      <option value={1}>Forward (+)</option>
                      <option value={-1}>Reverse (−)</option>
                    </select>
                  </Field>
                  <Field label="Mutability">
                    <select
                      value={feature.mutability}
                      onChange={(e) =>
                        setFeature({ ...feature, mutability: e.target.value })
                      }
                    >
                      {[
                        "LOCKED",
                        "SYNONYMOUS_ONLY",
                        "FREE_MUTATION",
                        "DESIGNABLE_RBS",
                      ].map((row) => (
                        <option key={row}>{row}</option>
                      ))}
                    </select>
                  </Field>
                  {feature.feature_type === "CDS" && (
                    <>
                      <Field label="Genetic code">
                        <select
                          value={feature.translation_table}
                          onChange={(e) =>
                            setFeature({
                              ...feature,
                              translation_table: Number(e.target.value),
                            })
                          }
                        >
                          {catalog?.genetic_codes.map((row) => (
                            <option value={row.id} key={row.id}>
                              {row.id} · {row.name}
                            </option>
                          ))}
                        </select>
                      </Field>
                      <Field label="CDS mode">
                        <select
                          value={feature.cds_mode}
                          onChange={(e) =>
                            setFeature({ ...feature, cds_mode: e.target.value })
                          }
                        >
                          <option>complete</option>
                          <option>fragment</option>
                        </select>
                      </Field>
                      {feature.cds_mode === "fragment" && (
                        <Field label="Leading frame offset">
                          <input
                            type="number"
                            min={0}
                            max={2}
                            value={feature.frame_offset}
                            onChange={(e) =>
                              setFeature({
                                ...feature,
                                frame_offset: Number(e.target.value),
                              })
                            }
                          />
                        </Field>
                      )}
                    </>
                  )}
                  {feature.mutability === "DESIGNABLE_RBS" && (
                    <Field label="Oriented RBS mask (N/ACGT)">
                      <input
                        value={feature.rbs_mask}
                        onChange={(e) =>
                          setFeature({
                            ...feature,
                            rbs_mask: e.target.value.toUpperCase(),
                          })
                        }
                        placeholder="NNAGGAGGNNN"
                      />
                    </Field>
                  )}
                </div>
                <div className="actions">
                  <button
                    className="primary"
                    disabled={busy || !feature.name}
                    onClick={() => saveAnnotation()}
                  >
                    Save annotation
                  </button>
                  {feature.id && (
                    <button
                      disabled={busy}
                      onClick={() => saveAnnotation(true)}
                    >
                      Delete annotation
                    </button>
                  )}
                  <button onClick={() => setFeature(emptyFeature)}>
                    New annotation
                  </button>
                </div>
                {help("locked")}
              </section>
              <section className="card">
                <h3>Hard motif constraints</h3>
                <div className="form-grid">
                  <Field label="Mode">
                    <select
                      value={motifMode}
                      onChange={(e) => setMotifMode(e.target.value)}
                    >
                      {["AVOID", "REQUIRED", "PRESERVE"].map((row) => (
                        <option key={row}>{row}</option>
                      ))}
                    </select>
                  </Field>
                  <Field label="Restriction enzyme">
                    <input
                      value={enzyme}
                      onChange={(e) => setEnzyme(e.target.value)}
                      placeholder="EcoRI"
                    />
                  </Field>
                  <Field label="IUPAC motif">
                    <input
                      value={motif}
                      onChange={(e) => setMotif(e.target.value.toUpperCase())}
                      placeholder="GAATTC"
                    />
                  </Field>
                </div>
                <p className="muted">
                  Constraints search both strands. PRESERVE binds existing
                  occurrences and their original bases.
                </p>
                <button
                  disabled={busy || (!motif && !enzyme)}
                  onClick={() =>
                    guard(async () => {
                      adopt(
                        await api<Project>("/motifs", {
                          project,
                          id: "motif-" + crypto.randomUUID(),
                          mode: motifMode,
                          motif: motif || null,
                          enzyme: enzyme || null,
                        }),
                      );
                    })
                  }
                >
                  Add hard constraint
                </button>
                {project.original.constraints.map((row) => (
                  <div className="constraint" key={row.id}>
                    <span>
                      <strong>{row.mode}</strong> {row.enzyme ?? row.motif}
                    </span>
                    <button
                      disabled={busy}
                      onClick={() =>
                        guard(async () => {
                          adopt(
                            await api<Project>("/motifs", {
                              project,
                              id: row.id,
                              delete: true,
                            }),
                          );
                        })
                      }
                    >
                      Remove
                    </button>
                  </div>
                ))}
              </section>
            </div>
          </>
        )}

        {project && page === "Design" && (
          <div className="grid two">
            <section className="card">
              <h3>Optimization settings</h3>
              <div className="form-grid">
                <Field label="Host codon reference">
                  <select
                    value={host}
                    onChange={(e) => changeHost(e.target.value)}
                  >
                    {catalog?.hosts.map((row) => (
                      <option value={row.id} key={row.id}>
                        {row.organism} · {row.strain}
                      </option>
                    ))}
                  </select>
                </Field>
                <Field label="Preference preset">
                  <select
                    value={preset}
                    onChange={(e) => setPresetChoice(e.target.value)}
                  >
                    {[
                      "Balanced",
                      "Expression priority",
                      "Synthesis priority",
                      "Custom",
                    ].map((row) => (
                      <option key={row}>{row}</option>
                    ))}
                  </select>
                </Field>
                <Field label="Custom codon reference JSON">
                  <input
                    type="file"
                    accept=".json"
                    onChange={(e) =>
                      guard(async () => {
                        if (e.target.files?.[0])
                          setCustomProfile(
                            JSON.parse(await e.target.files[0].text()),
                          );
                      })
                    }
                  />
                </Field>
                {customProfile && (
                  <p className="note">
                    Custom reference: {String(customProfile.id)} ·{" "}
                    {String(customProfile.version)}
                    <button onClick={() => setCustomProfile(null)}>
                      Use built-in reference
                    </button>
                  </p>
                )}
                <Field label="Random seed">
                  <input
                    type="number"
                    value={seed}
                    onChange={(e) => setSeed(Number(e.target.value))}
                  />
                </Field>
                <Field label="Search proposals">
                  <input
                    type="number"
                    min={0}
                    max={20000}
                    value={iterations}
                    onChange={(e) => setIterations(Number(e.target.value))}
                  />
                </Field>
                <Field label="Alternative designs">
                  <input
                    type="number"
                    min={1}
                    max={32}
                    value={alternatives}
                    onChange={(e) => setAlternatives(Number(e.target.value))}
                  />
                </Field>
                <Field label="Synthesis profile">
                  <select
                    value={policyID}
                    onChange={(e) => {
                      setPolicyID(e.target.value);
                      setLoadedPolicy(null);
                    }}
                  >
                    {catalog?.synthesis_profiles.map((row) => (
                      <option value={row.id} key={row.id}>
                        {row.id}
                      </option>
                    ))}
                    {loadedPolicy &&
                      !catalog?.synthesis_profiles.some(
                        (row) => row.id === loadedPolicy.id,
                      ) && (
                        <option value={loadedPolicy.id}>
                          {loadedPolicy.id}
                        </option>
                      )}
                  </select>
                </Field>
              </div>
              <p className="note">{policy?.approximation}</p>
              <label className="checkbox">
                <input
                  type="checkbox"
                  checked={advanced}
                  onChange={(e) => setAdvanced(e.target.checked)}
                />
                Advanced settings
              </label>
              {advanced && (
                <>
                  <div className="form-grid">
                    <Field label="Algorithm">
                      <select
                        value={algorithm}
                        onChange={(e) => setAlgorithm(e.target.value)}
                      >
                        <option>simulated_annealing</option>
                        <option>hill_climb</option>
                      </select>
                    </Field>
                    <Field label="Protected first coding bases">
                      <input
                        type="number"
                        min={0}
                        value={nTerminal}
                        onChange={(e) => setNTerminal(Number(e.target.value))}
                      />
                    </Field>
                    <Field label="Codon adaptation weight">
                      <input
                        type="number"
                        min={0}
                        step="0.1"
                        value={codonWeight}
                        onChange={(e) => {
                          setCodonWeight(Number(e.target.value));
                          setPreset("Custom");
                        }}
                      />
                    </Field>
                    <Field label="GC weight">
                      <input
                        type="number"
                        min={0}
                        value={gcWeight}
                        onChange={(e) => {
                          setGCWeight(Number(e.target.value));
                          setPreset("Custom");
                        }}
                      />
                    </Field>
                    <Field label="Synthesis risk weight">
                      <input
                        type="number"
                        min={0}
                        step="0.1"
                        value={synthesisWeight}
                        onChange={(e) => {
                          setSynthesisWeight(Number(e.target.value));
                          setPreset("Custom");
                        }}
                      />
                    </Field>
                    <Field label="Explicit GC target (%)">
                      <input
                        type="number"
                        min={0}
                        max={100}
                        value={gcTarget}
                        onChange={(e) => setGCTarget(e.target.value)}
                      />
                    </Field>
                    <Field label="Rare-codon weight">
                      <input
                        type="number"
                        min={0}
                        value={rareWeight}
                        onChange={(e) => {
                          setRareWeight(Number(e.target.value));
                          setPreset("Custom");
                        }}
                      />
                    </Field>
                    <Field label="Explicit rare threshold (0–1)">
                      <input
                        type="number"
                        min={0}
                        max={1}
                        step="0.01"
                        value={rareThreshold}
                        onChange={(e) => setRareThreshold(e.target.value)}
                      />
                    </Field>
                    <Field label="Homopolymer maximum">
                      <input
                        type="number"
                        min={1}
                        value={homopolymer}
                        onChange={(e) => setHomopolymer(e.target.value)}
                        placeholder="Unassessed"
                      />
                    </Field>
                    <Field label="Direct repeat minimum (bp)">
                      <input
                        type="number"
                        min={2}
                        value={repeatLength}
                        onChange={(e) => setRepeatLength(e.target.value)}
                        placeholder="Unassessed"
                      />
                    </Field>
                    <Field label="Inverted repeat minimum (bp)">
                      <input
                        type="number"
                        min={2}
                        value={invertedLength}
                        onChange={(e) => setInvertedLength(e.target.value)}
                        placeholder="Unassessed"
                      />
                    </Field>
                    <Field label="Global GC minimum (%)">
                      <input
                        type="number"
                        value={gcMin}
                        onChange={(e) => setGCMin(e.target.value)}
                      />
                    </Field>
                    <Field label="Global GC maximum (%)">
                      <input
                        type="number"
                        value={gcMax}
                        onChange={(e) => setGCMax(e.target.value)}
                      />
                    </Field>
                    <Field label="Entropy window (bp)">
                      <input
                        type="number"
                        value={entropyWindow}
                        onChange={(e) => setEntropyWindow(e.target.value)}
                      />
                    </Field>
                    <Field label="Minimum entropy (0–2 bits)">
                      <input
                        type="number"
                        step="0.1"
                        value={entropyMin}
                        onChange={(e) => setEntropyMin(e.target.value)}
                      />
                    </Field>
                  </div>
                  <p className="muted">
                    Weights are transparent engineering preferences. Blank
                    thresholds remain unassessed.
                  </p>
                </>
              )}
              <div className="actions">
                <button
                  className="primary"
                  disabled={busy}
                  onClick={() => run("codon")}
                >
                  Optimize codons
                </button>
                <button disabled={busy} onClick={() => run("synthesis")}>
                  Repair synthesis risks
                </button>
              </div>
              {help("codon_adaptation")}
              {help("initiation")}
              {help("rare_threshold")}
              {help("synthesis")}
            </section>
            <section className="card">
              <h3>Design result</h3>
              <p className="muted">
                All designs are checked against the original construct. Limited
                search may retain risks or miss a solution.
              </p>
              {payload?.alternatives && (
                <div className="alternatives">
                  {payload.alternatives.map((_, i) => (
                    <button
                      key={i}
                      disabled={busy || !job}
                      onClick={() =>
                        guard(async () => {
                          if (job)
                            setProject(
                              await api<Project>(
                                "/jobs/" + job.id + "/select",
                                { index: i },
                              ),
                            );
                        })
                      }
                    >
                      Use alternative {i + 1}
                    </button>
                  ))}
                </div>
              )}
              {payload && (
                <>
                  <p>
                    {payload.operation} · {payload.alternatives?.length ?? 0}{" "}
                    alternative designs
                  </p>
                  <ScientificResult payload={payload} />
                </>
              )}
              {analysis && (
                <>
                  <div className="metric-strip">
                    <div>
                      <small>Global GC</small>
                      <strong>{fmt((analysis.gc ?? 0) * 100, 1)}%</strong>
                    </div>
                    <div>
                      <small>Changed bases</small>
                      <strong>{analysis.changes.length}</strong>
                    </div>
                  </div>
                  {analysis.codon_metrics.map((row) => (
                    <p key={row.feature_id}>
                      {row.feature_id} · adaptation {fmt(row.adaptation_score)}{" "}
                      · CAI {fmt(row.cai)}
                    </p>
                  ))}
                </>
              )}
            </section>
          </div>
        )}

        {project && page === "RBS" && (
          <div className="grid two">
            <section className="card">
              <h3>Oriented initiation context</h3>
              <p className="muted">
                The backend extracts upstream and coding context on the
                annotated strand. Fixed mask bases remain unchanged.
              </p>
              <div className="form-grid">
                <Field label="CDS feature">
                  <select
                    value={selectedCDS}
                    onChange={(e) => setSelectedCDS(e.target.value)}
                  >
                    <option value="">Select CDS</option>
                    {project.original.features
                      .filter((row) => row.feature_type === "CDS")
                      .map((row) => (
                        <option value={row.id} key={row.id}>
                          {row.name}
                        </option>
                      ))}
                  </select>
                </Field>
                <Field label="RBS feature">
                  <select
                    value={selectedRBS}
                    onChange={(e) => setSelectedRBS(e.target.value)}
                  >
                    <option value="">Select RBS</option>
                    {project.original.features
                      .filter((row) => row.feature_type === "RBS")
                      .map((row) => (
                        <option value={row.id} key={row.id}>
                          {row.name}
                        </option>
                      ))}
                  </select>
                </Field>
                <Field label="Host">
                  <select
                    value={host}
                    onChange={(e) => changeHost(e.target.value)}
                  >
                    {catalog?.hosts.map((row) => (
                      <option value={row.id} key={row.id}>
                        {row.organism} · {row.strain}
                      </option>
                    ))}
                  </select>
                </Field>
                <Field label="Target relative TIR">
                  <input
                    type="number"
                    min={0.001}
                    value={target}
                    onChange={(e) => setTarget(Number(e.target.value))}
                  />
                </Field>
                <Field label="Joint synonymous CDS bases">
                  <input
                    type="number"
                    min={0}
                    value={joint}
                    onChange={(e) => setJoint(Number(e.target.value))}
                  />
                </Field>
                <Field label="Protected downstream bases">
                  <input
                    type="number"
                    min={3}
                    value={protectDownstream}
                    onChange={(e) =>
                      setProtectDownstream(Number(e.target.value))
                    }
                  />
                </Field>
                <Field label="RBS alternatives">
                  <input
                    type="number"
                    min={1}
                    max={32}
                    value={rbsAlternatives}
                    onChange={(e) => setRBSAlternatives(Number(e.target.value))}
                  />
                </Field>
                <Field label="Library targets (comma-separated)">
                  <input
                    value={targets}
                    onChange={(e) => setTargets(e.target.value)}
                  />
                </Field>
                <Field label="Random seed">
                  <input
                    type="number"
                    value={rbsSeed}
                    onChange={(e) => setRBSSeed(Number(e.target.value))}
                  />
                </Field>
                <Field label="RBS search proposals">
                  <input
                    type="number"
                    min={0}
                    max={5000}
                    value={rbsIterations}
                    onChange={(e) => setRBSIterations(Number(e.target.value))}
                  />
                </Field>
                <Field label="Upstream context (bp)">
                  <input
                    type="number"
                    min={1}
                    value={upstreamNT}
                    onChange={(e) => setUpstreamNT(Number(e.target.value))}
                  />
                </Field>
                <Field label="Downstream context (bp)">
                  <input
                    type="number"
                    min={3}
                    value={downstreamNT}
                    onChange={(e) => setDownstreamNT(Number(e.target.value))}
                  />
                </Field>
              </div>
              {(host !== "ecoli-mg1655" || uncalibrated) && (
                <div className="uncalibrated">
                  <p>
                    Host transfer is uncalibrated. Supply a sourced nine-base
                    anti-SD; no experimental PAO1 accuracy is claimed.
                  </p>
                  <Field label="Anti-SD (DNA)">
                    <input
                      value={antiSD}
                      onChange={(e) => setAntiSD(e.target.value.toUpperCase())}
                    />
                  </Field>
                  <Field label="Parameter source">
                    <input
                      value={antiSource}
                      onChange={(e) => setAntiSource(e.target.value)}
                    />
                  </Field>
                  <label className="checkbox">
                    <input
                      type="checkbox"
                      checked={uncalibrated}
                      onChange={(e) => setUncalibrated(e.target.checked)}
                    />
                    Use this explicitly uncalibrated host model
                  </label>
                </div>
              )}
              <div className="actions">
                <button
                  disabled={busy || !selectedCDS || !selectedRBS}
                  onClick={() => run("rbs-predict")}
                >
                  Predict current TIR
                </button>
                <button
                  className="primary"
                  disabled={busy || !selectedCDS || !selectedRBS}
                  onClick={() => run("rbs")}
                >
                  Design RBS
                </button>
                <button
                  disabled={busy || !selectedCDS || !selectedRBS}
                  onClick={() => run("rbs-library")}
                >
                  Build target library
                </button>
              </div>
              {!catalog?.rbs_available && (
                <p className="note">
                  Local thermodynamic backend is unavailable. Install the
                  optional RBS package to use these actions.
                </p>
              )}
              {help("tir")}
              {help("initiation")}
            </section>
            <section className="card">
              <h3>Model predictions</h3>
              {payload?.operation.startsWith("rbs") && payload.alternatives && (
                <div className="alternatives">
                  {payload.alternatives.map((_, i) => (
                    <button
                      key={i}
                      disabled={busy || !job}
                      onClick={() =>
                        guard(async () => {
                          if (job)
                            setProject(
                              await api<Project>(
                                "/jobs/" + job.id + "/select",
                                { index: i },
                              ),
                            );
                        })
                      }
                    >
                      Use RBS alternative {i + 1}
                    </button>
                  ))}
                </div>
              )}
              <ScientificResult
                payload={payload?.operation.startsWith("rbs") ? payload : null}
              />
            </section>
          </div>
        )}

        {project && page === "Assembly" && (
          <>
            <div className="grid two">
              <section className="card">
                <h3>Fragment ordering limits</h3>
                <div className="form-grid">
                  <Field label="Maximum ordered length (bp)">
                    <input
                      type="number"
                      value={maximum}
                      onChange={(e) => setMaximum(Number(e.target.value))}
                    />
                  </Field>
                  <Field label="Overlap length (bp)">
                    <input
                      type="number"
                      value={overlap}
                      onChange={(e) => setOverlap(Number(e.target.value))}
                    />
                  </Field>
                  <Field label="Minimum ordered length (bp)">
                    <input
                      type="number"
                      value={minimum}
                      onChange={(e) => setMinimum(Number(e.target.value))}
                    />
                  </Field>
                  <Field label="Preferred overlap Tm (°C, optional)">
                    <input
                      type="number"
                      value={preferredTm}
                      onChange={(e) => setPreferredTm(e.target.value)}
                    />
                  </Field>
                </div>
                <label className="checkbox">
                  <input
                    type="checkbox"
                    checked={avoidProtected}
                    onChange={(e) => setAvoidProtected(e.target.checked)}
                  />
                  Require overlaps outside protected DNA
                </label>
                <details>
                  <summary>Tm solution assumptions</summary>
                  <p>
                    Defaults are Biopython solution assumptions, not a measured
                    Gibson reaction mix.
                  </p>
                  <div className="form-grid">
                    <Field label="Overlap GC minimum (%)">
                      <input
                        type="number"
                        value={overlapGCMin}
                        onChange={(e) => setOverlapGCMin(e.target.value)}
                      />
                    </Field>
                    <Field label="Overlap GC maximum (%)">
                      <input
                        type="number"
                        value={overlapGCMax}
                        onChange={(e) => setOverlapGCMax(e.target.value)}
                      />
                    </Field>
                    <Field label="Preferred boundary window (bp)">
                      <input
                        type="number"
                        min={1}
                        value={boundaryWindow}
                        onChange={(e) =>
                          setBoundaryWindow(Number(e.target.value))
                        }
                      />
                    </Field>
                    <Field label="Strand 1 concentration (nM)">
                      <input
                        type="number"
                        value={strandConcentration}
                        onChange={(e) =>
                          setStrandConcentration(Number(e.target.value))
                        }
                      />
                    </Field>
                    <Field label="Strand 2 concentration (nM)">
                      <input
                        type="number"
                        value={strand2Concentration}
                        onChange={(e) =>
                          setStrand2Concentration(Number(e.target.value))
                        }
                      />
                    </Field>
                    <Field label="Sodium (mM)">
                      <input
                        type="number"
                        value={sodium}
                        onChange={(e) => setSodium(Number(e.target.value))}
                      />
                    </Field>
                    <Field label="Magnesium (mM)">
                      <input
                        type="number"
                        value={magnesium}
                        onChange={(e) => setMagnesium(Number(e.target.value))}
                      />
                    </Field>
                    <Field label="Potassium (mM)">
                      <input
                        type="number"
                        min={0}
                        value={potassium}
                        onChange={(e) => setPotassium(Number(e.target.value))}
                      />
                    </Field>
                    <Field label="Tris (mM)">
                      <input
                        type="number"
                        min={0}
                        value={tris}
                        onChange={(e) => setTris(Number(e.target.value))}
                      />
                    </Field>
                    <Field label="dNTPs (mM)">
                      <input
                        type="number"
                        min={0}
                        value={dntps}
                        onChange={(e) => setDNTPs(Number(e.target.value))}
                      />
                    </Field>
                    <Field label="Salt correction method">
                      <select
                        value={saltCorrection}
                        onChange={(e) =>
                          setSaltCorrection(Number(e.target.value))
                        }
                      >
                        <option value={5}>5 · entropy correction</option>
                        <option value={7}>7 · mixed-ion correction</option>
                      </select>
                    </Field>
                  </div>
                </details>
                <button
                  className="primary"
                  disabled={busy}
                  onClick={() => run("gibson")}
                >
                  Plan Gibson fragments
                </button>
                {help("overlap")}
                {help("tm")}
              </section>
              <section className="card">
                <h3>Assembly verification</h3>
                <p>
                  Ordered lengths include all shared bases. Plans must
                  reconstruct the final construct exactly, including circular
                  closure, with unique full overlaps.
                </p>
                {project.gibson ? (
                  <>
                    <div className="metric-strip">
                      <div>
                        <small>Fragments</small>
                        <strong>{project.gibson.fragments.length}</strong>
                      </div>
                      <div>
                        <small>Geometric minimum</small>
                        <strong>
                          {project.gibson.minimum_geometric_fragments}
                        </strong>
                      </div>
                    </div>
                    <button
                      disabled={busy || !analysis?.valid}
                      onClick={() => exportDesign("csv")}
                    >
                      Download ordering CSV
                    </button>
                    {project.gibson.warnings.map((row) => (
                      <p className="note" key={row}>
                        {row}
                      </p>
                    ))}
                  </>
                ) : (
                  <p className="muted">No validated fragment plan yet.</p>
                )}
                {payload?.operation === "gibson" && !project.gibson && (
                  <ScientificResult payload={payload} />
                )}
              </section>
            </div>
            {project.gibson && (
              <section className="card">
                <h3>Ordered fragments</h3>
                <div className="table-wrap">
                  <table>
                    <thead>
                      <tr>
                        <th>Fragment</th>
                        <th>Coordinates</th>
                        <th>Unique bp</th>
                        <th>Ordered bp</th>
                        <th>GC</th>
                        <th>Sequence</th>
                      </tr>
                    </thead>
                    <tbody>
                      {project.gibson.fragments.map((row) => (
                        <tr key={row.number}>
                          <td>{row.number}</td>
                          <td>
                            {row.start + 1}–{row.end}
                          </td>
                          <td>{row.unique_length}</td>
                          <td>{row.ordered_length}</td>
                          <td>{fmt(row.gc * 100, 1)}%</td>
                          <td>
                            <details>
                              <summary>Show DNA</summary>
                              <pre>{row.sequence}</pre>
                            </details>
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
                <h3>Overlap metrics</h3>
                <table>
                  <thead>
                    <tr>
                      <th>Start</th>
                      <th>Sequence</th>
                      <th>GC</th>
                      <th>Tm</th>
                      <th>Unique</th>
                    </tr>
                  </thead>
                  <tbody>
                    {project.gibson.overlaps.map((row) => (
                      <tr key={row.start}>
                        <td>{row.start + 1}</td>
                        <td className="mono">{row.sequence}</td>
                        <td>{fmt(row.gc * 100, 1)}%</td>
                        <td>{fmt(row.tm_c, 1)} °C</td>
                        <td>{row.unique ? "Yes" : "No"}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </section>
            )}
          </>
        )}

        {project && page === "Reports" && analysis && (
          <>
            <section className="card">
              <div className="section-heading">
                <h3>Independent final validation</h3>
                <span
                  className={"status " + (analysis.valid ? "pass" : "fail")}
                >
                  {analysis.valid ? "PASS" : "ERROR"}
                </span>
              </div>
              <div className="actions">
                {["fasta", "genbank", "json", "csv"].map((kind) => (
                  <button
                    key={kind}
                    disabled={
                      busy ||
                      !analysis.valid ||
                      (kind === "csv" && !project.gibson)
                    }
                    onClick={() => exportDesign(kind)}
                  >
                    Export {kind.toUpperCase()}
                  </button>
                ))}
              </div>
              <p className="muted">
                Final exports revalidate from scratch. Project files preserve
                settings, source DNA and history; FASTA contains DNA only.
              </p>
              <table>
                <thead>
                  <tr>
                    <th>Check</th>
                    <th>Status</th>
                    <th>Finding</th>
                  </tr>
                </thead>
                <tbody>
                  {analysis.validation.checks.map((row, i) => (
                    <tr key={i}>
                      <td>{row.code}</td>
                      <td>
                        <span
                          className={
                            "status " +
                            (row.status === "PASS"
                              ? "pass"
                              : row.status === "ERROR"
                                ? "fail"
                                : "warn")
                          }
                        >
                          {row.status}
                        </span>
                      </td>
                      <td>{row.message}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
              {help("validation")}
            </section>
            <div className="grid two">
              <section className="card">
                <h3>Sequence metrics</h3>
                <div className="metric-strip">
                  <div>
                    <small>Global GC</small>
                    <strong>{fmt((analysis.gc ?? 0) * 100, 1)}%</strong>
                  </div>
                  <div>
                    <small>Changed nucleotides</small>
                    <strong>{analysis.changes.length}</strong>
                  </div>
                </div>
                <Field label="Local GC window (bp)">
                  <input
                    type="number"
                    min={1}
                    value={gcWindow}
                    onChange={(e) => setGCWindow(Number(e.target.value))}
                  />
                </Field>
                <svg
                  className="gc-chart"
                  viewBox="0 0 500 120"
                  aria-label="Local GC chart"
                >
                  <line x1={0} y1={119} x2={500} y2={119} stroke="#d5dfda" />
                  <polyline
                    fill="none"
                    stroke="#2f7569"
                    strokeWidth={2}
                    points={analysis.local_gc
                      .filter(
                        (_, i) =>
                          i %
                            Math.max(
                              1,
                              Math.floor(analysis.local_gc.length / 500),
                            ) ===
                          0,
                      )
                      .map(
                        ([p, value]) =>
                          `${(p / Math.max(1, project.sequence.length - gcWindow)) * 500},${115 - value * 110}`,
                      )
                      .join(" ")}
                  />
                </svg>
                {analysis.metric_warning && (
                  <p className="note">{analysis.metric_warning}</p>
                )}
                {analysis.codon_metrics.map((row) => (
                  <div className="constraint" key={row.feature_id}>
                    <strong>{row.feature_id}</strong>
                    <span>
                      Adaptation {fmt(row.adaptation_score)} · CAI{" "}
                      {fmt(row.cai)}
                      <br />
                      Rare fraction {fmt(row.rare_fraction)} · clusters{" "}
                      {row.rare_clusters.length}
                    </span>
                  </div>
                ))}
                {help("gc")}
                {help("codon_adaptation")}
              </section>
              <section className="card">
                <h3>Protein integrity</h3>
                {analysis.features
                  .filter((row) => row.feature_type === "CDS")
                  .map((row) => (
                    <details key={row.id} open>
                      <summary>
                        {row.name} · {row.start_1based}–{row.end_1based} ·
                        strand {row.strand}
                      </summary>
                      {row.translation_error ? (
                        <p className="note">{row.translation_error}</p>
                      ) : (
                        <>
                          <p>Original</p>
                          <pre>{row.original_protein}</pre>
                          <p>Designed</p>
                          <pre>{row.designed_protein}</pre>
                        </>
                      )}
                    </details>
                  ))}
              </section>
            </div>
            <section className="card">
              <h3>Synthesis findings</h3>
              <table>
                <thead>
                  <tr>
                    <th>Criterion</th>
                    <th>Status</th>
                    <th>Observed</th>
                    <th>Finding</th>
                  </tr>
                </thead>
                <tbody>
                  {analysis.synthesis?.findings.slice(0, 150).map((row, i) => (
                    <tr key={i}>
                      <td>{row.metric}</td>
                      <td>{row.status}</td>
                      <td>{fmt(row.observed)}</td>
                      <td>{row.message}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
              <p className="muted">
                Showing up to 150 findings. The full report remains available in
                the run result.
              </p>
              {help("synthesis")}
            </section>
            <section className="card">
              <h3>Mutation log</h3>
              <details>
                <summary>Changed codons and encoded amino acids</summary>
                <table>
                  <thead>
                    <tr>
                      <th>CDS</th>
                      <th>Codon</th>
                      <th>Positions</th>
                      <th>Original</th>
                      <th>Designed</th>
                      <th>Amino acid</th>
                      <th>Category</th>
                    </tr>
                  </thead>
                  <tbody>
                    {analysis.codon_changes.slice(0, 200).map((row, i) => (
                      <tr key={i}>
                        <td>{row.feature_id}</td>
                        <td>{row.codon_index_1based}</td>
                        <td>{row.positions_1based.join(", ")}</td>
                        <td>{row.original}</td>
                        <td>{row.designed}</td>
                        <td>{row.amino_acid}</td>
                        <td>{row.category}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </details>
              <table>
                <thead>
                  <tr>
                    <th>Position (1-based)</th>
                    <th>Original</th>
                    <th>Designed</th>
                  </tr>
                </thead>
                <tbody>
                  {analysis.changes.slice(0, 200).map((row) => (
                    <tr key={row.position_1based}>
                      <td>{row.position_1based}</td>
                      <td>{row.original}</td>
                      <td>{row.designed}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
              <p className="muted">
                Showing up to 200 changes; source and final DNA are retained in
                the project.
              </p>
            </section>
          </>
        )}

        {page === "Reliability" && (
          <div className="grid two">
            <section className="card">
              <h3>Experimental data evaluation</h3>
              <p>
                Load a local benchmark CSV and its provenance manifest.
                Evaluation requires explicit units and held-out versus
                calibration labels.
              </p>
              <Field label="Benchmark CSV file">
                <input
                  type="file"
                  accept=".csv"
                  onChange={(e) =>
                    guard(async () => {
                      if (e.target.files?.[0])
                        setBenchmarkText(await e.target.files[0].text());
                    })
                  }
                />
              </Field>
              <Field label="Provenance manifest JSON">
                <textarea
                  value={benchmarkManifest}
                  onChange={(e) => setBenchmarkManifest(e.target.value)}
                  placeholder="Versioned data source, checksum, license, measurement units and split policy"
                />
              </Field>
              <Field label="Local manifest file">
                <input
                  type="file"
                  accept=".json"
                  onChange={(e) =>
                    guard(async () => {
                      if (e.target.files?.[0])
                        setBenchmarkManifest(await e.target.files[0].text());
                    })
                  }
                />
              </Field>
              <details>
                <summary>CSV format and study preparation</summary>
                <p>
                  Required mapped columns: id, sequence, group, host, condition,
                  split, measured, predicted. Splits are calibration, heldout or
                  excluded; replicates and related sequences must stay in one
                  group. The manifest records source, license, checksum,
                  measurement units and model version. Predictions must be
                  generated with that declared model; this evaluator reads
                  supplied predictions.
                </p>
                <p>
                  See docs/benchmark_workflow.md for the manifest schema and
                  public source links. Ordinal scores support ranking; fold
                  errors require a defensible common ratio scale or separate
                  calibration.
                </p>
              </details>
              <button
                className="primary"
                disabled={busy || !benchmarkText || !benchmarkManifest}
                onClick={() =>
                  guard(async () => {
                    setBenchmarkResult({
                      csv: benchmarkText,
                      manifest: benchmarkManifest,
                      report: await api<BenchmarkReport>(
                        "/benchmarks/evaluate",
                        {
                          csv: benchmarkText,
                          manifest: JSON.parse(benchmarkManifest),
                        },
                      ),
                    });
                  })
                }
              >
                Evaluate held-out predictions
              </button>
              <p className="note">
                The reliability study and final data choices will be planned
                together after the tool is built. Synthetic engineering tests do
                not establish experimental prediction accuracy.
              </p>
            </section>
            <section className="card">
              <h3>Quantitative report</h3>
              {currentBenchmark ? (
                <BenchmarkResult report={currentBenchmark} />
              ) : (
                <p className="muted">
                  No evaluation report for the current inputs. Load data and a
                  manifest, then evaluate.
                </p>
              )}
              <h3>Validation questions</h3>
              <p>
                Protein/lock integrity, ranking against measured expression,
                target TIR error, synthesis-risk repair and exact assembly
                reconstruction answer different questions. Calibration records
                must stay separate from independent validation.
              </p>
            </section>
          </div>
        )}
        <footer>
          Local scientific computation · Original-reference validation ·
          Reproducible seeds
        </footer>
      </main>
    </div>
  );
}
