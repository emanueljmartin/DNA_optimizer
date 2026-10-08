export type Json =
  string | number | boolean | null | Json[] | { [key: string]: Json };
export type BenchmarkReport = {
  report_version: string;
  manifest: Record<string, Json> & {
    dataset_id: string;
    dataset_version: string;
    model_version: string;
    evidence: string;
    scale: string;
    calibration: string;
  };
  heldout_records: number;
  calibration_records: number;
  excluded_records: number;
  independence_status: string;
  warnings: string[];
  groups: {
    host: string;
    condition: string;
    calibration_records: number;
    calibration_log_scale: number;
    metrics: {
      records: number;
      positive_pairs: number;
      nonpositive_pairs_excluded_from_logs: number;
      spearman: number | null;
      pearson_linear: number | null;
      pearson_log: number | null;
      median_fold_error: number | null;
      fractions_within_fold: Record<string, number> | null;
      fold_metrics_assessed: boolean;
    };
  }[];
};
export type Feature = {
  id: string;
  name: string;
  feature_type: string;
  start: number;
  end: number;
  strand: number;
  mutability: string;
  translation_table: number;
  cds_mode: string;
  frame_offset: number;
  rbs_mask: string | null;
  qualifiers: Record<string, string[]>;
};
export type Construct = {
  name: string;
  sequence: string;
  topology: string;
  features: Feature[];
  constraints: {
    id: string;
    mode: string;
    motif: string;
    enzyme: string | null;
    strand: string;
  }[];
  host_profile: { id: string } | null;
};
export type Project = {
  schema_version: number;
  original: Construct;
  sequence: string;
  initiation_locks: Json[];
  gibson: GibsonPlan | null;
  synthesis_policy: Policy | null;
  require_synthesis_policy: boolean;
  settings: Record<string, Json>;
  history: Json[];
  model_versions: Record<string, string>;
};
export type Policy = Record<string, Json> & {
  id: string;
  approximation: string;
};
export type Finding = {
  status: string;
  code: string;
  message: string;
  positions: number[];
  feature_id: string | null;
};
export type Validation = { checks: Finding[] };
export type SynthesisFinding = {
  status: string;
  metric: string;
  observed: number | null;
  message: string;
  positions: number[];
  penalty: number;
};
export type Analysis = {
  codon_changes: {
    feature_id: string;
    codon_index_1based: number;
    positions_1based: number[];
    original: string;
    designed: string;
    amino_acid: string;
    category: string;
  }[];
  valid: boolean;
  validation: Validation;
  gc: number | null;
  local_gc: [number, number][];
  codon_metrics: {
    feature_id: string;
    adaptation_score: number;
    cai: number | null;
    rare_fraction: number | null;
    rare_clusters: number[][];
  }[];
  features: (Feature & {
    start_1based: number;
    end_1based: number;
    original_protein?: string;
    designed_protein?: string;
    translation_error?: string;
  })[];
  permissions: { start: number; end: number; state: string }[];
  synthesis: { findings: SynthesisFinding[] } | null;
  metric_warning: string | null;
  changes: { position_1based: number; original: string; designed: string }[];
  motif_matches: {
    start: number;
    length: number;
    strand: number;
    positions: number[];
  }[];
};
export type Help = {
  title: string;
  definition: string;
  effect: string;
  limitations: string;
  reference: string;
};
export type Catalog = {
  hosts: {
    id: string;
    version: string;
    organism: string;
    strain: string;
    reference_kind: string;
    source_url: string;
    limitations: string[];
  }[];
  synthesis_profiles: Policy[];
  genetic_codes: { id: number; name: string }[];
  rbs_available: boolean;
  help: Record<string, Help>;
};
export type GibsonPlan = {
  status: string;
  minimum_geometric_fragments: number;
  warnings: string[];
  validation: Validation;
  settings: Record<string, Json>;
  overlaps: {
    start: number;
    sequence: string;
    gc: number;
    tm_c: number;
    unique: boolean;
  }[];
  fragments: {
    number: number;
    start: number;
    end: number;
    sequence: string;
    unique_length: number;
    ordered_length: number;
    gc: number;
    left_overlap: string;
    right_overlap: string;
  }[];
};
export type JobPayload = {
  operation: string;
  project: Project;
  result: Json;
  alternatives?: string[];
};
export type Job = {
  id: string;
  status: string;
  progress: number;
  stage: string;
  error: string | null;
  result: JobPayload | null;
};
