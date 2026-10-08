import { downloadText } from "./client";
import type { BenchmarkReport } from "./types";

const value = (number: number | null) =>
  number === null
    ? "Unassessed"
    : number.toLocaleString(undefined, { maximumFractionDigits: 3 });

export function BenchmarkResult({ report }: { report: BenchmarkReport }) {
  return (
    <div>
      <h4>
        {report.manifest.dataset_id} · {report.manifest.dataset_version}
      </h4>
      <p>Model: {report.manifest.model_version}</p>
      <p className="note">
        {report.manifest.evidence === "synthetic"
          ? "Synthetic arithmetic check · no experimental accuracy evidence"
          : "Experimental observations · review study provenance and model-training overlap"}
      </p>
      <div className="metric-strip">
        <div>
          <small>Held-out records</small>
          <strong>{report.heldout_records}</strong>
        </div>
        <div>
          <small>Calibration records</small>
          <strong>{report.calibration_records}</strong>
        </div>
        <div>
          <small>Excluded records</small>
          <strong>{report.excluded_records}</strong>
        </div>
      </div>
      <p>
        Model-training independence:{" "}
        {report.independence_status === "user_documented_disjoint"
          ? "Documented by the supplied manifest; verify its source"
          : "Unverified"}
        .
      </p>
      <ul>
        {report.warnings.map((warning) => (
          <li key={warning}>{warning}</li>
        ))}
      </ul>
      {report.groups.map((group) => (
        <div
          className="model-result"
          key={JSON.stringify([group.host, group.condition])}
        >
          <h4>
            {group.host} · {group.condition}
          </h4>
          <p>
            {group.metrics.records} held-out records;{" "}
            {group.metrics.positive_pairs} positive pairs;{" "}
            {group.metrics.nonpositive_pairs_excluded_from_logs} pairs excluded
            from log/fold metrics.
          </p>
          <table aria-label={`Metrics for ${group.host} · ${group.condition}`}>
            <tbody>
              <tr>
                <th scope="row">Spearman rank correlation</th>
                <td>{value(group.metrics.spearman)}</td>
              </tr>
              <tr>
                <th scope="row">Pearson linear correlation</th>
                <td>{value(group.metrics.pearson_linear)}</td>
              </tr>
              <tr>
                <th scope="row">Pearson log correlation</th>
                <td>{value(group.metrics.pearson_log)}</td>
              </tr>
              <tr>
                <th scope="row">Median symmetric fold error</th>
                <td>{value(group.metrics.median_fold_error)}</td>
              </tr>
              {Object.entries(group.metrics.fractions_within_fold ?? {}).map(
                ([fold, fraction]) => (
                  <tr key={fold}>
                    <th scope="row">Within {fold}-fold</th>
                    <td>{value(fraction * 100)}%</td>
                  </tr>
                ),
              )}
            </tbody>
          </table>
          {!group.metrics.fold_metrics_assessed && (
            <p className="muted">
              Fold errors are unassessed. They require ratio-scale measurements
              and comparable units or separate calibration.
            </p>
          )}
          <p>
            Calibration:{" "}
            {report.manifest.calibration === "none"
              ? "Not fitted"
              : `${group.calibration_records} records; log scale ${value(group.calibration_log_scale)}`}
            .
          </p>
        </div>
      ))}
      <button
        onClick={() =>
          downloadText(
            JSON.stringify(report, null, 2) + "\n",
            "benchmark-report.json",
            "application/json",
          )
        }
      >
        Download evaluation report
      </button>
      <details>
        <summary>Full report and provenance</summary>
        <pre>{JSON.stringify(report, null, 2)}</pre>
      </details>
    </div>
  );
}
