import { test, expect } from "@playwright/test";
import { createHash } from "node:crypto";
import { readFile } from "node:fs/promises";

const benchmarkCSV =
  "id,sequence,group,host,condition,split,measured,predicted\n" +
  "a,AAA,a,test,c,heldout,1,2\nb,AAC,b,test,c,heldout,2,4\nc,AAG,c,test,c,heldout,4,8\n";
const benchmarkManifest = (scale: "ratio" | "ordinal") => ({
  dataset_id: "synthetic-browser-arithmetic",
  dataset_version: "v1",
  source: "Synthetic browser fixture",
  license: "Project-generated synthetic observations",
  license_source: "frontend/tests/workbench.spec.ts",
  accessed_date: "2026-10-08",
  csv_sha256: createHash("sha256").update(benchmarkCSV).digest("hex"),
  evidence: "synthetic",
  measurement: "Synthetic value",
  measured_units: "test-unit",
  prediction: "Fictional values",
  predicted_units: "test-unit",
  model_version: "test-double-v1",
  split_policy: "Synthetic disjoint groups",
  group_definition: "One group per synthetic sequence",
  scale,
  comparable_units: scale === "ratio",
});

test("switching host clears the previous reference and sourced RBS parameters", async ({
  page,
}) => {
  const project = await (await page.request.get("/api/demo")).json();
  const catalog = await (await page.request.get("/api/catalog")).json();
  const context = await (
    await page.request.post("/api/rbs-context", {
      headers: { "X-DNA-Design-Client": "local-ui" },
      data: { project, cds_id: "cds", rbs_id: "rbs" },
    })
  ).json();
  project.settings = {
    host: "ecoli-mg1655",
    search: { seed: 7, iterations: 0 },
    codon_profile: catalog.hosts[0],
  };
  project.history = [
    { operation: "rbs-predict", settings: { rbs_context: context } },
  ];
  await page.goto("/");
  await page.getByLabel("Input format").selectOption("project");
  await page
    .getByLabel("Sequence or file contents")
    .fill(JSON.stringify(project));
  await page.getByRole("button", { name: "Load reference" }).click();
  await page.getByRole("button", { name: "03 Design" }).click();
  await expect(
    page.getByRole("button", { name: "Use built-in reference" }),
  ).toBeVisible();
  await page
    .getByLabel("Host codon reference")
    .selectOption("pseudomonas-pao1");
  await expect(
    page.getByRole("button", { name: "Use built-in reference" }),
  ).toHaveCount(0);
  const requested = page.waitForRequest(
    (request) =>
      request.url().endsWith("/api/jobs") && request.method() === "POST",
  );
  await page.getByRole("button", { name: "Optimize codons" }).click();
  expect((await requested).postDataJSON().codon_profile).toBeNull();
  await expect(
    page.getByText("Local run completed.", { exact: false }),
  ).toBeVisible();
  const downloaded = page.waitForEvent("download");
  await page.getByRole("button", { name: "Save project" }).click();
  const saved = JSON.parse(
    await readFile((await (await downloaded).path())!, "utf-8"),
  );
  expect(saved.project.settings.codon_profile.id).toBe("pseudomonas-pao1");
  await page.getByRole("button", { name: "04 RBS" }).click();
  await expect(page.getByLabel("Anti-SD (DNA)")).toHaveValue("");
  await expect(page.getByLabel("Parameter source")).toHaveValue("");
  await expect(
    page.getByRole("checkbox", {
      name: "Use this explicitly uncalibrated host model",
    }),
  ).not.toBeChecked();
});

for (const scale of ["ratio", "ordinal"] as const) {
  test(`local ${scale} evaluation displays evidence and downloads the exact report`, async ({
    page,
  }) => {
    await page.goto("/");
    await page.getByRole("button", { name: "07 Reliability" }).click();
    await page.getByLabel("Benchmark CSV file").setInputFiles({
      name: "synthetic.csv",
      mimeType: "text/csv",
      buffer: Buffer.from(benchmarkCSV),
    });
    await page
      .getByLabel("Provenance manifest JSON")
      .fill(JSON.stringify(benchmarkManifest(scale)));
    const evaluated = page.waitForResponse((response) =>
      response.url().endsWith("/api/benchmarks/evaluate"),
    );
    await page
      .getByRole("button", { name: "Evaluate held-out predictions" })
      .click();
    const report = await (await evaluated).json();
    await expect(
      page.getByText(
        "Synthetic arithmetic check · no experimental accuracy evidence",
      ),
    ).toBeVisible();
    await expect(
      page.getByText("Model-training independence: Unverified."),
    ).toBeVisible();
    const metrics = page.getByRole("table", { name: "Metrics for test · c" });
    await expect(
      metrics.getByRole("row", {
        name: "Spearman rank correlation 1",
        exact: true,
      }),
    ).toBeVisible();
    await expect(
      metrics.getByRole("row", {
        name: `Median symmetric fold error ${scale === "ratio" ? "2" : "Unassessed"}`,
        exact: true,
      }),
    ).toBeVisible();
    if (scale === "ratio")
      await expect(
        metrics.getByRole("row", { name: "Within 2-fold 100%", exact: true }),
      ).toBeVisible();
    const downloaded = page.waitForEvent("download");
    await page
      .getByRole("button", { name: "Download evaluation report" })
      .click();
    const file = await downloaded;
    expect(file.suggestedFilename()).toBe("benchmark-report.json");
    expect(JSON.parse(await readFile((await file.path())!, "utf-8"))).toEqual(
      report,
    );
    await page
      .getByLabel("Provenance manifest JSON")
      .fill(
        JSON.stringify({ ...benchmarkManifest(scale), dataset_version: "v2" }),
      );
    await expect(
      page.getByRole("button", { name: "Download evaluation report" }),
    ).toHaveCount(0);
    await expect(
      page.getByText("No evaluation report for the current inputs.", {
        exact: false,
      }),
    ).toBeVisible();
  });
}

test("an earlier evaluation response cannot label newly edited inputs", async ({
  page,
}) => {
  let release!: () => void;
  let started!: () => void;
  const pending = new Promise<void>((resolve) => {
    release = resolve;
  });
  const arrived = new Promise<void>((resolve) => {
    started = resolve;
  });
  await page.route("**/api/benchmarks/evaluate", async (route) => {
    const response = await route.fetch();
    started();
    await pending;
    await route.fulfill({ response });
  });
  await page.goto("/");
  await page.getByRole("button", { name: "07 Reliability" }).click();
  await page.getByLabel("Benchmark CSV file").setInputFiles({
    name: "synthetic.csv",
    mimeType: "text/csv",
    buffer: Buffer.from(benchmarkCSV),
  });
  await page
    .getByLabel("Provenance manifest JSON")
    .fill(JSON.stringify(benchmarkManifest("ratio")));
  await page
    .getByRole("button", { name: "Evaluate held-out predictions" })
    .click();
  await arrived;
  await page.getByLabel("Provenance manifest JSON").fill(
    JSON.stringify({
      ...benchmarkManifest("ratio"),
      dataset_version: "changed-during-evaluation",
    }),
  );
  release();
  await expect(
    page.getByRole("button", { name: "Evaluate held-out predictions" }),
  ).toBeEnabled();
  await expect(
    page.getByRole("button", { name: "Download evaluation report" }),
  ).toHaveCount(0);
  await expect(
    page.getByRole("table", { name: "Metrics for test · c" }),
  ).toHaveCount(0);
});

test("multi-stage projects restore earlier search and RBS controls", async ({
  page,
}) => {
  const project = await (await page.request.get("/api/demo")).json();
  const headers = { "X-DNA-Design-Client": "local-ui" };
  const context = await (
    await page.request.post("/api/rbs-context", {
      headers,
      data: {
        project,
        cds_id: "cds",
        rbs_id: "rbs",
        upstream_nt: 15,
        downstream_nt: 42,
      },
    })
  ).json();
  const search = {
    seed: 21,
    iterations: 4,
    alternatives: 1,
    n_terminal_lock_nt: 15,
    algorithm: "hill_climb",
    initial_temperature: 0.123,
    weights: {
      codon_adaptation: 0.25,
      global_gc: 0,
      rare_codons: 0,
      custom: 0,
    },
    gc_target: null,
    rare_threshold: null,
    allowed_positions: [58, 61],
    additional_locks: [
      {
        start: 0,
        end: 3,
        sequence: project.original.sequence.slice(0, 3),
        reason: "Synthetic replay lock",
      },
    ],
  };
  const rbsSettings = {
    seed: 39,
    iterations: 0,
    alternatives: 1,
    joint_cds_nt: 12,
    protect_downstream_nt: 45,
    additional_locks: [],
  };
  project.history = [
    {
      operation: "codon",
      settings: { host: "ecoli-mg1655", search, synthesis_weight: 2 },
    },
    {
      operation: "rbs",
      settings: {
        rbs_context: context,
        rbs_settings: rbsSettings,
        target_relative_tir: 777,
        library_targets: [111, 777],
      },
    },
  ];
  project.settings = {
    operation: "gibson",
    host: "ecoli-mg1655",
    gibson_settings: { maximum_ordered_length: 2500 },
  };
  await page.goto("/");
  await page.getByLabel("Input format").selectOption("project");
  await page
    .getByLabel("Sequence or file contents")
    .fill(JSON.stringify(project));
  await page.getByRole("button", { name: "Load reference" }).click();
  await page.getByRole("button", { name: "03 Design" }).click();
  await expect(page.getByLabel("Random seed")).toHaveValue("21");
  await page.getByLabel("Advanced settings").check();
  await expect(page.getByLabel("Synthesis risk weight")).toHaveValue("2");
  const codonRequest = page.waitForRequest(
    (request) =>
      request.url().endsWith("/api/jobs") && request.method() === "POST",
  );
  await page.getByRole("button", { name: "Optimize codons" }).click();
  expect((await codonRequest).postDataJSON().search).toEqual(search);
  await expect(
    page.getByText("Local run completed.", { exact: false }),
  ).toBeVisible();
  await page.getByRole("button", { name: "04 RBS" }).click();
  await expect(page.getByLabel("Random seed")).toHaveValue("39");
  await expect(page.getByLabel("Joint synonymous CDS bases")).toHaveValue("12");
  await expect(page.getByLabel("Upstream context (bp)")).toHaveValue("15");
  await expect(page.getByLabel("Downstream context (bp)")).toHaveValue("42");
  await expect(page.getByLabel("Target relative TIR")).toHaveValue("777");
  const rbsRequest = page.waitForRequest(
    (request) =>
      request.url().endsWith("/api/jobs") && request.method() === "POST",
  );
  await page.getByRole("button", { name: "Design RBS", exact: true }).click();
  const replay = (await rbsRequest).postDataJSON();
  expect(replay.rbs_context).toEqual(context);
  expect(replay.rbs_settings).toEqual(rbsSettings);
  await expect(
    page.getByText("Local run completed.", { exact: false }),
  ).toBeVisible();
  await expect(page.getByRole("alert")).toHaveCount(0);
});

test("loading another construct clears previous scientific controls", async ({
  page,
}) => {
  await page.goto("/");
  await page.getByRole("button", { name: "Load synthetic demo" }).click();
  await page.getByRole("button", { name: "03 Design" }).click();
  await page.getByLabel("Advanced settings").check();
  await page.getByLabel("Random seed").fill("99");
  await page.getByLabel("Explicit GC target (%)").fill("71");
  await page.getByLabel("Explicit rare threshold (0–1)").fill("0.13");
  await page.getByRole("button", { name: "04 RBS" }).click();
  await page.getByLabel("Joint synonymous CDS bases").fill("12");
  await page.getByLabel("Upstream context (bp)").fill("15");
  await page.getByRole("button", { name: "01 Input" }).click();
  await page.getByRole("button", { name: "Load synthetic demo" }).click();
  await page.getByRole("button", { name: "03 Design" }).click();
  await expect(page.getByLabel("Random seed")).toHaveValue("7");
  await expect(page.getByLabel("Explicit GC target (%)")).toHaveValue("");
  await expect(page.getByLabel("Explicit rare threshold (0–1)")).toHaveValue(
    "",
  );
  await page.getByRole("button", { name: "04 RBS" }).click();
  await expect(page.getByLabel("Joint synonymous CDS bases")).toHaveValue("0");
  await expect(page.getByLabel("Upstream context (bp)")).toHaveValue("35");
});

test("assembly replay retains all declared Tm solution assumptions", async ({
  page,
}) => {
  const project = await (await page.request.get("/api/demo")).json();
  const conditions = {
    strand1_nm: 17,
    strand2_nm: 29,
    sodium_mm: 40,
    potassium_mm: 5,
    tris_mm: 8,
    magnesium_mm: 2,
    dntps_mm: 0.2,
    salt_correction: 7,
  };
  const settings = {
    maximum_ordered_length: 2500,
    overlap_length: 30,
    minimum_ordered_length: 60,
    search_window: 100,
    maximum_fragments: 3,
    overlap_gc_range: null,
    preferred_tm_c: null,
    require_unique_overlaps: true,
    avoid_protected_overlaps: false,
    conditions,
    synthesis_profile: null,
  };
  project.settings = { operation: "gibson", gibson_settings: settings };
  await page.goto("/");
  await page.getByLabel("Input format").selectOption("project");
  await page
    .getByLabel("Sequence or file contents")
    .fill(JSON.stringify(project));
  await page.getByRole("button", { name: "Load reference" }).click();
  await page.getByRole("button", { name: "05 Assembly" }).click();
  const request = page.waitForRequest(
    (request) =>
      request.url().endsWith("/api/jobs") && request.method() === "POST",
  );
  await page.getByRole("button", { name: "Plan Gibson fragments" }).click();
  const replay = (await request).postDataJSON().gibson_settings;
  expect(replay.conditions).toEqual(conditions);
  expect(replay.maximum_fragments).toBe(3);
  await expect(
    page.getByRole("heading", { name: "Ordered fragments" }),
  ).toBeVisible();
});

test("saved projects restore search seed and selected custom synthesis policy", async ({
  page,
}) => {
  const project = await (await page.request.get("/api/demo")).json();
  const catalog = await (await page.request.get("/api/catalog")).json();
  project.settings = {
    host: "ecoli-mg1655",
    search: {
      seed: 21,
      iterations: 12,
      alternatives: 2,
      n_terminal_lock_nt: 15,
      algorithm: "hill_climb",
      weights: {
        codon_adaptation: 0.25,
        global_gc: 0,
        rare_codons: 0,
        custom: 0,
      },
      gc_target: null,
      rare_threshold: null,
    },
    synthesis_weight: 2,
  };
  project.synthesis_policy = {
    ...catalog.synthesis_profiles[0],
    id: "custom",
    version: "saved-v1",
    source: "Explicit synthetic engineering policy",
    homopolymer_limit: 3,
  };
  await page.goto("/");
  await page.getByLabel("Input format").selectOption("project");
  await page
    .getByLabel("Sequence or file contents")
    .fill(JSON.stringify(project));
  await page.getByRole("button", { name: "Load reference" }).click();
  await page.getByRole("button", { name: "03 Design" }).click();
  await expect(page.getByLabel("Random seed")).toHaveValue("21");
  await expect(page.getByLabel("Search proposals")).toHaveValue("12");
  await expect(page.getByLabel("Synthesis profile")).toHaveValue("custom");
  await page.getByLabel("Advanced settings").check();
  await expect(page.getByLabel("Protected first coding bases")).toHaveValue(
    "15",
  );
  await expect(page.getByLabel("Synthesis risk weight")).toHaveValue("2");
  await expect(page.getByRole("alert")).toHaveCount(0);
});

test("cancelling a real constrained search preserves the input DNA", async ({
  page,
}) => {
  await page.goto("/");
  await page
    .getByLabel("Sequence or file contents")
    .fill("ATG" + "GAA".repeat(1000) + "TAA");
  await page.getByRole("checkbox", { name: /Explicitly annotate/ }).check();
  await page.getByRole("button", { name: "Load reference" }).click();
  await page.getByRole("button", { name: "03 Design" }).click();
  await page.getByLabel("Search proposals").fill("20000");
  await page.getByRole("button", { name: "Optimize codons" }).click();
  await page.getByRole("button", { name: "Cancel run" }).click();
  await expect(
    page.getByText("Run cancelled. Your previous construct was preserved."),
  ).toBeVisible();
  await page.getByRole("button", { name: "06 Reports" }).click();
  await expect(page.getByText("Changed nucleotides")).toBeVisible();
  await expect(
    page.getByText("Changed nucleotides").locator("..").locator("strong"),
  ).toHaveText("0");
  await expect(page.getByRole("alert")).toHaveCount(0);
});

test("demo, coordinate selection and existing annotation editing stay local", async ({
  page,
}) => {
  const external: string[] = [];
  page.on("request", (request) => {
    if (!request.url().startsWith("http://127.0.0.1:8765/"))
      external.push(request.url());
  });
  await page.goto("/");
  await page.getByRole("button", { name: "Load synthetic demo" }).click();
  await expect(
    page.getByRole("heading", { name: "Annotate and protect" }),
  ).toBeVisible();
  await page
    .getByRole("button", { name: "Synthetic CDS", exact: true })
    .click();
  await expect(page.getByLabel("Start (1-based)")).toHaveValue("36");
  await expect(page.getByLabel("End (inclusive)")).toHaveValue("101");
  await page.getByLabel("Feature name").fill("Edited CDS");
  await page.getByRole("button", { name: "Save annotation" }).click();
  await expect(
    page.getByRole("button", { name: "Edited CDS", exact: true }),
  ).toBeVisible();
  await expect(page.getByRole("alert")).toHaveCount(0);
  await page.getByRole("button", { name: "Base 1 A", exact: true }).click();
  await page
    .getByRole("button", { name: "Base 5 A", exact: true })
    .click({ modifiers: ["Shift"] });
  await expect(page.getByLabel("Start (1-based)")).toHaveValue("1");
  await expect(page.getByLabel("End (inclusive)")).toHaveValue("5");
  expect(external).toEqual([]);
});

test("seeded codon job, validated export and draft download", async ({
  page,
}) => {
  await page.goto("/");
  await page.getByRole("button", { name: "Load synthetic demo" }).click();
  await page.getByRole("button", { name: "03 Design" }).click();
  await page.getByLabel("Search proposals").fill("20");
  await page.getByRole("button", { name: "Optimize codons" }).click();
  await expect(
    page.getByText("Local run completed.", { exact: false }),
  ).toBeVisible();
  await page.getByRole("button", { name: "06 Reports" }).click();
  await expect(
    page.getByRole("button", { name: "Export FASTA" }),
  ).toBeEnabled();
  const finalDownload = page.waitForEvent("download");
  await page.getByRole("button", { name: "Export FASTA" }).click();
  expect((await finalDownload).suggestedFilename()).toBe("dna-design.fasta");
  const draftDownload = page.waitForEvent("download");
  await page.getByRole("button", { name: "Save project" }).click();
  expect((await draftDownload).suggestedFilename()).toBe("dna-project.json");
});

test("single-fragment assembly and ordering CSV use fresh validation", async ({
  page,
}) => {
  await page.goto("/");
  await page.getByRole("button", { name: "Load synthetic demo" }).click();
  await page.getByRole("button", { name: "05 Assembly" }).click();
  await page.getByRole("button", { name: "Plan Gibson fragments" }).click();
  await expect(
    page.getByRole("heading", { name: "Ordered fragments" }),
  ).toBeVisible();
  const download = page.waitForEvent("download");
  await page.getByRole("button", { name: "Download ordering CSV" }).click();
  expect((await download).suggestedFilename()).toBe("dna-design.csv");
  await expect(page.getByRole("alert")).toHaveCount(0);
});

test("RBS prediction uses the real installed local backend", async ({
  page,
}) => {
  await page.goto("/");
  await page.getByRole("button", { name: "Load synthetic demo" }).click();
  await page.getByRole("button", { name: "04 RBS" }).click();
  await page.getByRole("button", { name: "Predict current TIR" }).click();
  await expect(
    page.getByText("Relative model TIR", { exact: true }),
  ).toBeVisible();
  await expect(page.getByRole("alert")).toHaveCount(0);
});

test("invalid imports are explained and mobile layout stays usable", async ({
  page,
}) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto("/");
  await page.getByLabel("Sequence or file contents").fill("ACGTXYZ");
  await page.getByRole("button", { name: "Load reference" }).click();
  await expect(page.getByRole("alert")).toBeVisible();
  await page.getByRole("button", { name: "Dismiss" }).click();
  await page.getByRole("button", { name: "Load synthetic demo" }).click();
  await expect(
    page.getByRole("heading", { name: "Annotate and protect" }),
  ).toBeVisible();
  expect(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= window.innerWidth,
    ),
  ).toBe(true);
});
