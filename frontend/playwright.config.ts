import { defineConfig } from "@playwright/test";

export default defineConfig({
  testDir: "./tests",
  timeout: 30_000,
  workers: 1,
  use: { baseURL: "http://127.0.0.1:8765", trace: "retain-on-failure" },
  webServer: {
    command: "../.venv/bin/dna-design serve --port 8765",
    url: "http://127.0.0.1:8765/api/health",
    reuseExistingServer: !process.env.CI,
  },
});
