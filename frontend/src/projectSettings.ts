import type { Json, Project } from "./types";

export function object(value: Json | undefined): Record<string, Json> {
  return value && typeof value === "object" && !Array.isArray(value)
    ? value
    : {};
}

// A project's most recent operation may belong to another stage. History holds
// the full requests, including controls that the browser does not expose.
export function latestSettings(
  project: Project,
  field: string,
): Record<string, Json> {
  const requests = [
    project.settings,
    ...project.history
      .slice()
      .reverse()
      .map((entry) => object(object(entry).settings)),
  ];
  return (
    requests.find(
      (request) => Object.keys(object(request[field])).length > 0,
    ) ?? {}
  );
}
