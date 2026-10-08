export async function api<T>(
  path: string,
  payload?: unknown,
  method?: string,
): Promise<T> {
  const response = await fetch("/api" + path, {
    method: method ?? (payload === undefined ? "GET" : "POST"),
    headers: {
      "Content-Type": "application/json",
      "X-DNA-Design-Client": "local-ui",
    },
    body: payload === undefined ? undefined : JSON.stringify(payload),
  });
  if (!response.ok) {
    const error = await response.json();
    throw new Error(
      typeof error.detail === "string"
        ? error.detail
        : JSON.stringify(error.detail),
    );
  }
  return response.json() as Promise<T>;
}
export async function download(
  path: string,
  payload: unknown,
  filename: string,
): Promise<void> {
  const response = await fetch("/api" + path, {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      "X-DNA-Design-Client": "local-ui",
    },
    body: JSON.stringify(payload),
  });
  if (!response.ok) {
    const error = await response.json();
    throw new Error(
      typeof error.detail === "string"
        ? error.detail
        : JSON.stringify(error.detail),
    );
  }
  downloadBlob(await response.blob(), filename);
}

export function downloadText(
  text: string,
  filename: string,
  type: string,
): void {
  downloadBlob(new Blob([text], { type }), filename);
}

function downloadBlob(blob: Blob, filename: string): void {
  const link = document.createElement("a");
  const url = URL.createObjectURL(blob);
  link.href = url;
  link.download = filename;
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
