import { FORMATS, convert, parseGenBank } from "./genbank.mjs";

const fileInput = document.querySelector("#file-input");
const dropZone = document.querySelector("#drop-zone");
const summary = document.querySelector("#file-summary");
const formatFieldset = document.querySelector("#format-fieldset");
const formatGrid = document.querySelector("#format-grid");
const convertButton = document.querySelector("#convert-button");
const status = document.querySelector("#status");
const results = document.querySelector("#results");
const resultList = document.querySelector("#result-list");
let sourceFile = null;
let records = null;
let objectUrls = [];

for (const format of FORMATS) {
  const option = document.createElement("div");
  option.className = "format-option";
  option.innerHTML = `<input id="format-${format.id}" type="checkbox" value="${format.id}"><label for="format-${format.id}"><strong>${format.name}</strong><small>${format.description}</small><span class="extension">${format.extension}</span></label>`;
  formatGrid.append(option);
}

dropZone.addEventListener("click", () => fileInput.click());
fileInput.addEventListener("change", () => loadFile(fileInput.files[0]));
document.querySelector("#remove-file").addEventListener("click", reset);
formatGrid.addEventListener("change", updateButton);
convertButton.addEventListener("click", generate);

for (const event of ["dragenter", "dragover"]) dropZone.addEventListener(event, dragEvent => {
  dragEvent.preventDefault();
  dropZone.classList.add("dragging");
});
for (const event of ["dragleave", "drop"]) dropZone.addEventListener(event, dragEvent => {
  dragEvent.preventDefault();
  dropZone.classList.remove("dragging");
});
dropZone.addEventListener("drop", event => loadFile(event.dataTransfer.files[0]));

async function loadFile(file) {
  if (!file) return;
  clearResults();
  setStatus("Reading and validating your file…");
  try {
    const text = await file.text();
    records = parseGenBank(text);
    sourceFile = file;
    dropZone.hidden = true;
    summary.hidden = false;
    summary.style.display = "flex";
    document.querySelector("#file-name").textContent = file.name;
    const features = records.reduce((total, record) => total + record.features.length, 0);
    document.querySelector("#file-meta").textContent = `${formatBytes(file.size)} · ${records.length} record${records.length === 1 ? "" : "s"} · ${features.toLocaleString()} features`;
    formatFieldset.disabled = false;
    setStatus("File parsed locally and ready to convert.");
    updateButton();
  } catch (error) {
    reset();
    setStatus(error.message || "This file could not be parsed as GenBank.", true);
  }
}

function generate() {
  clearResults();
  const selected = [...formatGrid.querySelectorAll("input:checked")].map(input => input.value);
  try {
    for (const formatId of selected) {
      const format = FORMATS.find(item => item.id === formatId);
      const content = convert(records, formatId);
      const filename = `${baseName(sourceFile.name)}${format.extension}`;
      const url = URL.createObjectURL(new Blob([content], { type: "text/plain;charset=utf-8" }));
      objectUrls.push(url);
      const item = document.createElement("div");
      item.className = "result-item";
      item.innerHTML = `<div><strong></strong><small></small></div><a>Download</a>`;
      item.querySelector("strong").textContent = filename;
      item.querySelector("small").textContent = `${format.name} · ${formatBytes(new Blob([content]).size)}`;
      const link = item.querySelector("a");
      link.href = url;
      link.download = filename;
      resultList.append(item);
    }
    results.hidden = false;
    results.scrollIntoView({ behavior: "smooth", block: "start" });
    setStatus(`${selected.length} output file${selected.length === 1 ? " is" : "s are"} ready.`);
  } catch (error) {
    setStatus(error.message || "Conversion failed.", true);
  }
}

function reset() {
  sourceFile = null;
  records = null;
  fileInput.value = "";
  summary.hidden = true;
  summary.style.display = "none";
  dropZone.hidden = false;
  formatFieldset.disabled = true;
  formatGrid.querySelectorAll("input").forEach(input => { input.checked = false; });
  convertButton.disabled = true;
  clearResults();
  setStatus("");
}

function clearResults() {
  objectUrls.forEach(URL.revokeObjectURL);
  objectUrls = [];
  resultList.innerHTML = "";
  results.hidden = true;
}

function updateButton() { convertButton.disabled = !records || !formatGrid.querySelector("input:checked"); }
function setStatus(message, error = false) { status.textContent = message; status.classList.toggle("error", error); }
function baseName(filename) { return filename.replace(/\.(genbank|gbff|gbk?|txt)$/i, "") || "genbank"; }
function formatBytes(bytes) {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 ** 2) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / 1024 ** 2).toFixed(1)} MB`;
}
