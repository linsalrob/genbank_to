const COMPLEMENT = { A: "T", T: "A", G: "C", C: "G", U: "A", R: "Y", Y: "R", S: "S", W: "W", K: "M", M: "K", B: "V", V: "B", D: "H", H: "D", N: "N" };

export const FORMATS = [
  { id: "fna", name: "Genome sequence", extension: ".fna", description: "Complete nucleotide records in FASTA format" },
  { id: "orfs", name: "Coding sequences", extension: ".orfs.fna", description: "Extracted CDS nucleotide sequences" },
  { id: "faa", name: "Proteins", extension: ".faa", description: "Translated CDS amino acid sequences" },
  { id: "functions", name: "Protein functions", extension: ".functions.tsv", description: "Sequence, protein ID, and product table" },
  { id: "gff3", name: "Genome annotations", extension: ".gff3", description: "Features and sequence in GFF3 format" },
  { id: "ptt", name: "Protein table", extension: ".ptt", description: "Legacy NCBI-style protein table" },
];

export function parseGenBank(text) {
  const blocks = text.replace(/\r\n?/g, "\n").split(/^\/\/\s*$/m);
  const records = blocks.map(parseRecord).filter(Boolean);
  if (!records.length) throw new Error("No GenBank records were found in this file.");
  return records;
}

function parseRecord(block) {
  if (!/^LOCUS\s+/m.test(block)) return null;
  const locus = block.match(/^LOCUS\s+(\S+)/m)?.[1] || "unknown";
  const accession = block.match(/^ACCESSION\s+(\S+)/m)?.[1];
  const version = block.match(/^VERSION\s+(\S+)/m)?.[1];
  const definition = collectHeader(block, "DEFINITION");
  const id = version || accession || locus;
  const originIndex = block.search(/^ORIGIN\b/m);
  if (originIndex < 0) throw new Error(`Record ${id} has no ORIGIN sequence.`);
  const featureMatch = block.match(/^FEATURES\s+Location\/Qualifiers\s*$/m);
  const featureText = featureMatch ? block.slice(featureMatch.index + featureMatch[0].length, originIndex) : "";
  const sequence = block.slice(originIndex).replace(/^ORIGIN[^\n]*\n?/, "").replace(/[^A-Za-z]/g, "").toUpperCase();
  return { id, locus, accession, definition, sequence, features: parseFeatures(featureText) };
}

function collectHeader(block, name) {
  const match = block.match(new RegExp(`^${name}\\s+(.+(?:\\n {12}.+)*)`, "m"));
  return match ? match[1].replace(/\n\s+/g, " ").trim() : "";
}

function parseFeatures(text) {
  const features = [];
  let feature = null;
  let qualifier = null;
  for (const line of text.split("\n")) {
    const start = line.match(/^ {5}(\S+)\s+(.+)$/);
    if (start) {
      feature = { type: start[1], location: start[2].trim(), qualifiers: {} };
      features.push(feature);
      qualifier = null;
      continue;
    }
    if (!feature) continue;
    const q = line.match(/^ {21}\/([^=\s]+)(?:=(.*))?$/);
    if (q) {
      qualifier = q[1];
      const value = q[2] === undefined ? "" : q[2];
      feature.qualifiers[qualifier] ||= [];
      feature.qualifiers[qualifier].push(value);
      continue;
    }
    const continuation = line.match(/^ {21}(.+)$/);
    if (continuation && qualifier) {
      const values = feature.qualifiers[qualifier];
      const separator = qualifier === "translation" ? "" : " ";
      values[values.length - 1] += separator + continuation[1].trim();
    } else if (continuation && !line.trim().startsWith("/")) {
      feature.location += continuation[1].trim();
    }
  }
  for (const feature of features) {
    for (const [key, values] of Object.entries(feature.qualifiers)) {
      feature.qualifiers[key] = values.map(value => unquote(value));
    }
    feature.parts = parseLocation(feature.location);
  }
  return features;
}

function unquote(value) {
  if (value.startsWith('"')) value = value.slice(1);
  if (value.endsWith('"')) value = value.slice(0, -1);
  return value.replace(/""/g, '"');
}

export function parseLocation(location) {
  const cleaned = location.replace(/\s/g, "");
  const reverse = cleaned.startsWith("complement(");
  const inner = unwrap(cleaned, "complement") ?? cleaned;
  const joined = unwrap(inner, "join") ?? unwrap(inner, "order") ?? inner;
  const parts = splitTopLevel(joined).map(part => {
    const partReverse = part.startsWith("complement(");
    const range = unwrap(part, "complement") ?? part;
    const numbers = [...range.matchAll(/\d+/g)].map(match => Number(match[0]));
    if (!numbers.length) throw new Error(`Unsupported feature location: ${location}`);
    return { start: numbers[0], end: numbers.at(-1), reverse: reverse !== partReverse };
  });
  return reverse ? parts.reverse() : parts;
}

function unwrap(value, fn) {
  const prefix = `${fn}(`;
  return value.startsWith(prefix) && value.endsWith(")") ? value.slice(prefix.length, -1) : null;
}

function splitTopLevel(value) {
  const parts = [];
  let depth = 0;
  let start = 0;
  for (let index = 0; index < value.length; index += 1) {
    if (value[index] === "(") depth += 1;
    if (value[index] === ")") depth -= 1;
    if (value[index] === "," && depth === 0) {
      parts.push(value.slice(start, index));
      start = index + 1;
    }
  }
  parts.push(value.slice(start));
  return parts;
}

export function convert(records, format) {
  const converters = { fna: genomeFasta, orfs: orfFasta, faa: proteinFasta, functions: functionsTsv, gff3, ptt };
  if (!converters[format]) throw new Error(`Unknown output format: ${format}`);
  return converters[format](records);
}

function genomeFasta(records) {
  return records.map(record => fasta(record.id, record.sequence, record.definition)).join("");
}

function orfFasta(records) {
  return codingFeatures(records).map(({ record, feature, index }) => fasta(featureId(feature, index), featureSequence(record, feature))).join("");
}

function proteinFasta(records) {
  return codingFeatures(records).map(({ record, feature, index }) => {
    const translation = first(feature, "translation") || translate(featureSequence(record, feature), Number(first(feature, "codon_start") || 1) - 1);
    return fasta(featureId(feature, index), translation.replace(/\*$/, ""), first(feature, "product"));
  }).join("");
}

function functionsTsv(records) {
  const rows = ["sequence_id\tprotein_id\tproduct"];
  for (const { record, feature, index } of codingFeatures(records)) rows.push([record.id, featureId(feature, index), first(feature, "product") || ""].join("\t"));
  return `${rows.join("\n")}\n`;
}

function gff3(records) {
  const lines = ["##gff-version 3"];
  for (const record of records) {
    lines.push(`##sequence-region ${escapeGff(record.id)} 1 ${record.sequence.length}`);
    record.features.forEach((feature, index) => {
      const attributes = { ...feature.qualifiers };
      if (feature.parts.length > 1 && !attributes.ID) attributes.ID = [featureId(feature, index + 1, record.id)];
      feature.parts.forEach((part, partIndex) => {
        const strand = part.reverse ? "-" : "+";
        const phase = feature.type === "CDS" ? cdsPhase(feature, partIndex) : ".";
        lines.push([record.id, "genbank_to", feature.type, part.start, part.end, ".", strand, phase, formatAttributes(attributes)].join("\t"));
      });
    });
  }
  lines.push("##FASTA", genomeFasta(records).trimEnd());
  return `${lines.join("\n")}\n`;
}

function ptt(records) {
  const rows = [];
  for (const record of records) {
    const cds = record.features.filter(feature => feature.type === "CDS");
    rows.push(`${record.definition || record.id} - 1..${record.sequence.length}`, `${cds.length} proteins`, "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct");
    cds.forEach((feature, index) => {
      const start = Math.min(...feature.parts.map(part => part.start));
      const end = Math.max(...feature.parts.map(part => part.end));
      const aaLength = Math.floor(featureSequence(record, feature).length / 3);
      rows.push([`${start}..${end}`, feature.parts[0]?.reverse ? "-" : "+", aaLength, featureId(feature, index + 1), first(feature, "gene") || "-", first(feature, "locus_tag") || "-", "-", "-", first(feature, "product") || "-"].join("\t"));
    });
  }
  return `${rows.join("\n")}\n`;
}

function codingFeatures(records) {
  const result = [];
  for (const record of records) record.features.filter(feature => feature.type === "CDS").forEach((feature, index) => result.push({ record, feature, index: index + 1 }));
  return result;
}

function featureSequence(record, feature) {
  return feature.parts.map(part => {
    const sequence = record.sequence.slice(part.start - 1, part.end);
    return part.reverse ? reverseComplement(sequence) : sequence;
  }).join("");
}

function reverseComplement(sequence) {
  return [...sequence].reverse().map(base => COMPLEMENT[base] || "N").join("");
}

function featureId(feature, index, recordId = "feature") {
  return first(feature, "protein_id") || first(feature, "locus_tag") || first(feature, "ID") || `${recordId}:${feature.type}:${index}`;
}

function first(feature, key) { return feature.qualifiers[key]?.[0]; }

function fasta(id, sequence, description = "") {
  const header = description ? `>${id} ${description}` : `>${id}`;
  const wrapped = sequence.match(/.{1,60}/g)?.join("\n") || "";
  return `${header}\n${wrapped}\n`;
}

function translate(sequence, offset) {
  const table = { TTT:"F",TTC:"F",TTA:"L",TTG:"L",TCT:"S",TCC:"S",TCA:"S",TCG:"S",TAT:"Y",TAC:"Y",TAA:"*",TAG:"*",TGT:"C",TGC:"C",TGA:"*",TGG:"W",CTT:"L",CTC:"L",CTA:"L",CTG:"L",CCT:"P",CCC:"P",CCA:"P",CCG:"P",CAT:"H",CAC:"H",CAA:"Q",CAG:"Q",CGT:"R",CGC:"R",CGA:"R",CGG:"R",ATT:"I",ATC:"I",ATA:"I",ATG:"M",ACT:"T",ACC:"T",ACA:"T",ACG:"T",AAT:"N",AAC:"N",AAA:"K",AAG:"K",AGT:"S",AGC:"S",AGA:"R",AGG:"R",GTT:"V",GTC:"V",GTA:"V",GTG:"V",GCT:"A",GCC:"A",GCA:"A",GCG:"A",GAT:"D",GAC:"D",GAA:"E",GAG:"E",GGT:"G",GGC:"G",GGA:"G",GGG:"G" };
  let protein = "";
  for (let index = offset; index + 2 < sequence.length; index += 3) protein += table[sequence.slice(index, index + 3)] || "X";
  return protein;
}

function cdsPhase(feature, partIndex) {
  const initial = Number(first(feature, "codon_start") || 1) - 1;
  if (partIndex === 0) return String(initial);
  const previousBases = feature.parts.slice(0, partIndex).reduce((total, part) => total + part.end - part.start + 1, -initial);
  return String((3 - (previousBases % 3)) % 3);
}

function escapeGff(value) { return encodeURIComponent(String(value)).replace(/%20/g, " ").replace(/%3A/gi, ":"); }
function formatAttributes(attributes) {
  const entries = Object.entries(attributes).sort(([a], [b]) => a.localeCompare(b));
  if (!entries.length) return ".";
  return entries.map(([key, values]) => `${escapeGff(key)}=${values.map(escapeGff).join(",")}`).join(";");
}
