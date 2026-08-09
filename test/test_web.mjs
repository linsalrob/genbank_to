import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";
import { convert, parseGenBank, parseLocation } from "../web/genbank.mjs";

test("parses the repository GenBank fixture", async () => {
  const source = await readFile(new URL("./NC_001417.gbk", import.meta.url), "utf8");
  const records = parseGenBank(source);

  assert.equal(records.length, 1);
  assert.equal(records[0].id, "NC_001417.2");
  assert.equal(records[0].sequence.length, 3569);
  assert.ok(records[0].features.some(feature => feature.type === "CDS"));
});

test("extracts joined and reverse-complement locations", () => {
  assert.deepEqual(parseLocation("join(1..3,8..10)"), [
    { start: 1, end: 3, reverse: false },
    { start: 8, end: 10, reverse: false },
  ]);
  assert.deepEqual(parseLocation("complement(join(1..3,8..10))"), [
    { start: 8, end: 10, reverse: true },
    { start: 1, end: 3, reverse: true },
  ]);
});

test("generates each browser output format", async () => {
  const source = await readFile(new URL("./NC_001417.gbk", import.meta.url), "utf8");
  const records = parseGenBank(source);

  for (const format of ["fna", "orfs", "faa", "functions", "gff3", "ptt"]) {
    const output = convert(records, format);
    assert.ok(output.length > 100, `${format} output should not be empty`);
  }
  assert.match(convert(records, "gff3"), /^##gff-version 3/);
  assert.match(convert(records, "functions"), /^sequence_id\tprotein_id\tproduct/);
});
