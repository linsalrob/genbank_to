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

test("rejects remote feature locations", () => {
  assert.throws(
    () => parseLocation("J00194.1:100..202"),
    /Remote feature locations are not supported/,
  );
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
  assert.match(convert(records, "ptt"), /130\.\.1311\t\+\t393\tYP_009640124\.1/);
});

test("uses annotated translation tables for un-translated CDS features", () => {
  const table4 = record("ATGTGATAA", feature(1, 9, { protein_id: ["table4"], transl_table: ["4"] }));
  const table11 = record("GTGAAATAG", feature(1, 9, { protein_id: ["table11"], transl_table: ["11"] }));

  assert.match(convert([table4], "faa"), />table4\nMW\n/);
  assert.match(convert([table11], "faa"), />table11\nMK\n/);
});

test("skips pseudo CDS features from protein-producing outputs", () => {
  const pseudo = feature(1, 9, { protein_id: ["pseudo1"], pseudo: [""] });
  const coding = feature(10, 18, { protein_id: ["coding1"], translation: ["MK"] });
  const records = [record("ATGAAATAAATGAAATAA", pseudo, coding)];

  for (const format of ["orfs", "faa", "functions", "ptt"]) {
    const output = convert(records, format);
    assert.doesNotMatch(output, /pseudo1/);
    assert.match(output, /coding1/);
  }
  assert.match(convert(records, "gff3"), /pseudo=/);
});

test("does not count a terminal stop codon in PTT protein lengths", () => {
  const records = [record("ATGAAATAA", feature(1, 9, { protein_id: ["protein1"], translation: ["MK"] }))];
  assert.match(convert(records, "ptt"), /1\.\.9\t\+\t2\tprotein1/);
});

function feature(start, end, qualifiers = {}) {
  return { type: "CDS", location: `${start}..${end}`, qualifiers, parts: [{ start, end, reverse: false }] };
}

function record(sequence, ...features) {
  return { id: "record1", definition: "Test record", sequence, features };
}
