import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";


const projectRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const outputsDir = path.join(projectRoot, "outputs");
const dataPath = path.join(outputsDir, "Table1_workbook_data.json");
const outputPath = path.join(outputsDir, "Table1_participant_characteristics.xlsx");
const previewDir = path.join(outputsDir, "table1_workbook_previews");

const payload = JSON.parse(await fs.readFile(dataPath, "utf8"));
const workbook = Workbook.create();

const colors = {
  navy: "#1F4E78",
  blue: "#D9EAF7",
  pale: "#F4F7FA",
  border: "#B7C4CE",
  text: "#1F2933",
  muted: "#536471",
  white: "#FFFFFF",
};

function writeTableSheet({ name, title, subtitle, matrix, notes, tabColor }) {
  const sheet = workbook.worksheets.add(name);
  sheet.showGridLines = false;
  sheet.tabColor = tabColor;

  const colCount = matrix[0].length;
  const lastCol = String.fromCharCode(64 + colCount);
  sheet.mergeCells(`A1:${lastCol}1`);
  sheet.getRange("A1").values = [[title]];
  sheet.getRange(`A1:${lastCol}1`).format = {
    fill: colors.navy,
    font: { bold: true, color: colors.white, size: 15, name: "Arial" },
    verticalAlignment: "center",
  };
  sheet.getRange(`A1:${lastCol}1`).format.rowHeight = 30;

  sheet.mergeCells(`A2:${lastCol}2`);
  sheet.getRange("A2").values = [[subtitle]];
  sheet.getRange(`A2:${lastCol}2`).format = {
    font: { italic: true, color: colors.muted, size: 9, name: "Arial" },
    verticalAlignment: "center",
    wrapText: true,
  };
  sheet.getRange(`A2:${lastCol}2`).format.rowHeight = 28;

  const startRow = 4;
  const endRow = startRow + matrix.length - 1;
  const tableRange = sheet.getRange(`A${startRow}:${lastCol}${endRow}`);
  tableRange.values = matrix;
  tableRange.format = {
    font: { color: colors.text, size: 9.5, name: "Arial" },
    verticalAlignment: "center",
    wrapText: true,
    borders: {
      insideHorizontal: { style: "thin", color: colors.border },
      bottom: { style: "thin", color: colors.border },
    },
  };

  const header = sheet.getRange(`A${startRow}:${lastCol}${startRow}`);
  header.format = {
    fill: colors.navy,
    font: { bold: true, color: colors.white, size: 9.5, name: "Arial" },
    verticalAlignment: "center",
    wrapText: true,
    borders: { preset: "all", style: "thin", color: colors.white },
  };
  header.format.rowHeight = 35;

  for (let r = 1; r < matrix.length; r += 1) {
    const excelRow = startRow + r;
    const label = String(matrix[r][0] ?? "");
    const rowRange = sheet.getRange(`A${excelRow}:${lastCol}${excelRow}`);
    if (label === "Clinical diagnosis" || label === "T reference") {
      rowRange.format.fill = colors.blue;
      rowRange.format.font = { bold: true, color: colors.navy, size: 9.5, name: "Arial" };
      rowRange.format.rowHeight = 22;
    } else if (r % 2 === 0) {
      rowRange.format.fill = colors.pale;
    }
  }

  sheet.getRange(`A${startRow + 1}:A${endRow}`).format.font = {
    color: colors.text,
    size: 9.5,
    name: "Arial",
  };
  sheet.getRange(`B${startRow + 1}:${lastCol}${endRow}`).format.horizontalAlignment = "center";
  sheet.getRange("A:A").format.columnWidth = 34;
  sheet.getRange(`B:${lastCol}`).format.columnWidth = 23;

  const noteStart = endRow + 2;
  sheet.getRange(`A${noteStart}`).values = [["Notes"]];
  sheet.getRange(`A${noteStart}`).format.font = {
    bold: true,
    color: colors.navy,
    size: 10,
    name: "Arial",
  };
  notes.forEach((note, i) => {
    const row = noteStart + 1 + i;
    sheet.mergeCells(`A${row}:${lastCol}${row}`);
    sheet.getRange(`A${row}`).values = [[`${i + 1}. ${note}`]];
    sheet.getRange(`A${row}:${lastCol}${row}`).format = {
      font: { color: colors.muted, size: 8.5, name: "Arial" },
      wrapText: true,
      verticalAlignment: "top",
    };
    sheet.getRange(`A${row}:${lastCol}${row}`).format.rowHeight = 26;
  });

  sheet.freezePanes.freezeRows(startRow);
  sheet.freezePanes.freezeColumns(1);
  return { sheet, endRow: noteStart + notes.length };
}

const table1Notes = [
  ...payload.metadata.notes,
  "MMSE and MoCA-P cells include the available sample size; missing values were not imputed.",
  "AD/AD-like combines AD, ADlike, and recorded mixed diagnoses containing an AD/AD-like component; all remaining diagnoses are grouped as other clinical diagnoses.",
];

writeTableSheet({
  name: "Table 1",
  title: "Table 1. Participant characteristics by A/T group",
  subtitle: "Overall cohort, classified A/T groups, and participants without an available T reference",
  matrix: payload.table1,
  notes: table1Notes,
  tabColor: colors.navy,
});

writeTableSheet({
  name: "Biomarkers",
  title: "Supplementary table. Absolute biomarker concentrations by A/T group",
  subtitle: "Values are median [IQR]; units must be confirmed from the laboratory report before submission",
  matrix: payload.biomarkers,
  notes: [
    "The raw workbook column MTBR-Tau243 is displayed scientifically as eMTBR-tau243, based on manufacturer confirmation.",
    "Ratios are omitted from this descriptive table because the absolute concentrations are the primary measurements.",
    "No missing biomarker values were identified in the source workbook.",
  ],
  tabColor: "#5B9BD5",
});

workbook.recalculate();

const check1 = await workbook.inspect({
  kind: "table",
  sheetId: "Table 1",
  range: "A1:F28",
  include: "values,formulas",
  tableMaxRows: 28,
  tableMaxCols: 6,
  maxChars: 12000,
});
console.log(check1.ndjson);

const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log(errors.ndjson);

await fs.mkdir(outputsDir, { recursive: true });
await fs.mkdir(previewDir, { recursive: true });
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);

for (const sheetName of ["Table 1", "Biomarkers"]) {
  const preview = await workbook.render({
    sheetName,
    autoCrop: "all",
    scale: 1.5,
    format: "png",
  });
  const bytes = new Uint8Array(await preview.arrayBuffer());
  await fs.writeFile(path.join(previewDir, `${sheetName.replaceAll(" ", "_")}.png`), bytes);
}

console.log(JSON.stringify({ outputPath, previewDir }, null, 2));
