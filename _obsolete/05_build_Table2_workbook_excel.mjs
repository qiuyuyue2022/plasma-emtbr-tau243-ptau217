import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

process.on("uncaughtException", (error) => {
  console.error(`Workbook build failed: ${error.name}: ${error.message}`);
  process.exit(1);
});
process.on("unhandledRejection", (error) => {
  console.error(`Workbook build failed: ${error?.name ?? "Error"}: ${error?.message ?? error}`);
  process.exit(1);
});

const scriptDir = path.dirname(fileURLToPath(import.meta.url));
const outputDir = process.env.BIOMARKER_OUTPUT_DIR ?? path.join(scriptDir, "..", "outputs");
const outputPath = path.join(outputDir, "Table2_key_AT_performance_workbook.xlsx");
const records = JSON.parse(await fs.readFile(`${outputDir}/biomarker_performance_results.json`, "utf8"));
const sampleSize = (analysis) => {
  const row = records.find((record) => record.Analysis === analysis);
  if (!row) throw new Error(`Missing analysis in results file: ${analysis}`);
  return row.N;
};

const workbook = Workbook.create();
const fontFamily = "Arial";
const headerFill = "#1F4E78";
const sectionFill = "#D9EAF7";
const lightBlue = "#EDF4FA";

const sheets = [
  { name: "A prediction", analysis: "A status", title: "Biomarker performance for A status", note: "A+ versus A−; all participants" },
  { name: "A non-control", analysis: "A status, non-control cohort", title: "A prediction after excluding healthy controls", note: "A+ versus A− in the non-control clinical cohort" },
  { name: "T prediction", analysis: "T status", title: "Biomarker performance for T status", note: "T+ versus T−; Tau PET/CSF" },
  { name: "T Tau PET", analysis: "T status, Tau PET subset", title: "T prediction among participants with Tau PET results", note: "T+ versus T−; Tau PET" },
  { name: "T in A+", analysis: "T status, A+ subset", title: "T prediction among A+ participants", note: "T+ versus T− in the A+ cohort" },
];

for (const spec of sheets) {
  const sheet = workbook.worksheets.add(spec.name);
  const rows = records.filter((row) => row.Analysis === spec.analysis).sort((a, b) => b.AUC - a.AUC);
  if (!rows.length) throw new Error(`Missing analysis in results file: ${spec.analysis}`);
  const first = rows[0];

  sheet.showGridLines = false;
  sheet.getRange("A2:J2").format.borders = { preset: "outside", style: "thin", color: "#9FBAD0" };
  sheet.getRange("A2").values = [[spec.title]];
  sheet.getRange("A2").format.font = { name: fontFamily, size: 14, bold: true, color: "#17365D" };
  sheet.getRange("A3").values = [[`${spec.note}; n=${first.N} (${first["Positive n"]} positive, ${first["Negative n"]} negative)`]];
  sheet.getRange("A3").format.font = { name: fontFamily, size: 9, italic: true, color: "#595959" };

  const headers = [[
    "Biomarker", "AUC", "AUC 95% CI", "CV accuracy", "CV balanced accuracy",
    "CV sensitivity", "CV specificity", "Cutoff direction", "Cutoff", "N"
  ]];
  sheet.getRange("A5:J5").values = headers;
  sheet.getRange("A5:J5").format = {
    fill: headerFill,
    font: { name: fontFamily, size: 9, bold: true, color: "#FFFFFF" },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    wrapText: true,
    borders: { preset: "all", style: "thin", color: "#FFFFFF" },
  };

  const values = rows.map((row) => [
    row.Biomarker,
    row.AUC,
    `${row["AUC 95% CI lower"].toFixed(3)}–${row["AUC 95% CI upper"].toFixed(3)}`,
    row["CV accuracy"],
    row["CV balanced accuracy"],
    row["CV sensitivity"],
    row["CV specificity"],
    row["Cutoff direction"],
    row.Cutoff,
    row.N,
  ]);
  sheet.getRange(`A6:J${5 + values.length}`).values = values;
  const body = sheet.getRange(`A6:J${5 + values.length}`);
  body.format.font = { name: fontFamily, size: 9, color: "#222222" };
  body.format.verticalAlignment = "center";
  body.format.borders = { insideHorizontal: { style: "thin", color: "#D9E2F3" } };
  sheet.getRange(`A6:A${5 + values.length}`).format.font = { name: fontFamily, size: 9, bold: true, color: "#17365D" };
  sheet.getRange(`B6:B${5 + values.length}`).format.numberFormat = "0.000";
  sheet.getRange(`D6:G${5 + values.length}`).format.numberFormat = "0.0%";
  sheet.getRange(`H6:H${5 + values.length}`).format.horizontalAlignment = "center";
  sheet.getRange(`I6:I${5 + values.length}`).format.numberFormat = "0.000";
  sheet.getRange(`J6:J${5 + values.length}`).format.numberFormat = "0";

  const maxAuc = Math.max(...rows.map((row) => row.AUC));
  const maxAccuracy = Math.max(...rows.map((row) => row["CV accuracy"]));
  const maxBalanced = Math.max(...rows.map((row) => row["CV balanced accuracy"]));
  rows.forEach((row, index) => {
    const excelRow = 6 + index;
    if (row.AUC === maxAuc) sheet.getRange(`B${excelRow}:C${excelRow}`).format.font = { name: fontFamily, size: 9, bold: true, color: "#17365D" };
    if (row["CV accuracy"] === maxAccuracy) sheet.getRange(`D${excelRow}`).format.font = { name: fontFamily, size: 9, bold: true, color: "#17365D" };
    if (row["CV balanced accuracy"] === maxBalanced) sheet.getRange(`E${excelRow}`).format.font = { name: fontFamily, size: 9, bold: true, color: "#17365D" };
  });

  const notesRow = 7 + values.length;
  sheet.getRange(`A${notesRow}:J${notesRow}`).values = [[
    "Notes", "AUC CI: 3,000 subject-level bootstrap resamples.", "", "CV metrics: mean of 50 repeated stratified 5-fold validations.", "", "", "", "Cutoff: full-sample Youden index; exploratory, not externally validated.", "", ""
  ]];
  sheet.getRange(`A${notesRow}:J${notesRow}`).format = {
    fill: lightBlue,
    font: { name: fontFamily, size: 8, italic: true, color: "#44546A" },
    verticalAlignment: "center",
    wrapText: true,
    borders: { preset: "outside", style: "thin", color: "#B4C7E7" },
  };

  sheet.getRange(`A1:J${notesRow}`).format.verticalAlignment = "center";
  sheet.getRange("A:A").format.columnWidth = 18;
  sheet.getRange("B:B").format.columnWidth = 10;
  sheet.getRange("C:C").format.columnWidth = 17;
  sheet.getRange("D:G").format.columnWidth = 15;
  sheet.getRange("H:H").format.columnWidth = 13;
  sheet.getRange("I:I").format.columnWidth = 11;
  sheet.getRange("J:J").format.columnWidth = 8;
  sheet.getRange("5:5").format.rowHeight = 30;
  sheet.getRange(`${notesRow}:${notesRow}`).format.rowHeight = 56;
  sheet.freezePanes.freezeRows(5);
}

const methods = workbook.worksheets.add("Methods");
methods.showGridLines = false;
methods.getRange("A2:F2").format.borders = { preset: "outside", style: "thin", color: "#9FBAD0" };
methods.getRange("A2").values = [["Analysis definitions and statistical methods"]];
methods.getRange("A2").format.font = { name: fontFamily, size: 14, bold: true, color: "#17365D" };
methods.getRange("A4:B4").values = [["Item", "Definition"]];
methods.getRange("A4:B4").format = {
  fill: headerFill,
  font: { name: fontFamily, size: 9, bold: true, color: "#FFFFFF" },
  horizontalAlignment: "center",
  verticalAlignment: "center",
};
methods.getRange("A5:B14").values = [
  ["A analysis", `A+ versus A− in all ${sampleSize("A status")} participants.`],
  ["A non-control analysis", "A+ versus A− after excluding records with DIAG = 健康对照. The outcome is A status, not the diagnosis label itself."],
  ["T analysis", `T+ versus T− in ${sampleSize("T status")} participants assessed by Tau PET or CSF.`],
  ["Tau PET subset", `T+ versus T− among ${sampleSize("T status, Tau PET subset")} participants assessed by Tau PET.`],
  ["T analysis in A+", `T+ versus T− among ${sampleSize("T status, A+ subset")} A+ participants.`],
  ["AUC", "Single-marker ROC AUC, oriented so that AUC ≥ 0.5; direction is shown with each cutoff."],
  ["AUC 95% CI", "Percentile interval from 3,000 nonparametric subject-level bootstrap resamples of the same analysis cohort."],
  ["Cutoff", "Full-sample threshold maximizing the Youden index (sensitivity + specificity − 1). Exploratory only."],
  ["Cross-validation", "Stratified 5-fold validation repeated 50 times. Direction and Youden cutoff were selected in each training fold and evaluated in its held-out fold."],
  ["Reported CV metrics", "Mean test-fold accuracy, balanced accuracy, sensitivity, and specificity across 250 held-out folds."],
];
methods.getRange("A5:B14").format.font = { name: fontFamily, size: 9, color: "#222222" };
methods.getRange("A5:B14").format.verticalAlignment = "top";
methods.getRange("A5:B14").format.wrapText = true;
methods.getRange("A5:B14").format.borders = { insideHorizontal: { style: "thin", color: "#D9E2F3" } };
methods.getRange("A5:A14").format = { fill: sectionFill, font: { name: fontFamily, size: 9, bold: true, color: "#17365D" } };
methods.getRange("A:A").format.columnWidth = 22;
methods.getRange("B:B").format.columnWidth = 82;
methods.getRange("5:14").format.rowHeight = 32;

workbook.recalculate();

const check = await workbook.inspect({
  kind: "table",
  sheetId: "A prediction",
  range: "A2:J14",
  include: "values,formulas",
  tableMaxRows: 20,
  tableMaxCols: 12,
});
console.log(check.ndjson);

const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log(errors.ndjson);

await fs.mkdir(outputDir, { recursive: true });
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
