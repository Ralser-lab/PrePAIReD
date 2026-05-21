#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path
import os, platform

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    Image,
    ListFlowable,
    ListItem,
    KeepTogether,
)

# ================================================================
# Input helpers
# ================================================================

def find_single_params_file(workdir: Path) -> Path:
    matches = sorted(workdir.glob("*_Parameters.json"))
    if not matches:
        raise FileNotFoundError(f"No '*_Parameters.json' found in {workdir.resolve()}")
    if len(matches) > 1:
        raise RuntimeError(
            "Multiple '*_Parameters.json' files found:\n"
            + "\n".join(f"  - {p.name}" for p in matches)
        )
    return matches[0]

def infer_title_from_params_file(params_path: Path) -> str:
    suffix = "_Parameters.json"
    name = params_path.name
    return name[: -len(suffix)] if name.endswith(suffix) else params_path.stem

def infer_step_from_title(title: str) -> str:
    m = re.match(r"^\d+_(.+)$", title.strip())
    return (m.group(1) if m else title).strip()

def load_params_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError(f"Expected JSON object in {path}, got {type(data)}")
    return data

def extract_tagged_warnings_errors(log_path: Path) -> Tuple[List[str], List[str]]:
    warnings: List[str] = []
    errors: List[str] = []
    seen_w, seen_e = set(), set()

    def clean_msg(s: str) -> str:
        s = s.strip()
        s = re.sub(r"^:\s*", "", s)  # remove leading ":" artifacts
        return s

    with log_path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n").strip()
            if not line:
                continue

            if "###ERROR###" in line:
                msg = clean_msg(line.replace("###ERROR###", "").strip())
                if msg and msg not in seen_e:
                    errors.append(msg)
                    seen_e.add(msg)

            if "###WARNING###" in line:
                msg = clean_msg(line.replace("###WARNING###", "").strip())
                if msg and msg not in seen_w:
                    warnings.append(msg)
                    seen_w.add(msg)

    return warnings, errors

# ================================================================
# Documentation TSV
# ================================================================

@dataclass(frozen=True)
class MethodDoc:
    language: str
    step: str
    flag_variable: str
    flag_value: str
    short_description: str
    reference: str
    package: str

def _norm_header(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "").strip())

def find_doc_tsv(workdir: Path, preferred: Optional[str]) -> Optional[Path]:
    if preferred:
        p = Path(preferred)
        if not p.is_absolute():
            p = workdir / p
        if p.exists():
            return p

    candidates = sorted(workdir.glob("*Documentation*.tsv"))
    return candidates[0] if candidates else None

def load_method_docs(tsv_path: Path) -> List[MethodDoc]:
    docs: List[MethodDoc] = []
    with tsv_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            return docs

        header_map = {_norm_header(h): h for h in reader.fieldnames if h is not None}

        def get(row: Dict[str, str], want: str) -> str:
            key = header_map.get(want, "")
            return (row.get(key, "") or "").strip()

        for row in reader:
            docs.append(
                MethodDoc(
                    language=get(row, "Language"),
                    step=get(row, "Step"),
                    flag_variable=get(row, "Flag_Variable"),
                    flag_value=get(row, "Flag_Value"),
                    short_description=get(row, "Short Description"),
                    reference=get(row, "Reference (Method)"),
                    package=get(row, "Package"),
                )
            )
    return docs

def normalize_references(ref_field: str) -> List[str]:
    """
    Robust reference parsing:
    - ignore '-', empty
    - split on commas and newlines
    - preserve DOIs/URLs
    """
    if not ref_field:
        return []
    s = ref_field.strip()
    if not s or s == "-":
        return []

    # split on comma OR newline
    parts = [p.strip() for p in re.split(r"[,\n\r]+", s) if p.strip()]
    out: List[str] = []
    seen = set()
    for p in parts:
        if not p or p == "-":
            continue
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out

# ================================================================
# Figures + captions
# ================================================================

def load_figure_captions_tsv(tsv_path: Path) -> Dict[str, Dict[str, str]]:
    """Load figure metadata keyed by Figure ID (e.g., 'F01').

    Returns:
      meta[FigureID] = {"title": Title, "description": Description}
    """
    meta: Dict[str, Dict[str, str]] = {}
    if not tsv_path.exists():
        return meta

    with tsv_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            return meta

        # Build a case-insensitive header lookup
        header_ci = {(_norm_header(h).lower()): h for h in reader.fieldnames if h}

        def pick(*cands: str) -> str:
            for c in cands:
                h = header_ci.get(_norm_header(c).lower())
                if h:
                    return h
            return ""

        id_col = pick("Figure ID", "Figure_ID", "FigureID", "ID", "Figure")
        title_col = pick("Title", "Figure Title", "Figure_Title")
        desc_col = pick("Description", "DEscription", "Caption", "Figure Description")

        # fallback: first col = id, second col = title, third col = desc (best effort)
        fallback_cols = [h for h in reader.fieldnames if h]

        for row in reader:
            if id_col:
                fig_id = (row.get(id_col, "") or "").strip()
            else:
                fig_id = (row.get(fallback_cols[0], "") or "").strip() if fallback_cols else ""

            if not fig_id:
                continue

            if title_col:
                title = (row.get(title_col, "") or "").strip()
            else:
                title = (row.get(fallback_cols[1], "") or "").strip() if len(fallback_cols) > 1 else ""

            if desc_col:
                desc = (row.get(desc_col, "") or "").strip()
            else:
                desc = (row.get(fallback_cols[2], "") or "").strip() if len(fallback_cols) > 2 else ""

            meta[fig_id] = {"title": title, "description": desc}

    return meta
    
def list_figures_sorted(figures_dir: Path) -> List[Path]:
    if not figures_dir.exists() or not figures_dir.is_dir():
        return []

    figs = sorted(figures_dir.glob("*.png"))

    def sort_key(p: Path) -> Tuple[int, str]:
        stem = p.stem
        m = re.search(r"(\d+)", stem)  # numeric part after removing F-like prefix
        n = int(m.group(1)) if m else 10**9
        return (n, stem)

    return sorted(figs, key=sort_key)

def load_figure_log_captions(log_path: Path) -> Dict[str, str]:
    """
    Parse report_folder/Figure_log.txt:

      ###FIGURE###: report_folder/F01.png
      ###CAPTION###: some caption ...

    Returns:
      captions["F01"] = "some caption ..."
    """
    captions: Dict[str, str] = {}
    if not log_path.exists():
        return captions

    current_fig_id: Optional[str] = None

    with log_path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            if line.startswith("###FIGURE###:"):
                # example: ###FIGURE###: report_folder/F01.png
                fig_path = line.split(":", 1)[1].strip()
                current_fig_id = Path(fig_path).stem  # -> "F01"
                continue

            if line.startswith("###CAPTION###:") and current_fig_id:
                cap = line.split(":", 1)[1].strip()
                captions[current_fig_id] = cap
                current_fig_id = None
                continue

    return captions
# ================================================================
# Matching logic (improved)
# ================================================================

def param_value_contains_flag_value(param_val: Any, flag_value: str) -> bool:
    if param_val is None:
        return False

    fv = str(flag_value).strip()
    if not fv:
        return False

    if isinstance(param_val, str):
        s = param_val.strip()
        if s == fv:
            return True
        parts = [p.strip() for p in re.split(r"[;,]", s) if p.strip()]
        return fv in parts

    if isinstance(param_val, (list, tuple, set)):
        return any(param_value_contains_flag_value(x, fv) for x in param_val)

    if isinstance(param_val, dict):
        return any(param_value_contains_flag_value(x, fv) for x in param_val.values())

    return str(param_val).strip() == fv

# ================================================================
# PDF styling helpers
# ================================================================

def escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def build_title_banner(title: str, subtitle: str, styles) -> Table:
    title_style = ParagraphStyle(
        "BannerTitle",
        parent=styles["Heading1"],
        fontSize=18,
        leading=22,
        spaceAfter=2,
        textColor=colors.HexColor("#0F172A"),
    )
    sub_style = ParagraphStyle(
        "BannerSub",
        parent=styles["BodyText"],
        fontSize=9.5,
        leading=12,
        textColor=colors.HexColor("#475569"),
    )

    data = [
        [Paragraph(escape(title), title_style)],
        [Paragraph(escape(subtitle), sub_style)],
    ]
    t = Table(data, colWidths=[17.5 * cm])
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#F3F6FB")),
                ("BOX", (0, 0), (-1, -1), 0.8, colors.HexColor("#D7E3F4")),
                ("LEFTPADDING", (0, 0), (-1, -1), 10),
                ("RIGHTPADDING", (0, 0), (-1, -1), 10),
                ("TOPPADDING", (0, 0), (-1, -1), 10),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 10),
            ]
        )
    )
    return t

def section_heading(text: str, number: int, styles) -> Paragraph:
    h = ParagraphStyle(
        f"SecH{number}",
        parent=styles["Heading2"],
        fontSize=12.5,
        leading=15,
        textColor=colors.HexColor("#111827"),
        spaceBefore=10,
        spaceAfter=6,
    )
    return Paragraph(f"{number}. {escape(text)}", h)

def build_params_table(params: Dict[str, Any]) -> Table:
    base = getSampleStyleSheet()
    header_style = ParagraphStyle("P9B", parent=base["BodyText"], fontName="Helvetica-Bold", fontSize=9, leading=11)
    cell_style = ParagraphStyle("P88", parent=base["BodyText"], fontSize=8.8, leading=11)

    rows = [
        [Paragraph("Parameter", header_style), Paragraph("Value", header_style)]
    ]

    for k in sorted(params, key=lambda x: str(x).lower()):
        v = params[k]
        if isinstance(v, (dict, list)):
            v = json.dumps(v, ensure_ascii=False)
        v_str = "" if v is None else str(v)
        rows.append([Paragraph(escape(str(k)), cell_style), Paragraph(escape(v_str), cell_style)])

    table = Table(rows, colWidths=[6.5 * cm, 11.0 * cm])
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EEF2F7")),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#D1D5DB")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#FBFCFE")]),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return table

def method_card(display_name: str, docrow: MethodDoc, styles) -> Table:
    title = Paragraph(
        escape(display_name),
        ParagraphStyle(
            "MethodName",
            parent=styles["Heading3"],
            fontSize=11,
            leading=13,
            textColor=colors.HexColor("#0F172A"),
            spaceAfter=4,
        ),
    )

    desc = Paragraph(
        escape(docrow.short_description or "No description available."),
        ParagraphStyle(
            "MethodDesc",
            parent=styles["BodyText"],
            fontSize=9.5,
            leading=12.5,
            textColor=colors.HexColor("#111827"),
        ),
    )

    meta_parts: List[str] = []
    if docrow.language:
        meta_parts.append(f"<b>Language:</b> {escape(docrow.language)}")
    if docrow.package:
        meta_parts.append(f"<b>Package:</b> {escape(docrow.package)}")

    meta = Paragraph(
        " &nbsp;&nbsp;•&nbsp;&nbsp; ".join(meta_parts) if meta_parts else "",
        ParagraphStyle(
            "MethodMeta",
            parent=styles["BodyText"],
            fontSize=9,
            leading=11,
            textColor=colors.HexColor("#374151"),
        ),
    )

    data = [[title], [desc], [Spacer(1, 6)], [meta]]
    t = Table(data, colWidths=[17.5 * cm])
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#FFFFFF")),
                ("BOX", (0, 0), (-1, -1), 0.8, colors.HexColor("#E5E7EB")),
                ("LEFTPADDING", (0, 0), (-1, -1), 10),
                ("RIGHTPADDING", (0, 0), (-1, -1), 10),
                ("TOPPADDING", (0, 0), (-1, -1), 8),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
            ]
        )
    )
    return t

def callout(title: str, items: List[str], kind: str, styles) -> Table:
    if kind == "error":
        bg = colors.HexColor("#FEF2F2")
        border = colors.HexColor("#FCA5A5")
        head = colors.HexColor("#991B1B")
    elif kind == "warning":
        bg = colors.HexColor("#FFFBEB")
        border = colors.HexColor("#FCD34D")
        head = colors.HexColor("#92400E")
    else:
        bg = colors.HexColor("#EFF6FF")
        border = colors.HexColor("#93C5FD")
        head = colors.HexColor("#1D4ED8")

    head_style = ParagraphStyle(
        f"CalloutHead_{kind}",
        parent=styles["BodyText"],
        fontName="Helvetica-Bold",
        fontSize=9.5,
        leading=12,
        textColor=head,
        spaceAfter=4,
    )
    body_style = ParagraphStyle(
        f"CalloutBody_{kind}",
        parent=styles["BodyText"],
        fontSize=9,
        leading=11.5,
        textColor=colors.HexColor("#111827"),
    )

    cleaned = [re.sub(r"^:\s*", "", x).strip() for x in items if x.strip()]
    if not cleaned:
        cleaned = ["None"]

    lst = ListFlowable(
        [ListItem(Paragraph(escape(x), body_style)) for x in cleaned],
        bulletType="bullet",
        leftIndent=16,
    )

    data = [[Paragraph(escape(title), head_style)], [lst]]
    t = Table(data, colWidths=[17.5 * cm])
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), bg),
                ("BOX", (0, 0), (-1, -1), 0.9, border),
                ("LEFTPADDING", (0, 0), (-1, -1), 10),
                ("RIGHTPADDING", (0, 0), (-1, -1), 10),
                ("TOPPADDING", (0, 0), (-1, -1), 8),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
            ]
        )
    )
    return t

# ================================================================
# Per-module PDF chunk (NO PageBreak here)
# ================================================================

def build_module_pdf(
    title: str,
    params: Dict[str, Any],
    warnings: List[str],
    errors: List[str],
    method_sections: List[Tuple[str, MethodDoc]],
    method_references: List[str],
    out_pdf: Path,
    figures_dir: Optional[Path] = None,
    figure_captions: Optional[Dict[str, Dict[str, str]]] = None,
) -> None:
    styles = getSampleStyleSheet()

    doc = SimpleDocTemplate(
        str(out_pdf),
        pagesize=A4,
        rightMargin=2 * cm,
        leftMargin=2 * cm,
        topMargin=1.6 * cm,
        bottomMargin=1.6 * cm,
        title=title,
    )

    step = infer_step_from_title(title)

    story: List[Any] = []
    story.append(build_title_banner(title=title, subtitle=f"Step: {step}", styles=styles))
    story.append(Spacer(1, 10))

    # 1) Parameters
    story.append(section_heading("Parameters", 1, styles))
    story.append(build_params_table(params))

    # 2) Methods
    sec_no = 2
    if method_sections:
        story.append(section_heading("Methods", sec_no, styles))
        for display_name, docrow in method_sections:
            story.append(KeepTogether([method_card(display_name, docrow, styles), Spacer(1, 8)]))

        # References once, at end of Methods
        refs: List[str] = []
        seen = set()
        for r in method_references:
            r = (r or "").strip()
            if not r or r == "-":
                continue
            if r not in seen:
                seen.add(r)
                refs.append(r)

        if refs:
            ref_title = ParagraphStyle(
                "RefTitle",
                parent=styles["Heading3"],
                fontSize=10.5,
                leading=13,
                textColor=colors.HexColor("#111827"),
                spaceBefore=6,
                spaceAfter=4,
            )
            
            story.append(Paragraph(f"{sec_no}.1 References", ref_title))

            ref_mono = ParagraphStyle(
                "RefMono",
                parent=styles["BodyText"],
                fontName="Courier",
                fontSize=8.3,
                leading=11,
                textColor=colors.HexColor("#111827"),
            )

            # Make URLs/DOIs clickable
            def as_link(s: str) -> str:
                s = s.strip()
                if not s:
                    return ""
                # If it's a bare DOI, convert to https://doi.org/...
                if s.lower().startswith("10."):
                    url = f"https://doi.org/{s}"
                    return f'<a href="{escape(url)}">{escape(s)}</a>'
                # If it's already a URL
                if s.lower().startswith("http://") or s.lower().startswith("https://"):
                    return f'<a href="{escape(s)}">{escape(s)}</a>'
                return escape(s)

            ref_rows = []
            for i, r in enumerate(refs, start=1):
                ref_rows.append([
                    Paragraph(f"<b>[{i}]</b>", ParagraphStyle(
                        "RefIdx",
                        parent=styles["BodyText"],
                        fontSize=9,
                        leading=11,
                        textColor=colors.HexColor("#334155"),
                    )),
                    Paragraph(as_link(r), ref_mono),
                ])

            ref_table = Table(ref_rows, colWidths=[1.0 * cm, 16.5 * cm])
            ref_table.setStyle(TableStyle([
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 2),
                ("RIGHTPADDING", (0, 0), (-1, -1), 2),
                ("TOPPADDING", (0, 0), (-1, -1), 2),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            ]))

            # Wrap the references in a subtle box (consistent with your style)
            ref_box = Table([[ref_table]], colWidths=[17.5 * cm])
            ref_box.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#F8FAFC")),
                ("BOX", (0, 0), (-1, -1), 0.8, colors.HexColor("#E2E8F0")),
                ("LEFTPADDING", (0, 0), (-1, -1), 10),
                ("RIGHTPADDING", (0, 0), (-1, -1), 10),
                ("TOPPADDING", (0, 0), (-1, -1), 8),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
            ]))

            story.append(ref_box)

        sec_no += 1

    # Evaluation (Figures)  <-- renamed from Results
    fig_caps = figure_captions or {}
    figs: List[Path] = list_figures_sorted(figures_dir) if figures_dir else []
    if figs:
        story.append(section_heading("Evaluation", sec_no, styles))

        sub_style = ParagraphStyle(
            "ResultsSubhead",
            parent=styles["Heading3"],
            fontSize=10.5,
            leading=13,
            textColor=colors.HexColor("#111827"),
            spaceBefore=6,
            spaceAfter=6,
        )

        cap_style = ParagraphStyle(
            "FigureCaption",
            parent=styles["BodyText"],
            fontSize=9,
            leading=11.5,
            textColor=colors.HexColor("#111827"),
            spaceBefore=4,
            spaceAfter=10,
        )

        max_w = 17.5 * cm
        max_h = 11.5 * cm

        for i, img_path in enumerate(figs, start=1):
            fig_id = img_path.stem  # e.g. "F01"

            meta = fig_caps.get(fig_id, {})
            title_txt = (meta.get("title", "") or "").strip()
            desc = (meta.get("description", "") or "").strip()

            # subsection title = Title column (fallback to fig_id if missing)
            subsection_title = title_txt if title_txt else fig_id
            story.append(Paragraph(f"{sec_no}.{i} {escape(subsection_title)}", sub_style))

            img = Image(str(img_path))
            # Scale image to fit within max_w x max_h while preserving aspect ratio
            try:
                iw, ih = float(img.imageWidth), float(img.imageHeight)
                if iw > 0 and ih > 0:
                    scale = min(max_w / iw, max_h / ih, 1.0)
                    img.drawWidth = iw * scale
                    img.drawHeight = ih * scale
            except Exception:
                pass

            # caption = ONLY description (no "F01:" prefix)
            caption_txt = escape(desc) if desc else ""

            story.append(KeepTogether([img, Paragraph(caption_txt, cap_style)]))

        sec_no += 1

    # Errors & Warnings  <-- comes AFTER Results
    story.append(section_heading("Errors & Warnings", sec_no, styles))
    if errors:
        story.append(callout(f"Errors ({len(errors)})", errors, "error", styles))
        story.append(Spacer(1, 8))
    if warnings:
        story.append(callout(f"Warnings ({len(warnings)})", warnings, "warning", styles))
    if not errors and not warnings:
        story.append(callout("Status", ["No tagged warnings or errors found."], "info", styles))

    doc.build(story)

# ================================================================
# PDF merge
# ================================================================

def merge_pdfs(old: Path, new: Path, out: Path) -> None:
    try:
        from pypdf import PdfReader, PdfWriter
    except ImportError:
        from PyPDF2 import PdfReader, PdfWriter

    writer = PdfWriter()
    for src in (old, new):
        reader = PdfReader(str(src))
        for page in reader.pages:
            writer.add_page(page)

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wb") as f:
        writer.write(f)

# ================================================================
# main
# ================================================================

def main() -> None:

    print("platform.system():", platform.system())
    print("cwd:", os.getcwd())
    print("__file__:", Path(__file__).resolve())
    print("script_dir:", Path(__file__).resolve().parent)
    
    ap = argparse.ArgumentParser(description="Pipeline PDF report step")
    ap.add_argument("--log", default="log.txt", help="Path to log.txt")
    ap.add_argument("--in_pdf", default=None, help="Existing report to append")
    ap.add_argument("--out", required=True, help="Output report.pdf")
    ap.add_argument("--doc_tsv", default="Method_description.tsv", help="TSV documentation table path")
    args = ap.parse_args()

    workdir = Path(".").resolve()
    params_path = find_single_params_file(workdir)
    title = infer_title_from_params_file(params_path)
    step = infer_step_from_title(title).upper()

    log_path = Path(args.log).resolve()
    if not log_path.exists():
        raise FileNotFoundError(f"log file not found: {log_path}")

    params = load_params_json(params_path)
    warnings, errors = extract_tagged_warnings_errors(log_path)

    # Figures (optional)
    figures_dir = workdir / "report_folder"
    figure_caps_path = workdir / "Figures_Captions.tsv"
    figure_captions = load_figure_captions_tsv(figure_caps_path)

    # NEW: captions (Description) come from report_folder/Figure_log.txt
    figure_log_path = figures_dir / "Figure_log.txt"
    log_captions = load_figure_log_captions(figure_log_path)

    # overwrite description with log captions
    for fig_id, cap in log_captions.items():
        if fig_id not in figure_captions:
            figure_captions[fig_id] = {"title": "", "description": ""}
        figure_captions[fig_id]["description"] = cap
    	
    # load documentation
    doc_path = find_doc_tsv(workdir, args.doc_tsv)
    docs: List[MethodDoc] = load_method_docs(doc_path) if doc_path else []
    docs_for_step = [d for d in docs if (d.step or "").strip().upper() == step]

    # Match methods by Flag_Variable -> params key, and Flag_Value within param value
    method_sections: List[Tuple[str, MethodDoc]] = []
    references_ordered: List[str] = []
    refs_seen = set()
    uniq_seen = set()  # de-dup by (matched_key, flag_val)

    for d in docs_for_step:
        flag_vars_raw = (d.flag_variable or "").strip()
        flag_val = (d.flag_value or "").strip()
        if not flag_vars_raw or not flag_val:
            continue

        flag_vars = [x.strip() for x in flag_vars_raw.split(";") if x.strip()]

        matched_key: Optional[str] = None
        for key in flag_vars:
            if key in params and param_value_contains_flag_value(params[key], flag_val):
                matched_key = key
                break

        if not matched_key:
            continue

        uniq = (matched_key, flag_val)
        if uniq in uniq_seen:
            continue
        uniq_seen.add(uniq)

        display = f"{matched_key} = {flag_val}"
        method_sections.append((display, d))

        # FIX: collect references from matched rows
        for r in normalize_references(d.reference):
            if r not in refs_seen:
                refs_seen.add(r)
                references_ordered.append(r)

    with tempfile.TemporaryDirectory() as td:
        chunk = Path(td) / "module.pdf"

        build_module_pdf(
            title=title,
            params=params,
            warnings=warnings,
            errors=errors,
            method_sections=method_sections,
            method_references=references_ordered,
            figures_dir=figures_dir,
            figure_captions=figure_captions,
            out_pdf=chunk,
        )

        in_pdf = Path(args.in_pdf).resolve() if args.in_pdf else None
        out_pdf = Path(args.out).resolve()

        if in_pdf and in_pdf.exists() and in_pdf.stat().st_size > 0:
            merge_pdfs(in_pdf, chunk, out_pdf)
        else:
            out_pdf.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(chunk, out_pdf)

    print(f"[OK] Report written: {Path(args.out).resolve()}")

if __name__ == "__main__":
    main()
