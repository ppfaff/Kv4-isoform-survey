"""Add the Kv4 isoform-survey conclusions write-up to the lab database's documents table as a
PROVISIONAL summary document (paulp, 2026-10-29). Mirrors analysis/add_key_findings_document.py
in the HHModels repo and uses the same LabDB API.

Run once on the Windows machine that holds the lab database, from this folder, with the hh env:

    D:\\anaconda3\\envs\\hh\\python.exe add_kv4_survey_document.py

(or double-click run_add_kv4_survey_doc.bat, which does that and shows the log).

Idempotent: if a document with the same title already exists, it saves a new VERSION of it
(update_document) instead of adding a duplicate. Only run it when the markdown has changed: a re-run with
unchanged content still creates a version carrying this script's fixed change_note. Link entries use the
LabDB schema {"label", "target"} (lab_forms._doc_render_links reads "target"). The body is KV4_SURVEY_CONCLUSIONS.md from
this folder (the markdown export of the Claude doc); figures are linked to the GitHub repo.
"""
import os, sys

HH_DIR = os.environ.get("HHMODELS_DIR", r"D:\Dropbox\Programming\githome\HHModels\Claude_code_gpu_dev")
sys.path.insert(0, HH_DIR)
import lab_db  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MD = os.path.join(HERE, "KV4_SURVEY_CONCLUSIONS.md")

TITLE = "Kv4 auxiliary-subunit isoform survey: conclusions and reliability (V1 + hippocampus scRNA-seq)"
SUMMARY = ("Survey of Kv4 (Kcnd1-3), KChIP (Kcnip1-4) and DPP (Dpp6/Dpp10) splice-variant usage across "
           "~11,300 single neurons: Tasic 2016 V1 (1,677 cells), Allen SSv4 hippocampal formation "
           "(8,131 cells) and an Allen SSv4 VISp chemistry control (1,511 cells). Each conclusion is "
           "rated reliable / tentative / flagged. Key results: interneurons vs principal cells use "
           "disjoint KChIP programs; somatic Kv4-null cells (DG 18%, Vip 27%, Pvalb 20%) still "
           "transcribe Kcnd2 (post-transcriptional; REVISED 2026-10-30: not dendritic localization, see compartment document); hippocampus-"
           "specific Sst and Lamp5 Lhx6 types have the KChIP genes off; near-projecting pyramidal "
           "cells carry a ~92% Kv4.2 pore with an interneuron-style KChIP1 complement and are the "
           "pyramidal class most exposed to Kv4.2-V404M; Vip/Sncg/Lamp5 interneurons (Kv4.3-based) "
           "are least exposed. Read-chemistry effects between the two V1 datasets are quantified. "
           "Includes proposed experiments (Kcnd2 exon/intron smFISH, NP-cell validation, V404M "
           "physiology with internal controls).")
TAGS = ["Kv4.2", "V404M", "KChIP", "DPP6", "DPP10", "splice variants", "scRNA-seq", "SMART-seq",
        "hippocampus", "V1", "near-projecting", "dentate gyrus", "provisional"]
LINKS = [
    {"label": "Live write-up (Claude doc, with embedded figures)",
     "target": "https://claude.ai/artifact/BdsYCYmxtyRcaqM968Bj9Z"},
    {"label": "Code, tables and figures (GitHub)",
     "target": "https://github.com/ppfaff/Kv4-isoform-survey"},
    {"label": "Tasic 2016 data (GSE71585)",
     "target": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585"},
    {"label": "Allen SSv4 mouse cortex + hippocampus (Yao 2021)",
     "target": "https://brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq"},
]
SESSIONS = []   # no electrophysiology recording sessions are associated with this document


def body_from_file():
    txt = open(MD, encoding="utf-8").read()
    return ("_Provisional summary generated from the single-cell isoform survey; the live version "
            "with embedded figures is the Claude doc linked below, and the code, tables and figures "
            "are in the GitHub repository. Figures and statistics are as of 2026-10-29._\n\n"
            + txt.strip() + "\n")


def main():
    db = lab_db.LabDB(archive_path=r"D:\physiology\abf_files", mode="production", analyst="claude")
    body = body_from_file()
    docs = db.list_documents()
    same = docs[docs["title"] == TITLE] if len(docs) else docs
    if len(same):
        doc_id = int(same.iloc[0]["id"])
        ver = db.update_document(doc_id, change_note="Oct 30 revision note added: Kv4-null cells have low Kv4 mRNA, not dendritic mRNA (see compartment-analysis document)",
                                 body=body, summary=SUMMARY, tags=TAGS, links=LINKS,
                                 sessions=SESSIONS, status="provisional", doc_type="summary")
        print("updated document %d -> version %d" % (doc_id, ver))
    else:
        doc_id = db.create_document(TITLE, body=body, summary=SUMMARY, status="provisional",
                                    doc_type="summary", tags=TAGS, links=LINKS, sessions=SESSIONS)
        print("created document %d (provisional)" % doc_id)
    d = db.get_document(doc_id)
    print("title:", d["title"])
    print("status:", d["status"], "| type:", d["doc_type"], "| body chars:", len(d["body"] or ""))


if __name__ == "__main__":
    main()
