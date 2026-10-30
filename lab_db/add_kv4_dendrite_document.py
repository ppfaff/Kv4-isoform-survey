"""Add the Kv4 dendrite-vs-soma compartment write-up to the lab database's documents table as a
PROVISIONAL summary document (paulp, 2026-10-30). Mirrors analysis/add_key_findings_document.py
in the HHModels repo and uses the same LabDB API.

Run once on the Windows machine that holds the lab database, from this folder, with the hh env:

    D:\\anaconda3\\envs\\hh\\python.exe add_kv4_dendrite_document.py

(or double-click run_add_kv4_dendrite_doc.bat, which does that and shows the log).

Idempotent: if a document with the same title already exists, it saves a new VERSION of it
(update_document) instead of adding a duplicate. Only run it when the markdown has changed: a re-run with
unchanged content still creates a version carrying this script's fixed change_note. Link entries use the
LabDB schema {"label", "target"} (lab_forms._doc_render_links reads "target"). The body is KV4_DENDRITE_COMPARTMENT.md from
this folder; figures are linked to the GitHub repo.
"""
import os, sys

HH_DIR = os.environ.get("HHMODELS_DIR", r"D:\Dropbox\Programming\githome\HHModels\Claude_code_gpu_dev")
sys.path.insert(0, HH_DIR)
import lab_db  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MD = os.path.join(HERE, "KV4_DENDRITE_COMPARTMENT.md")

TITLE = "Kv4 complex mRNAs are somatic: dendrite-vs-soma compartment analysis (Farris 2019 laminar, Middleton 2019 single-neuron)"
SUMMARY = ("Tests whether Kv4-complex mRNAs (Kcnd2, Kcnip4, Kcnip2, Dpp6 and the other five genes) are transported "
           "into dendrites, using laser-microdissected cell-body vs dendritic layers of adult mouse CA1/CA2/CA3/DG "
           "(Farris 2019, 24 libraries) and 16 paired soma/dendrite single cultured neurons (Middleton 2019). "
           "Kcnd2, Kcnip4, Dpp6 and Kcnip2 are depleted 2-4x in dendritic layers while Camk2a/Map2/Shank1 are "
           "enriched 3-9x, and the neuropil Kcnd2 signal is 52% intronic (nuclear pre-mRNA from resident "
           "interneurons/glia). Conclusion: Kv4.2 mRNA is somatic in vivo; the survey's somatic Kv4-null cells "
           "(DG 18%) have low Kv4 mRNA, not dendritic mRNA, and show no Kcnd1/Kcnd3 compensation - supersedes the "
           "dendritic-localization reading in the Oct 29 document. Two variants survive cell-type-mixing controls "
           "at the splice-junction level: Kv4.3L exon inclusion 0.17 -> 0.58 (p=0.001) and the Kcnip3 cys-rich "
           "acceptor 0.19 -> 0.04 (p=0.002; neuropil Kcnip3 is ~96% KChIP3a). Both are rated tentative; concerns "
           "(laminar interneurons, nuclear retention, single in vivo dataset) and eight follow-up analyses are listed. "
           "Also records that salmon underestimates Kv4.3L ~4-fold because only the Kv4.3S model carries the 3'UTR. "
           "Read-level follow-up on the 446 Kv4-null DG cells: baseline Kcnd2 transcription present, burst state nearly absent (1% vs 45%), "
           "no mature mRNA at matched pre-mRNA; null cells are the low-IEG, Calb1-low quiet granule cells. Lab dFosB staining shows the "
           "V404M dentate granule layer is chronically hyperactive; going-forward plan: dFosB-positive fraction, re-analysis of existing "
           "littermate bulk RNA-seq (variants, allele-specific Kcnd2, exon:intron, dFosB:FosB), then smFISH.")
TAGS = ["Kv4.2", "Kcnd2", "deltaFosB", "dentate gyrus hyperactivity", "transcriptional bursting", "dendritic mRNA", "mRNA localization", "neuropil", "KChIP3", "DREAM", "Kv4.3L",
        "splice variants", "laminar RNA-seq", "LCM", "hippocampus", "dentate gyrus", "V404M", "provisional"]
LINKS = [
    {"label": "Code, tables and figures (GitHub)",
     "target": "https://github.com/ppfaff/Kv4-isoform-survey"},
    {"label": "Companion document: isoform survey conclusions (Claude doc)",
     "target": "https://claude.ai/artifact/BdsYCYmxtyRcaqM968Bj9Z"},
    {"label": "Farris et al. 2019 data (GSE116342)",
     "target": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116342"},
    {"label": "Middleton et al. 2019 data (GSE115480)",
     "target": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115480"},
    {"label": "Middleton et al. 2019, BMC Biology",
     "target": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6344992"},
]
SESSIONS = []   # no electrophysiology recording sessions are associated with this document


def body_from_file():
    txt = open(MD, encoding="utf-8").read()
    return ("_Provisional summary of the dendrite-vs-soma compartment analysis; code, tables and figures "
            "are in the GitHub repository linked below. Supersedes the dendritic-localization reading in the "
            "Oct 29 survey document. Figures and statistics are as of 2026-10-30._\n\n"
            + txt.strip() + "\n")


def main():
    db = lab_db.LabDB(archive_path=r"D:\physiology\abf_files", mode="production", analyst="claude")
    body = body_from_file()
    docs = db.list_documents()
    same = docs[docs["title"] == TITLE] if len(docs) else docs
    if len(same):
        doc_id = int(same.iloc[0]["id"])
        ver = db.update_document(doc_id, change_note="added activity-state result, V404M dFosB observation and the going-forward plan (bulk littermate RNA-seq re-analysis, dFosB counts, smFISH)",
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
