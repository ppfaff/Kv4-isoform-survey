"""Add the V404M dentate-gyrus activity-program write-up to the lab database's documents table as a
PROVISIONAL summary document (paulp, 2026-11-09). Same pattern as add_kv4_dendrite_document.py (LabDB API
from the HHModels repo). Run on the Windows machine holding lab_data.db, from this folder, with the hh env:

    D:\\anaconda3\\envs\\hh\\python.exe add_kv4_dentate_document.py

(or double-click run_add_kv4_dentate_doc.bat). Idempotent: an existing document with the same title gets a
new version instead of a duplicate. Links use the LabDB schema {"label", "target"}.
"""
import os, sys

HH_DIR = os.environ.get("HHMODELS_DIR", r"D:\Dropbox\Programming\githome\HHModels\Claude_code_gpu_dev")
sys.path.insert(0, HH_DIR)
import lab_db  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MD = os.path.join(HERE, "KV4_V404M_DENTATE_PROGRAM.md")

TITLE = "V404M/+ hippocampus: chronic activation of the Kv4.2-dependent principal cells (DG and CA1) - bulk RNA-seq program, dFosB/c-Fos histology, Viaat-Cre conditional"
SUMMARY = ("In the 3-month littermate bulk RNA-seq, V404M/+ hippocampus shows the transcriptional signature of chronically "
           "hyperactive dentate granule cells: Penk up (+1.0 to +1.2 log2), Pdyn down (-0.9), Npy up, Drd1 up, Ecel1 up, Cck down, "
           "Igfbpl1 and Rprm up, acute IEGs flat or lower; Penk, Npy and Drd1 separate all 6 mutants from all 7 wild types "
           "(exact p 0.001). Hippocampus-specific, robust to tissue-composition covariates, and the transcriptional counterpart of "
           "the dFosB immunohistochemistry (DG 2.4x, CA1 2.7x, CA2/CA3 flat, retrosplenial cortex 1.6x) and the qPCR c-Fos decrease (cortex ~0.3x, hippocampus 0.7x). The field pattern tracks Kv4.2 expression by cell type (DG, CA1 high; CA3 low), favouring cell-autonomous activation over propagation from the dentate. Kcnip3 (DREAM) is up with Pdyn down, consistent with DREAM repression of the Pdyn DRE. "
           "Viaat-Cre x conditional V404M (GABAergic neurons only) shows frontal spike-wave discharges without dentate dFosB, so the "
           "dentate activation requires V404M in glutamatergic cells - most simply a cell-autonomous loss of the granule-cell A-current gate. "
           "Records the Nov 9 finding that the bulk blocks carried genotype-correlated contamination (olfactory bulb in cortex, choroid plexus "
           "in hippocampus) that withdraws the cortical astrogliosis and leaves the hippocampal astrogliosis unconfirmed. Proposes granule-cell "
           "(Prox1/Dock10-Cre) and Emx1-Cre conditionals, RNAscope cell assignment, granule-cell A-current recordings and paired dentate/frontal EEG.")
TAGS = ["V404M", "Kv4.2", "Kcnd2", "dentate gyrus", "granule cells", "deltaFosB", "Penk", "Pdyn", "Npy", "Drd1", "Kcnip3", "DREAM",
        "spike-wave discharge", "Viaat-Cre", "conditional", "bulk RNA-seq", "tissue composition", "hippocampus", "provisional"]
LINKS = [
    {"label": "Code, tables and figures (GitHub)", "target": "https://github.com/ppfaff/Kv4-isoform-survey"},
    {"label": "Figure: v404m_dg_program.png", "target": "https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/v404m_dg_program.png"},
    {"label": "Figure: v404m_admixture.png", "target": "https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/v404m_admixture.png"},
    {"label": "Companion: V404M bulk RNA-seq and heterogeneity write-up (revised Nov 9)", "target": "https://github.com/ppfaff/Kv4-isoform-survey/blob/main/lab_db/KV4_V404M_BULK_AND_HETEROGENEITY.md"},
    {"label": "Companion: isoform survey conclusions (Claude doc)", "target": "https://claude.ai/artifact/BdsYCYmxtyRcaqM968Bj9Z"},
]
SESSIONS = []


def body_from_file():
    txt = open(MD, encoding="utf-8").read()
    return ("_Provisional summary; code, tables and figures are in the GitHub repository linked below. "
            "Figures and statistics as of 2026-11-09._\n\n" + txt.strip() + "\n")


def main():
    db = lab_db.LabDB(archive_path=r"D:\physiology\abf_files", mode="production", analyst="claude")
    body = body_from_file()
    docs = db.list_documents()
    same = docs[docs["title"] == TITLE] if len(docs) else docs
    if len(same):
        doc_id = int(same.iloc[0]["id"])
        ver = db.update_document(doc_id, change_note="revised dentate program write-up (2026-11-09)",
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
