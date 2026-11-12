"""Add the Kv4 thalamocortical-loop write-up to the lab database's documents table as a
PROVISIONAL summary document (paulp, 2026-11-09). Same pattern as add_kv4_dendrite_document.py (LabDB API
from the HHModels repo). Run on the Windows machine holding lab_data.db, from this folder, with the hh env:

    D:\\anaconda3\\envs\\hh\\python.exe add_kv4_thalamus_document.py

(or double-click run_add_kv4_thalamus_doc.bat). Idempotent: an existing document with the same title gets a
new version instead of a duplicate. Links use the LabDB schema {"label", "target"}.
"""
import os, sys

HH_DIR = os.environ.get("HHMODELS_DIR", r"D:\Dropbox\Programming\githome\HHModels\Claude_code_gpu_dev")
sys.path.insert(0, HH_DIR)
import lab_db  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MD = os.path.join(HERE, "KV4_THALAMUS_LOOP.md")

TITLE = "The Kv4 complex across the thalamocortical loop: which thalamic neurons depend on Kv4.2 and with which auxiliary subunits (ABC Atlas 10x / THALMANAC + ThalamoSeq full-length variants)"
SUMMARY = ("ABC Atlas 10x thalamus (183,553 neurons, THALMANAC annotations), calibrated for gene length against the ThalamoSeq full-length relay-neuron TPM. "
           "The reticular nucleus is the most Kv4.2-dependent thalamic population: Kcnd3 is 10-30x lower than in any relay group (calibrated Kv4.2 share 0.71, Kv4.3 0.17), "
           "with the cortical-interneuron auxiliary complement (KChIP1 76% + KChIP2 15%, DPP10 81%, subunits in large excess). The Spp1+ core TRN cells of "
           "Hartley 2024 are the extreme case (Kcnd3 0.5% of Kcnd2, KChIP1-rich, KChIP4-poor, strongest T-type expression). Relay neurons co-express "
           "Kcnd2/Kcnd3 with a first-order (Kv4.3-balanced; LGd 0.31, VPM/VPL 0.30) to higher-order/midline (Kcnd2-dominant; LP 0.70, PVT 0.69) gradient, "
           "partners KChIP3/KChIP4 + DPP6. Viaat-Cre covers 93% of RT cells and ~1% of relay cells. TRN is the natural substrate of the Viaat-Cre spike-wave "
           "phenotype: a Kv4.2-only A-current with no Kv4.3 buffer. Relay-neuron N-terminal variants from ThalamoSeq (120 full-length pools, 22 nuclei) "
           "are being quantified; TRN variants need new full-length libraries. Proposes EEG in Vglut2-cKI, Spp1-Cre/Ecel1-Cre x cKI, sorted TRN SMART-seq, "
           "TRN A-current recordings, and KChIP1/DPP10 structural models.")
TAGS = ["Kv4.2", "Kcnd2", "Kcnd3", "thalamus", "reticular nucleus", "TRN", "thalamocortical loop", "spike-wave discharge", "KChIP1", "KChIP2", "DPP10", "KChIP3", "KChIP4", "DPP6", "ABC Atlas", "THALMANAC", "ThalamoSeq", "Spp1", "Ecel1", "Viaat-Cre", "V404M", "provisional"]
LINKS = [
    {"label": "Code, tables and figures (GitHub)", "target": "https://github.com/ppfaff/Kv4-isoform-survey"},
    {"label": "Figure: th10x_kv4_calibrated.png", "target": "https://github.com/ppfaff/Kv4-isoform-survey/blob/main/figures/th10x_kv4_calibrated.png"},
    {"label": "Thalamus README (methods, run steps, tier-2 plan)", "target": "https://github.com/ppfaff/Kv4-isoform-survey/blob/main/thalamus/README.md"},
    {"label": "Companion: V404M hippocampal activation / conditional dissociation write-up", "target": "https://github.com/ppfaff/Kv4-isoform-survey/blob/main/lab_db/KV4_V404M_DENTATE_PROGRAM.md"},
    {"label": "THALMANAC (Allen Institute for Neural Dynamics)", "target": "https://github.com/AllenNeuralDynamics/thalmanac"},
]
SESSIONS = []


def body_from_file():
    txt = open(MD, encoding="utf-8").read()
    return ("_Provisional summary; code, tables and figures are in the GitHub repository linked below. "
            "Figures and statistics as of 2026-11-10; ThalamoSeq relay-neuron variant section added 2026-11-12._\n\n" + txt.strip() + "\n")


def main():
    db = lab_db.LabDB(archive_path=r"D:\physiology\abf_files", mode="production", analyst="claude")
    body = body_from_file()
    docs = db.list_documents()
    same = docs[docs["title"] == TITLE] if len(docs) else docs
    if len(same):
        doc_id = int(same.iloc[0]["id"])
        ver = db.update_document(doc_id, change_note="revised thalamus loop write-up (2026-11-10)",
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
