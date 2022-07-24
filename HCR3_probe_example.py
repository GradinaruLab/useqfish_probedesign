import probe_design
import os

RESULT_PATH = os.path.join(os.getcwd(), "HCR3_probe_design_files")
if not os.path.isdir(RESULT_PATH):
    os.mkdir(RESULT_PATH)

# inputs
email = "mjjang@caltech.edu"
gene_id = "NM_001317043.1"
gene_name = "vegfa"
gene_synonym = []
hairpin_id = 5
sequence=""

# other parameters
prb_length = 20
gc_range = [40, 60]
prb_space = 5
dg_thresh = -9

# run probe design
resultdf = probe_design.designHCR3Probes(gene_id=gene_id, gene_name=gene_name, 
                        hairpin_id=hairpin_id, 
                        db=os.path.join("db/rat/mratbn7.2"),
                        result_path=RESULT_PATH,
                        prb_length=prb_length,
                        gc_range=gc_range,
                        prb_space=prb_space,
                        dg_thresh=dg_thresh,
                        to_excel=True,
                        email=email)