import probe_design
import os

RESULT_PATH = os.path.join(os.getcwd(), "probe_design_files")

gene_id = "NM_182993.2"
gene_name = "slc17a7"
gene_synonym = []
hairpin_id = 2

prb_length = 20
gc_range = [40, 60]
hits_thresh = 5
num_offtar_thresh = 0
prb_space = 1
dg_thresh = -9

probe_design.designProbes(gene_id=gene_id, gene_name=gene_name, 
                        hairpin_id=hairpin_id, 
                        db=os.path.join(os.getcwd(), "db/mouse/mouse_refseq_rna"),
                        result_path=RESULT_PATH,
                        prb_length=prb_length,
                        gc_range=gc_range,
                        hits_thresh=hits_thresh,
                        num_offtar_thresh=num_offtar_thresh,
                        prb_space=prb_space,
                        dg_thresh=dg_thresh,
                        to_excel=True)