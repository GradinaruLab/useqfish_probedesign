import unittest
import os
import probe_design
from api import ProbeDesignHCR3 as PDHCR3

class TestProbeDesigner(unittest.TestCase):

    def test_HCR3_probes(self):
        RESULT_PATH = os.path.join(os.getcwd(), "HCR3_probe_design_files")
        if not os.path.isdir(RESULT_PATH):
            os.mkdir(RESULT_PATH)

        # inputs
        gene_id = "NM_182993.2"
        gene_name = "slc17a7"
        gene_synonym = []
        hairpin_id = 2

        # other parameters
        prb_length = 20
        gc_range = [40, 60]
        prb_space = 1
        dg_thresh = -9

        # run probe design
        modified_code_df = probe_design.designHCR3Probes(gene_id=gene_id, gene_name=gene_name, 
                                hairpin_id=hairpin_id, 
                                db=os.path.join(os.getcwd(), "db/mouse/mouse_refseq_rna"),
                                result_path=RESULT_PATH,
                                prb_length=prb_length,
                                gc_range=gc_range,
                                prb_space=prb_space,
                                dg_thresh=dg_thresh,
                                to_excel=False)
        
        # run probe design with original code
        original_code_df = PDHCR3.designHCR3Probes()

        self.assertTrue(original_code_df.equals(modified_code_df))

if __name__ == '__main__':
    unittest.main()