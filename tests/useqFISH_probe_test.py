import unittest
import os
import probe_design
from api import ProbeDesign as PD

class TestProbeDesigner(unittest.TestCase):

    def test_HCR3_probes(self):
        RESULT_PATH = os.path.join(os.getcwd(), "useqFISH_probe_design_files")
        if not os.path.isdir(RESULT_PATH):
            os.mkdir(RESULT_PATH)

        email = "mjjang@caltech.edu"
        gene_id = "NM_009891.2"
        gene_name = "chat"
        gene_host = "mus musculus"
        sequence = ""

        barcode_num = 44

        ## parameters
        prb_length = 20
        gc_range = [40, 60]
        prb_space = 1
        dg_thresh = -9

        primer_end = "TAATGTTATCTT"
        padlock_start = "ACATTA"
        padlock_end = "AAGATA"
        spacer1 = "attta"
        spacer2 = "atta"

        # run probe design
        modified_code_df = probe_design.designuseqFISHProbes(gene_id=gene_id, 
                                gene_name=gene_name,
                                gene_host=gene_host,
                                email=email,
                                sequence=sequence,
                                barcode_num=barcode_num,
                                barcode_path="db/mouse/barcodes/barcodes.xlsx",
                                db=os.path.join(os.getcwd(), "db/mouse/mouse_refseq_rna"),
                                result_path=RESULT_PATH,
                                prb_length=prb_length,
                                gc_range=gc_range,
                                prb_space=prb_space,
                                dg_thresh=dg_thresh,
                                primer_end=primer_end,
                                padlock_start=padlock_start,
                                padlock_end=padlock_end,
                                spacer1=spacer1,
                                spacer2=spacer2,
                                to_excel=False)
        
        # run probe design with original code
        original_code_df = PD.designuseqFISHProbes()

        self.assertTrue(original_code_df.equals(modified_code_df))

if __name__ == '__main__':
    unittest.main()