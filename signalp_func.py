#!/home/liuxuanming/miniconda3/envs/ProSeqAnnotation/bin/python
# -*- coding: utf-8 -*- 
#===============================================================================
# Data      : 20250312
# Author    : Xuanming
# Annotation: This script is used to detect the signal peptide of the sequence.
#===============================================================================
import os
from typing import Dict
import tempfile
import subprocess
#===============================================================================
#===============================================================================
class Signalp:

    def __init__(self, sequence: str):
        self._sequence = sequence
    
    def run_signalp(self) -> Dict[str, str]:
        """run the signalp for searching signal peptide.
        """
        with tempfile.TemporaryDirectory() as signalp_tmp_dir:
            with tempfile.NamedTemporaryFile() as input_fasta_file:
                with open(input_fasta_file.name, 'w') as f:
                    f.write(f">signalp_query_sequence\n{self._sequence}")
                signalp_cmd = f"signalp -fasta {input_fasta_file.name} -tmp {signalp_tmp_dir}"
                # The application of signalp5 will generate the output file named {input_file_name}_summary.signalp5. Cannot substitute the output with  the temporary file.
                signalp_run = subprocess.run(signalp_cmd.split(" "), capture_output=True, text=True)
                if signalp_run.returncode == 0:
                    signalp_result_filename = os.path.join(os.getcwd(), f"{input_fasta_file.name.split('/')[-1]}_summary.signalp5")
                    print(f"\n>>> signalp processed successfully.")
                else:
                    raise Exception(f"Error: signalp processed failed.")
        signalp_result = Signalp.get_signalp_result(signalp_result_filename)
        return signalp_result
    
    @staticmethod
    def get_signalp_result(signalp_result_filename: str) -> Dict[str, str]:
        """parse a signalp result
        """
        with open(signalp_result_filename, 'r') as f:
            content = f.readlines()
        if len(content) > 2:  # a normal result
            sigp_record = content[2].strip().split('\t')
            if len(sigp_record) >= 5:
                sigp_flag = 'Y'  # the quality of signalp result
                sigp_idx = int(sigp_record[4].split()[2].split('.')[0].split('-')[0])  # extract the end index os signalp region
            else:
                sigp_flag = 'N'
                sigp_idx = 'None'
        signalp_record = {
            'signalp_flag': sigp_flag, 
            'signalp_idx': sigp_idx
        }
        os.system(f"rm {signalp_result_filename}")
        return signalp_record

# test
if __name__ == '__main__':
    signalp_record = Signalp("DVQLVQSGAEVKKPGASVKVSCKASGYTFTRYTMHWVRQAPGQGLEWIGYINPSRGYTNYADSVKGRFTITTDKSTSTAYMELSSLRSEDTATYYCARYYDDHYCLDYWGQGTTVTVSSGGGGSDIVLTQSPATLSLSPGERATLSCRASQSVSYMNWYQQKPGKAPKRWIYDTSKVASGVPARFSGSGSGTDYSLTINSLEAEDAATYYCQQWSSNPLTFGGGTKVEIKGGGGSRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC").run_signalp()  # no signalp
    print(signalp_record)

    signalp_record = Signalp("MGWTLVFLFLLSVTAGVHSQEQLEESGGDLVKPGASLTLTCTASGFSFSSSYWICWVRQAPGKGLEWIACIYAGGIGSTYYASWAKGRFTISKTSSTTVTLQMTSLTAADTATYFCAGDLPGGGYYTLTRLDLWGQGTLVTVSS").run_signalp()  # signalp

    print(signalp_record)