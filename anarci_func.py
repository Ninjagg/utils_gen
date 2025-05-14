# -*- coding: utf-8 -*-
# ===============================================================================
# Data      : 20250312
# Author    : Xuanming
# Annotation: This script is used to search the V region of antibodies.
# ===============================================================================

# ====================================================
# Load packages
# ====================================================
import os
import tempfile
import subprocess
from typing import Dict, List

# ===============================================================================
# ===============================================================================
class Anarci:
    """The class for runing the ANARCI for searching the V region.
    If there's no v region found, the return value will be `[]`.
    If there are v regions found, the return value will be `List[Dict[str, str]]`
    """

    def __init__(self, sequence: str):
        self._sequence = sequence

    def run_anarci(self) -> List[Dict[str, str]]:
        """Run the ANARCI for searching the V region. 
        """
        with tempfile.NamedTemporaryFile() as input_fasta_file:
            with open(input_fasta_file.name, 'w') as f:
                f.write(f">anarci_query_sequence\n{self._sequence}")
            output_file = f"{input_fasta_file.name}_anarci_result.txt"
            anarci_cmd = f"ANARCI -s k -i {input_fasta_file.name} -o {output_file}"
            anarci_run = subprocess.run(anarci_cmd.split(" "), capture_output=True, text=True)
            if anarci_run.returncode == 0:
                print(f"\n>>> ANARCI processed successfully.")
            else:
                raise Exception(f"Error: ANARCI processed failed.")
        anarci_result = Anarci.get_anarci_result(output_file)
        return anarci_result

    @staticmethod
    def get_anarci_result(output_file) -> List[Dict[str, str]]:
        anarci_records_list = []
        with open(output_file, "r") as f:
            content = f.readlines()
        if len(content) == 2:
            pass
        else:
            num_v_region = int(content[2].strip().split()[-1])
            # =======================================================================
            for i in range (len(content)-1):
                line = content[i]
                if line.strip() == '# Most significant HMM hit':
                    v_kabat_sequence = []
                    v_kabat_idx = []
                    v_region_tmp = {}
                    v_region_infor = content[i+2].strip().split('|')
                    species = v_region_infor[1]
                    chain_type = v_region_infor[2]
                    v_start_idx = int(v_region_infor[5])
                    v_end_idx = int(v_region_infor[6])+1
                if line.strip().split()[0] != '#' and\
                line.strip().split('|')[0] != '#' and\
                line.strip().split()[0] != '//':
                    v_kabat_sequence.append(line.strip().split()[-1])
                    if len(line.strip().split()) == 3:
                        v_kabat_idx.append(line.strip().split()[1])
                    if len(line.strip().split()) == 4:
                        v_kabat_idx.append(line.strip().split()[1] + line.strip().split()[2])
                if content[i].strip().split()[0] != '#' and\
                content[i].strip().split('|')[0] != '#' and\
                (content[i+1].strip().split()[0] == '#' or content[i+1].strip().split()[0] == '//'):
                    v_region_tmp = {'num_v_region': num_v_region,
                                    'species' : species,
                                    'chain_type' : chain_type,
                                    'v_start_idx': v_start_idx,
                                    'v_end_idx': v_end_idx,
                                    'v_kabat_sequence': v_kabat_sequence,
                                    'v_kabat_idx' : v_kabat_idx}
                    anarci_records_list.append(v_region_tmp)
        # ===========================================================================
        os.system(f"rm {output_file}")
        return anarci_records_list

# test
if __name__ == "__main__":
    a = Anarci(
        "MGWTLVFLFLLSVTAGVHSQEQLEESGGDLVKPGASLTLTCTASGFSFSSSYWICWVRQAPGKGLEWIACIYAGGIGSTYYASWAKGRFTISKTSSTTVTLQMTSLTAADTATYFCAGDLPGGGYYTLTRLDLWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
    ).run_anarci()
    print(a)