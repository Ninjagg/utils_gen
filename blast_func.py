# -*- coding: utf-8 -*- 
#===============================================================================
# Data      : 20250312
# Author    : Xuanming
# Annotation: This script is used to blast the sequence of the C region and V region.
#===============================================================================

# ====================================================
# Load packages
# ====================================================
import tempfile
import subprocess
from typing import Dict

import Bio
import operator
from Bio.Blast import NCBIXML as nx
#===============================================================================
#===============================================================================
class Blast:
    """The base class of blast
    """
    def __init__(self, sequence: str):
        self._sequence = sequence
    
    def run_blast(self):
        """Waiting to override in the subclass.
        """
        pass

    @staticmethod
    def get_blast_result(blast_result_filename: str):
        """parse a blast result
        """
        blast_record = {}
        with open(blast_result_filename) as blast_result:
            for b in nx.parse(blast_result):
                lobj = BlastRecordObj(b)
                if not lobj.nohit:
                    blast_record.update(lobj.toDict())
                else:
                    blast_record['nohit'] = True
        return blast_record


class Blastp(Blast):
    """The object was need to be aligned by blastp.
    """
    def __init__(self, sequence: str):
        super().__init__(sequence)

    def run_blast(
            self,
            database_path="./Database/cregion/cregion",
            paras="-outfmt 5 -task blastp-short -num_threads 12"
    ) -> Dict[str, str]:
        with tempfile.NamedTemporaryFile() as input_fasta_file:
            with open(input_fasta_file.name, 'w') as f:
                f.write(f">blastp_query_sequence\n{self._sequence}")
            with tempfile.NamedTemporaryFile(suffix=".xml") as output_file:
                blastp_cmd = f"blastp -query {input_fasta_file.name} -db {database_path} -out {output_file.name} {paras}"
                blastp_run = subprocess.run(blastp_cmd.split(" "), capture_output=True, text=True)
                if blastp_run.returncode == 0:
                    print("\n>>> blastp processed successfully.")
                else:
                    raise Exception(f"Error: blastp processed failed.")
                blastp_result = Blast.get_blast_result(output_file.name)
        return blastp_result


class TBlastn(Blast):
    """The object was need to be aligned by tblasn.
    """
    def __init__(self, sequence: str):
        super().__init__(sequence)
    
    def run_blast(
            self,
            database_path="./Database/vregion/vregion",
            paras="-outfmt 5 -num_threads 12 -evalue 1000"
    ) -> Dict[str, str]:
        with tempfile.NamedTemporaryFile() as input_fasta_file:
            with open(input_fasta_file.name, 'w') as f:
                f.write(f">tblastn_query_sequence\n{self._sequence}")
            with tempfile.NamedTemporaryFile(suffix=".xml") as output_file:
                tblastn_cmd = f"tblastn -query {input_fasta_file.name} -db {database_path} -out {output_file.name} {paras}"
                tblastn_run = subprocess.run(tblastn_cmd.split(" "), capture_output=True, text=True)
                if tblastn_run.returncode == 0:
                    print("\n>>> tblastn processed successfully.")
                else:
                    raise Exception(f"Error: tblastn processed failed.")
                tblastn_result = Blast.get_blast_result(output_file.name)
        return tblastn_result


class BlastRecordObj(object):
    """The object was used to parse records of blast.
    """
    def __init__(self, blast_record):
        if isinstance(blast_record, Bio.Blast.Record.Blast):
            if blast_record.alignments != []:
                self.nohit = False
                self.query_name = blast_record.query 
                score_list = []
                for ind, rec in enumerate(blast_record.alignments):
                    temp_hsp = rec.hsps[0]
                    temp_align_length = float(temp_hsp.align_length)
                    temp_identities = float(temp_hsp.identities)
                    temp_score = float(temp_hsp.score)
                    temp_identities_ratio = float(temp_identities)/float(temp_align_length)
                    result_score = temp_score*temp_identities_ratio*temp_identities_ratio
                    score_list.append([ind, result_score])
                score_list = sorted(score_list, key=operator.itemgetter(1), reverse=True)
                max_ind = score_list[0][0]
                #===============================================================
                local_description = blast_record.descriptions[max_ind]
                self.bits = local_description.bits
                self.e = local_description.e
                self.num_alignments = local_description.num_alignments
                self.score = local_description.score
                self.title = local_description.title
                #===============================================================
                local_alignment = blast_record.alignments[max_ind]
                self.accession = local_alignment.accession
                self.hit_def = local_alignment.hit_def
                self.hit_id = local_alignment.hit_id
                self.length = local_alignment.length
                #===============================================================
                local_hsp = local_alignment.hsps[0]
                self.expect = local_hsp.expect
                self.identities = local_hsp.identities
                self.positives = local_hsp.positives
                self.gaps = local_hsp.gaps
                self.strand = local_hsp.strand
                self.frame = local_hsp.frame
                self.query = local_hsp.query
                self.query_start = local_hsp.query_start
                self.query_end = local_hsp.query_end
                self.match = local_hsp.match
                self.sbjct = local_hsp.sbjct
                self.sbjct_start = local_hsp.sbjct_start
                self.sbjct_end = local_hsp.sbjct_end
                self.align_length = local_hsp.align_length
                self.identities_ratio = float(self.identities) / float(self.align_length)  # the ratio of aligned query_seq (do NOT include '+', '-' special string) and aligned sbjct_seq  e.g. 12 / 31
                self.sequence_similarity = float(self.align_length) / float(self.length)  # the ratio of aligned query_seq (include '+', '-' special string) and the length of sbjct_seq (entire length)  e.g. 31 / 106
            else:
                self.nohit = True
    #===========================================================================
    #===========================================================================
    def toDict(self):
        return {'nohit': self.nohit,
                'query_name': str(self.query_name),
                'bits': str(self.bits),
                'e': str(self.e),
                'num_alignments': str(self.num_alignments),
                'score': str(self.score),
                'title': str(self.title),
                'accession': str(self.accession),
                'hit_def': str(self.hit_def),
                'hit_id': str(self.hit_id),
                'length': str(self.length),
                'expect': str(self.expect),
                'identities': str(self.identities),
                'positives': str(self.positives),
                'gaps': str(self.gaps),
                'strand': str(self.gaps),
                'frame': str(self.frame),
                'query': str(self.query),
                'query_start': str(self.query_start),
                'query_end': str(self.query_end),
                'match': str(self.match),
                'sbjct': str(self.sbjct),
                'sbjct_start': str(self.sbjct_start),
                'sbjct_end': str(self.sbjct_end),
                'align_length': str(self.align_length),}


if __name__ == '__main__':
    # test for partial blast
    a = Blastp("EDQVTQSPEALRLQEGESSSLNCSYTSRMLRGLFWYRQDPGKGPEFLFTLYSAGEEKEKERLKATLTKKESFLHITAPKPEDSATYLCAVQADSWPSYALNFGKGTSLLVTPYIQNPDPAVYQLRDSKSSDKKVCLFTDFDSQTNVSQSKDSDVYITDKCVLDDPSEDFKSNSAVAWSNKPDFACANAFNNSIIPEDTFFPSPESSC").run_blast()

    # test for blastp
    a = Blastp("MGWTLVFLFLLSVTAGVHSQEQLEESGGDLVKPGASLTLTCTASGFSFSSSYWICWVRQAPGKGLEWIACIYAGGIGSTYYASWAKGRFTISKTSSTTVTLQMTSLTAADTATYFCAGDLPGGGYYTLTRLDLWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK").run_blast()  # a sequence with signal peptide + V region + C region
    print(a)

    # test for tblastn
    b = TBlastn("MGWTLVFLFLLSVTAGVHSQEQLEESGGDLVKPGASLTLTCTASGFSFSSSYWICWVRQAPGKGLEWIACIYAGGIGSTYYASWAKGRFTISKTSSTTVTLQMTSLTAADTATYFCAGDLPGGGYYTLTRLDLWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK").run_blast()
    