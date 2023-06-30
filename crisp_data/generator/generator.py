import json
import os.path
import random
import time
import CrispTestDataGenerator_3.crisp_data.generator.sequences_generator as sg
import CrispTestDataGenerator_3.crisp_data.generator.read_quality as rq


class Generator:

    """
    The class to generate sequences.

    The generation happens based on the content of the config.json, which contains the array of required
    sequences and their desired properties.
    """

    def __init__(self):
        self.to_parse = []
        with open("config.json") as config:
            self.json_data = json.loads(config.read())
            self.to_parse = self.json_data["mutationsArray"]
            self.output_folder = self.json_data["outputFolder"]
            self.time_format = self.json_data["timeFormat"]

    @staticmethod
    def join_dict(map_to_join: dict) -> str:
        """

        :param map_to_join: is the map where keys are mutation types,
        and the values are lists with the letters mutated during the CRISP and their positions
        :return: the string denoting the type of expected mutation resulting from the CRISP:
        *Example:*

        :S:C-55&C-60&A-73
        where S denotes the type of mutation
        C-55 is the resulting letter C on the position 55
        C-60 is the resulting letter C on the position 60
        A-73 is the resulting letter A on the position 73
        """
        res = ""
        for item in map_to_join.keys():
            res += ":" + item + ":" + "&".join(["-".join(sub_array) for sub_array in map_to_join[item]])
        return res

    def create_file(self, sequence_description: dict) -> None:
        """

        :param sequence_description: the JSON dictionary from the array stored in the mutationsArray in the config.json
        This JSON dictionary contains the following information:
        - initLength -- the initial length of the sequence which the CRISP should modify
        - isCoding -- the flag specifying if the sequence is the coding one (true) or not (false)
        - mutationSite -- the array with mutation sites
        - mutationLength -- the array with mutation lengths (could be negative for deletions, zero for SNPs,
          and positive for insertions)
        - direction -- the array with mutation directions (true denotes that from the corresponding position
          mutation is applied in the right direction, false denotes that the mutation is applied in the
          left directio)
        - quality -- the percent of mutated sequences (must be 5-fold)
        - count -- the number of reads
        - sizeScale -- the range in percents of the initial sequence length, where to bottom value means the minimal
          possible length of the resulting sequence, and the top value means the maximal possible length of the resulting
          sequence
        - scaleDirection -- the value which takes -1, 0, 1 values, where
          -1 means that the sequence length modification happens to the left side from the modification site
          0 means that the sequence length modification does not happen
          1 means that the sequence length modification happens to the right side from the modification site
        - reader -- the sequence reader name (MiSeq, MiniSeq, HiSeq2500, HiSeqXTen, NextSeq500, NextSeq550, NovaSeq6000)
        :return: none

        The method creates four files as output:
        - test_sequences.fasta -- the set of sequences which are produced by the reader as the result of CRISP outcome
          analysis
        - fastq_output.fastq -- the set of sequences along with their read quality
        - expected_analysis.csv -- what the clustering pipeline should produce as the outcome of CRISP results reading
          analysis
        - length_scale_report.csv -- the read sequences lengths data

        The above mentioned files are time stamped as per the timeFormat field value from the config.json.
        The formatter string follows the rules from the strftime method of the time.py module.
        """
        s = os.path.sep
        s = s if s == "/" else "\\"
        time_stamp = time.strftime(self.time_format, time.localtime()) + "_" \
                     + str(random.randint(a=100, b=999))
        file_name = self.output_folder + s + "test_sequences.fasta"
        fasta_q_name = self.output_folder + s + "fastq_output.fastq"
        csv_output_name = self.output_folder + s + time_stamp + "_expected_analysis.csv"
        csv_report_name = self.output_folder + s + time_stamp + "_length_scale_report.csv"
        seq_gen = sg.SequencesGenerator(sequence_description)
        init_seq, ref_seq, mutation, m_map = seq_gen.generate_coding()
        corr_seq, orig_arr, c_orig, report = seq_gen.corrupt_generated(init_seq, ref_seq)
        res_array = [init_seq,
                     ref_seq,
                     seq_gen.reverse_complement(ref_seq),
                     self.join_dict(m_map),
                     "-".join([str(x) for x in sequence_description["mutationSite"]]),
                     str(c_orig)]
        with open(file_name, "w") as raw_data:
            raw_data.write(">ref\n")
            raw_data.write(init_seq + "\n")
            raw_data.write(">alt\n")
            raw_data.write(ref_seq + "\n")
        with open(fasta_q_name, "w") as fastq:
            for seq in corr_seq:
                qual = rq.ReadQuality(seq)
                read_id = random.randint(0, 10000000)
                fastq.write("@SEQ_ID" + str(read_id) + "\n")
                fastq.write(seq + "\n")
                fastq.write("+\n")
                fastq.write(qual.randomize_quality() + "\n")
        with open(csv_output_name, "w") as csv_output:
            csv_output.write("Init,Reference,ReverseComplem,Mutation,Positions,# untouched\n")
            csv_output.write(",".join(res_array))
        with open(csv_report_name, "w") as csv_report:
            csv_report.write("Original,Mutated,OutputSeq,Type,Left,Right,Modification\n")
            for i in range(len(corr_seq)):
                curr_val = ",".join([orig_arr[i],
                                     ref_seq,
                                     corr_seq[i],
                                     report[i][0],
                                     str(report[i][1]),
                                     str(report[i][2]),
                                     report[i][3]]
                                    ) + "\n"
                csv_report.write(curr_val)

    def batch_create_file(self):
        """

        :return: none

        This method goes over the array of nested JSONs within the mutationsArray field of the config.json
        and creates files as per the specification in the current array element.
        """
        for item in self.to_parse:
            self.create_file(item)






