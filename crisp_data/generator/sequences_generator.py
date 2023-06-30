import random
import CrispTestDataGenerator_3.crisp_data.generator.reader_error as re


class SequencesGenerator:
    """

    This class is used to create random sequence and its CRISP modifications coming from the specified reader.
    """
    letters = ['A', "T", 'G', 'C']

    l_map = {
        "A": ["T", "G", "C"],
        "T": ["A", "G", "C"],
        "G": ["T", "A", "C"],
        "C": ["T", "A", "G"]
    }

    complementary_nucleotides = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    start = ['TGA']
    ends = ['TAG', 'TAA', 'TGA']

    def __init__(self, sequences_description):
        self.description = sequences_description
        self.positions = self.description["mutationSite"]
        self.length = self.description["initLength"]
        self.mutation_lens = self.description["mutationLength"]
        self.mutation_types = []
        self.deletions = [True if mutation < 0 else False for mutation in self.mutation_lens]
        self.is_coding = self.description["isCoding"]
        self.directions = self.description["direction"]
        self.quality = self.description["quality"]
        if self.quality % 5 != 0:
            raise "Quality must be 5-fold."
        self.count = self.description["count"]
        self.size_scale = self.description["sizeScale"]
        self.size_direction = self.description["scaleDirection"]
        self.reader = self.description["reader"]

    def cumulate_shifts(self, ind: int) -> int:
        """

        :param ind: index of the position where a modification happens
        :return: returns cumulated modification length (i.e. sum of insertions/deletions/SNPs lengths up to the index)
        """
        return sum(self.mutation_lens[:ind])

    def generate_coding(self) -> tuple:
        """

        :return:
        - init_sequence -- the random DNA sequence of the length specified in the initLength field of the config.json
          Depending on the isCoding fla value it can be coding or non-coding
        - final_sequence -- the sequence resulting from the CRISP
        - mutation -- the list of mutation positions introduced by the CRISP
        - m_map -- the dictionary with mutation types along with their positions and letters introduced by each mutation

        Mutation types are described in the mutationSite array, mutationLength array, and direction array in the
        mutationsArray nested JSONs stored in the config.json (see the create_file method of the Generator class for
        more details).
        """
        init_sequence = self.start[0] if self.is_coding else ""
        mutation = [""] * len(self.mutation_lens)
        m_map = {}
        seek_range = self.generate_range()
        init_sequence += self.generate_random_seq(seek_range)
        self.mutation_types = [abs(m_len) if m_len != 0 else 1 for m_len in self.mutation_lens]
        init_sequence = self.strip_inner_stop_codons(init_sequence) if self.is_coding else init_sequence
        left_sequence = ""
        right_sequence = ""
        count = 0
        cur_sequence = init_sequence
        for i in range(len(self.mutation_lens)):
            mutation[i] = self.get_mutation(self.mutation_lens[i], self.positions[i], init_sequence)
            m_type = self.mutation_lens[i]
            m_pos = self.positions[i]
            if m_type < 0:
                self.update_map(m_map, "D", [mutation[i], str(m_pos)])
            elif m_type == 0:
                self.update_map(m_map, "S", [mutation[i], str(m_pos)])
            else:
                self.update_map(m_map, "I", [mutation[i], str(m_pos)])
        for i in range(len(self.positions)):
            pos = self.positions[i]
            pos = pos + self.cumulate_shifts(i)
            m_len = self.mutation_lens[i]
            direction = self.directions[i]
            left_sequence, right_sequence, mutation = self.mutate(pos,
                                                                  m_len,
                                                                  direction,
                                                                  cur_sequence,
                                                                  mutation,
                                                                  count)
            cur_sequence = left_sequence + right_sequence
            count += 1
        final_sequence = left_sequence + right_sequence
        return init_sequence, final_sequence, mutation, m_map

    def get_mutation(self, m_len: int, pos: int, init_sequence: str) -> str:
        """

        :param m_len: the length of the mutation
        :param pos: the position of the mutation
        :param init_sequence: the initial sequence to which mutation happens
        :return: the mutation letters on the position with the defined length
        """
        res = ""
        if m_len > 0:
            for i in range(m_len):
                res += self.letters[random.randint(a=0, b=3)]
        elif m_len == 0:
            res = self.random_mutate(init_sequence[pos - 1:pos])
        return res

    @staticmethod
    def update_map(m_map: dict, letter: str, elem: list) -> None:
        """

        :param m_map: the dictionary with mutations where key is the mutation type, and value is the list of
        corresponding mutations
        :param letter: the mutation type
        :param elem: the mutation letters in the list
        :return: the method does not return anything, and it takes the mutation dictionary by reference and updates it.
        """
        if letter in m_map.keys():
            cur_array = m_map[letter]
            cur_array.append(elem)
            m_map[letter] = cur_array
        else:
            m_map[letter] = [elem]

    def corrupt_generated(self, orig_sequence: str, sequence: str) -> tuple:
        """

        :param orig_sequence: the initial sequence which needs to be corrupted
        :param sequence: the sequence resulting from CRISP
        :return:
        - s_array -- array with either original sequence, or CRISPed sequence, or their corrupted versions
        - orig_array -- array with CRISPed sequences
        - c_int -- number of non-modified sequences
        - report -- the array of arrays where description of changes to the sequence is provided, i.e.
          ["corrupted", left, right, mod] -- meaning that the sequence has been corrupted, with starting position
          being left, ending position being right, and mod being the type of corruption (deletion/addition/no change)

        The level of corruption is defined by the quality field in the mutationsArray nested JSON stored in the
        config.json. The higher quality level is the less corruption happens.
        """
        if self.quality > 100:
            raise "Quality cannot exceed 100%"
        c_count = int(self.count * (100 - self.quality) / 100)
        c_int = random.randint(0, c_count // 2)
        r_count = int(self.count * self.quality / 200)
        s_array = []
        orig_array = []
        report = {}
        for i in range(0, r_count):
            seq = self.reverse_complement(sequence)
            seq, d = re.introduce_reader_error(seq, reader=self.get_reader())
            orig_array.append(seq)
            report[i] = ["reverse complem", 0, 0, "no change"]
            s_array.append(seq)
        for i in range(r_count, c_int + r_count):
            seq, d = re.introduce_reader_error(orig_sequence, reader=self.get_reader())
            orig_array.append(seq)
            report[i] = ["original", 0, 0, "no change"]
            s_array.append(seq)
        for i in range(c_int + r_count, c_count + r_count):
            seq = self.corrupt(sequence)
            seq, d = re.introduce_reader_error(seq, reader=self.get_reader())
            orig_array.append(seq)
            seq, left, right, mod = self.modify_length(seq)
            report[i] = ["corrupted", left, right, mod]
            s_array.append(seq)
        for i in range(c_count + r_count, self.count):
            report[i] = ["mutated", 0, 0, "no change"]
            seq, d = re.introduce_reader_error(sequence, reader=self.get_reader())
            s_array.append(seq)
            orig_array.append(seq)
        return s_array, orig_array, c_int, report

    def corrupt(self, sequence: str) -> str:
        """

        :param sequence: the sequence to be corrupted
        :return: the corrupted sequence
        The method uses random number of mutations and then introduces random SNPs to the original sequence.
        Maximal number of mutations is half of the sequence length.
        """
        init_len = len(sequence)
        max_mutation = init_len // 2
        c_count = random.randint(1, max_mutation)
        mutations = []
        res = sequence
        for i in range(c_count):
            mutations.append(random.randint(0, init_len - 1))
        for i in range(c_count):
            left = res[0:mutations[i]]
            mutation = res[mutations[i]:mutations[i] + 1]
            right = res[mutations[i] + 1:len(sequence)]
            mutation = self.random_mutate(mutation)
            res = left + mutation + right
        return res

    def random_mutate(self, letter: str) -> str:
        """

        :param letter: initial letter to be replaced by a random SNP
        :return:
        The random letter replacing the original letter, it cannot be the same letter as the original one.
        """
        if letter not in self.l_map.keys():
            return ""
        return self.l_map[letter][random.randint(a=0, b=2)]

    def mutate(self, pos: int, m_len: int, direction: bool, init_sequence: str, mutation: dict, count: int) -> tuple:
        """

        :param pos: position of the mutation
        :param m_len: length of the mutation
        :param direction: the direction of the mutation: True is the right-sided mutation,
        False is the left-sided mutation
        :param init_sequence: sequence to be mutated
        :param mutation: the dictionary with mutations
        :param count: key of the required mutation which is the number of the mutation
        :return:
        - left_sequence -- the left side of the original sequence before the mutation site
        - right_sequence -- the right side of the original sequence after the mutation site
        - mutation -- the mutated sequence
        """
        if m_len > 0:
            left_sequence = init_sequence[0:pos]
            right_sequence = init_sequence[pos:len(init_sequence)]
            left_sequence = left_sequence + mutation[count]
        elif m_len < 0:
            cut = self.mutation_lens[count] if direction else 0
            cut_r = 0 if direction else -self.mutation_lens[count]
            left_sequence = init_sequence[0:pos + cut]
            right_sequence = init_sequence[pos + cut_r:len(init_sequence)]
            mutation[count] = init_sequence[pos + cut:pos] if direction \
                else init_sequence[pos:pos + cut_r]
        else:
            left_sequence = init_sequence[0:pos - 1]
            right_sequence = init_sequence[pos:len(init_sequence)]
            left_sequence = left_sequence + mutation[count]
        return left_sequence, right_sequence, mutation

    def generate_range(self) -> range:
        """

        :return:
        Returns the range within which the random sequence to be generated.
        """
        if self.length % 3 != 0 and self.is_coding:
            raise "Cannot generate coding sequence not being 3-fold"
        if self.is_coding:
            return range(3, self.length - 3)
        else:
            return range(0, self.length)

    def strip_inner_stop_codons(self, init_sequence: str) -> str:
        """

        :param init_sequence: the coding sequence where randomly introduced stop codons should be replaced by AA codons
        :return:
        The method returns the coding sequence where inner stop codons are replaced by AA codons.
        """
        strip_array = []
        codon_len = 3
        total = len(init_sequence) // 3
        count = 0
        for i in range(total):
            strip_array.append(init_sequence[count:count + codon_len])
            count += codon_len
        for i in range(len(strip_array)):
            if strip_array[i] in SequencesGenerator.ends:
                strip_array[i] = strip_array[i][::-1]
        init_sequence = "".join(strip_array) + self.ends[random.randint(a=0, b=2)]
        return init_sequence

    def reverse_complement(self, sequence: str) -> str:
        """

        :param sequence: initial sequence
        :return:
        The method returns the reverse compliment of the initial sequence.
        """
        return "".join([self.complementary_nucleotides[x] for x in sequence[::-1]])

    def modify_length(self, sequence: str) -> tuple:
        """

        :param sequence: the sequence which length should be modified
        The range of modification is stored in the sizeScale field of the config.json
        :return:
        - res -- the resulting sequence
        - left -- the left side before the modification position
        - right -- the right side after the modification position
        - modification -- type of modification (addition/deletion/no change)
        """
        up_limit = self.size_scale[1] * len(sequence) // 100
        down_limit = self.size_scale[0] * len(sequence) // 100
        limit = random.randint(down_limit, up_limit)
        if limit < len(sequence):
            res, left, right = self.cut(sequence, limit)
            modification = "cut"
        elif limit > len(sequence):
            res, left, right = self.add(sequence, limit)
            modification = "add"
        else:
            left = 0
            right = len(sequence)
            res = sequence
            modification = "no change"
        return res, left, right, modification

    def cut(self, sequence: str, limit: int) -> tuple:
        """

        :param sequence: the sequence to which deletion is applied
        :param limit: the parameter defining the limit of the sequence length depending on the scaleDirection value from
        the config.json (-1 means that it is the bottom limit, 1 means that it is the top limit, 0 means that sequence
        is cut from both beginning and end).
        :return:
        The cut sequence, the start position of the cut, and the end position of the cut.
        """
        seq_len = len(sequence)
        shift = seq_len - limit
        five_percent = seq_len // 20
        if self.size_direction == -1:
            left = shift
            right = seq_len
        elif self.size_direction == 0:
            left = random.randint(min(five_percent, shift - five_percent), max(five_percent, shift - five_percent))
            right = seq_len + left - shift
        else:
            left = 0
            right = seq_len - shift
        return sequence[left:right], left, right

    def generate_random_seq(self, length: range) -> str:
        """

        :param length: the range of the sequence length
        :return:
        The random DNA sequence within the specified range.
        """
        return "".join([self.letters[random.randint(0, 3)] for _ in length])

    def add(self, sequence: str, limit: int) -> tuple:
        """

        :param sequence: the sequence to which insertion is applied
        :param limit: the parameter defining the limit of the sequence length depending on the scaleDirection value from
        the config.json (-1 means that it is the bottom limit, 1 means that it is the top limit, 0 means that sequence
        can be enlarged from both beginning and end).
        :return:
        The sequence with insertion, the start position of the insertion, and the end position of the insertion.
        """
        seq_len = len(sequence)
        shift = limit - seq_len
        five_percent = seq_len // 20
        if self.size_direction == -1:
            left = shift
            right = seq_len + left
            res = self.generate_random_seq(range(left)) + sequence
        elif self.size_direction == 0:
            left = random.randint(min(five_percent, shift - five_percent), max(five_percent, shift - five_percent))
            right = shift - left
            res = self.generate_random_seq(range(left)) + sequence + self.generate_random_seq(range(right))
            right = left + seq_len
        else:
            left = seq_len
            right = shift
            res = sequence + self.generate_random_seq(range(right))
            right += seq_len
        return res, left, right

    def get_reader(self) -> re.Reader:
        """

        :return: returns the Reader class instance depending on the name specified in the reader field in the
        config.json
        """
        try:
            return list(filter(lambda x: x.name == self.reader, list(re.ReadersEnum)))[0].value
        except:
            return re.Readers.MiSeq
