import random


class ReadQuality:

    """
    This is the class generating a random read quality string.
    """

    quality_str = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"

    def __init__(self, sequence: str):
        self.sequence = sequence

    def randomize_quality(self) -> str:
        """

        :return:
        This method returns the random string describing the read quality of the sequence added to the
        instance of the class.
        """
        res = ""
        for i in range(len(self.sequence)):
            ind = random.randint(0, len(self.quality_str) - 1)
            res += self.quality_str[ind]
        return res

