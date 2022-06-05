from error_handling import log_mssg, premature_exit


class Variant:
    def __init__(self, position1: int, length1: int, position2: int = None, length2: int = None):
        """
        DUP: pos1, len1, pos2=pos1+len1, len2=NULL
        INV: pos1, len1, pos2=NULL, len2=-len1
        DEL: pos1, len1, pos2=NULL, len2=NULL
        INS: pos1, len1, pos2, len2=NULL
        TRA: pos1, len1, pos2, len2
        :param position1: first position of variant. This will be the origin of the variant.
        :param length1: Length of the variant.
        :param position2: second position of the variant. This will be, for example,
        where the duplication startsor where the translocation was inserted
        :param length2: Mainly used to determine the length of the other segment that was swapped in the case of
        translocation. A negative number represents and inversion.
        """

        self.position1 = position1
        self.length1 = length1
        self.position2 = position2
        self.length2 = length2
        # There isn't a situation where we want a length2 == 0. Maybe length2 should get replaced by a None, but for now
        # we'll try this check
        if self.length2 == 0:
            log_mssg(f"Invalid variant detected, length2 should never = 0. length2 = {self.length2}", 'critical')
            premature_exit(1)

    def is_translocation(self):
        return self.length2 > 0

    def is_inversion(self):
        return self.length2 < 0

    def is_deletion(self):
        return not self.position2 and not self.length2

    def is_ins_or_dup(self):
        return self.position2 and not self.length2

    def is_duplication(self):
        return self.is_ins_or_dup() and self.position2 == self.position1 + self.length1

    def is_insertion(self):
        return self.is_ins_or_dup() and not self.is_duplication()

