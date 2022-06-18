from Bio.Seq import Seq
from Bio import SeqRecord


class TSamRecord:
    def __init__(self, qname: str, flag: int, rname: str, pos: int, mapq: int, ref_rec: SeqRecord,
                 rnext: str, pnext: int, tlen: int, segment: Seq, qual: str):
        """
        Each parameter correpsonds to the official SAM format field indicated

        :param qname: QNAME - query name
        :param flag: FLAEG - combination bitwise flags
        :param rname: RNAME - reference sequence name of the alignment, must match an SQ header
        :param pos: POS - 1-based leftmost mapping position. This will be stored as 0-based to make calculations easier
        :param mapq: MAPQ - mapping quality of the read, integer
        :param ref_rec: Unique to this, this is the SeqRecord object of the reference sequence.
         Will use to calculate CIGAR
        :param rnext: RNEXT - Ref name of the mate/next read
        :param pnext: PNEXT - position of the mate/next read (stored as 0-based)
        :param tlen: TLEN - Observerd template length (fragment length)
        :param segment: SEQ - The sequence from the read, called segment to avoid confusion with Biopython Seq
        :param qual: QUAL - ASCII of Phred-based quality score (quality + 33)
        """

        self.qname = qname
        self.flag = flag
        self.rname = rname
        self.pos = pos
        self.mapq = mapq
        self.ref_rec = ref_rec
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = tlen
        self.segment = segment
        self.segment_length = len(segment)
        self.qual = qual

    def get_reference_sequence(self):
        return self.ref_rec[self.pos: self.pos + self.segment_length]

    def get_sam_pos(self):
        return self.pos + 1

    def get_pnext_pos(self):
        return self.pnext + 1

    def generate_cigar(self):
        pass

