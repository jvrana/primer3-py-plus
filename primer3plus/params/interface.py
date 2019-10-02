class ParameterDescriptor:
    """A data descriptor that sets params key."""

    def __init__(self, name):
        self.name = name

    def __get__(self, obj, objtype):
        return obj.params[self.name]

    def __set__(self, obj, val):
        obj.params[self.name] = val


class ParameterInterface:
    PRIMER_DNA_CONC = ParameterDescriptor("PRIMER_DNA_CONC")
    PRIMER_MASK_KMERLIST_PATH = ParameterDescriptor("PRIMER_MASK_KMERLIST_PATH")
    PRIMER_PAIR_WT_PRODUCT_SIZE_LT = ParameterDescriptor(
        "PRIMER_PAIR_WT_PRODUCT_SIZE_LT"
    )
    PRIMER_DNTP_CONC = ParameterDescriptor("PRIMER_DNTP_CONC")
    PRIMER_MASK_KMERLIST_PREFIX = ParameterDescriptor("PRIMER_MASK_KMERLIST_PREFIX")
    PRIMER_PAIR_WT_PRODUCT_TM_GT = ParameterDescriptor("PRIMER_PAIR_WT_PRODUCT_TM_GT")
    PRIMER_EXPLAIN_FLAG = ParameterDescriptor("PRIMER_EXPLAIN_FLAG")
    PRIMER_MASK_TEMPLATE = ParameterDescriptor("PRIMER_MASK_TEMPLATE")
    PRIMER_PAIR_WT_PRODUCT_TM_LT = ParameterDescriptor("PRIMER_PAIR_WT_PRODUCT_TM_LT")
    PRIMER_FIRST_BASE_INDEX = ParameterDescriptor("PRIMER_FIRST_BASE_INDEX")
    PRIMER_MAX_END_GC = ParameterDescriptor("PRIMER_MAX_END_GC")
    PRIMER_PAIR_WT_PR_PENALTY = ParameterDescriptor("PRIMER_PAIR_WT_PR_PENALTY")
    PRIMER_GC_CLAMP = ParameterDescriptor("PRIMER_GC_CLAMP")
    PRIMER_MAX_END_STABILITY = ParameterDescriptor("PRIMER_MAX_END_STABILITY")
    PRIMER_PAIR_WT_TEMPLATE_MISPRIMING = ParameterDescriptor(
        "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING"
    )
    PRIMER_INSIDE_PENALTY = ParameterDescriptor("PRIMER_INSIDE_PENALTY")
    PRIMER_MAX_GC = ParameterDescriptor("PRIMER_MAX_GC")
    PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH = ParameterDescriptor(
        "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH"
    )
    PRIMER_INTERNAL_DNA_CONC = ParameterDescriptor("PRIMER_INTERNAL_DNA_CONC")
    PRIMER_MAX_HAIRPIN_TH = ParameterDescriptor("PRIMER_MAX_HAIRPIN_TH")
    PRIMER_PICK_ANYWAY = ParameterDescriptor("PRIMER_PICK_ANYWAY")
    PRIMER_INTERNAL_DNTP_CONC = ParameterDescriptor("PRIMER_INTERNAL_DNTP_CONC")
    PRIMER_MAX_LIBRARY_MISPRIMING = ParameterDescriptor("PRIMER_MAX_LIBRARY_MISPRIMING")
    PRIMER_PICK_INTERNAL_OLIGO = ParameterDescriptor("PRIMER_PICK_INTERNAL_OLIGO")
    PRIMER_INTERNAL_MAX_GC = ParameterDescriptor("PRIMER_INTERNAL_MAX_GC")
    PRIMER_MAX_NS_ACCEPTED = ParameterDescriptor("PRIMER_MAX_NS_ACCEPTED")
    PRIMER_PICK_LEFT_PRIMER = ParameterDescriptor("PRIMER_PICK_LEFT_PRIMER")
    PRIMER_INTERNAL_MAX_HAIRPIN_TH = ParameterDescriptor(
        "PRIMER_INTERNAL_MAX_HAIRPIN_TH"
    )
    PRIMER_MAX_POLY_X = ParameterDescriptor("PRIMER_MAX_POLY_X")
    PRIMER_PICK_RIGHT_PRIMER = ParameterDescriptor("PRIMER_PICK_RIGHT_PRIMER")
    PRIMER_INTERNAL_MAX_LIBRARY_MISHYB = ParameterDescriptor(
        "PRIMER_INTERNAL_MAX_LIBRARY_MISHYB"
    )
    PRIMER_MAX_SELF_ANY = ParameterDescriptor("PRIMER_MAX_SELF_ANY")
    PRIMER_PRODUCT_MAX_TM = ParameterDescriptor("PRIMER_PRODUCT_MAX_TM")
    PRIMER_INTERNAL_MAX_NS_ACCEPTED = ParameterDescriptor(
        "PRIMER_INTERNAL_MAX_NS_ACCEPTED"
    )
    PRIMER_MAX_SELF_ANY_TH = ParameterDescriptor("PRIMER_MAX_SELF_ANY_TH")
    PRIMER_PRODUCT_MIN_TM = ParameterDescriptor("PRIMER_PRODUCT_MIN_TM")
    PRIMER_INTERNAL_MAX_POLY_X = ParameterDescriptor("PRIMER_INTERNAL_MAX_POLY_X")
    PRIMER_MAX_SELF_END = ParameterDescriptor("PRIMER_MAX_SELF_END")
    PRIMER_PRODUCT_OPT_SIZE = ParameterDescriptor("PRIMER_PRODUCT_OPT_SIZE")
    PRIMER_INTERNAL_MAX_SELF_ANY = ParameterDescriptor("PRIMER_INTERNAL_MAX_SELF_ANY")
    PRIMER_MAX_SELF_END_TH = ParameterDescriptor("PRIMER_MAX_SELF_END_TH")
    PRIMER_PRODUCT_OPT_TM = ParameterDescriptor("PRIMER_PRODUCT_OPT_TM")
    PRIMER_INTERNAL_MAX_SELF_ANY_TH = ParameterDescriptor(
        "PRIMER_INTERNAL_MAX_SELF_ANY_TH"
    )
    PRIMER_MAX_SIZE = ParameterDescriptor("PRIMER_MAX_SIZE")
    PRIMER_PRODUCT_SIZE_RANGE = ParameterDescriptor("PRIMER_PRODUCT_SIZE_RANGE")
    PRIMER_INTERNAL_MAX_SELF_END = ParameterDescriptor("PRIMER_INTERNAL_MAX_SELF_END")
    PRIMER_MAX_TEMPLATE_MISPRIMING = ParameterDescriptor(
        "PRIMER_MAX_TEMPLATE_MISPRIMING"
    )
    PRIMER_QUALITY_RANGE_MAX = ParameterDescriptor("PRIMER_QUALITY_RANGE_MAX")
    PRIMER_INTERNAL_MAX_SELF_END_TH = ParameterDescriptor(
        "PRIMER_INTERNAL_MAX_SELF_END_TH"
    )
    PRIMER_MAX_TEMPLATE_MISPRIMING_TH = ParameterDescriptor(
        "PRIMER_MAX_TEMPLATE_MISPRIMING_TH"
    )
    PRIMER_QUALITY_RANGE_MIN = ParameterDescriptor("PRIMER_QUALITY_RANGE_MIN")
    PRIMER_INTERNAL_MAX_SIZE = ParameterDescriptor("PRIMER_INTERNAL_MAX_SIZE")
    PRIMER_MAX_TM = ParameterDescriptor("PRIMER_MAX_TM")
    PRIMER_SALT_CORRECTIONS = ParameterDescriptor("PRIMER_SALT_CORRECTIONS")
    PRIMER_INTERNAL_MAX_TM = ParameterDescriptor("PRIMER_INTERNAL_MAX_TM")
    PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION = ParameterDescriptor(
        "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION"
    )
    PRIMER_SALT_DIVALENT = ParameterDescriptor("PRIMER_SALT_DIVALENT")
    PRIMER_INTERNAL_MIN_GC = ParameterDescriptor("PRIMER_INTERNAL_MIN_GC")
    PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION = ParameterDescriptor(
        "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION"
    )
    PRIMER_SALT_MONOVALENT = ParameterDescriptor("PRIMER_SALT_MONOVALENT")
    PRIMER_INTERNAL_MIN_QUALITY = ParameterDescriptor("PRIMER_INTERNAL_MIN_QUALITY")
    PRIMER_MIN_END_QUALITY = ParameterDescriptor("PRIMER_MIN_END_QUALITY")
    PRIMER_SECONDARY_STRUCTURE_ALIGNMENT = ParameterDescriptor(
        "PRIMER_SECONDARY_STRUCTURE_ALIGNMENT"
    )
    PRIMER_INTERNAL_MIN_SIZE = ParameterDescriptor("PRIMER_INTERNAL_MIN_SIZE")
    PRIMER_MIN_GC = ParameterDescriptor("PRIMER_MIN_GC")
    PRIMER_SEQUENCING_ACCURACY = ParameterDescriptor("PRIMER_SEQUENCING_ACCURACY")
    PRIMER_INTERNAL_MIN_TM = ParameterDescriptor("PRIMER_INTERNAL_MIN_TM")
    PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE = ParameterDescriptor(
        "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE"
    )
    PRIMER_SEQUENCING_INTERVAL = ParameterDescriptor("PRIMER_SEQUENCING_INTERVAL")
    PRIMER_INTERNAL_MISHYB_LIBRARY = ParameterDescriptor(
        "PRIMER_INTERNAL_MISHYB_LIBRARY"
    )
    PRIMER_MIN_QUALITY = ParameterDescriptor("PRIMER_MIN_QUALITY")
    PRIMER_SEQUENCING_LEAD = ParameterDescriptor("PRIMER_SEQUENCING_LEAD")
    PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME = ParameterDescriptor(
        "PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME"
    )
    PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE = ParameterDescriptor(
        "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE"
    )
    PRIMER_SEQUENCING_SPACING = ParameterDescriptor("PRIMER_SEQUENCING_SPACING")
    PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME = ParameterDescriptor(
        "PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME"
    )
    PRIMER_MIN_SIZE = ParameterDescriptor("PRIMER_MIN_SIZE")
    PRIMER_TASK = ParameterDescriptor("PRIMER_TASK")
    PRIMER_INTERNAL_OPT_GC_PERCENT = ParameterDescriptor(
        "PRIMER_INTERNAL_OPT_GC_PERCENT"
    )
    PRIMER_MIN_THREE_PRIME_DISTANCE = ParameterDescriptor(
        "PRIMER_MIN_THREE_PRIME_DISTANCE"
    )
    PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = ParameterDescriptor(
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT"
    )
    PRIMER_INTERNAL_OPT_SIZE = ParameterDescriptor("PRIMER_INTERNAL_OPT_SIZE")
    PRIMER_MIN_TM = ParameterDescriptor("PRIMER_MIN_TM")
    PRIMER_THERMODYNAMIC_PARAMETERS_PATH = ParameterDescriptor(
        "PRIMER_THERMODYNAMIC_PARAMETERS_PATH"
    )
    PRIMER_INTERNAL_OPT_TM = ParameterDescriptor("PRIMER_INTERNAL_OPT_TM")
    PRIMER_MISPRIMING_LIBRARY = ParameterDescriptor("PRIMER_MISPRIMING_LIBRARY")
    PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT = ParameterDescriptor(
        "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT"
    )
    PRIMER_INTERNAL_SALT_DIVALENT = ParameterDescriptor("PRIMER_INTERNAL_SALT_DIVALENT")
    PRIMER_MUST_MATCH_FIVE_PRIME = ParameterDescriptor("PRIMER_MUST_MATCH_FIVE_PRIME")
    PRIMER_TM_FORMULA = ParameterDescriptor("PRIMER_TM_FORMULA")
    PRIMER_INTERNAL_SALT_MONOVALENT = ParameterDescriptor(
        "PRIMER_INTERNAL_SALT_MONOVALENT"
    )
    PRIMER_MUST_MATCH_THREE_PRIME = ParameterDescriptor("PRIMER_MUST_MATCH_THREE_PRIME")
    PRIMER_WT_END_QUAL = ParameterDescriptor("PRIMER_WT_END_QUAL")
    PRIMER_INTERNAL_WT_END_QUAL = ParameterDescriptor("PRIMER_INTERNAL_WT_END_QUAL")
    PRIMER_NUM_RETURN = ParameterDescriptor("PRIMER_NUM_RETURN")
    PRIMER_WT_END_STABILITY = ParameterDescriptor("PRIMER_WT_END_STABILITY")
    PRIMER_INTERNAL_WT_GC_PERCENT_GT = ParameterDescriptor(
        "PRIMER_INTERNAL_WT_GC_PERCENT_GT"
    )
    PRIMER_OPT_GC_PERCENT = ParameterDescriptor("PRIMER_OPT_GC_PERCENT")
    PRIMER_WT_GC_PERCENT_GT = ParameterDescriptor("PRIMER_WT_GC_PERCENT_GT")
    PRIMER_INTERNAL_WT_GC_PERCENT_LT = ParameterDescriptor(
        "PRIMER_INTERNAL_WT_GC_PERCENT_LT"
    )
    PRIMER_OPT_SIZE = ParameterDescriptor("PRIMER_OPT_SIZE")
    PRIMER_WT_GC_PERCENT_LT = ParameterDescriptor("PRIMER_WT_GC_PERCENT_LT")
    PRIMER_INTERNAL_WT_HAIRPIN_TH = ParameterDescriptor("PRIMER_INTERNAL_WT_HAIRPIN_TH")
    PRIMER_OPT_TM = ParameterDescriptor("PRIMER_OPT_TM")
    PRIMER_WT_HAIRPIN_TH = ParameterDescriptor("PRIMER_WT_HAIRPIN_TH")
    PRIMER_INTERNAL_WT_LIBRARY_MISHYB = ParameterDescriptor(
        "PRIMER_INTERNAL_WT_LIBRARY_MISHYB"
    )
    PRIMER_OUTSIDE_PENALTY = ParameterDescriptor("PRIMER_OUTSIDE_PENALTY")
    PRIMER_WT_LIBRARY_MISPRIMING = ParameterDescriptor("PRIMER_WT_LIBRARY_MISPRIMING")
    PRIMER_INTERNAL_WT_NUM_NS = ParameterDescriptor("PRIMER_INTERNAL_WT_NUM_NS")
    PRIMER_PAIR_MAX_COMPL_ANY = ParameterDescriptor("PRIMER_PAIR_MAX_COMPL_ANY")
    PRIMER_WT_MASK_FAILURE_RATE = ParameterDescriptor("PRIMER_WT_MASK_FAILURE_RATE")
    PRIMER_INTERNAL_WT_SELF_ANY = ParameterDescriptor("PRIMER_INTERNAL_WT_SELF_ANY")
    PRIMER_PAIR_MAX_COMPL_ANY_TH = ParameterDescriptor("PRIMER_PAIR_MAX_COMPL_ANY_TH")
    PRIMER_WT_NUM_NS = ParameterDescriptor("PRIMER_WT_NUM_NS")
    PRIMER_INTERNAL_WT_SELF_ANY_TH = ParameterDescriptor(
        "PRIMER_INTERNAL_WT_SELF_ANY_TH"
    )
    PRIMER_PAIR_MAX_COMPL_END = ParameterDescriptor("PRIMER_PAIR_MAX_COMPL_END")
    PRIMER_WT_POS_PENALTY = ParameterDescriptor("PRIMER_WT_POS_PENALTY")
    PRIMER_INTERNAL_WT_SELF_END = ParameterDescriptor("PRIMER_INTERNAL_WT_SELF_END")
    PRIMER_PAIR_MAX_COMPL_END_TH = ParameterDescriptor("PRIMER_PAIR_MAX_COMPL_END_TH")
    PRIMER_WT_SELF_ANY = ParameterDescriptor("PRIMER_WT_SELF_ANY")
    PRIMER_INTERNAL_WT_SELF_END_TH = ParameterDescriptor(
        "PRIMER_INTERNAL_WT_SELF_END_TH"
    )
    PRIMER_PAIR_MAX_DIFF_TM = ParameterDescriptor("PRIMER_PAIR_MAX_DIFF_TM")
    PRIMER_WT_SELF_ANY_TH = ParameterDescriptor("PRIMER_WT_SELF_ANY_TH")
    PRIMER_INTERNAL_WT_SEQ_QUAL = ParameterDescriptor("PRIMER_INTERNAL_WT_SEQ_QUAL")
    PRIMER_PAIR_MAX_LIBRARY_MISPRIMING = ParameterDescriptor(
        "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING"
    )
    PRIMER_WT_SELF_END = ParameterDescriptor("PRIMER_WT_SELF_END")
    PRIMER_INTERNAL_WT_SIZE_GT = ParameterDescriptor("PRIMER_INTERNAL_WT_SIZE_GT")
    PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING = ParameterDescriptor(
        "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING"
    )
    PRIMER_WT_SELF_END_TH = ParameterDescriptor("PRIMER_WT_SELF_END_TH")
    PRIMER_INTERNAL_WT_SIZE_LT = ParameterDescriptor("PRIMER_INTERNAL_WT_SIZE_LT")
    PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH = ParameterDescriptor(
        "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH"
    )
    PRIMER_WT_SEQ_QUAL = ParameterDescriptor("PRIMER_WT_SEQ_QUAL")
    PRIMER_INTERNAL_WT_TM_GT = ParameterDescriptor("PRIMER_INTERNAL_WT_TM_GT")
    PRIMER_PAIR_WT_COMPL_ANY = ParameterDescriptor("PRIMER_PAIR_WT_COMPL_ANY")
    PRIMER_WT_SIZE_GT = ParameterDescriptor("PRIMER_WT_SIZE_GT")
    PRIMER_INTERNAL_WT_TM_LT = ParameterDescriptor("PRIMER_INTERNAL_WT_TM_LT")
    PRIMER_PAIR_WT_COMPL_ANY_TH = ParameterDescriptor("PRIMER_PAIR_WT_COMPL_ANY_TH")
    PRIMER_WT_SIZE_LT = ParameterDescriptor("PRIMER_WT_SIZE_LT")
    PRIMER_LIBERAL_BASE = ParameterDescriptor("PRIMER_LIBERAL_BASE")
    PRIMER_PAIR_WT_COMPL_END = ParameterDescriptor("PRIMER_PAIR_WT_COMPL_END")
    PRIMER_WT_TEMPLATE_MISPRIMING = ParameterDescriptor("PRIMER_WT_TEMPLATE_MISPRIMING")
    PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS = ParameterDescriptor(
        "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS"
    )
    PRIMER_PAIR_WT_COMPL_END_TH = ParameterDescriptor("PRIMER_PAIR_WT_COMPL_END_TH")
    PRIMER_WT_TEMPLATE_MISPRIMING_TH = ParameterDescriptor(
        "PRIMER_WT_TEMPLATE_MISPRIMING_TH"
    )
    PRIMER_LOWERCASE_MASKING = ParameterDescriptor("PRIMER_LOWERCASE_MASKING")
    PRIMER_PAIR_WT_DIFF_TM = ParameterDescriptor("PRIMER_PAIR_WT_DIFF_TM")
    PRIMER_WT_TM_GT = ParameterDescriptor("PRIMER_WT_TM_GT")
    PRIMER_MASK_3P_DIRECTION = ParameterDescriptor("PRIMER_MASK_3P_DIRECTION")
    PRIMER_PAIR_WT_IO_PENALTY = ParameterDescriptor("PRIMER_PAIR_WT_IO_PENALTY")
    PRIMER_WT_TM_LT = ParameterDescriptor("PRIMER_WT_TM_LT")
    PRIMER_MASK_5P_DIRECTION = ParameterDescriptor("PRIMER_MASK_5P_DIRECTION")
    PRIMER_PAIR_WT_LIBRARY_MISPRIMING = ParameterDescriptor(
        "PRIMER_PAIR_WT_LIBRARY_MISPRIMING"
    )
    PRIMER_MASK_FAILURE_RATE = ParameterDescriptor("PRIMER_MASK_FAILURE_RATE")
    PRIMER_PAIR_WT_PRODUCT_SIZE_GT = ParameterDescriptor(
        "PRIMER_PAIR_WT_PRODUCT_SIZE_GT"
    )
    SEQUENCE_EXCLUDED_REGION = ParameterDescriptor("SEQUENCE_EXCLUDED_REGION")
    SEQUENCE_INCLUDED_REGION = ParameterDescriptor("SEQUENCE_INCLUDED_REGION")
    SEQUENCE_PRIMER_REVCOMP = ParameterDescriptor("SEQUENCE_PRIMER_REVCOMP")
    SEQUENCE_FORCE_LEFT_END = ParameterDescriptor("SEQUENCE_FORCE_LEFT_END")
    SEQUENCE_INTERNAL_EXCLUDED_REGION = ParameterDescriptor(
        "SEQUENCE_INTERNAL_EXCLUDED_REGION"
    )
    SEQUENCE_QUALITY = ParameterDescriptor("SEQUENCE_QUALITY")
    SEQUENCE_FORCE_LEFT_START = ParameterDescriptor("SEQUENCE_FORCE_LEFT_START")
    SEQUENCE_INTERNAL_OLIGO = ParameterDescriptor("SEQUENCE_INTERNAL_OLIGO")
    SEQUENCE_START_CODON_POSITION = ParameterDescriptor("SEQUENCE_START_CODON_POSITION")
    SEQUENCE_FORCE_RIGHT_END = ParameterDescriptor("SEQUENCE_FORCE_RIGHT_END")
    SEQUENCE_OVERLAP_JUNCTION_LIST = ParameterDescriptor(
        "SEQUENCE_OVERLAP_JUNCTION_LIST"
    )
    SEQUENCE_TARGET = ParameterDescriptor("SEQUENCE_TARGET")
    SEQUENCE_FORCE_RIGHT_START = ParameterDescriptor("SEQUENCE_FORCE_RIGHT_START")
    SEQUENCE_PRIMER = ParameterDescriptor("SEQUENCE_PRIMER")
    SEQUENCE_TEMPLATE = ParameterDescriptor("SEQUENCE_TEMPLATE")
    SEQUENCE_ID = ParameterDescriptor("SEQUENCE_ID")
    SEQUENCE_PRIMER_PAIR_OK_REGION_LIST = ParameterDescriptor(
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"
    )

    def __init__(self, params):
        self.params = params


class HasParameters:
    """A data accessor that exposes parameters to attribute level.

    e.g.  ::      P = HasParameters(); P.SEQUENCE_ID = "MY SEQ"
    """

    def __get__(self, obj, objtype):
        return ParameterInterface(obj.params)


# def print_and_check_parameters():
#
#     for p in _expected_opts:
#         print("{name} = {cls}(\"{name}\")".format(name=p, cls=ParameterDescriptor.__name__))
