s = """
SEQUENCE_ID={}
SEQUENCE_TEMPLATE={}
SEQUENCE_TARGET=37,21
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=18
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=21
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_PRODUCT_SIZE_RANGE=75-100
P3_FILE_FLAG=1
SEQUENCE_INTERNAL_EXCLUDED_REGION=37,21
PRIMER_EXPLAIN_FLAG=1
=
""".strip()
import random


def random_seq():

    seq = ""
    for l in range(1000):
        i = random.randint(0, 3)
        seq += "AGTC"[i]
    return seq


with open("example2.io", "w") as f:
    for i in range(25):
        f.write(s.format("example" + str(i), random_seq()))
        f.write("\n")
