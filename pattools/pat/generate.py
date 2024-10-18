import pysam
from pysam.libcalignedsegment import AlignedSegment
from pattools.io import Output, Open
from pattools.pat.pat import PatStep
from pattools.log import logger

def bam2pat(input_file, output):
    samfile = pysam.AlignmentFile(input_file, "rb")
    # print(samfile.header)
    # for k,v in samfile.header.items():
    #     print(1111111111111111111)
    #     print(k,v)

    for read in samfile.fetch(): # AlignedSegment
        print(read.reference_name, read.reference_id)


    samfile.close()

def compress_and_tabix_pat(input_file, output):
    logger.info(f"Compress and tabixing {input_file}")
    with Open(input_file) as infile, Output(output, file_format='pat', bgzip=True) as of:
        for line in infile:
            of.write(line)


def merge_pat(pat_files, output):
    pats = []
    top = []
    for pat_file in pat_files:
        pat = PatStep(pat_file)
        pats.append(pat)
        top.append(pat.read_item())
    with Output(output, file_format='pat', bgzip=True) as of:
        while True:
            if all(x is None for x in top):
                break
            cpg_least = 9999999999999
            for item in top:
                if item:
                    if item.cpg < cpg_least:
                        cpg_least = item.cpg
            pat_item = None
            for i, item in enumerate(top):
                if item:
                    if item.cpg == cpg_least:
                        if pat_item is None:
                            pat_item = item
                        else:
                            pat_item += item
                        top[i] = pats[i].read_item()
            for line in pat_item.get_lines():
                of.write(line + '\n')
