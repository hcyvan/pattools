from pattools.io import Output, Open
from pattools.log import logger


def compress_and_tabix_pat(input_file, output):
    logger.info(f"Compress and tabixing {input_file}")
    with Open(input_file) as infile, Output(output, file_format='pat', bgzip=True) as of:
        for line in infile:
            of.write(line)
