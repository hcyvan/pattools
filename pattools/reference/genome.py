import shutil
from pathlib import Path
from collections import OrderedDict
from pattools.io import Open, Output
from pattools.log import logger

def genome_build(fasta, output, chromosomes=None):
    logger.info("Start building genome...")
    output = Path(output)
    fasta = Path(fasta)
    tmp = output / 'tmp'
    if not output.exists():
        output.mkdir()
    if not tmp.exists():
        tmp.mkdir()
    fasta_new = output / fasta.name
    out_cpg_bed = output / 'CpG.bed.gz'
    if chromosomes is None:
        shutil.copy(fasta, fasta_new)
    else:
        logger.info(f"generate genome include chromosomes {chromosomes}")
        chrom_map=OrderedDict(zip(chromosomes.split(','), [tmp / f'{fasta.name}.{x}' for x in chromosomes.split(',')]))
        chrom_fo=None
        with Open(fasta) as fin:
            for line in fin:
                if line.startswith('>'):
                    chrom = line.split()[0][1:]
                    chrom_fo = None
                    if chrom in chrom_map:
                        chrom_fo = Output(chrom_map[chrom])
                        chrom_fo.write(line)
                else:
                    if chrom_fo is not None:
                        chrom_fo.write(line)
        if chrom_fo is not None:
            chrom_fo.close()
        with Output(fasta_new, file_format=None, bgzip=True) as fo:
            for k,v in chrom_map.items():
                with Open(v) as fi:
                    for line in fi:
                        fo.write(line)
    logger.info(f"build Cpg.bed.gz")
    generate_genome_cpg_map(fasta_new, out_cpg_bed)
    shutil.rmtree(tmp)


def generate_genome_cpg_map(fasta, output):
    pre = ""
    genome_index = 0
    cpg_index = 1
    chromosomes_array = []
    with Open(fasta) as fin, Output(output, bgzip=True) as fo:
        for line in fin:
            line = line.strip()
            if line.startswith('>'):
                chrom = line.split()[0][1:]
                chromosomes_array.append(chrom)
                pre = ""
                genome_index = 0
            else:
                for i, p in enumerate(line):
                    if (pre == 'C' or pre == 'c') and (p == 'G' or p == 'g'):
                        fo.write(f"{chromosomes_array[-1]}\t{genome_index}\t{cpg_index}\n")
                        cpg_index += 1
                    genome_index += 1
                    pre = p
