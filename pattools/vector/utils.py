def parse_file_list(file_list):
    input_files = []
    groups = []
    samples = []
    with open(file_list, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            input_files.append(line[0])
            groups.append(line[1])
            samples.append(line[2])
    return input_files, groups, samples
