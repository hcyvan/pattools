def parse_file_list(file_list, target_groups=None):
    input_files = []
    groups = []
    samples = []
    if target_groups is not None:
        target_groups = target_groups.split(',')
    with open(file_list, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if target_groups is not None:
                if line[1] in target_groups:
                    input_files.append(line[0])
                    groups.append(line[1])
                    samples.append(line[2])
            else:
                input_files.append(line[0])
                groups.append(line[1])
                samples.append(line[2])
    return input_files, groups, samples
