def is_gzip_file(filepath):
    with open(filepath, 'rb') as file:
        file_header = file.read(2)
    return file_header == b'\x1f\x8b'
