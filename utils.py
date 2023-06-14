import re

def extract_number_from_filename(file_path):
    pattern = r'\d+'
    match = re.search(pattern, file_path)
    if match:
        return int(match.group())
    else:
        return 0