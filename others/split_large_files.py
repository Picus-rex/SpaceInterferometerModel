import os, argparse

MAX_SIZE_MB = 50
MAX_SIZE_BYTES = MAX_SIZE_MB * 1024 * 1024

def split_large_txt_file(filepath):
    base, ext = os.path.splitext(filepath)
    part_num = 1
    current_size = 0
    output_lines = []

    def write_part(lines, num):
        out_path = f"{base}.{num}{ext}"
        with open(out_path, 'w', encoding='utf-8') as f:
            f.writelines(lines)
        print(f"Created: {out_path} ({len(lines)} rows)")

    with open(filepath, 'r', encoding='utf-8', errors='ignore') as file:
        for line in file:
            line_size = len(line.encode('utf-8'))
            if current_size + line_size > MAX_SIZE_BYTES:
                write_part(output_lines, part_num)
                part_num += 1
                output_lines = []
                current_size = 0
            output_lines.append(line)
            current_size += line_size

    # Scrive l'ultima parte se c'è qualcosa rimasto
    if output_lines:
        write_part(output_lines, part_num)


def split_large_txt_files_in_directory(root_dir):
    for foldername, subfolders, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.lower().endswith('.txt'):
                filepath = os.path.join(foldername, filename)
                try:
                    size = os.path.getsize(filepath)
                    if size > MAX_SIZE_BYTES:
                        print(f"Big file found: {filepath} ({size / (1024*1024):.2f} MB)")
                        split_large_txt_file(filepath)
                except Exception as e:
                    print(f"Error on {filepath}: {e}")



def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Split files larger than 50 MB.')
    parser.add_argument('folder', type=str, help='Path to the folder of functions.')

    args = parser.parse_args()

    # Call the function with the provided arguments
    split_large_txt_files_in_directory(args.folder)


if __name__ == '__main__':
    main()
