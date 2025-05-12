import os, argparse

def search_string_in_files(directory, search_string):
    for root, dirs, files in os.walk(directory):
        for filename in files:
            filepath = os.path.join(root, filename)
            try:
                with open(filepath, 'r', encoding='utf-8', errors='ignore') as file:
                    for i, line in enumerate(file, 1):
                        if search_string in line:
                            print(f"Found: {filepath}, riga {i}: {line.strip()}")
            except Exception as e:
                print(f"Error reading {filepath}: {e}")


def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Check for specific strings in functions.')
    parser.add_argument('folder', type=str, help='Path to the folder of functions.')
    parser.add_argument('string', type=str, help='String to search for.')

    args = parser.parse_args()

    # Call the function with the provided arguments
    search_string_in_files(args.folder, args.string)



if __name__ == '__main__':
    main()
