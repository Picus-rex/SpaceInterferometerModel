import argparse

def add_headers_to_file(input_file, output_file, interval):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    modified_lines = []
    header_count = 1

    for i in range(0, len(lines), interval):
        # Add the header for the current group
        modified_lines.append(f"== {header_count}")
        header_count += 1

        # Add the specified number of lines, transforming them into CSV format
        for line in lines[i:i + interval]:
            # Split the line into elements and join them with commas
            csv_line = ','.join(line.strip().split())
            modified_lines.append(csv_line)

    # Write the modified lines to the output file
    with open(output_file, 'w') as file:
        for line in modified_lines:
            file.write(line + '\n')

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Add headers to a text file at specified intervals.')
    parser.add_argument('input_file', type=str, help='The input text file to read from.')
    parser.add_argument('output_file', type=str, help='The output text file to write to.')
    parser.add_argument('interval', type=int, help='The number of lines between headers.')

    args = parser.parse_args()

    # Call the function with the provided arguments
    add_headers_to_file(args.input_file, args.output_file, args.interval)

if __name__ == '__main__':
    main()
