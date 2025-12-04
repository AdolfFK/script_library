def read_fasta(fasta_file):
    """
    从FASTA文件中读取序列并返回一个字典，其中键为ID，值为序列。
    """
    fasta_dict = {}
    current_id = ""
    with open(fasta_file, "r") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:]
                fasta_dict[current_id] = ""
            else:
                fasta_dict[current_id] += line
    return fasta_dict

def write_sorted_fasta(fasta_dict, output_file):
    """
    将按ID排序的FASTA序列写入输出文件。
    """
    with open(output_file, "w") as outfile:
        for seq_id in sorted(fasta_dict.keys()):
            outfile.write(f">{seq_id}\n")
            sequence = fasta_dict[seq_id]
            outfile.write(sequence + "\n")
            # for i in range(0, len(sequence), 60):
            #     outfile.write(sequence[i:i+60] + "\n")

def main():
    import sys

    if len(sys.argv) != 3 or sys.argv[1] == "-h":
        print("Usage: python sort_fasta_by_id.py input.fasta output_sorted.fasta")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_sorted_fasta = sys.argv[2]

    fasta_dict = read_fasta(input_fasta)
    write_sorted_fasta(fasta_dict, output_sorted_fasta)

if __name__ == "__main__":
    main()