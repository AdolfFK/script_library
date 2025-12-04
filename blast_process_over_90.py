import sys

def print_usage():
    usage = """
===================== process_blast_over_90.py =====================
                            Version 1.0

Usage:   python process_blast_over_90.py <blast_output_file> <qualified_fasta_file>
Example: python process_blast_over_90.py RSV-ref-A-F.blast.out RSV-ref-A-F.blast.extraction.fasta
Description:
- This program is used to filter blast results.
- It takes a blast output file and a qualified fasta file as input, and outputs the filtered fasta sequences.
"""
    print(usage)

def process_blast(blast_file, fasta_file):
    blast_output = blast_file
    fasta_output = fasta_file

    with open(blast_output, 'r') as blast_in, open(fasta_output, 'w') as fasta_out:
        n = 1
        lenths = 0
        for line in blast_in:
            line = line.strip()
            temp = line.split('\t')
            if n == 1:
                lenths = float(temp[3])
            linelen = float(temp[3])
            lencutoff = lenths * 0.9
            similarity = float(temp[2])
            if similarity >= 0.9 and linelen >= lencutoff:
                fasta_out.write(f">{temp[1]}\n{temp[6]}\n")
            n += 1

def main():
    if len(sys.argv) != 3 or sys.argv[1] == "-h":
        print_usage()
        sys.exit(1)

    blast_file, fasta_file = sys.argv[1], sys.argv[2]
    process_blast(blast_file, fasta_file)

if __name__ == "__main__":
    main()
