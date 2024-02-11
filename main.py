import numpy as np


def read_fasta(file):  # read fasta file format
    with open(file, "r") as file:
        lines_without_header = file.readlines()[
                               1:]  # skips the first line, returns variable as a list where each line is a item in
        # the list
        lines_without_header = ''.join(
            lines_without_header)  # joins the list of string with a '' nothing as a seperator. Continuous string.
        FASTA = lines_without_header.strip()  # remove white space
        print(f"Here is your FASTA sequence:\n{FASTA}")
        return FASTA.upper()  # convert all characters to uppercase for consistency

# while True: #exception handling
#     try:
#         file_1 = input("Write the  name of the FASTA file you wish to open to extract the nucleotide/amino acid sequence: ")
#         sequence_1 = read_fasta(file=file_1)
#         break
#     except FileNotFoundError:
#         print(f"File '{file_1}' not found. Please try again.")

#
# while True:
#     try:
#         file_2 = input("Write the  name of the FASTA file you wish to open to extract the nucleotide/amino acid sequence: ")
#         sequence_2 = read_fasta(file=file_2)
#         break
#     except FileNotFoundError:
#         print(f"File '{file_2}' not found. Please try again.")
#

# #test sequences (WORKS on small sequences, big file has formatting issues due to length of it)
sequence_1 = "GCATGCG"
sequence_2 = "GATTACA"

# create matrices:
main_matrix = np.zeros((len(sequence_1) + 1, len(sequence_2) + 1))
match_checker_matrix = np.zeros((len(sequence_1), len(sequence_2)))

# Provide the scores for matches, mismatches and gaps.
match_score = input("Please provide the match score based on the appropriate scoring matrices: ")
mismatch_score = input("Please provide the mismatch score based on the appropriate scoring matrices: ")
gap_score = input("Please provide the gap score based on the appropriate scoring matrices: ")

# Upload the scoring matrix file to be used for assigning scores for matching nucleotides or amino acids
# while True:
#     try:
#         scoring_file = input("Please provide the file name/file path to be uploaded for assigning scores: ")
#         # Read and parse the substitution matrix file
#         matrix_file = open(scoring_file, "r")
#         lines = matrix_file.read().splitlines()
#         matrix_file.close()
#         break
#     except FileNotFoundError:
#         print(f"File '{scoring_file}' not found. Please try again.")




# Initialize a dictionary to store the substitution matrix scoring file uploaded from user
# substitution_matrix = {}
# header = lines[0].split("\t")[1:]
#
# for line in lines[1:]:
#     parts = line.split("\t")
#     row_label = parts[0]
#     values = [int(val) for val in parts[1:] if val != '']
#     substitution_matrix[row_label] = dict(zip(header, values))
#
# print(substitution_matrix)

# Ask user to input gap score
while True:
    try:
        gap_score = int(input("Please provide the gap score for scoring purposes: "))
        break
    except ValueError:
        print("Invalid input. Please enter an integer and try again.")

# Fill the match checker matrix according to match or mismatch
for i in range(len(sequence_1)):
    for j in range(len(sequence_2)):
        char1 = sequence_1[i]
        char2 = sequence_2[j]

        if char1 in substitution_matrix and char2 in substitution_matrix[char1]:
            match_checker_matrix[i][j] = substitution_matrix[char1][char2]

# print(match_checker_matrix)

# Filling up the matrix using Needleman_Wunsch algorithm
# STEP 1: Initialization (first row and first column)
for i in range(len(sequence_1) + 1):
    main_matrix[i][0] = i * gap_score
for j in range(len(sequence_2) + 1):
    main_matrix[0][j] = j * gap_score

# STEP 2: Matrix Filling
for i in range(1, len(sequence_1) + 1):  # starts from 1 because doesn't use the gap_score columns
    for j in range(1, len(sequence_2) + 1):
        main_matrix[i][j] = max(main_matrix[i - 1][j - 1] + match_checker_matrix[i - 1][j - 1],
                                # diagonal/checks for matches
                                main_matrix[i - 1][j] + gap_score,  # left
                                main_matrix[i][j - 1] + gap_score)  # top

print(main_matrix)

# STEP 3: Traceback
aligned_1 = ""
aligned_2 = ""

ti = len(sequence_1)
tj = len(sequence_2)

while ti > 0 or tj > 0:
    if (ti > 0 and tj > 0 and main_matrix[ti][tj] == main_matrix[ti - 1][tj - 1] + match_checker_matrix[ti - 1][
        tj - 1]):
        aligned_1 = sequence_1[ti - 1] + aligned_1
        aligned_2 = sequence_2[tj - 1] + aligned_2

        ti = ti - 1
        tj = tj - 1

    elif ti > 0 and main_matrix[ti][tj] == main_matrix[ti - 1][tj] + gap_score:
        aligned_1 = sequence_1[ti - 1] + aligned_1
        aligned_2 = "-" + aligned_2
        ti = ti - 1

    else:
        aligned_1 = "-" + aligned_1
        aligned_2 = sequence_2[tj - 1] + aligned_2
        tj = tj - 1

print(aligned_1)
print(aligned_2)


# print(len(aligned_1))
# print(len(aligned_2))


# assumes both sequences are of same length (post-alignment)
def print_align(seq1, seq2, length):
    """

    Returns the output of the alignment in a cleaner format for when sequences are long lengths.

    """

    result = ""

    while len(seq1) > 0:
        chunk1 = seq1[:length]  # divides the sequence into predetermined chunks of chosen length
        chunk2 = seq2[:length]
        match_line = ''.join('|' if c1 == c2 else ' ' for c1, c2 in zip(chunk1, chunk2))  # adds a line between
        # matching pairs of the sequences

        result += f"seq1: {chunk1}\n"
        result += f"      {match_line}\n"
        result += f"seq2: {chunk2}\n"
        result += f"\n"

        seq1 = seq1[length:]
        seq2 = seq2[length:]

    return result


optimized = print_align(seq1=aligned_1, seq2=aligned_2, length=38)  # call print_align with chosen length
print(optimized)

# Save optimized alignment to a text file
file_path_alignment = "optimized_alignment.txt"

with open(file_path_alignment, 'w') as text_file:  # use with as this closes the file automatically
    text_file.write(optimized)
