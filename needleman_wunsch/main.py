import sys
import copy
import blosum
from printy import raw_format
from tabulate import tabulate

if len(sys.argv) < 4:
    print("USAGE: python3 main.py <SEQ1> <SEQ2> <gap>")
    exit(1)

matrix = blosum.BLOSUM(62)

seq_col = sys.argv[1].upper()
seq_row = sys.argv[2].upper()
gap = int(sys.argv[3])

# Should be a command-line flag:
show_steps = True

for char in [*seq_row, *seq_col]:
    if char not in matrix:
        print(f"Character {char} not a valid protein base")
        exit(1)

# Initialize matrix with numbers:
result_matrix = [
        [None, None, *list(seq_row)],
        [None] + [i * -gap for i in range(len(seq_col) + 2)],
    ]
for i, char in enumerate(seq_col):
    empty_row = [char, (i + 1) * -gap] + [None for _ in range(len(seq_col) + 1)]
    result_matrix.append(empty_row)

# Initialize matrix with directions we've taken so we can go back along them later:
direction_matrix = [
        [set(), set(), *list(seq_row)],
        [None] + [{'→'} for _ in range(len(seq_col) + 2)],
    ]
for i, char in enumerate(seq_col):
    empty_row = [char, {'↓'}] + [set() for _ in range(len(seq_col) + 1)]
    direction_matrix.append(empty_row)

# Special case: initial 0 is the starting point of all paths:
direction_matrix[1][1] = {'→', '↓', '↘'}

if show_steps:
    print("Initial state:")
    formatted_matrix = copy.deepcopy(result_matrix)
    for col_index in range(1, len(seq_col) + 2):
        formatted_matrix[col_index][1] = raw_format(str(result_matrix[col_index][1]), 'r')
    for row_index in range(1, len(seq_row) + 2):
        formatted_matrix[1][row_index] = raw_format(str(result_matrix[1][row_index]), 'r')
    print(tabulate(formatted_matrix, tablefmt="grid"))
    input("Press Enter to continue")

step_count = 1

# Loop over the content of the matrix while skipping row/column 0 with the
# sequences and row/column 1 with the initial gap values:
for row_index in range(2, len(seq_row) + 2):
    for col_index in range(2, len(seq_col) + 2):
        char_row = seq_row[row_index - 2]
        char_col = seq_col[col_index - 2]

        value = int(matrix[char_row][char_col])

        diag_value = result_matrix[col_index - 1][row_index - 1]
        left_value = result_matrix[col_index    ][row_index - 1]
        top_value  = result_matrix[col_index - 1][row_index    ]

        diag_candidate = diag_value + value
        left_candidate = left_value - gap
        top_candidate  = top_value - gap

        best_candidate = max([diag_candidate, top_candidate, left_candidate])
        result_matrix[col_index][row_index] = best_candidate

        if best_candidate == diag_candidate:
            direction_matrix[col_index - 1][row_index - 1].add('↘')
        elif best_candidate == left_candidate:
            direction_matrix[col_index    ][row_index - 1].add('→')
        else:
            direction_matrix[col_index - 1][row_index    ].add('↓')

        if show_steps:
            print(f"\nStep {step_count}")
            print(f"Diagonal candidate: {diag_candidate} = {diag_value} + {value} (value)")
            print(f"Left candidate:     {left_candidate} = {left_value} - {gap} (gap)")
            print(f"Top candidate:      {top_candidate} = {top_value} - {gap} (gap)")
            print("--------------------")
            print(f"Best:               {best_candidate}")

            formatted_matrix = copy.deepcopy(result_matrix)
            formatted_matrix[col_index - 1][row_index - 1] = raw_format(result_matrix[col_index - 1][row_index - 1], 'r')
            formatted_matrix[col_index - 1][row_index    ] = raw_format(result_matrix[col_index - 1][row_index    ], 'r')
            formatted_matrix[col_index    ][row_index - 1] = raw_format(result_matrix[col_index    ][row_index - 1], 'r')
            print(tabulate(formatted_matrix, tablefmt="grid"))
            input("Press Enter to continue")
            step_count += 1

print("Filled in:")
print(tabulate(result_matrix, tablefmt="grid"))
print("Paths taken:")
print(tabulate(direction_matrix, tablefmt="grid"))

col_index = len(seq_col) + 1
row_index = len(seq_row) + 1

aligned_row = ""
aligned_col = ""

while row_index > 1 or col_index > 1:
    formatted_matrix = copy.deepcopy(result_matrix)
    best_candidate = None

    if row_index > 1 and col_index > 1 and '↘' in direction_matrix[col_index - 1][row_index - 1]:
        best_candidate = result_matrix[col_index - 1][row_index - 1]
        formatted_matrix[col_index - 1][row_index - 1] = raw_format(result_matrix[col_index - 1][row_index - 1], 'r')
        print(f"\nCame from the diagonal: {best_candidate}")
        aligned_row = seq_row[row_index - 2] + aligned_row
        aligned_col = seq_col[col_index - 2] + aligned_col
        row_index -= 1
        col_index -= 1

    if row_index > 1 and best_candidate is None and '→' in direction_matrix[col_index][row_index - 1]:
        best_candidate = result_matrix[col_index][row_index - 1]
        formatted_matrix[col_index][row_index - 1] = raw_format(result_matrix[col_index][row_index - 1], 'r')
        print(f"\nCame from the left: {best_candidate}")
        aligned_row = seq_row[row_index - 2] + aligned_row
        aligned_col = '-' + aligned_col
        row_index -= 1

    if col_index > 1 and best_candidate is None and '↓' in direction_matrix[col_index - 1][row_index]:
        best_candidate = result_matrix[col_index - 1][row_index]
        formatted_matrix[col_index - 1][row_index] = raw_format(result_matrix[col_index - 1][row_index], 'r')
        print(f"\nCame from the top: {best_candidate}")
        aligned_col = seq_col[col_index - 2] + aligned_col
        aligned_row = '-' + aligned_row
        col_index -= 1

    if best_candidate is None:
        print(tabulate(direction_matrix, tablefmt="grid"))
        raise Exception("Couldn't find a candidate to go in reverse")

    if show_steps:
        print(tabulate(formatted_matrix, tablefmt="grid"))
        print("Result so far:")
        print(aligned_row)
        print(aligned_col)
        input("Press Enter to continue")

print("\n\nFinal result:")
print(aligned_row)
print(aligned_col)
