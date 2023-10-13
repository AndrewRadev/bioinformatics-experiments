import sys
import copy
import blosum
from termcolor import colored
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

# Initialize matrix with directions taken:
direction_matrix = [
        [[], [], *list(seq_row)],
        [None] + [['→'] for _ in range(len(seq_col) + 2)],
    ]
for i, char in enumerate(seq_col):
    empty_row = [char, ['↓']] + [[] for _ in range(len(seq_col) + 1)]
    direction_matrix.append(empty_row)

# Special case: initial 0
direction_matrix[1][1] = ['→', '↓', '↘']

if show_steps:
    print("Initial state:")
    formatted_matrix = copy.deepcopy(result_matrix)
    for col_index in range(1, len(seq_col) + 2):
        formatted_matrix[col_index][1] = colored(str(result_matrix[col_index][1]), "red")
    for row_index in range(1, len(seq_row) + 2):
        formatted_matrix[1][row_index] = colored(str(result_matrix[1][row_index]), "red")
    print(tabulate(formatted_matrix, tablefmt="grid"))
    input("Press Enter to continue")

step_count = 1

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
            direction_matrix[col_index - 1][row_index - 1].append('↘')
        elif best_candidate == left_candidate:
            direction_matrix[col_index    ][row_index - 1].append('→')
        else:
            direction_matrix[col_index - 1][row_index    ].append('↓')

        if show_steps:
            print(f"\nStep {step_count}")
            print(f"Diagonal candidate: {diag_candidate} = {diag_value} + {value} (value)")
            print(f"Left candidate:     {left_candidate} = {left_value} - {gap} (gap)")
            print(f"Top candidate:      {top_candidate} = {top_value} - {gap} (gap)")
            print("--------------------")
            print(f"Best:               {best_candidate}")

            formatted_matrix = copy.deepcopy(result_matrix)
            formatted_matrix[col_index - 1][row_index - 1] = colored(result_matrix[col_index - 1][row_index - 1], "red")
            formatted_matrix[col_index - 1][row_index    ] = colored(result_matrix[col_index - 1][row_index    ], "red")
            formatted_matrix[col_index    ][row_index - 1] = colored(result_matrix[col_index    ][row_index - 1], "red")
            print(tabulate(formatted_matrix, tablefmt="grid"))
            input("Press Enter to continue")
            step_count += 1

print("Filled in:")
print(tabulate(result_matrix, tablefmt="grid"))

col_index = len(seq_col) + 1
row_index = len(seq_row) + 1

aligned_row = ""
aligned_col = ""

while row_index > 1 or col_index > 1:
    formatted_matrix = copy.deepcopy(result_matrix)
    best_candidate = None

    if row_index > 1 and col_index > 1 and direction_matrix[col_index - 1][row_index - 1].count('↘') > 0:
        best_candidate = result_matrix[col_index - 1][row_index - 1]
        formatted_matrix[col_index - 1][row_index - 1] = colored(result_matrix[col_index - 1][row_index - 1], "red")
        print(f"Came from the diagonal: {best_candidate}")
        aligned_row += seq_row[row_index - 2]
        aligned_col += seq_col[col_index - 2]
        row_index -= 1
        col_index -= 1

    if row_index > 1 and best_candidate is None and direction_matrix[col_index][row_index - 1].count('→') > 0:
        best_candidate = result_matrix[col_index][row_index - 1]
        formatted_matrix[col_index][row_index - 1] = colored(result_matrix[col_index][row_index - 1], "red")
        print(f"Came from the left: {best_candidate}")
        aligned_row += seq_row[row_index - 2]
        aligned_col += '-'
        row_index -= 1

    if col_index > 1 and best_candidate is None and direction_matrix[col_index - 1][row_index].count('↓') > 0:
        best_candidate = result_matrix[col_index - 1][row_index]
        formatted_matrix[col_index - 1][row_index] = colored(result_matrix[col_index - 1][row_index], "red")
        print(f"Came from the top: {best_candidate}")
        aligned_col += seq_col[col_index - 2]
        aligned_row += '-'
        col_index -= 1

    if best_candidate is None:
        print(tabulate(direction_matrix, tablefmt="grid"))
        raise Exception("Couldn't find a candidate to go in reverse")

    if show_steps:
        print(tabulate(formatted_matrix, tablefmt="grid"))
        print("Result so far:")
        print(aligned_row[::-1])
        print(aligned_col[::-1])
        input("Press Enter to continue")

print("Final result:")
print(aligned_row[::-1])
print(aligned_col[::-1])
