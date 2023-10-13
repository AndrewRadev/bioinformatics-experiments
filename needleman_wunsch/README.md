## Installation

Install dependencies

``` .sh-session
$ pip install -r requirements.txt
```

Or install them manually: blosum, tabulate, termcolor

## Usage

Run `python3 main.py` to get usage instructions:

```
USAGE: python3 main.py <SEQ1> <SEQ2> <gap>
```

The two sequences are passed along as the first two arguments, the final one is the gap penalty as a positive number.

The output should look like this:

| Forward step | Backward traversal | Final result |
| --- | --- | --- |
| ![Forwards](https://github.com/AndrewRadev/bioinformatics-experiments/assets/124255/02db6d06-cd05-4d1e-80b9-108b4cbd48c2) | ![Backwards](https://github.com/AndrewRadev/bioinformatics-experiments/assets/124255/ac0e4c9a-cf4e-45f0-8c03-4f7e4b5a916a) | ![Final](https://github.com/AndrewRadev/bioinformatics-experiments/assets/124255/a03e4e7d-4ff3-4ef5-bc41-0c22d01b09d9) |
