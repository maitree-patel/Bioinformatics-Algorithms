The pairwiseAlignment.py script contains a function to calculate an optimal pairwise alignment score between two nucleotide or amino acid sequences. It uses dynami programming to calculate the alignment matrix and do a trace back to find the optimal alignment score.

The algorithm follows 3 main step:
1. Initialization of the matrix
2. Main iteration to fill the matrix using match, mismatch and gap scores
3. Termination and traceback to find optimal alignment score
