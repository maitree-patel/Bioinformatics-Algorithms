# Dynamic Programming for pairwise sequence alignment
# Main goal: to find the optimal alignment out of all possible alignments

# The algorithm will follow three steps
# 1.Initialization
# 2.Main iteration
# 3.Termination and traceback to find optimal solution

# Loading libraries
import numpy as np
import pandas as pd

# pairwiseAlign: a custom function to perform pariwise alignment of two protein or DNA sequences
# Parameters include:
    # seqX: first sequence
    # seqY: second sequence
    # match: score for a base or amino acid match between the two sequences
    # mismatch: mismatch score for the alignment
    # gap: penalty score for a gap in the alignment
def pairwiseAlign(seqX, seqY, match, mismatch, gap):
    # 1.Initialization
    # Setting up an empty matrix
    rows = len(seqX)+1
    cols = len(seqY)+1

    matrix = np.empty((rows, cols))
    rowNames = [nt for nt in seqX]
    rowNames.insert(0, '*')

    colNames = [nt for nt in seqY]
    colNames.insert(0, '*')

    matrixDf = pd.DataFrame(matrix, 
                            index=rowNames, 
                            columns=colNames)

    # Filling the top left cell with 0
    matrixDf.iloc[0,0] = 0

    # Filling the first row with 0s
    for i in range(1,len(seqX)+1):
        matrixDf.iloc[i,0] = 0

    # Filling the first column with mismatch rows
    for j in range(1,len(seqY)+1):
        matrixDf.iloc[0,j] = 0

    globalMax = 0
    count = 0

    # 2.Main Iteration
    for i in range(1,len(seqX)+1):
        for j in range(1,len(seqY)+1):
            # Calculating diagonal, horizontal and vertical for each cell
            # Case 1: Diagonal
            if rowNames[i]==colNames[j]:
                diagonal = matrixDf.iloc[i-1,j-1]+match
            else:
                diagonal = matrixDf.iloc[i-1,j-1]+mismatch
            
            # Case2 2: Horizontal
            horizontal = matrixDf.iloc[i,j-1]+gap

            # Case 3: Vertical
            vertical = matrixDf.iloc[i-1,j]+gap

            # Calculating maximum among the 3 cases and a case 4 i.e. 0
            maxScore = max(diagonal, horizontal, vertical, 0)
            
            if maxScore > globalMax:
                globalMax = maxScore
            
            # Setting maximum score as the value for the cell
            matrixDf.iloc[i,j] = maxScore

    print(matrixDf)

    print(f"The optimal local alignment score is {globalMax}.")

    print(f"The optimal score of {globalMax} occurs in {(matrixDf == 14).sum().sum()} times.")

# Applying the function to two test sequences
seqA = "CTTAAGTCAAT"
seqB = "CCTAGCT"

matchScore = 5
mismatchScore = -3
gapPenalty = -4

pairwiseAlign(seqA, seqB, match=matchScore, mismatch=mismatchScore, gap=gapPenalty)

