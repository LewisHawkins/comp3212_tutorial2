# Tutorial 2, Question 2: Extension of Needleman-Wunsch, implements Smith-Waterman
# The blosum50 matrix is supplied below
blosum50 = [
    [5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0],
    [-2, 7, -1, -2, -4, 1, 0, -3, 0, -4, -3, 3, -2, -3, -3, -1, -1, -3, -1, -3],
    [-1, -1, 7, 2, -2, 0, 0, 0, 1, -3, -4, 0, -2, -4, -2, 1, 0, -4, -2, -3],
    [-2, -2, 2, 8, -4, 0, 2, -1, -1, -4, -4, -1, -4, -5, -1, 0, -1, -5, -3, -4],
    [-1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1],
    [-1, 1, 0, 0, -3, 7, 2, -2, 1, -3, -2, 2, 0, -4, -1, 0, -1, -1, -1, -3],
    [-1, 0, 0, 2, -3, 2, 6, -3, 0, -4, -3, 1, -2, -3, -1, -1, -1, -3, -2, -3],
    [0, -3, 0, -1, -3, -2, -3, 8, -2, -4, -4, -2, -3, -4, -2, 0, -2, -3, -3, -4],
    [-2, 0, 1, -1, -3, 1, 0, -2, 10, -4, -3, 0, -1, -1, -2, -1, -2, -3, 2, -4],
    [-1, -4, -3, -4, -2, -3, -4, -4, -4, 5, 2, -3, 2, 0, -3, -3, -1, -3, -1, 4],
    [-2, -3, -4, -4, -2, -2, -3, -4, -3, 2, 5, -3, 3, 1, -4, -3, -1, -2, -1, 1],
    [-1, 3, 0, -1, -3, 2, 1, -2, 0, -3, -3, 6, -2, -4, -1, 0, -1, -3, -2, -3],
    [-1, -2, -2, -4, -2, 0, -2, -3, -1, 2, 3, -2, 7, 0, -3, -2, -1, -1, 0, 1],
    [-3, -3, -4, -5, -2, -4, -3, -4, -1, 0, 1, -4, 0, 8, -4, -3, -2, 1, 4, -1],
    [-1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3],
    [1, -1, 1, 0, -1, 0, -1, 0, -1, -3, -3, 0, -2, -3, -1, 5, 2, -4, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 2, 5, -3, -2, 0],
    [-3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1, 1, -4, -4, -3, 15, 2, -3],
    [-2, -1, -2, -3, -3, -1, -2, -3, 2, -1, -1, -2, 0, 4, -3, -2, -2, 2, 8, -1],
    [0, -3, -3, -4, -1, -3, -3, -4, -4, 4, 1, -3, 1, -1, -3, -2, 0, -3, -1, 5]]


# A helper method to retrieve the substitution cost of protein1 into protein2 from the blosum50 matrix
def getSubstitutionCost(protein1, protein2):
    # Convert the letters into their index in the matrix
    indexDic = {"A": 0, "R": 1, "N": 2, "D": 3, "C": 4,
                "Q": 5, "E": 6, "G": 7, "H": 8, "I": 9,
                "L": 10, "K": 11, "M": 12, "F": 13, "P": 14,
                "S": 15, "T": 16, "W": 17, "Y": 18, "V": 19}
    # ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    protein1index = indexDic[protein1]
    protein2index = indexDic[protein2]
    # Return the number found in that cell in the matrix
    return blosum50[protein2index][protein1index]


# Forwards algorithm, build up the matrix F(i, j)
def forwards(stringA, stringB):
    # Create an empty matrix full of zeroes, with the dimension of the two words + 1
    fMatrix = [[0 for i in range(len(stringA) + 1)] for j in range(len(stringB) + 1)]
    # Also create an empty matrix to record the arrows
    dMatrix = [["  " for i in range(len(stringA) + 1)] for j in range(len(stringB) + 1)]

    # Assume linear insertion/deletion cost of 8
    d = 8
    # Complete the matrix iteratively
    for i in range(len(fMatrix)):
        for j in range(len(fMatrix[0])):
            # Handle the first row
            if i == 0:
                fMatrix[i][j] = -1 * j * d
                dMatrix[i][j] = "--"
            # Handle the first column
            elif j == 0:
                fMatrix[i][j] = -1 * i * d
                dMatrix[i][j] = "| "
            else:
                # Compare the three different options
                option1 = fMatrix[i-1][j-1] + getSubstitutionCost(stringA[j-1], stringB[i-1])
                option2 = fMatrix[i-1][j] - d
                option3 = fMatrix[i][j-1] - d
                # Smith-Waterman now has a fourth option, to ignore the match
                option4 = 0

                # Insert the maximum of the three into the matrix
                maxValue = max(option1, option2, option3, option4)
                fMatrix[i][j] = maxValue

                # Also insert the direction into the matrix for the arrows
                if maxValue == option1:
                    dMatrix[i][j] = "\\"
                elif maxValue == option2 or j == 0:
                    dMatrix[i][j] = "| "
                elif maxValue == option3 or i == 0:
                    dMatrix[i][j] = "--"
                else:
                    # Ignore the match
                    dMatrix[i][j] = "  "

    dMatrix[0][0] = "  "

    # Forwards algorithm is now complete, return the matrix
    return fMatrix, dMatrix


# Backward algorithm, find the minimum cost / edit distance
def backwards(fMatrix, dMatrix, stringA, stringB):
    alignedStringA = ""
    alignedStringB = ""

    # Find the cell with the greatest value, and start there
    greatestValue = 0
    i = 0
    j = 0
    for p in range(len(fMatrix)):
        for q in range(len(fMatrix[0])):
            if fMatrix[p][q] > greatestValue:
                greatestValue = fMatrix[p][q]
                i = p
                j = q

    # Then follow the arrows until reaching a 0
    # If the cell is diagonal ( \\ ) then match
    # If vertical arrow ( | ) then gap in B
    # If horizontal arrow  ( -- ) then gap in A

    while i + j > 0:
        if dMatrix[i][j] == "\\":
            # Match the strings, by adding a character to the end of the output for both strings
            alignedStringA = stringA[j-1] + alignedStringA
            alignedStringB = stringB[i-1] + alignedStringB
            # Move to the next cell diagonally
            i -= 1
            j -= 1

        elif dMatrix[i][j] == "| ":
            # Represent a gap in B by a dash ( - ) in string A
            alignedStringA = "-" + alignedStringA
            alignedStringB = stringB[i-1] + alignedStringB
            # Move to the next cell vertically
            i -= 1

        elif dMatrix[i][j] == "--":
            # Represent a gap in A by a dash ( - ) in string B
            alignedStringB = "-" + alignedStringB
            alignedStringA = stringA[j-1] + alignedStringA
            # Move to the next cell horizontally
            j -= 1

        else:
            # Match is ignored
            i = 0
            j = 0

    return alignedStringA, alignedStringB


# Run the full algorithm
def runNeedlemanWunsch():
    # stringA = "HEAGAWGHEE"
    # stringB = "PAWHEAE"
    stringA = "MQNSHSGVNQLGGVFVNGRPLPDSTRQKIVELAHSGARPCDISRILQVSNGCVSKILGRY"
    stringB = "TDDECHSGVNQLGGVFVGGRPLPDSTRQKIVELAHSGARPCDISRI"
    print("Input String A:   " + stringA)
    print("Input String B:   " + stringB)
    print()

    fMatrix, dMatrix = forwards(stringA, stringB)
    # print("Forwards:")
    # for row in fMatrix:
    #     print(row)
    # print()
    # for row in dMatrix:
    #     print(row)

    alignedStringA, alignedStringB = backwards(fMatrix, dMatrix, stringA, stringB)
    print("Aligned string A: " + alignedStringA)
    print("Aligned string B: " + alignedStringB)


# Run an example
if __name__ == "__main__":
    runNeedlemanWunsch()
