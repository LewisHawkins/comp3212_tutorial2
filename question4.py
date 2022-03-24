# Tutorial 2, Question 4: Generating CG rich regions with an HMM
import random

# Build up a sequence
sequence = ""

# First there is a coinflip to determine the starting state, "AT" or "CG"
coinflip = random.randint(0, 1)
if coinflip == 0:
    state = "AG"
else:
    state = "CG"

# Loops n times to generate a sequence of n length
n = 100


# Emission probabilities
# A function to return a random codon when the HMM is in the "AT" state
def getNextCodonAT():
    r = random.randint(1, 10000)
    if r <= 2698:           # p(A) = 0.2698
        return "A"
    elif 2698 < r <= 5935:  # Difference of 3237, to reflect p(T) = 0.3237
        return "T"
    elif 5935 < r <= 8015:  # Difference of 2080, to reflect p(C) = 0.2080
        return "C"
    else:                   # The remaining probability is equal to 0.1985
        return "G"


# A function to return a random codon when the HMM is in the "CG" state
def getNextCodonCG():
    r = random.randint(1, 10000)
    if r <= 2459:           # p(A) = 0.2459
        return "A"
    elif 2459 < r <= 4538:  # Difference of 2079, to reflect p(T) = 0.2079
        return "T"
    elif 4538 < r <= 7016:  # Difference of 2478, to reflect p(C) = 0.2478
        return "C"
    else:                   # The remaining probability is equal to 0.2984
        return "G"


# Transition probability
# A function to determine the next state, given the current state of the HMM
def getNextState(currentState):
    r = random.randint(1, 10000)
    if currentState == "AT":
        if r <= 9998:
            # The state stays in "AT" 9998 times out of 10000
            return "AT"
        else:
            return "CG"
    else:
        if r <= 9997:
            # The state stays in "CG" 9997 times out of 10000
            return "CG"
        else:
            return "AT"


if __name__ == "__main__":
    for i in range(n):
        # First, generate the next codon to add to the sequence based on the probabilities
        if state == "AT":
            sequence += getNextCodonAT()
        else:
            sequence += getNextCodonCG()

        # There is a then a probability that the HMM switches to the next state
        state = getNextState(state)

    print(sequence)
