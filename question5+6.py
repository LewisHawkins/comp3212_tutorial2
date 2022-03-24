# Tutorial 2, Questions 5 and 6: Writing and testing the Viterbi algorithm, finding the most likely "CG" regions
import math

# Get the genome from the file, and concatenate it into a single string
o = open("phaseLambda.fasta", "r")
genome = "".join(o.readlines()).replace("\n", "")
# genome = "GGGCGGCGACCT. . . "
# print(genome)

# Create the various matrices needed for the viterbi algorithm
# The transition matrix, giving the probabilities of moving from one state to another
transitionMatrix = [[0.9998, 0.0002],
                    [0.0003, 0.9997]]
# The emission matrix, giving the probabilities of each codon being emitted in each state
emissionMatrix = [[0.2698, 0.3237, 0.2080, 0.1985],
                  [0.2459, 0.2079, 0.2478, 0.2984]]
# The initial matrix, giving the probabilities of beginning in each of the states
initialMatrix = [0.5, 0.5]


def getTransitionProbability(fromState, toState):
    if fromState == "AT":
        if toState == "AT":
            return transitionMatrix[0][0]  # AT -> AT
        else:
            return transitionMatrix[0][1]  # AT -> CG
    else:
        if toState == "AT":
            return transitionMatrix[1][0]  # CG -> AT
        else:
            return transitionMatrix[1][1]  # CG -> CG


def getEmissionProbability(state, emission):
    if state == "AT":
        if emission == "A":
            return emissionMatrix[0][0]
        elif emission == "T":
            return emissionMatrix[0][1]
        elif emission == "C":
            return emissionMatrix[0][2]
        else:  # "G"
            return emissionMatrix[0][3]
    else:  # "CG" state
        if emission == "A":
            return emissionMatrix[1][0]
        elif emission == "T":
            return emissionMatrix[1][1]
        elif emission == "C":
            return emissionMatrix[1][2]
        else:  # "G"
            return emissionMatrix[1][3]


def getInitialProbability(state):
    # For the given example, there's an equal chance of starting in either of the two states
    return 0.5


def viterbi(genome):
    # The two lists to store the calculated probabilities as the algorithms works through the genome
    viterbiATprobs = []
    viterbiCGprobs = []
    blueArrows = []

    # The algorithm must first work out the probability of beginning in either state.
    # This is equal to the probability of starting in that state, multiplied by the probability of
    # the first codon in the genome (in this example, G) being emitted given that the HMM was in that state
    nextATstate = math.log(getInitialProbability("AT") * getEmissionProbability("AT", genome[0]))
    viterbiATprobs.append(nextATstate)
    nextCGstate = math.log(getInitialProbability("CG") * getEmissionProbability("CG", genome[0]))
    viterbiCGprobs.append(nextCGstate)

    # The algorithm then works its way forward through the genome
    for nextCodon in genome[1:]:
        # Calculate the 4 paths
        ATtoAT = nextATstate + math.log(getTransitionProbability("AT", "AT"))
        CGtoAT = nextCGstate + math.log(getTransitionProbability("CG", "AT")) 

        ATtoCG = nextATstate + math.log(getTransitionProbability("AT", "CG"))
        CGtoCG = nextCGstate + math.log(getTransitionProbability("CG", "CG"))

        # Add the emission probability to the greater (max) of the two paths
        nextATstate = math.log(getEmissionProbability("AT", nextCodon)) + max(ATtoAT, CGtoAT)
        nextCGstate = math.log(getEmissionProbability("CG", nextCodon)) + max(ATtoCG, CGtoCG)

        # Add the results to the storage
        viterbiATprobs.append(nextATstate)
        viterbiCGprobs.append(nextCGstate)

        # Also keep the direction of these better paths now that they've been found (the blue arrows on the slides)
        # This is stored as a list of pairs, the first pair being:
        # (the better of the two initial states to go to v1(AT) from, the better state to go to v1(CG) from)
        newArrowPair = ["", ""]
        if ATtoAT > CGtoAT:
            newArrowPair[0] = "AT"
        else:
            newArrowPair[0] = "CG"
        if ATtoCG > CGtoCG:
            newArrowPair[1] = "AT"
        else:
            newArrowPair[1] = "CG"
        blueArrows.append(newArrowPair)

    # Once the loop has gone forward through the genome, we can compare the two end values
    # and use the blue arrows to trace back through the codon to the state
    listOfStates = []

    if viterbiATprobs[-1] > viterbiCGprobs[-1]:
        tracebackState = "AT"
        listOfStates.append("AT")
    else:
        tracebackState = "CG"
        listOfStates.append("CG")
 
    # Reverse the blue arrows list
    blueArrows.reverse()

    for arrowPath in blueArrows:
        if tracebackState == "AT":
            listOfStates.append(arrowPath[0])
            tracebackState = arrowPath[0]
        else:
            listOfStates.append(arrowPath[1])
            tracebackState = arrowPath[1]

    # Reverse the list so that it lines up with the genome, forwards
    listOfStates.reverse()

    return listOfStates


if __name__ == "__main__":
    # Run the viterbi algorithm on the supplied genome
    likelyStates = viterbi(genome)
    # Output the full sequence of states
    # print(likelyStates)
    # Output the length of sequence of states
    print("Total states: " + str(len(likelyStates)))
    # Output how many AT / CG states there are
    # print("Total AT states: " + str(likelyStates.count("AT")))
    # print("Total CG states: " + str(likelyStates.count("CG")))

    # Calculate the length of each of the regions
    sequenceCount = 0
    previousState = "AT"
    regions = []
    for state in likelyStates:
        sequenceCount += 1
        if state != previousState:
            if previousState != "":
                regions.append((previousState, sequenceCount))
                previousState = state
                sequenceCount = 0
    regions.append((previousState, sequenceCount))

    print("The order and length of the alternating regions is shown below: ")
    print(regions)
