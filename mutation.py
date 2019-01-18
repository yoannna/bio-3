print("Hello")

newLine = "\n"
# resultFile = open("res.txt", "w")
resultFile = open("res.txt", "a")

from enum import Enum
class Letter(Enum):
    A = "A"
    T = "T"
    G = "G"
    C = "C"
    none = "_"

class Mutation:
    def __init__(self, previous: Letter, current: Letter):
        self.previous: Letter = previous
        self.current: Letter = current

class Coordinate:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Pointer(Enum):
    Diagonal = Coordinate(-1, -1)
    Up = Coordinate(0, -1)
    Left = Coordinate(-1, 0)
    none = Coordinate(0, 0)

class MatrixField:
    def __init__(self, value = 0):
        self.value = value
        self.pointsTo = []
    def __eq__(self, other):
        return self.value == other.value
    def __lt__(self, other):
        return self.value < other.value
    def __gt__(self, other):
        return self.value > other.value


class Reads:
    def __init__(self):
        self.sequences = []

#region Read refererence
referenceFile = open("reference.fa")
referenceFile.readline()
reference = []
referenceStr = ""
referenceLine = referenceFile.readline()
while referenceLine != "":
    # remove new line
    referenceLine = referenceLine[:-1]
    referenceStr += referenceLine
    for letter in referenceLine:
        reference.append(Letter(letter))
    referenceLine = referenceFile.readline()
#endregion

def ReadReadsFile(fileName):
    #1 some code
    #2 sequence
    #3 +
    #4 quality values
    readsFile = open(fileName)
    readsFile.readline() # line 1
    sequenceLine = readsFile.readline() # line 2
    reads = Reads()
    read = True

    while read:
        sequenceLine = sequenceLine[:-1]
        if sequenceLine not in referenceStr:
            sequence = []
            for letter in sequenceLine:
                sequence.append(Letter(letter))
            reads.sequences.append(sequence)

        readsFile.readline() # line 3
        if readsFile.readline() == "": # line 4
            read = False
        else:
            readsFile.readline() # line 1
            sequenceLine = readsFile.readline() # line 2

    return reads

reads16 = ReadReadsFile("16_reads.fastq")
# reads18 = ReadReadsFile("18_reads.fastq")

from typing import Sequence
def GenerateMatrix(reference, sequence):
    matrix = []
    for i in range(len(sequence) + 1):
        matrix.append([])
        for j in range(len(reference) + 1):
            matrix[i].append(MatrixField())

    for i in range(len(reference)):
        y = i + 1
        for j in range(len(sequence)):
            x = j + 1
            matchScore = 2 if reference[i] == sequence[j] else -1
            diagonalScore = matrix[x-1][y-1].value + matchScore
            upScore = matrix[x][y-1].value - 2
            leftScore = matrix[x-1][y].value - 2

            maxValue = max([diagonalScore, upScore, leftScore])

            field = MatrixField(maxValue)
            if diagonalScore == maxValue:
                field.pointsTo.append(Pointer.Diagonal)
            if upScore == maxValue:
                field.pointsTo.append(Pointer.Left)
            if leftScore == maxValue:
                field.pointsTo.append(Pointer.Up)
            matrix[x][y] = field

    return matrix

def FindPath(matrix: Sequence[Sequence[MatrixField]], reference, sequence):
    sequenceR1 = []
    sequenceR2 = []
    matches = []
    currentY = len(matrix) - 1
    bottomLine = matrix[currentY]
    currentX = len(matrix[0]) - 1
    currentMax = bottomLine[currentX]

    for i in range(len(matrix[0]) - 1):
        if bottomLine[i] > currentMax:
            currentX = i
            currentMax = bottomLine[i]

    currentCoord = Coordinate(currentX, currentY)

    matchingStarted = False
    mutations = []
    while currentCoord.x != 0 and currentCoord.y != 0:
        currentField: MatrixField = matrix[currentCoord.y][currentCoord.x]
        currentPointer: Pointer = currentField.pointsTo[0]
        refLetter = reference[currentCoord.x - 1]
        seqLetter = sequence[currentCoord.y - 1]
        if currentPointer == Pointer.Diagonal:
            sequenceR1.append(refLetter)
            sequenceR2.append(seqLetter)
            isMatch = refLetter == seqLetter
            matches.append(isMatch)
            matchingStarted = True
            if not isMatch:
                mutations.append(Mutation(refLetter, seqLetter))
              
        elif currentPointer == Pointer.Left:
            matchingStarted = True
            sequenceR1.append(refLetter)
            sequenceR2.append(Letter.none)
            matches.append(False)
            mutations.append(Mutation(refLetter, Letter.none))

        elif currentPointer == Pointer.Up:
            sequenceR1.append(Letter.none)
            sequenceR2.append(seqLetter)
            matches.append(False)
            if matchingStarted:
                mutations.append(Mutation(Letter.none, seqLetter))

        currentCoord.x += currentPointer.value.x
        currentCoord.y += currentPointer.value.y

    length = len(sequenceR1)
    for i in range(length):
        value = sequenceR1[length - i - 1].value
        resultFile.write(value)
    resultFile.write(newLine)
    
    score = 0
    for i in range(length):
        if matches[length - i - 1]:
            score += 1
            value = "|"
        else: value = " "
        resultFile.write(value)
    resultFile.write(newLine)

    for i in range(length):
        value = sequenceR2[length - i - 1].value
        resultFile.write(value)
    resultFile.write(newLine)
    resultFile.write("Matches: " + str(score))
    resultFile.write(newLine)
    
    resultFile.write("Mutations: " + str(len(mutations)))
    resultFile.write(newLine)
    for mutation in mutations:
        resultFile.write(mutation.previous.value + " -> " + mutation.current.value)
        resultFile.write(newLine)
    resultFile.write(newLine)
    resultFile.write(newLine)
    return mutations

# sequence1 = [Letter.G, Letter.A , Letter.A, Letter.T, Letter.T, Letter.C, Letter.A, Letter.G, Letter.T, Letter.T, Letter.A, Letter.A, Letter.A]
# sequence2 = [Letter.G, Letter.G, Letter.A, Letter.T, Letter.C, Letter.G, Letter.A, Letter.C]
# matrix = GenerateMatrix(sequence2, sequence1)
# FindPath(matrix, sequence2, sequence1)

ref = [Letter.G, Letter.A , Letter.A, Letter.T, Letter.T, Letter.C, Letter.A, Letter.G, Letter.T, Letter.T, Letter.A, Letter.A, Letter.A]
seq = [Letter.G, Letter.A, Letter.T, Letter.T, Letter.C]
matrix = GenerateMatrix(ref, seq)
FindPath(matrix, ref, seq)

# seq = [Letter.G, Letter.A, Letter.T]
# ref = [Letter.G, Letter.A, Letter.T, Letter.T, Letter.T]
# matrix = GenerateMatrix(ref, seq)
# FindPath(matrix, ref, seq)

for i in range(1):
    matrix2 = GenerateMatrix(reference, reads16.sequences[i])
    FindPath(matrix2, reference, reads16.sequences[i])

# from Bio.Blast import NCBIWWW
# result_handle = NCBIWWW.qblast("blastn", "nt", referenceStr, megablast=True)
# res = result_handle.read()
# f = open("blast.xml", "w")
# f.write(res)
# f.close()

print("Bye")