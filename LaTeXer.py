import copy

def formatLatexMatrix(matrix):
	result = "$\n\\begin{bmatrix}\n"
	for row in matrix:
		result += "\t"
		for entry in row:
			result += str(entry) + " & "
		result += "\\\\" + "\n"
	result += "\\end{bmatrix}\n$\\\\\\\\"
	return result

def switchRows(matrix, r1, r2):
	r1 -= 1
	r2 -= 1
	assert(r1 < len(matrix) and r2 < len(matrix))
	temp = []
	for value in matrix[r1]:
		temp.append(value)
	matrix[r1] = matrix[r2]
	matrix[r2] = temp

		
def addRows(matrix, resRow, row1ToAdd, scalarForRow1, row2ToAdd, scalarForRow2):
	resRow -= 1
	row1ToAdd -=1
	row2ToAdd -= 1
	assert(resRow < len(matrix) and row1ToAdd < len(matrix) and row2ToAdd < len(matrix))
	sumRow1 = []
	sumRow2 = []
	for entry in matrix[row1ToAdd]:
		sumRow1.append(entry*scalarForRow1)
	for entry in matrix[row2ToAdd]:
		sumRow2.append(entry*scalarForRow2)
	sumRow = []
	assert(len(sumRow1) == len(sumRow2))
	for i in range(len(sumRow1)):
		sumRow.append(sumRow1[i] + sumRow2[i])
	matrix[resRow] = sumRow


class ElementaryRowOperation:
	def __init__(self, type, r1, r2, resRow, s1, s2):
		self.type = type
		self.r1 = r1
		self.r2 = r2
		self.resRow = resRow
		self.s1 = s1
		self.s2 = s2
	
	def perform(self, M):
		if self.type == "SWAP":
			switchRows(M, self.r1 , self.r2)
		else:
			addRows(M, self.resRow, self.r1 , self.s1 , self.r2 , self.s2 )
	
	def getOperationLatexCode(self):
		if self.type == "SWAP":
			return "swap $R_{" + str(self.r1) + "}$ with $R_{" + str(self.r2) + "}$\\\\\\\\\n"
		else:
			return "let $R_{" + str(self.resRow) + "new} = " + str(self.s1) + "R_{" + str(self.r1) + "} + " + str(self.s2) + "R_{" + str(self.r2) + "}$\\\\\\\\\n"


class Operations:
	def __init__(self, matrix):
		self.matrix = matrix
		self.alteredMatrix = copy.deepcopy(matrix)
		self.operations = []
		self.generatedLatex = ""
	
	def addOperation(self, op):
		self.operations.append(op)
		op.perform(self.alteredMatrix)
	
	
	def carryOutOpsAndGetLaTex(self):
		self.generatedLatex += "start: \n"
		self.generatedLatex += formatLatexMatrix(self.matrix)
		matrixForPrinting = copy.deepcopy(self.matrix)
		for op in self.operations:
			self.generatedLatex += "do: " + op.getOperationLatexCode()
			op.perform(matrixForPrinting)
			self.generatedLatex += "we get: " + formatLatexMatrix(matrixForPrinting)
		return self.generatedLatex
	
	#return operations that undo all the operations that this operations object has
	def getInverseOperations(self):
		inverseOperations = Operations(self.alteredMatrix)
		for op in self.operations[::-1]:
			if op.type == "SWAP":
				inverseOperations.addOperation(copy.deepcopy(op))
			elif  op.resRow == op.r1:
				inverseOperations.addOperation(ElementaryRowOperation("ADD", op.resRow, op.r2, op.resRow, 1, -op.s2))
				inverseOperations.addOperation(ElementaryRowOperation("ADD", op.resRow, 1, op.resRow, 1/op.s1, 0))
			else:
				assert(op.resRow == op.r2)
				inverseOperations.addOperation(ElementaryRowOperation("ADD", op.resRow, op.r1, op.resRow, 1, -op.s1))
				inverseOperations.addOperation(ElementaryRowOperation("ADD", op.resRow, 1, op.resRow, 1/op.s2, 0))	

		return inverseOperations


#return the operations to put matrix M in row echelon form (and put the matrix in row echelon form)
def reducedRowEchelon(M):
	highestRowWithPossiblePivot = 0
	res = Operations(M)
	pivotRows = []
	#A) ROW ECHELON
	#return a list of elementary row operations to:
	#	-1. find a possible pivot that lies underneath highestRowWithPossiblePivot (in col)
	#	-2. swap the row that contains the pivot with highestRowWithPossiblePivot (in col)
	#	-3. use the pivot to set zeros in all rows below highestRowwithPossiblePivot (in col)
	#	-4. set highestRowWithPossiblePivot to one more than the previous (If step one succeeded), otherwise just don't change anything
	for col in range(len(M[0])):
		for rowIndex in range(len(M)):
			if rowIndex >= highestRowWithPossiblePivot and res.alteredMatrix[rowIndex][col] != 0: #found a pivot
				pivotRows.append(rowIndex)
				if rowIndex != highestRowWithPossiblePivot:
					res.addOperation(ElementaryRowOperation("SWAP", rowIndex + 1, highestRowWithPossiblePivot + 1, 0, 0, 0)) #2
				rowToWipeCol = highestRowWithPossiblePivot + 1
				while rowToWipeCol < len(M):
					valueAtPivot = res.alteredMatrix[highestRowWithPossiblePivot][col]
					valueAtWipePoint = res.alteredMatrix[rowToWipeCol][col]
					res.addOperation(ElementaryRowOperation("ADD", rowToWipeCol + 1, highestRowWithPossiblePivot + 1, rowToWipeCol + 1, -valueAtPivot, valueAtWipePoint))
					rowToWipeCol += 1
				highestRowWithPossiblePivot += 1
				break	
	#B) REDUCED ROW ECHELON
	#	-1. divide the pivot row in the col by the pivot value
	#	-2. for all values above the pivot, subtract by the pivot to cancel those terms
	for col in range(len(M[0])):
		if col < len(pivotRows):
			pivotRow = pivotRows[col]
			pivotValue = res.alteredMatrix[pivotRow][col]
			if pivotValue != 0:
				res.addOperation(ElementaryRowOperation("ADD", pivotRow + 1, 1, pivotRow + 1, 1/pivotValue, 0))
				pivotValue = 1
				goBackUpRow = pivotRow - 1
				while goBackUpRow >= 0:
					if res.alteredMatrix[goBackUpRow][col] != 0:
						numberToGetRidOf = res.alteredMatrix[goBackUpRow][col]
						res.addOperation(ElementaryRowOperation("ADD", goBackUpRow + 1, pivotRow + 1, goBackUpRow + 1, pivotValue, -numberToGetRidOf))
					goBackUpRow -= 1			
	return res



#USAGE:

foo = reducedRowEchelon([[-1, 2/3, 0, 0, 0, -1], [1/3, -2/3, 1/3, 0, 0, -1], [0, 1/3, -2/3, 1/3, 0, -1], [0, 0, 2/3, -1, 1/3, -1], [0, 0, 0, 1, -1, -1]]) #enter your matrix here!
foo.carryOutOpsAndGetLaTex()
print(foo.generatedLatex)

