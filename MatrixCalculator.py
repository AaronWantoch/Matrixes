import math
from copy import deepcopy

class MatrixCalculator:
    @staticmethod
    def calculatJacobi(A, b):
        iterations=0
        U = MatrixCalculator.getUpperMatrix(A)
        L = MatrixCalculator.getLowerMatrix(A)
        D = MatrixCalculator.getDiagonalMatrix(A)
        inversedD = MatrixCalculator.inverseDiagonal(D)

        factor1 = MatrixCalculator.multiplyMatrixByMatrix(inversedD, L+U)
        factor2 = MatrixCalculator.multiplyMatrixByVector(inversedD, b)

        r = [1]*len(A)
        res = MatrixCalculator.multiplyMatrixByVector(A, r) + -1*b
        norm = MatrixCalculator.vectorNorm(res)
        while norm > math.pow(10, -9):
            r = MatrixCalculator.multiplyMatrixByVector(factor1, r) + factor2
            res = MatrixCalculator.multiplyMatrixByVector(A, r) + -1*b
            norm = MatrixCalculator.vectorNorm(res)

            iterations+=1
        return  r, iterations

    def calculateGaussSeidle(A, b):
        iterations=0
        U = MatrixCalculator.getUpperMatrix(A)
        L = MatrixCalculator.getLowerMatrix(A)
        D = MatrixCalculator.getDiagonalMatrix(A)
        inversedD = MatrixCalculator.inverseDiagonal(D)

        factor1 = MatrixCalculator.multiplyMatrixByMatrix(inversedD, L+U)
        factor2 = MatrixCalculator.multiplyMatrixByVector(inversedD, b)

        r = [1]*len(A)
        res = MatrixCalculator.multiplyMatrixByVector(A, r) + -1*b
        norm = MatrixCalculator.vectorNorm(res)
        while norm > math.pow(10, -9):
            r = MatrixCalculator.multiplyMatrixByVector(factor1, r) + factor2
            res = MatrixCalculator.multiplyMatrixByVector(A, r) + -1*b
            norm = MatrixCalculator.vectorNorm(res)

            iterations+=1
        return  r, iterations



    #-----------------UTILS--------------------
    @staticmethod
    def getUpperMatrix(A):
        ret = deepcopy(A)
        for i in range(len(ret)):
            for j in range(len(ret[i])):
                if j <= i:
                    ret[i][j] = 0
        return ret

    @staticmethod
    def getLowerMatrix(A):
        ret = deepcopy(A)
        for i in range(len(ret)):
            for j in range(len(ret[i])):
                if j >= i:
                    ret[i][j] = 0
        return ret

    @staticmethod
    def getDiagonalMatrix(A):
        ret = deepcopy(A)
        for i in range(len(ret)):
            for j in range(len(ret[i])):
                if j != i:
                    ret[i][j] = 0
        return ret

    @staticmethod
    def eliminate(r1, r2, col, target=0):
        factor = (r2[col] - target) / r1[col]
        for i in range(len(r2)):
            r2[i] -= factor * r1[i]

    @staticmethod
    def gaussElimination(A):
        for i in range(len(A)):
            #If 0 is on diagonal switch rows
            if A[i][i] == 0:
                for j in range(i + 1, len(A)):
                    if A[i][j] != 0:
                        A[i], A[j] = A[j], A[i]
                        break

            for j in range(i + 1, len(A)):
                MatrixCalculator.eliminate(A[i], A[j], i)
        for i in range(len(A) - 1, -1, -1):
            for j in range(i - 1, -1, -1):
                MatrixCalculator.eliminate(A[i], A[j], i)
        for i in range(len(A)):
            MatrixCalculator.eliminate(A[i], A[i], i, 1) #Set the left ones to one
        return A

    @staticmethod
    def inverse(A):
        tmp = [[] for _ in A]
        for i, row in enumerate(A):
            assert len(row) == len(A)
            tmp[i].extend(row + [0] * i + [1] + [0] * (len(A) - i - 1))
        MatrixCalculator.gaussElimination(tmp)
        ret = []
        for i in range(len(tmp)):
            ret.append(tmp[i][len(tmp[i]) // 2:])
        return ret

    @staticmethod
    def inverseDiagonal(A):
        ret = deepcopy(A)
        for i in range(len(A)):
            ret[i][i]=1/ret[i][i]
        return ret

    @staticmethod
    def multiplyMatrixByMatrix(A,B):
        ret = []
        for i in range(len(A)):
            row=A[i][:]
            column=[i[0] for i in B]

            ret.append(MatrixCalculator.multiplyVectorByVector(row, column))
        return ret

    @staticmethod
    def multiplyMatrixByVector(A, b):
        ret = []
        for i in range(len(A)):
            current = 0
            for j in range(len(A[i])):
                current += A[i][j]*b[i]
            ret.append(current)
        return ret

    @staticmethod
    def multiplyVectorByVector(a, b):
        ret = []
        for i in range(len(a)):
            ret.append(a[i]*b[i])
        return ret

    @staticmethod
    def vectorNorm(a):
        ret = 0
        for i in range(len(a)):
            ret += math.pow(a[i], 2)
        return math.sqrt(ret)

