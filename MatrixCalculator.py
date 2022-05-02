import math
import numpy as np
from copy import deepcopy

class MatrixCalculator:
    @staticmethod
    def calculateJacobi(A, b):
        iterations=0
        U = MatrixCalculator.getUpperMatrix(A)
        L = MatrixCalculator.getLowerMatrix(A)
        D = MatrixCalculator.getDiagonalMatrix(A)
        inversedD = MatrixCalculator.inverseDiagonal(D)
        minusInversedD = MatrixCalculator.multiplyMatrixByNumber(inversedD, -1)
        LUsum = MatrixCalculator.sumMatrixes(L, U)

        factor1 = MatrixCalculator.multiplyMatrixByMatrix(
            minusInversedD,
            LUsum)
        factor2 = MatrixCalculator.multiplyMatrixByVector(inversedD, b)

        minusb = MatrixCalculator.multiplyVectorByNumber(b, -1)
        r = [1]*len(A)

        res = MatrixCalculator.sumVectors(
            MatrixCalculator.multiplyMatrixByVector(A, r),
            minusb) #Ar-b
        norm = MatrixCalculator.vectorNorm(res)

        while norm > math.pow(10, -9):
            r = MatrixCalculator.sumVectors(
                MatrixCalculator.multiplyMatrixByVector(factor1, r),
                factor2)
            res = MatrixCalculator.sumVectors(
                  MatrixCalculator.multiplyMatrixByVector(A, r),
                  minusb)
            norm = MatrixCalculator.vectorNorm(res)
            iterations += 1
        return r, iterations

    @staticmethod
    def calculateGaussSeidle(A, b):
        iterations = 0
        U = MatrixCalculator.getUpperMatrix(A)
        L = MatrixCalculator.getLowerMatrix(A)
        D = MatrixCalculator.getDiagonalMatrix(A)
        inversedDLSum = MatrixCalculator.inverse(MatrixCalculator.sumMatrixes(D, L))

        factor1 = MatrixCalculator.multiplyMatrixByNumber(inversedDLSum, -1)
        factor2 = MatrixCalculator.multiplyMatrixByVector(inversedDLSum, b)

        minusb = MatrixCalculator.multiplyVectorByNumber(b, -1)
        r = [1] * len(A)

        res = MatrixCalculator.sumVectors(
            MatrixCalculator.multiplyMatrixByVector(A, r),
            minusb)  # Ar-b
        norm = MatrixCalculator.vectorNorm(res)

        while norm > math.pow(10, -9):
            r = MatrixCalculator.sumVectors(
                MatrixCalculator.multiplyMatrixByVector(factor1, MatrixCalculator.multiplyMatrixByVector(U, r)),
                factor2)
            res = MatrixCalculator.sumVectors(
                MatrixCalculator.multiplyMatrixByVector(A, r),
                minusb)
            norm = MatrixCalculator.vectorNorm(res)
            iterations += 1
        return r, iterations



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
        if factor == 0.0:
            return
        for i in range(len(r2)):
            r2[i] -= factor * r1[i]

    @staticmethod
    def gaussElimination(A):
        for i in range(len(A)):
            # If 0 is on diagonal switch rows
            if A[i][i] == 0:
                for j in range(i + 1, min(i + 3, len(A))): #only check 2 next cause rest are 0
                    if A[j][i] != 0:
                        A[i], A[j] = A[j], A[i]
                        break

            for j in range(i + 1, min(i + 3, len(A))):#only eliminate 2 next cause rest are 0
                MatrixCalculator.eliminate(A[i], A[j], i)
        for i in range(len(A) - 1, -1, -1):
            for j in range(i - 1, max(i - 3, -1), -1):
                MatrixCalculator.eliminate(A[i], A[j], i)
        for i in range(len(A)):
            MatrixCalculator.eliminate(A[i], A[i], i, 1)  # Set the diagonal ones to one
        return A

    @staticmethod
    def inverse(A):
        tmp = [[] for _ in A]
        for i, row in enumerate(A):
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
    def inverseSymethric(A):
        ret = deepcopy(A)
        for i in range(len(A)):
            ret[i][i] = 1 / ret[i][i]
        return ret

    @staticmethod
    def multiplyMatrixByMatrix(A,B):
        ret = []
        for i in range(len(A)):
            current=[]
            row = A[i]
            for j in range(len(B)):
                column = [k[j] for k in B]
                current.append(MatrixCalculator.multiplyVectorByVector(row, column))
            ret.append(current)
        return ret

    @staticmethod
    def multiplyVectorByVector(a, b):
        ret = 0
        for i in range(len(a)):
            ret += a[i]*b[i]
        return ret

    @staticmethod
    def multiplyMatrixByVector(A, b):
        ret = []
        for i in range(len(A)):
            current = 0
            for j in range(len(A[i])):
                current += A[j][i]*b[i]
            ret.append(current)
        return ret

    @staticmethod
    def vectorNorm(a):
        ret = 0
        for i in range(len(a)):
            ret += math.pow(a[i], 2)
        return math.sqrt(ret)

    @staticmethod
    def sumMatrixes(A, B):
        ret = []
        for i in range(len(A)):
            row = []
            for j in range(len(B)):
                row.append(A[i][j]+B[i][j])
            ret.append(row)
        return ret

    @staticmethod
    def sumVectors(A, B):
        ret = []
        for i in range(len(A)):
            ret.append(A[i] + B[i])
        return ret

    @staticmethod
    def multiplyVectorByNumber(a, x):
        return [element * x for element in a]

    @staticmethod
    def multiplyMatrixByNumber(A, x):
        ret = []
        for i in range(len(A)):
            ret.append([element * x for element in A[i]])
        return ret

    @staticmethod
    def forwardSubstitiution(A, b):
        ret = []
        for i in range(len(A)):
            sum=0
            for j in range(i):
                sum+=A[i][j] * ret[j]
            ret.append((b[i]-sum)/A[i][i])
        return ret