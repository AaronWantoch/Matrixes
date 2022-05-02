import math

class DataGenerator:
    @staticmethod
    def generateA(N, e, f):
        a1 = 5 + e
        a2 = a3 = -1
        A = [[0 for x in range(N)] for y in range(N)]
        for i in range(N):
            A[i][i] = a1
            if i+1<N:
                A[i + 1][i] = a2
                A[i][i + 1] = a2
            if i+2<N:
                A[i + 2][i] = a3
                A[i][i + 2] = a3

        b = []
        for i in range(N):
            b.append(math.sin(i*(f+1)))
        return A,b
