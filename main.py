from DataGenerator import DataGenerator
from MatrixCalculator import MatrixCalculator
# Index number: 184296
if __name__ == '__main__':
    c = 9
    d = 6
    N = 9*c*d
    A, b = DataGenerator.generateA(N, 2, 4)
    #_, iterationsJacobi = MatrixCalculator.calculateJacobi(A, b)
    #print("Using Jacobi method it was solved in", iterationsJacobi, "iterations")
    _, iterationsGaussSeidle = MatrixCalculator.calculateGaussSeidle(A, b)
    print("Using Gauss-Seidle method it was solved in", iterationsGaussSeidle, "iterations")