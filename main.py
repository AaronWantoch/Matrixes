import time

from DataGenerator import DataGenerator
from MatrixCalculator import MatrixCalculator
# Index number: 184296
if __name__ == '__main__':
    c = 9
    d = 6
    N = 9*c*d
    A, b = DataGenerator.generateA(N, 5+2, -1, -1, 4)

    startJacobi = time.time()
    _, iterationsJacobi = MatrixCalculator.calculateJacobi(A, b)
    endJacobi = time.time()
    print("Using Jacobi method it was solved in", iterationsJacobi, "iterations and lasted", endJacobi-startJacobi,
          "seconds")

    startGauss = time.time()
    _, iterationsGaussSeidle = MatrixCalculator.calculateGaussSeidle(A, b)
    endGauss = time.time()
    print("Using Gauss-Seidle method it was solved in", iterationsGaussSeidle, "iterations and lasted",
          endGauss-startGauss, "seconds")

    
