from DataGenerator import DataGenerator
from MatrixCalculator import MatrixCalculator
# Index number: 184296
if __name__ == '__main__':
    c = 9
    d = 6
    N=9*c*d
    A, b = DataGenerator.generateA(N,2, 4)
    _, iterations = MatrixCalculator.calculateJacobi(A, b)
    print("Using Jacob method it was solved in", iterations, "iterations")