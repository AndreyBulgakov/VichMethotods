package sem6.task3;

import sem6.task2.Task2;

/**
 * Created by andrey on 31.03.16.
 */
public class LinearSystem {

    double[][] matrix;
    double[] vector;
    int n ;

    public LinearSystem(double[][] matrix, double[] vector) {
        this.matrix = matrix;
        this.vector = vector;
        this.n = matrix.length;
    }

    public void minus(int fromLine, int secondLine){
        for (int i = 0; i < n + 1; i++) {
            matrix[fromLine][i] = matrix[fromLine][i] - matrix[secondLine][i];
        }
    }
    public void minusWithCoefficient(int fromLine, int secondString){
        for (int i = 0; i < n + 1; i++) {
            double coefficient = matrix[fromLine][secondString]/ matrix[secondString][secondString];
            matrix[fromLine][i] = matrix[fromLine][i] - matrix[secondString][i]*coefficient;
        }
    }

    /**
     * Возвращает диагональную матрицу
     * на главной диагонали которой стоят элементы основной матрицы
     * */
    public double[][] getD(){
        double[][] D = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if(i == j)
                    D[i][j] = matrix[i][j];
            }
        }

        return D;
    }

    /**
     * Возвращает cd = D^-1 * b
     * */
    public double[] getCd(double eps) {
        double[][] reverse = Task2.reverseMatrix(getD(), eps);
        double[] cd = Task2.multiplyMatrixAndVector(reverse, vector);
        return cd;
    }
}
