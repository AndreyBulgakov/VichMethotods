package sem6.task2;

import sem6.task1.Task1;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Точные методы решения систем линейных уравнений
 */
public class Task2 {
    public static void main(String[] args) {
        double[] x;
        double eps = 0.0001;
//        double[][] A = {
//                {3.0, 2.0, -5.0},
//                {2.0, -1.0, 3.0},
//                {1.0, 2.0, -1.0}
//        };
//        double[] b = {-1.0, 13.0, 9.0};
        double[][] A = {
                {6.6352e-06 , -7.8296e-03, 4.27056},
                {6.0352e-03 , -0.78296  , 1.50056},
                {0.88352    , 0.81704   , 0.95056}
        };
        double[] b = {3.70946, 1.58087, 1.88704};
        System.out.println("Матрица A|b");
        printAb(A,b);

        //Обратнаяматрица
        System.out.println("Обратная матрица");
        double[][] reverse = reverseMatrix(A,eps);
        Task1.printMatrix(reverse);

        //Классический метод
        System.out.println("\nГаус классический");
        x = gausClassic(A, b,eps);
        System.out.println(Arrays.toString(x));
        System.out.println("Невязка = "+neviazka(A,x,b));

        //Гаус с выбором
        System.out.println("\nГаус с выбором");
        x = gausWithChoiceRow(A, b, eps);
        System.out.println(Arrays.toString(x));
        System.out.println("Невязка = "+neviazka(A,x,b));

        //Обратная
        System.out.println("\nРешение с использованием обратной");
        x = multiplyMatrixAndVector(reverse, b);
        System.out.println(Arrays.toString(x));
        System.out.println("Невязка = "+neviazka(A,x,b));

        //Определитель
        double determinant = determinant(A);
        System.out.println("\nDeterminant = " + determinant);


    }

    /**
     * Обычный метод Гауса
     * Схема единственного деления
     * */
    public static double[] gausClassic(double[][] matrixA, double[] resultB, double eps){
        //Копируем матрицу
        double[][] A = Arrays.stream(matrixA).map(double[]::clone).toArray(double[][]::new);
        double[] b = resultB.clone();
        int n = A.length;
        double[] x = new double[n];

        //Прямой ход
        for (int k = 0; k < n ; k++) {
            //Вывод сообщения в случае погрешности
            if (A[k][k] < eps)
                System.out.printf("A[%d][%d] < eps\n",k,k);

            double Akktmp = A[k][k];

            for (int j = k; j < n ; j++) {
                A[k][j] = A[k][j] / Akktmp;
            }
            b[k] = b[k] / Akktmp;

            for (int i = k + 1; i < n ; i++) {
                double Aijtmp = A[i][k];
                for (int j = k; j < n ; j++) {
                    A[i][j] = A[i][j] - A[k][j]*Aijtmp;
                }
                b[i] = b[i]-b[k]*Aijtmp;
            }
        }

        //Обратный ход
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j]*x[j];
            }
            x[i] = b[i]-sum;
        }
        return x;
    }

    /**
     * Метод Гауса с выбором главного элемента
     * Выбор главного элемента по столбцу
     * */
    public static double[] gausWithChoiceRow(double[][] matrixA, double[] resultB, double eps){
        //Копируем матрицу
        double[][] A = Arrays.stream(matrixA).map(double[]::clone).toArray(double[][]::new);
        double[] b = resultB.clone();

        int n = A.length;
        double[] x = new double[n];

        for (int k = 0; k < n ; k++) {
            //Вывод сообщения в случае погрешности
            if (Math.abs(A[k][k]) < eps) {
//                System.out.printf("A[%d][%d] < eps\n", k, k);
                int p = findMaxAbs(A, k, n);
                swapRows(A,p,k);
                swapRows(b,p,k);
            }

            double Akktmp = A[k][k];

            for (int j = k; j < n ; j++) {
                A[k][j] = A[k][j] / Akktmp;
            }
            b[k] = b[k] / Akktmp;

            for (int i = k + 1; i < n ; i++) {
                double Aijtmp = A[i][k];
                for (int j = k; j < n ; j++) {
                    A[i][j] = A[i][j] - A[k][j]*Aijtmp;
                }
                b[i] = b[i]-b[k]*Aijtmp;
            }
        }

        //Обратный ход
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j]*x[j];
            }
            x[i] = b[i]-sum;
        }

        return x;
    }

    /**
     * Ищет максимальный по модулю элемент в столбце расположенный ниже k
     * */
    private static int findMaxAbs(double[][] A, int k, int n){
        int p = k;
        double max = Math.abs(A[k][k]);
        for (int i = k; i < n; i++) {
            if (Math.abs(A[i][k]) > max){
                max = Math.abs(A[i][k]);
                p = i;
            }
        }
        return p;
    }

    /**
     * Меняет местами строки
     * */
    private static void swapRows(double[][] A, int p, int k){
        double[] tmp = A[p];
        A[p] = A[k];
        A[k] = tmp;
    }
    private static void swapRows(double[] b, int p, int k){
        double tmp = b[p];
        b[p] = b[k];
        b[k] = tmp;
    }

    /**
     * Считает определитель матрицы
     * */
    public static double determinant(double[][] matrix){
        double determinant = 1.0;

        //Копируем матрицу
        double[][] A = Arrays.stream(matrix).map(double[]::clone).toArray(double[][]::new);
        int n = A.length;

        //Прямой ход
        for (int k = 0; k < n ; k++) {
            double Akktmp = A[k][k];
            determinant *= Akktmp;
            for (int j = k; j < n ; j++) {
                A[k][j] = A[k][j] / Akktmp;
            }

            for (int i = k + 1; i < n ; i++) {
                double Aijtmp = A[i][k];
                for (int j = k; j < n ; j++) {
                    A[i][j] = A[i][j] - A[k][j]*Aijtmp;
                }
            }
        }
        return determinant;
    }

    /**
     * Обратная матрица
     * */
    public static double[][] reverseMatrix(double[][] matrix, double eps){
        int n = matrix.length;
        double[][] reverse = new double[n][n];
//        double[][] res1 = IntStream.range(0, 3)
//                .mapToObj(i -> IntStream.range(0, 3)
//                        .mapToDouble(j -> i==j ? 1.0 : 0.0)
//                        .toArray()).toArray(double[][]::new);
        //Задаем единичную матрицу
        double[][] E = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                E[i][j] = i == j ? 1.0 : 0.0;
        //Считаем решение матрицы метоодом Гауса
        // используя в качестве стобца b строку единичной матрицы
        //Получаем транспонированную матрицу
        double[][] matrixT = new double[n][n];
        for (int i = 0; i < n; i++)
            matrixT[i] = gausWithChoiceRow(matrix, E[i], eps);

        //Переворачиваем транспонированную и получаем A^-1
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                reverse[i][j] = matrixT[j][i];

        return reverse;
    }

    /**
     * Умножает матрицу на вектор
     * Используется для решения системы используя обратную матрицу
     * */
    public static double[] multiplyMatrixAndVector(double[][] matrixA, double[] vector){
        //Копируем матрицу
        double[][] A = Arrays.stream(matrixA).map(double[]::clone).toArray(double[][]::new);
        double[] b = vector.clone();

        int n = A.length;
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            double xi = 0.0;
            for (int j = 0; j < n; j++) {
                xi += A[i][j]*b[j];
            }
            x[i] = xi;
        }
        return x;
    }

    /**
     * Невязка
     * */
    public static double neviazka(double[][] matrixA, double[] solutionX, double[] resultB) {
        double neviazka = 0.0;
        //Копируем матрицу
        double[][] A = Arrays.stream(matrixA).map(double[]::clone).toArray(double[][]::new);
        double[] b = resultB.clone();
        double[] x = solutionX.clone();
        int n = A.length;

        //Считаем новые b
        double[] solutionB = multiplyMatrixAndVector(A, x);

        //Считаем вектор невязки вычитая из имеющегося b полученный b1
        //После перепножения A на x
        double[] vectorNeviazkiAbs = IntStream.range(0, n).mapToDouble(i -> Math.abs(b[i]-solutionB[i])).toArray();

        //Считаем норму(значения внутри уже по модулю лежат
        neviazka = Arrays.stream(vectorNeviazkiAbs).sum();

        return neviazka;
    }

    /**
     * Печатае матрицу A|b принимая на вход A и b
     * */
    public static void printAb(double[][] A, double[] b){
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A.length; j++) {
                if (j != A.length - 1)
                    System.out.printf("|%12.9f", A[i][j]);
                else
                    System.out.printf("|%12.9f|%12.9f|\n", A[i][j], b[i]);
            }
        }
    }
}

