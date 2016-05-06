package sem6.task4;

import sem6.task2.Task2;
import sem6.task3.Task3;
import util.Matrix;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.Arrays;

public class Task4 {

    public static void main(String[] args) {
        double[][] A = {
                {-0.881923, -0.046444, 0.334218},
                {-0.046444, 0.560226, 0.010752},
                {0.334218, 0.010752, -0.883417}
        };
        System.out.println("Матрица А");
        printA(A);
        double eps = 10e-7;
//        L1 = Arrays.stream(stepennoyMethod(A)).reduce((i,j) -> Double.max(Math.abs(i), Math.abs(j))).orElse(0);
        System.out.println("----------------Степенной метод-------------");
        double L1 = stepennoyMethod(A, eps);
        System.out.println(L1);
        System.out.println("----------------Скалярный метод-------------");
        double L1scalar = scalarMethod(A, eps);

        //Обратная итерация
        System.out.println("----------------Метод обратной итерации-------------");
        PrintStream out = System.out;
        System.setOut(new PrintStream(new ByteArrayOutputStream()));

        double[][] A_1 = Task2.reverseMatrix(A, eps);
        double Ln = stepennoyMethod(A_1, eps);
        System.setOut(out);
        System.out.println("Ln = "+1.0/Ln);

        //Противоположенная граница
        System.out.println("----------------Противоположенная граница спектра-------------");
        System.setOut(new PrintStream(new ByteArrayOutputStream()));

        double opositeSpectreL = getOpositeSpectreL(A, eps);
        System.setOut(out);
        System.out.println("L_ = "+opositeSpectreL);

        System.out.println("Потерянное значение: " + missingLambdaInA3(A,L1,opositeSpectreL));

    }

    /**
     * Степенной метод
     * */
    public static double stepennoyMethod(double[][] A, double eps){
        int n = A.length;
        double[] Y0 = new double[n];

        //Set all values to 1
        Arrays.fill(Y0, 1.0);

        boolean isAnyLHigherThenEps = true;
        int k = 0;

        // Инициализация
        double[] Yk_1 = new double[0];
        double[] Yk = Y0;
        double[] Lk_1 = new double[n];
        double[] Lk = new double[n];

        while(isAnyLHigherThenEps){
            k++;
            // По формулам
            Yk_1 = Yk;
            Lk_1 = Lk;

            // Теперь заолним новыми значениями
            Yk = Task2.multiplyMatrixAndVector(A, Yk_1);
            Lk = new double[n];

            // Создаем массив собственных чисел на kой иттерации.
            // (не использую лямбды т.к. переменные не final)
            isAnyLHigherThenEps = false;
            Lk[0] = Yk[0]/Yk_1[0];

            //Проверка точности
            if (Math.abs(Lk[0]-Lk_1[0]) > eps)
                isAnyLHigherThenEps = true;
        }
        //Нормируем собсственный вектор
        double normYk = getMaxAbs(Yk);
        for (int i = 0; i < Yk.length; i++) {
            Yk[i] = Yk[i] / normYk;
        }

        // Вывод всего
        // Вывод собственного вектора
        // Можно нормаировать
        // В данном случае если в качетве Y0 взять вектор 0.1
        System.out.println("Eps = "+ eps);
        System.out.println("k = "+ k);
        System.out.println("L1 = " + Lk[0]);
        System.out.print("Собственный вектор (степенной):");
        System.out.println(Arrays.toString(Yk));
        System.out.println("Невязка = " + nevayzka(A, Yk, Lk[0]));

        return Lk[0];
    }

    /**
     * Скалярный метод
     * все кроме основной формулы аналогично предыдущему
     * Работает только для симметричной матрицы
     * если нужно больше возьмите
     * Xk+1 = A*Xk
     * Yk+1 = AT*Yk
     * Lk = (Xk+1, Yk+1)
     *      (Xk, Yk+1)
     * http://www.life-prog.ru/1_17667_algoritm-vichisleniya-ocherednogo-t---go-sobstvennogo-znacheniya-i.html
     * http://www.slideshare.net/tm_ssau/ss-12537172
     * */
    public static double scalarMethod(double[][] A ,double eps){

        int n = A.length;
        double[] Y0 = new double[n];

        //Set all values to 1
        Arrays.fill(Y0, 1.0);

        boolean isAnyLHigherThenEps = true;
        int k = 0;

        // Инициализация

        //Xk будут перемножаться с обычной матрицей
        double[] Xk_1 = new double[A.length];
        double[] Xk = Y0;

        //Масив Собственных значений
        double[] Lk_1 = new double[n];
        double[] Lk = new double[n];

        while(isAnyLHigherThenEps){
            k++;
            // По формулам
            Xk_1 = Xk;

            Lk_1 = Lk;

            // Теперь заолним новыми значениями
            Xk = Task2.multiplyMatrixAndVector(A, Xk_1);

            Lk = new double[n];

            // Создаем массив собственных чисел на kой иттерации.
            // (не использую лямбды т.к. переменные не final)
            isAnyLHigherThenEps = false;

            //Основная формула
            Lk[0] = scalarMul(Xk, Xk_1) / scalarMul(Xk_1, Xk_1);

            if (Math.abs(Lk[0]-Lk_1[0]) > eps)
                isAnyLHigherThenEps = true;

        }
        //Нормируем собсственный вектор
        double normYk = getMaxAbs(Xk);
        for (int i = 0; i < Xk.length; i++) {
            Xk[i] = Xk[i] / normYk;
        }
        // Вывод всего
        // Вывод собственного вектора
        // Можно нормаировать
        // В данном случае если в качетве Y0 взять вектор 0.1
        System.out.println("Eps = "+ eps);
        System.out.println("k = "+ k);
        System.out.println("L1 = " + Lk[0]);
        System.out.print("Собственный вектор (скалярный):");
        System.out.println(Arrays.toString(Xk));
        System.out.println("Невязка = " + nevayzka(A, Xk, Lk[0]));
        return Lk[0];
    }

    /**
     * Противоположенная граница спектра
     * По формулам
     * */
    public static double getOpositeSpectreL(double[][] matrix, double eps){
        double lA = stepennoyMethod(matrix, eps);
        double[][] B = Arrays.stream(matrix).map(double[]::clone).toArray(double[][]::new);
        int n = B.length;
        for (int i = 0; i < n; i++) {
            B[i][i] -= lA;
        }
        double lB = stepennoyMethod(B, eps);
        return  lA + lB;
    }

    /**
     * Недостаюзее собственное значение A
     * Матрица должна быть невырожденной размерности 3 и Ln != L*
     * */
    public static double missingLambdaInA3(double[][] A, double l1, double opositL){
        if (A.length == 3 && Task2.determinant(A) != 0){
            return (A[0][0]+A[1][1]+A[2][2]) - (l1 + opositL);
        }
        throw  new IllegalArgumentException("Матрица вырожденная или n != 3");
    }

    /**
     * Максимальный по модулю L
     * */
    public static double getMaxAbs(double[] x){
        // Берем максимальный по модулю.
        double L1 = 0;
        for (double arr : x) {
            if (Math.abs(arr) > Math.abs(L1))
                L1 = arr;
        }
        return L1;
    }

    /**
     * Печатае матрицу A|b принимая на вход A и b
     * */
    public static void printA(double[][] A){
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A.length; j++) {
                if (j != A.length - 1)
                    System.out.printf("|%12.9f", A[i][j]);
                else
                    System.out.printf("|%12.9f|\n", A[i][j]);
            }
        }
    }

    /**
     * Невязка
     * */
    public static double nevayzka(double[][] A, double[] X, double L){
        int n = A.length;
        //Задаем единичную матрицу
        double[][] E = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                E[i][j] = i == j ? L : 0.0;
        double[][] A_lE = Matrix.subtract(A, E);
        double[] res = Matrix.multiply(X, A_lE);
        double nevyazka = Task3.vectorNorm(res);
        return nevyazka;
    }
    /**
     * Перемножение векторов
     * */
    public static double scalarMul(double[] v1, double[] v2){
        double res = 0;
        for (int i = 0; i < v1.length; i++) {
            res += v1[i]*v2[i];
        }
        return res;
    }
}
