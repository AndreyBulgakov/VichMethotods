package sem6.task3;

import sem6.task2.Task2;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Решение систем иттерационными методами
 * Тут поступила инфа что априорная и аппостериорная должны на первом шаге совпадать
 * Здесь не сопадают =)
 */
public class Task3 {
    static int k_iter = 0;
    public static void main(String[] args) {
        double eps = 1e-32;
//        double eps = 1e-10;
//        double[][] A = {
//                {3.0, 2.0, -5.0},
//                {2.0, -1.0, 3.0},
//                {1.0, 2.0, -1.0}
//        };
//        double[] b = {-1.0, 13.0, 9.0};
//        double[][] A = {
//                {2.18432 , 0.37200 , 0.19955},
//                {0.37200 , 3.13980 , 0.58107},
//                {0.19955 , 0.58107 , 4.84558}
//        };
//        double[] b = {2.15308, 6.42023, 6.14104};
//

//        Димас
        double[][] A = {
                {2.19587, 0.34563, 0.17809},
                {0.34563, 3.16088, 0.55443},
                {0.17809, 0.55443, 4.89781}
        };
        double[] b = {2.16764, 6.52980, 6.29389};

        LinearSystem Axb =  new LinearSystem(A,b);

        double[] x0 = new double[A.length];

        //Точное решение методом гауса
        double[] x = Task2.gausClassic(A,b,eps);
        System.out.println("Решение методом гауса : " + Arrays.toString(x));

        //Круги Гершгорина
        ArrayList<Range> circles = gershgorinCircle(A);
        System.out.println("Круги гершгорина: " + circles.toString());

        LinearSystem Hxg =  transformToHgMatrix(A,b);
        double[][] H = Hxg.matrix;
        double[] g = Hxg.vector;
        double normH = matrixNorm(H);

        //Априорная оценка k
        int k = apriorK(Axb,x0,eps);

        //Итерационный и метод Зейделя
        System.out.println("k(eps) = "+k);
        double[] xiter      = iterationMethod(Axb,x0, x, eps);
        double[] xzeidel    = zeidelMethod(Axb, x0, x, k_iter);

        System.out.println("Фактическая погрешность = " + vectorNorm(minus(x, xzeidel)));
        System.out.println("Отклонение от итерационного метода = " + vectorNorm(minus(xzeidel, xiter)));
    }

    /**
     * Круги Гершгорина
     * Ищет промежутки в которых лажат собственные числа
     * см. http://sa.technolog.edu.ru/files%5Cchumakov%5CSobstvennye%20znachenija.pdf
     * */
    public static ArrayList<Range> gershgorinCircle(double[][] matrix){
        double[][] A = Arrays.stream(matrix).map(double[]::clone).toArray(double[][]::new);
        int n = A.length;
        ArrayList<Range> ranges = new ArrayList<>(n);

        //Ищем ri = sum(|aij|) i != j j=0..n-1
        double[] r = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j)
                    r[i] += Math.abs(A[i][j]);
            }
        }

        //Собственное значение лежит в промежутке [aii-ri; aii+ri]
        //Заполняем промежутки
        for (int i = 0; i < n; i++) {
            double begin = A[i][i] - r[i];
            double end = A[i][i] + r[i];
            ranges.add(new Range(begin,end));
        }

        return ranges;
    }

    /**
     * Преобразует матрицу вида Ax=b
     * К виду x = Bx + c
     * Где D(в индексах B и с) диагональная матрица у которой
     * на главной диагонали находятся диагональные элементы матрицы A
     * */
    public static LinearSystem transformToHgMatrix(double[][] matrixA, double[] resultb){
        //Копируем матрицу
        double[][] A = Arrays.stream(matrixA).map(double[]::clone).toArray(double[][]::new);
        double[] b = resultb.clone();
        int n = A.length;
        double[][] H = new double[n][n];
        double[] g = new double[n];

        //Считаем H
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j)
                    H[i][j] = 0;
                else
                    H[i][j] = -(A[i][j] / A[i][i]);
            }
        }
        //Считаем g
        for (int i = 0; i < n; i++) {
            g[i] = b[i]/A[i][i];
        }
        return new LinearSystem(H,g);
    }

    /**
     * Норма ||A|| infinity
     * */
    public static double matrixNorm(double[][] matrix){
        return Arrays.stream(matrix.clone()).map(arr -> Arrays.stream(arr.clone()).map(Math::abs).sum()).max(Double::compare).get();
    }

    /**
     * Норма вектора
     * max|xi|
     * */
    public static double vectorNorm(double[] vector){
        return Arrays.stream(vector.clone()).map(Math::abs).max().getAsDouble();
    }

    /**
     * Ищет априорную оценку того k, при котором ||x(*)-xK||< e
     * */
    public static int apriorK(LinearSystem Axb, double[] x0, double eps){
        LinearSystem Hxg = transformToHgMatrix(Axb.matrix, Axb.vector);
        //Копируем матрицу
        double[][] H = Arrays.stream(Hxg.matrix).map(double[]::clone).toArray(double[][]::new);
        double[] g = Hxg.vector.clone();
        int n = H.length;
        double normH = matrixNorm(H);

//        //Формулы H*b+c;
        double[] cd = Axb.getCd(eps);

        //H*xk
        double[] xk = Task2.multiplyMatrixAndVector(H,x0);
        //Прибавляем с
        xk = plus(xk,cd);

        //Ищем xk-x0
        double[] xk_x0 = minus(xk,x0);
        //Считаем норму вектора xk-x0

        //По формулам
        int k = 0;
        double normx0 = vectorNorm(x0);
        double normg = vectorNorm(g);
        //x^(k-1)
        double[] xk_1 = x0.clone();
        double aprior = Math.pow(normH,k)*normx0+(Math.pow(normH,k)/(1-normH))*normg;
        double apost = (normH/(1-normH))*vectorNorm(minus(xk, xk_1));

        double t = (eps * (1 - normH)) / aprior;

        // TODO: 04.04.16 long if double
        k = (int) Math.round(Math.log(t)/Math.log(normH));
        return k;
    }

    /**
     * Метод итераций
     * */
    public static double[] iterationMethod(LinearSystem Axb, double[] x0, double[] solutionX, double eps){
        LinearSystem Hxg = transformToHgMatrix(Axb.matrix, Axb.vector);
        //Копируем матрицу
        double[][] H = Arrays.stream(Hxg.matrix).map(double[]::clone).toArray(double[][]::new);
        double[] g = Hxg.vector.clone();
        int n = H.length;
        double normH = matrixNorm(H);

//        //Формулы H*b+c;
        double[] cd = Axb.getCd(eps);

        //H*xk
        double[] xk = Task2.multiplyMatrixAndVector(H,x0);
        //Прибавляем с
        xk = plus(xk,cd);

        //Ищем xk-x0
        double[] xk_x0 = minus(xk,x0);

        //Считаем норму вектора xk-x0
        //По формулам
        int k = 0;
        double normx0 = vectorNorm(x0);
        double normg = vectorNorm(g);
        double[] xk_1 = x0.clone();

        //Апрорная и апостериорная оценка по формулам
        double aprior = Math.pow(normH,k)*normx0+(Math.pow(normH,k)/(1-normH))*normg;
        double apost = (normH/(1-normH))*vectorNorm(minus(xk, xk_1));


        String formatter = "Count = %d, \n x_k = %s, \n Error = %9.8e, \n Apriority = %9.8e, \n Aposterior = %9.8e \n\n";
        System.out.printf(formatter, k, Arrays.toString(xk), vectorNorm(minus(xk,solutionX)), aprior, apost);

        while(vectorNorm(minus(xk,xk_1)) > eps){
            //Пересчитываем x^(k+1) = Hx^(k)+c
            xk_1 = xk.clone();
            xk = Task2.multiplyMatrixAndVector(H, xk);
            xk = plus(xk, cd);
            k++;

            //Пересчитыаем априорную и апостериорную формулы те же
            aprior = Math.pow(normH,k)*normx0+(Math.pow(normH,k)/(1-normH))*normg;
            apost = (normH/(1-normH))*vectorNorm(minus(xk, xk_1));

            //Вывод
            System.out.printf(formatter, k, Arrays.toString(xk), vectorNorm(minus(xk,solutionX)), aprior, apost);
        }

        k_iter = k;
        return xk;
    }

    /**
     * Метод Зейделя
     * */
    public static double[] zeidelMethod(LinearSystem Axb, double[] x0, double[] solutionX, int k_iter){
        LinearSystem Hxg = transformToHgMatrix(Axb.matrix, Axb.vector);
        //Копируем матрицу
        double[][] H = Arrays.stream(Hxg.matrix).map(double[]::clone).toArray(double[][]::new);
        double[] g = Hxg.vector.clone();
        int n = H.length;

        double[] xk = x0.clone();

        int k = 0;
        //По формулам
        while (true){
            k++;
            //x^(k-1)
            double[] xk_1 = xk.clone();

            //Расчетная формула
            for (int i = 0; i < n; i++) {
                double tmp = 0;

                for (int j = 0; j < i; j++)
                    tmp += H[i][j] * xk[j];

                for (int j = i; j < n; j++)
                    tmp += H[i][j] * xk_1[j];
//                xk[i] = (b[i] - tmp) / A[i][i];
                xk[i] = tmp + g[i];
            }

            if (k==k_iter) break;
        }
        System.out.println(k);
        return xk;
    }
    /**
     * Склаыдвает два вектора поэлементно
     * */
    public static double[] plus(double[] a, double[]b){
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++)
            result[i] = a[i]+b[i];
        return result;
    }

    /**
     * Вычитает два вектора поэлементно
     * */
    public static double[] minus(double[] a, double[]b){
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++)
            result[i] = a[i]-b[i];
        return result;
    }

}
