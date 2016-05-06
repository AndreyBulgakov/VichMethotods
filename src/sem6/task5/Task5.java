package sem6.task5;

import sem6.task2.Task2;

import java.util.Arrays;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.IntStream;

public class Task5 {

    public static void main(String[] args) {
        Function<Double, Double> f = x -> 2 - x + x*x;
        BiFunction<Double, Double, Double> K = (x,y) -> Math.sinh(x*y);
        int n = 3;
        int N = 10;
        double a = 0;
        double b = 1;
        double[] x = {-0.7745966692, 0, 0.7745966692};
        double[] A = {0.5555555556, 0.8888888889, 0.5555555556};

        for (int i = 0; i < x.length; i++) {
            x[i] = ((b-a)/2)*x[i] + (b+a)/2;
            A[i] = ((b-a)/2)*A[i];
        }

        Function<Double, Double> un = MMK(K, f, n, x, A);
        double[] res = solveGrid(un, a, b, N);

        System.out.println("n = "+ n);
        System.out.println("N = "+ N);
        System.out.println("Узлы");
        System.out.println(Arrays.toString(x));
        System.out.println("Коэффициенты");
        System.out.println(Arrays.toString(A));
        System.out.println(Arrays.toString(res));

        n = 6;
        N = 10;
        a = 0;
        b = 1;
        x = new double[]{-0.9324695142, -0.6612093865, -0.2386191861, 0.2386191861, 0.6612093865,  0.9324695142, };
        A = new double[]{0.1713244924, 0.3607615730, 0.4679139346, 0.4679139346, 0.3607615730, 0.1713244924};

        for (int i = 0; i < x.length; i++) {
            x[i] = ((b-a)/2)*x[i] + (b+a)/2;
            A[i] = ((b-a)/2)*A[i];
        }

        un = MMK(K, f, n, x, A);
        double[] res2 = solveGrid(un, a, b, N);
        System.out.println("--------------------------------");
        System.out.println("n = "+ n);
        System.out.println("N = "+ N);
        System.out.println("Узлы");
        System.out.println(Arrays.toString(x));
        System.out.println("Коэффициенты");
        System.out.println(Arrays.toString(A));
        System.out.println(Arrays.toString(res2));

        System.out.println("Delta = " + delta(res, res2));
    }
    /**
     * Считает значения приближенного решения в точках равномерной на [a;b] сетки
     * */
    public static double[] solveGrid(Function<Double, Double> un, double a, double b, int n){
        double h = (b - a)/n;
        double[] x = IntStream
                .range(0,n)
                .mapToDouble(k -> un.apply(a + k * h))
                .toArray();
        return x;
    }

    /**
     * Метод мехханических квадратур
     * Возвращает un из формул
     * */
    public static Function<Double, Double> MMK(BiFunction<Double, Double, Double> K, Function<Double, Double> f,
                           int n, double[] x, double[] A){
        //Заполняем G
        double[][] G = new double[n][n];
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                if (k==j)
                    G[k][k] = 1 + A[k] * K.apply(x[k],x[k]);
                else
                    G[k][j] = A[j] * K.apply(x[k],x[k]);
            }
        }

        //Заполняем значения в узлах
        double[] F = new double[n];
        for (int i = 0; i < n; i++) {
            F[i] = f.apply(x[i]);
        }

        //Решаем методом Гауса
        double[] u = Task2.gausWithChoiceRow(G,F,0.00001);

//        System.out.println(Arrays.toString(u));

        Function<Double, Double> sum = y -> IntStream
                .range(0, n)
                .mapToDouble(j -> A[j]*K.apply(y, x[j])*u[j])
                .sum();
        Function<Double, Double> un = y -> f.apply(y) - sum.apply(y);
        return un;
    }
    /**
     * Считает
     * delta = max | v1[i]-v2[i] | i=0..n
     * */
    public static double delta(double[] v1, double[]v2){
        int n = v1.length;
        double[] res = new double[n];
        for (int i = 0; i < n; i++) {
            res[i] = Math.abs(v1[i]-v2[i]);
        }
        double delta = Arrays.stream(res).max().getAsDouble();
        return delta;
    }
}
