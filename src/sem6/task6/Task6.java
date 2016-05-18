package sem6.task6;

import java.util.Arrays;
import java.util.function.Function;
import java.util.stream.IntStream;

import kotlin.jvm.functions.Function3;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import sem5.Task5.Task5;
import sem6.task2.Task2;
import org.apache.commons.math3.analysis.polynomials.PolynomialsUtils;
import util.Matrix;

/**
 * Метод Ритца
 * Метод Колокаций
 */
public class Task6 {
    public static void main(String[] args) {
//        int n = 3;
        int start = 3;
        int end = 15;
        double[][] res = new double[end-start+1][6];
        /////////////////////////////////////////////////////////////////////////
        for (int n = start; n <= end; n++) {

            //Формируем ОДУ с краевыми условиями
            Function<Double, Double> p = x -> 1.0 / (2.0 * x + 3.0);
            Function<Double, Double> r = x -> 2.0 - (x * x) / 2.0;
            Function<Double, Double> f = x -> 1.0 + x;
            double a1 = 1.0;
            double a2 = 0.0;
            double b1 = 0.0;
            double b2 = 1.0;
            Ly_f Ly_f = new Ly_f(p, r, f, a1, a2, b1, b2);

            // Заполняем wj
            // И их производные dwj по формуле дифференцирования
            // Полиномов якоби
            // Или через апач
            Function<Double, Double>[] w = new Function[n];
            Function<Double, Double>[] dw = new Function[n];
            Function<Double, Double>[] ddw = new Function[n];

            for (int i = 0; i < n; i++) {
                int PolyDegree = i - 2;
                if (i==0){
                    w[i] = x -> x*x-2*x-3;
                    dw[i] = x -> 2*x-2;
                    ddw[i] = x -> 2.0;

//                    w[i] = x -> 1.0+x;
//                    dw[i] = x -> 1.0;
//                    ddw[i] = x -> 0.0;
                }
                else if (i==1){

                    w[i] = x -> x*x*x-3*x-2;
                    dw[i] = x -> 3*x*x-3;
                    ddw[i] = x -> 6.0*x;

                }
                else {
//                    w[i] = x -> (1.0 - x * x) * PolynomialsUtils.createJacobiPolynomial(PolyDegree, 1, 1).value(x);
//                    dw[i] = x -> -2.0 * (PolyDegree + 1) * PolynomialsUtils.createJacobiPolynomial(PolyDegree + 1, 0, 0).value(x);
//                    ddw[i] = x -> -2.0 * (PolyDegree + 1)* (PolyDegree + 1 + 0 + 1) / 2.0 * PolynomialsUtils.createJacobiPolynomial(PolyDegree, 1, 1).value(x);
//
                    w[i] = x -> Math.pow((1.0 - x * x),2) * PolynomialsUtils.createJacobiPolynomial(PolyDegree, 2, 2)
                            .value(x);
                    dw[i] = x -> -2.0 * (PolyDegree + 1)* (1.0 - x * x) * PolynomialsUtils.createJacobiPolynomial(PolyDegree + 1, 1,1)
                            .value(x);
                    ddw[i] = x -> 4*(PolyDegree + 1) * (PolyDegree + 2 )* PolynomialsUtils.createJacobiPolynomial(PolyDegree + 2, 0, 0)
                            .value(x);
                }
            }

            double a = -0.5;
            double b = 0.5;
            double a_b = (a + b) / 2.0;
            //Получаем yn методом Ритца
            Function<Double, Double> y = ritzMethod(Ly_f, w);
            //Получаем yn методом Колакаций
            Function<Double, Double> y2 = kolockMethod(Ly_f, w, dw, ddw);

            /////////////////////////////////////////////////////////////////////////////
            //Надо для вывода
            res[n - start][0] = y.apply(a);
            res[n - start][1] = y.apply(a_b);
            res[n - start][2] = y.apply(b);
            res[n - start][3] = y2.apply(a);
            res[n - start][4] = y2.apply(a_b);
            res[n - start][5] = y2.apply(b);
        }
        System.out.println("|n|____________rn(x)__________|__________kn(x)___________|");
        for (int i = 0; i < res.length; i++) {
                System.out.printf("|%1d|%8.7f|%8.7f|%8.7f|%8.7f|%8.7f|%8.7f|\n",
                        i+start, res[i][0], res[i][1], res[i][2]
                        , res[i][3], res[i][4], res[i][5]);
        }
    }
    /**
     * Метод колокации
     * */
    public static Function<Double, Double> kolockMethod(Ly_f Ly_f,
                                                        Function<Double, Double>[] w,
                                                        Function<Double, Double>[] dw,
                                                        Function<Double, Double>[] ddw){
        int n = w.length;

        double[] t = new double[n];
        for (int i = 0; i < n; i++) {
            t[i] = Math.cos(((2 * i) * Math.PI) / (2 * n));
        }
        Arrays.sort(t);

        double[][] A = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = Ly_f.Ly(w[j], dw[j], ddw[j]).apply(t[i]);
//                A[i][j] = Ly_f.Ly(w[j]).apply(t[i]);
            }
        }

        double[] f = new double[n];
        for (int i = 0; i < n; i++) {
            f[i] = Ly_f.f.apply(t[i]);
        }
//        double[] c = Task2.gausWithChoiceRow(A, f, 1e-10);
        double[][] D = Task2.reverseMatrix(A,1e-10);

        double[] c = Matrix.multiply(D, f);

        Function<Double, Double> y = x -> IntStream
                .range(0, n-1)
                .mapToDouble(i -> c[i] * w[i].apply(x))
                .sum();

        return y;
    }


    /**
     * Метод Ритца
     * */
    public static Function<Double, Double> ritzMethod(Ly_f Ly_f, Function<Double, Double>[] w){
        int n = w.length;
        //Получаем cj
        double[][] A = new double[n][n];
        double[] f = new double[n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = y_z(Ly_f, w[j], w[i]);
            }
            f[i] = scalar(Ly_f.f, w[i]);
        }
        double[] c = Task2.gausWithChoiceRow(A, f, 1e-10);
        //Получаем yn
        Function<Double, Double> y = x -> IntStream
                .range(0, n)
                .mapToDouble(i -> c[i] * w[i].apply(x))
                .sum();
//        System.out.println("///////////////////////////////////////////////");
//        System.out.println("n = "+n);
//        System.out.println("-----Расширенная матрица системы-------");
//        Matrix.print(A, f);
//        System.out.println("-----Коэффисциенты сj-------");
//        System.out.println(Arrays.toString(c));
        return y;
    }
    /**
     * [y, z] = Int(py'z' + ryz)dx + Ql + Qr
     * */
    public static double y_z(Ly_f Ly_f, Function<Double, Double> y, Function<Double, Double> z){
        double res;
        // Параметры формулы
        // Количество узлов, промежуток интегрирования, h
        int m = 1000;
        double a = -1.0;
        double b = 1.0;
        double h = (b - a) / m;

        Function<Double, Double> dy = derivate(y, h);
        Function<Double, Double> dz = derivate(z, h);
        Function<Double, Double> p = Ly_f.p;
        Function<Double, Double> r = Ly_f.r;
        double Ql = getQl(Ly_f, y, z);
        double Qr = getQr(Ly_f, y, z);

        Function<Double, Double> integrated = x ->
                p.apply(x) * dy.apply(x) * dz.apply(x) + r.apply(x) * y.apply(x) * z.apply(x);

        res = Task5.simpson(integrated, a, b, h, m) + Ql + Qr;
        return res;
    }

    /**
     *      / 0 if I or II
     * Ql =|  a1 * p(-1)*y(-1)*z(-1) if III
     *      \ a2
     * */
    public static double getQl(Ly_f Ly_f, Function<Double, Double> y, Function<Double, Double> z){
        double a1 = Ly_f.a1;
        double a2 = Ly_f.a2;
        Function<Double, Double> p = Ly_f.p;

        if (Ly_f.taskClass == TaskClass.I || Ly_f.taskClass == TaskClass.II)
            return 0;
        else
            return (a1/a2) * p.apply(-1.0) * y.apply(-1.0) * z.apply(-1.0);

    }

    /**
     *      / 0 if I or II
     * Ql =|  b1 * p(1)*y(1)*z(1) if III
     *      \ b2
     * */
    public static double getQr(Ly_f Ly_f, Function<Double, Double> y, Function<Double, Double> z){
        double b1 = Ly_f.b1;
        double b2 = Ly_f.b2;
        Function<Double, Double> p = Ly_f.p;
        if (Ly_f.taskClass == TaskClass.I || Ly_f.taskClass == TaskClass.II)
            return 0;
        else
            return (b1/b2) * p.apply(1.0) * y.apply(1.0) * z.apply(1.0);

    }
    /**
     * Производная по формуле производных
     * */
    public static Function<Double, Double> derivate(Function<Double, Double> f, double h) {
        Function<Double,Double> res =  x -> (3 * f.apply(x) - 4 * f.apply(x - h) + f.apply(x - 2*h)) / (2*h);
        return res;
    }
    /**
     * Скалярное произведение
     * */
    public static double scalar(Function<Double, Double> f , Function<Double, Double> wi){
        int m = 1000;
        double a = -1.0;
        double b = 1.0;
        double h = (b - a) / m;

        Function<Double, Double> integrate = x -> f.apply(x) * wi.apply(x);
        return Task5.simpson(integrate, a, b, h, m);
//        return simpson2(integrate, a, b, m);

    }

    public static double simpson2(Function<Double, Double> f, double a, double b, int m){
        double h = (b - a) / m;
        double res = 0.0;
        double x1 = a;
        double x2 = x1 + h;
        for (int k = 0; k < m; k++)//k in 0.. m - 1)
        {
            res += (h / 6.0) * (f.apply(x1) + 4.0 * f.apply((x1 + x2) / 2.0) + f.apply(x2));
            x1 += h;
            x2 += h;
        }
        return res;
    }

}
