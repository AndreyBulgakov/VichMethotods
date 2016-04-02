package sem5.Task7;

import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * Created by andrey on 08.12.15.
 */
public class Task7 {

    public static void main(String[] args) throws IOException {
        //Решение
        Function<Double, Double> f = x -> 3/(2*Math.exp(3*x)+1);
        //N ая производная
        BiFunction<Double, Double, Double> d1y = (x, y) -> -3*y+y*y;
        BiFunction<Double, Double, Double> d2y = (x, y) -> (2*y-3)*d1y.apply(x,y);
        BiFunction<Double, Double, Double> d3y = (x, y) -> (2*y-3)*d2y.apply(x,y)+2*d1y.apply(x,y)*d1y.apply(x,y);
        BiFunction<Double, Double, Double> d4y = (x, y) -> (2*y-3)*d3y.apply(x,y)+6*d1y.apply(x,y)*d2y.apply(x,y);

//        Function<Double, Double> f = x -> Math.exp(Math.sin(x) - 2.0*x);
//        BiFunction<Double, Double, Double> d1y = (x, y) -> -1.0 * y * (2.0 - Math.cos(x));
//        BiFunction<Double, Double, Double> d2y = (x, y) -> d1y.apply(x, y)*(Math.cos(x) - 2.0) - y * Math.sin(x);
//        BiFunction<Double, Double, Double> d3y = (x, y) -> ( Math.cos(x) - 2.0 ) * d2y.apply(x, y) -2.0 * Math.sin(x) * d1y.apply(x, y) -y * Math.cos(x);
//        BiFunction<Double, Double, Double> d4y = (x, y) -> (Math.cos(x) - 2)*d3y.apply(x, y) - 3.0 * Math.sin(x) * d2y.apply(x, y) -3.0*Math.cos(x)*d1y.apply(x, y) +y*Math.sin(x);
//
        //Задача коши
        double x0 = 0;
        double y0 = 1;
        //N из условия а k=-2..N
        int n = 10;
        //Узлы будем считать равноотстоящими
        double h = 0.1;

        //Вывод таблицы настоящих значений
        System.out.println("|    x     |     y    |");
        System.out.println("|__________|__________|");
        for (int i = -2; i <= n; i++) {
            System.out.printf("|%10.5f|%10.5f|",x0+h*i,f.apply(x0+h*i));
            System.out.println();
        }
        //Считаем тейлора. Аргрументы начиная с -2
        System.out.println("------Тэйлор----------");
        System.out.println("Погрешность:");
        System.out.println("|    x     |     y    |Погрешность|");
        Map<Double,Double> tailorResult = tailor(d1y,d2y,d3y,d4y,x0,f.apply(x0),h);
//        Map<Double,Double> tailorResult = tailor(d1y,d2y,d3y,d4y,x0-3*h,f.apply(x0-3*h),h);
        for (Map.Entry<Double, Double> entry : tailorResult.entrySet()) {
            System.out.printf("|%10.5f|%10.5f|",entry.getKey(),entry.getValue());
            System.out.printf("%11.5f|",Math.abs(entry.getValue()-f.apply(entry.getKey())));
            System.out.println();
        }
        System.out.println("-----Эйлер----------");
        System.out.println("Погрешность: O(h^2)");
        System.out.println("|    x     |     y    |Погрешность|");
        Map<Double,Double> eilerResult = eiler(d1y,x0,f.apply(x0),h,10);
        for (Map.Entry<Double, Double> entry : eilerResult.entrySet()) {
            System.out.printf("|%10.5f|%10.5f|",entry.getKey(),entry.getValue());
            System.out.printf("%11.5f|",Math.abs(entry.getValue()-f.apply(entry.getKey())));
            System.out.println();
        }

        System.out.println("-----ЭйлерMod----------");
        System.out.println("Погрешность:O(h^2)");
        System.out.println("|    x     |     y    |Погрешность|");
        Map<Double,Double> eilerModResult = eilerMod(d1y,x0,f.apply(x0),h,10);
        for (Map.Entry<Double, Double> entry : eilerModResult.entrySet()) {
            System.out.printf("|%10.5f|%10.5f|",entry.getKey(),entry.getValue());
            System.out.printf("%11.5f|",Math.abs(entry.getValue()-f.apply(entry.getKey())));
            System.out.println();
        }
        System.out.println("-----ЭйлерКоши----------");
        System.out.println("Погрешность:");
        System.out.println("|    x     |     y    |Погрешность|");
        Map<Double,Double> eilerKoshi = eilerKoshi(d1y,x0,f.apply(x0),h,10);
        for (Map.Entry<Double, Double> entry : eilerKoshi.entrySet()) {
            System.out.printf("|%10.5f|%10.5f|",entry.getKey(),entry.getValue());
            System.out.printf("%11.5f|",Math.abs(entry.getValue()-f.apply(entry.getKey())));
            System.out.println();
        }

        System.out.println("-----Рунге-Кут----------");
        System.out.println("Погрешность:O(h^4)");
        System.out.println("|    x     |     y    |Погрешность|");
        Map<Double,Double> rungeKutResult = rungeKut(d1y,x0,f.apply(x0),h,10);
        for (Map.Entry<Double, Double> entry : rungeKutResult.entrySet()) {
            System.out.printf("|%10.5f|%10.5f|",entry.getKey(),entry.getValue());
            System.out.printf("%11.5f|",Math.abs(entry.getValue()-f.apply(entry.getKey())));
            System.out.println();
        }
        System.out.println("-----Адамса----------");
        System.out.println("Погрешность:O(h^k+2)");
        System.out.println("|    x     |     y    |Погрешность|");
        Map<Double,Double> adamsResult = adams(d1y,d2y,d3y,d4y,x0-3*h,f.apply(x0-3*h),h,13);
//        Map<Double,Double> adamsResult = adams(d1y,x0-3*h,f.apply(x0-3*h),h,13);
        for (Map.Entry<Double, Double> entry : adamsResult.entrySet()) {
            System.out.printf("|%10.5f|%10.5f|",entry.getKey(),entry.getValue());
            System.out.printf("%11.5f|",Math.abs(entry.getValue()-f.apply(entry.getKey())));
            System.out.println();
        }



    }

    public static Map<Double, Double> tailor(BiFunction<Double,Double,Double> d1y, BiFunction<Double,Double,Double> d2y,
                                             BiFunction<Double,Double,Double> d3y, BiFunction<Double,Double,Double> d4y,
                                             double x0, double y0, double h) {
//        double xm = x0;
//        double ym = y0;
        Map<Double,Double> map = new TreeMap<>();
        for (int i = -2; i <= 10; i++) {
            double xm = x0 + i*h;
//            ym = ym + h                               *  d1y.apply(xm,ym)
//                    + h * h                     *  d2y.apply(xm,ym) / 2
//                    + h * h * h           *  d3y.apply(xm,ym) / 6
//                    + h * h * h * h *  d4y.apply(xm,ym) / 24;
//            ym = ym + h                               *  d1y.apply(xm,ym)
//                    + h * h                     *  d2y.apply(xm,ym) * Math.pow(xm-x0,2) / 2.0
//                    + h * h * h           *  d3y.apply(xm,ym) * Math.pow(xm-x0,3)/ 6.0
//                    + h * h * h * h *  d4y.apply(xm,ym) * Math.pow(xm-x0,4)/ 24.0;
            double y = y0 +
                    + Math.pow((xm-x0), 1) *  d1y.apply(x0,y0)
                    + Math.pow((xm-x0), 2) *  d2y.apply(x0,y0) / 2.0
                    + Math.pow((xm-x0), 3) *  d3y.apply(x0,y0) / 6.0
                    + Math.pow((xm-x0), 4) *  d4y.apply(x0,y0) / 24.0;
//            xm += h;
            map.put(xm,y);
        }
        return map;
    }

    public static Map<Double, Double> eiler(BiFunction<Double,Double,Double> d1y, double x0, double y0, double h, int n) {
        double xm = x0;
        double ym = y0;
        Map<Double,Double> map = new TreeMap<>();
        for (int i = 1; i <= n; i++) {
            ym = ym + h * d1y.apply(xm,ym);
            xm += h;
            map.put(xm,ym);
        }
        return map;
    }

    public static Map<Double, Double> eilerMod(BiFunction<Double,Double,Double> d1y, double x0, double y0, double h, int n) {
        double xm = x0;
        double ym = y0;
        Map<Double,Double> map = new TreeMap<>();
        for (int i = 1; i <= n; i++) {
//            ym = ym + (h/2)*(d1y.apply(xm,ym)+d1y.apply(xm,ym+h*d1y.apply(xm,ym)));
            ym = ym + h * d1y.apply(xm + (h / 2), ym + (h / 2) * d1y.apply(xm, ym));
            xm += h;
            map.put(xm,ym);
        }
        return map;
    }

    //Суть
    //http://www.uchites.ru/files/nummethod_book_chapter4-1.pdf
    // ~yk+1 = yk + h * f(xk,yk)
    // yk+1 = yk + h * (f(xk,yk)+f(xk+1,~yk+1) / 2
    // xk+1 = xk + h
    // Не нашёл информации про погрешность
    public static Map<Double, Double> eilerKoshi(BiFunction<Double, Double, Double> f, double x0, double y0, double h, int N) {
        Map<Double, Double> map = new TreeMap<>();
        double ym = y0;
        double xm = x0;
        for (int i = 0; i <= N; ++i) {
            double _ym = ym + h*f.apply(xm, ym);
            ym = ym + h*(f.apply(xm, ym) + f.apply(xm, _ym))/2;
            xm += h;
            map.put(xm, ym);
        }
        return map;
    }

    public static Map<Double, Double> rungeKut(BiFunction<Double,Double,Double> d1y, double x0, double y0, double h, int n) {
        double xm = x0;
        double ym = y0;
        Map<Double,Double> map = new TreeMap<>();
        for (int i = 1; i <= n; i++) {
            double k1 = h * d1y.apply(xm,ym);
            double k2 = h * d1y.apply(xm + (h/2), ym + (k1/2));
            double k3 = h * d1y.apply(xm + (h/2), ym + (k2/2));
            double k4 = h * d1y.apply(xm + ( h ), ym + ( k3 ));
            ym = ym + (k1 + 2*k2 + 2*k3 + k4) / 6;
            xm += h;
            map.put(xm,ym);
        }
        return map;
    }


    public static Map<Double, Double> adams(BiFunction<Double,Double,Double> d1y, BiFunction<Double,Double,Double> d2y,BiFunction<Double,Double,Double> d3y,BiFunction<Double,Double,Double> d4y,double x0, double y0, double h, int n) {
        double xm = x0;
        double ym = y0;
        Map<Double,Double> map = new TreeMap<>();
        Map<Double,Double> rungeKut = rungeKut(d1y,x0,y0,h,n);


//        Map<Double,Double> rungeKut = new TreeMap<>();
//        for (int i = 1; i <= 13; i++) {
//            ym = ym + h             *  d1y.apply(xm,ym)
//                    + h * h         *  d2y.apply(xm,ym) / 2
//                    + h * h * h     *  d3y.apply(xm,ym) / 6
//                    + h * h * h * h *  d4y.apply(xm,ym) / 24;
//            xm += h;
//            rungeKut.put(xm,ym);
//        }

//        for (Map.Entry<Double, Double> aDouble : rungeKut.entrySet()) {
//            System.out.printf("%5.5f|%5.5f",aDouble.getKey(),aDouble.getValue());
//            System.out.println();
//        }
        double[] xj = rungeKut.keySet().stream().mapToDouble(Double::doubleValue).toArray();
        double[] yj = rungeKut.values().stream().mapToDouble(Double::doubleValue).toArray();
        double[][] qm = endDiff(xj,yj,(x,y) -> h * d1y.apply(x,y),n-1);
        xm += 5*h;
        for (int i = 4; i < 12; i++) {
            ym = yj[i] +       qm[0][i]
                    + (        qm[1][i-1] / 2.0)
                    + (5.0 *   qm[2][i-2] / 12.0)
                    + (3.0 *   qm[3][i-3] / 8.0)
                    + (251.0 * qm[4][i-4] / 720.0);
            map.put(xj[i]+h,ym);
        }
        return map;
    }

    public static double[][] endDiff(double[] xj, double[] yj, BiFunction<Double,Double, Double> f, int m) {
        double[][] result = new double[m+1][];
        result[0] = new double[m+1];
        /*
         * Заполняем массив значений от xi
         */
        for (int j = 0; j < m+1; j++) {
            result[0][j] = f.apply(xj[j],yj[j]);
        }
        /*
        * Строим таблицу конечных разностей
        * */
        for (int i = 1; i <= m; i++) {
            result[i] = new double[m - i + 1];
            for (int j = 0; j < (m - i + 1); j++) {
                result[i][j] = (result[i - 1][j + 1] - result[i - 1][j]);
            }
        }
        return result;
    }
}
