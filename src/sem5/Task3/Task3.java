package sem5.Task3;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.function.Function;

/**
 * Created by andrey on 23.11.15.
 */
public class Task3 {
    public static void main(String[] args) throws IOException {

        Function<Double,Double> f = (x) -> {return 2 * Math.sin(x) - x / 2;};
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        System.out.print("Ведите m :");
        int m = Integer.valueOf(reader.readLine());
        System.out.print("Введите a:");
        double a = Double.valueOf(reader.readLine());
        System.out.print("Введите b:");
        double b = Double.valueOf(reader.readLine());
        System.out.println("-------Считаю массив узлов-------");
        double[] points = createPoints(m,a,b);
        System.out.println("-------Строю таблицу узлов и значений-------");
        printTable(m,points,f);
        int n;
        do{
            System.out.print("Введите степень многочлена n<=m:");
            n = Integer.valueOf(reader.readLine());
        }
        while(n>m);
        double h = (b - a) / m;
        do {
            System.out.println("Введите точку из начала, конца или середины");
            System.out.printf("\nПромежутки: [%.2f;%.2f] или [%.2f;%.2f] или [%.2f;%.2f]",
                    a, a + h, a + h + h, b - h - h, b - h, b);
            double x = Double.parseDouble(reader.readLine());
            double[][] endDiff = endDiff(points, f, m);
            if (x >= a && x <= a + h) {
                System.out.println("-------Вычисляю по началу Ньютоном-------");
                double res = newtonBeg(n, x, points, f.apply(points[0]), endDiff, h);
                System.out.println("P(x)=" + res);
                System.out.println("f(x)=" + f.apply(x));
                System.out.println("|P(x)-f(x)|=" + Math.abs(res - f.apply(x)));

            }
            if (x >= a + h + h && x <= b - h - h) {
                System.out.println("-------Вычисляю по середине Гаусом-------");
                double res = gausMid(n, x, points, f.apply(points[0]), endDiff, h);
                System.out.println("P(x)=" + res);
                System.out.println("f(x)=" + f.apply(x));
                System.out.println("|P(x)-f(x)|=" + Math.abs(res - f.apply(x)));
            }
            if (x >= b - h && x <= b) {
                System.out.println("-------Вычисляю конец таблицы Ньютоном-------");
                double res = newtonEnd(n, x, points, f.apply(points[0]), endDiff, h);
                System.out.println("P(x)=" + res);
                System.out.println("f(x)=" + f.apply(x));
                System.out.println("|P(x)-f(x)|=" + Math.abs(res - f.apply(x)));
            }
            System.out.println("Вы хотите выйти y/n:");
        }
        while (!reader.readLine().equals("y"));
    }

    /*
    * Метод Ньютона для начала таблицы
    * Прямая формула интерполирования Ньютона
    * Как на вики
    * */
    static double newtonBeg(int n, double x, double[] points, double y, double[][]table, double h){
        double q = (x - points[0])/h;
        double res = 0;
        double den;
        for (int i = 0; i < n; i++) {
            den=1;
            for (int j = 1; j < i+1; j++) {
                den *= (q-j+1);
            }
            //Отличие от формулы для конца только здесь
            res += table[i][0]*den/factorial(i);
        }
        return res;
    }
    /*
    * Метод Ньютона для конца таблицы
    * Прямая формула интерполирования Ньютона
    * Как на вики
    * */
    static double newtonEnd(int n, double x, double[] points, double y, double[][]table, double h){
        double q = (x - points[points.length-1])/h;
        double res = 0;
        double den;
        for (int i = 0; i < n; i++) {
            den=1;
            for (int j = 1; j < i+1; j++) {
                den *= (q+j-1);
            }
            res += table[i][points.length-1-i]*den/factorial(i);
        }
        return res;
    }

    /*
    * Метод Гауса для середины таблицы
    * */
    static double gausMid(int n, double x, double[] points, double y, double[][]table, double h){
        int l = 0;
        for (int i = 0; i < points.length - 1; i++) {
            if(points[i]<=x && x<=points[i+1]) {
                l = i;
                break;
            }
        }
        double q = (x - points[l])/h;
        double res = 0;
        for (int i = 1; i <= n; i++) {
            double den = 1;
            for (int j = 0; j <= i-1; j++) {
                den *= (q + Math.pow(-1,j)*Math.floor((double)(j+1)/2)) ;
            }
            try {
                res += table[i][l-(int)Math.floor((double)(i/2))] * den / factorial(i);
            }catch (Exception e){}
        }
        return res+table[0][l];
    }

    /*
    * Создаем массив точек.
    * */
    static double[] createPoints(int m, double a, double b){
        double[] points = new double[m + 1];
        for (int j = 0; j <= m; j++)
        {
            points[j] = a + (j * (b - a)) / m;
        }
        return points;
    }

    /*
    * Печатаем таблицу.
    * */
    static void printTable(int m, double[] points, Function<Double, Double> f){
        System.out.println(" i|___xi__|__y___|");
        for (int i = 0; i < m+1; i++) {
            System.out.println(String.format("%2d|%5.5f|%5.5f",i,points[i],f.apply(points[i])));
        }
    }

    /*
    * Возвращаем матрицу конечных разностей
    * */
    public static double[][] endDiff(double[] nodes, Function<Double, Double> f, int m) {
        double[][] result = new double[m+1][];
        result[0] = new double[m+1];
        /*
         * Заполняем массив значений от xi
         */
        for (int j = 0; j < m+1; j++) {
            result[0][j] = f.apply(nodes[j]);
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
    static long factorial(long i) {
        if(i == 0) {
            return 1;
        }
        else
            return i * factorial(i-1);
    }
}

