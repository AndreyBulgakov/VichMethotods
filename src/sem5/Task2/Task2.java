package sem5.Task2;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;

public class Task2 {

    public static void main(String[] args) throws IOException {

        Function<Double,Double> f = (x) -> {return Math.log(1+x)-Math.exp(x);};
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
        System.out.print("Введите точку интерполиования:");
        double startPoint = Double.valueOf(reader.readLine());

        System.out.println(String.format("Значение f(%f)=%f", startPoint, f.apply(startPoint)));

        System.out.println("----Сортирую по дальности расположения узлов от точки----");
        points = sort(points,startPoint);
        //Создаем массив значений
        double[] values = new double[points.length];
        for (int i = 0; i < values.length; i++) {
            values[i] = f.apply(points[i]);
        }
        int n;
        do{
            System.out.print("Введите степень многочлена n<=m:");
            n = Integer.valueOf(reader.readLine());
        }
        while(n>m);
        System.out.println(String.format("Значение f(%f)=%f", startPoint, f.apply(startPoint)));
        System.out.println("Начинаю вычисление методом Ньютона");
        System.out.print("Получено значение:");
        double Newton = newton(n,startPoint,points,values);
        System.out.println(Newton);
        System.out.println("Абсолютная погрешность:"+Math.abs(f.apply(startPoint)-Newton));
        System.out.println("Начинаю вычисление методом Лагранжа");
        System.out.print("Получено значение:");
        double Lagranz = lagranz(n,startPoint,points,values);
        System.out.println(Lagranz);
        System.out.println("Абсолютная погрешность:"+Math.abs(f.apply(startPoint)-Lagranz));

        //Test
        /*double[] points = createPoints(15,0.4,0.9);
        Function<Double,Double> f = (x) -> {return Math.sin(x) + x*x;};
        printTable(15,points,f);
        points = sort(points,0.75);
        //Создаем массив значений
        double[] values = new double[points.length];
        for (int i = 0; i < values.length; i++) {
            values[i] = f.apply(points[i]);
        }
        int n = 7;
        newton(n, 0.75, points,values);
        lagranz(n, 0.75, points, values);
*/
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
    * Сортировка по дальности удаления от точки
    * */
    static double[] sort(double[] points, double polpoint){
        Set<Double> sorted = new TreeSet<>((o1, o2) -> Math.abs(polpoint-o1)<Math.abs(polpoint-o2)?-1:1);
        Arrays.stream(points).forEach(value -> sorted.add(value));
        return sorted.stream().mapToDouble(Double::doubleValue).toArray();
    }

    /*
    * Разделенная разность. Тут пока не используется.
    * */

    static double divdiff(int n, double[] x, double[]y){
        double f = 0;
        for (int j = 0; j < n; j++) {
            double den = 1;
            for (int i = 0; i < n; i++) {
                if(j!=i) den *= x[j]-x[i];
            }
            f += y[j]/den;
        }
        return f;
    }

    /*
    * Метод ньютона.
    * На вход талица иксов и соответсвуюящая вычисленная таблица значений.
    * Строим таблицу как тут http://pmpu.ru/vf4/interpolation
    * [xj,xk]= yj-yk
    *          xj-xk
    * Крайние игрики на самые крайние значения иксов
    *
    * */
    public static double newton(int n, double x0, double[] x, double[] y){
        double fx0 = 0;
        double[][] table = new double[n][];
        table[0] = new double[n];
        for (int i = 0; i < n; i++) {
            table[0][i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
        }
        for (int i = 1; i < n; i++) {
            table[i] = new double[n-i];
            for (int j = 0; j < (n - i); j++) {
                table[i][j] = (table[i-1][j+1]-table[i-1][j])/(x[j+i+1]-x[j]);
            }
        }

        for (int i = 0; i < n; i++) {
            double den = 1;
            for (int j = 0; j < i+1; j++) {
                den *= x0-x[j];
            }
            fx0 += table[i][0]*den;
        }

        return fx0+y[0];
    }
    /*
    * Тут можно посмотреть
    * http://aco.ifmo.ru/el_books/numerical_methods/lectures/glava3.html
    * Z(0..n)yi*П(1..n, i!=j)x-xj
    *                        xi-xj
    * */
    static double lagranz(int n, double x0, double[] x, double[] y){
        double res = 0;
        for (int i = 0; i < n; i++) {
            double li = 1;
            for (int j = 0; j < n; j++) {
                if (j!=i) li *= (x0-x[j])/(x[i]-x[j]);
            }
            res +=y[i]*li;
        }
        return res;
    }

}
