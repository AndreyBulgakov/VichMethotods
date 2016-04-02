package sem5.Task4;

import sem5.Task1.Task1;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.function.Function;

/**
 * Created by andrey on 23.11.15.
 */
public class Task4 {
    public static void main(String[] args) throws IOException {
        Function<Double,Double> f = (x) -> {return Math.log(1+x)-Math.exp(x);};//Math.sin(x) - x / 2;};
        Function<Double,Double> fdx = (x) -> {return 1/(x+1)-Math.exp(x);};//Math.sin(x) - x / 2;};
        Function<Double,Double> fd2x = (x) -> {return -Math.exp(x)-1/((x+1)*(x+1));};//Math.sin(x) - x / 2;};

        //Ввод параметров
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        System.out.print("Ведите m :");
        int m = Integer.valueOf(reader.readLine());
        System.out.print("Введите a:");
        double a = Double.valueOf(reader.readLine());
        System.out.print("Введите b:");
        double b = Double.valueOf(reader.readLine());
        double h = (b - a) / m;

        //Узлы
        System.out.println("-------Считаю массив узлов-------");
        double[] points = createPoints(m,a,b);
        double[] pointsOrdered = createPoints(m,a,b);
        //Значения
        double[] valuesOrdered = new double[pointsOrdered.length];
        for (int i = 0; i < pointsOrdered.length; i++) {
            valuesOrdered[i]=f.apply(pointsOrdered[i]);
        }

        System.out.println("-------Строю таблицу узлов и значений-------");
        double[] values = printTable(m,points,f);
        int n;
        do{
            System.out.print("Введите степень многочлена n<=m:");
            n = Integer.valueOf(reader.readLine());
        }while(n>m);

//        System.out.println("-------Производные-----------");
//        derivation1(pointsOrdered,valuesOrdered,n,m,h,f,fdx,fd2x);

        do {
            //Выберите точку интеполяции (y)
            System.out.print("Введите точку интерполирования");
            double y = Double.parseDouble(reader.readLine());
            System.out.println("----Сортирую по дальности расположения узлов от точки----");
            values = sort(values,points,y);
            System.out.println("-------Первый способ-------");
            double res1 = newton(n,y,values,points);
            System.out.println("X=Qn(F) : "+res1);
            System.out.println("Величина невязки:"+ Math.abs(f.apply(res1)-y));
            System.out.println("-------Второй способ-------");
            final int nf = n;
            final double[] finalValues = values;

            Function<Double,Double> func = (x -> newton(nf,x,points, finalValues)-y);
//            System.out.println("Введите eps");
            double eps = 0.00000001;//Double.parseDouble(reader.readLine());
            System.out.println("eps="+eps);
            List<Task1.Point> roots = Task1.createPoints(func,a,b,h);
            for (Task1.Point root : roots) {
                bisection(func,root.getX(), root.getY(),eps);
            }

            System.out.println("-------Производные-----------");
            derivation1(pointsOrdered,valuesOrdered,n,m,h,f,fdx,fd2x);
            System.out.println("Вы хотите выйти y/n:");
        }
        while (!reader.readLine().equals("y"));
    }

    /*
    * Дифференцирование
    * */
    public static void derivation1(double points[], double[] values, int n, int m, double h, Function<Double,Double> fx,Function<Double,Double> fdx, Function<Double,Double> fd2x){
        System.out.println("__________________________________________________________________________________________");
        System.out.println("| i|    xi    |   f(xi)  | f'(xi)ЧД | |f'(xi)T-f(xi)ЧД| | f''(xi)ЧД ||f''(xi)T-f''(xi)ЧД||");
        System.out.println("|__|__________|__________|__________|___________________|___________|____________________|");
        String formatter = "|%2d|%10.5f|%10.5f|%10.5f|%19.5f|%11.5f|%20.5f|";
        Derivation derivation = new Derivation(values,h,m);
        for (int i = 0; i <= m; i++) {
            double res = derivation.derivate1(i);
            double res2 = derivation.derivate2(i);
            System.out.printf(formatter,
                    i,
                    points[i],
                    values[i],
                    res,
                    Math.abs(fdx.apply(points[i])-res),
                    res2,
                    Math.abs(fd2x.apply(points[i])-res2)
            );
            System.out.println();
        }
        System.out.println("__________________________________________________________________________________________");
    }

    /*
    * Тут все просто.
    * Делим отрезок на два и проверяем погрешность.
    * Выбираем нужный отрезок.
    * */
    public static void bisection(Function<Double,Double> f, double a, double b, double eps){
        double xl = a;
        double xr = b;
        double xd = b-a;
        int n = 0;
        while (Math.abs(f.apply(xl))>eps || Math.abs(f.apply(xr))>eps){
            n++;
            xd = xd/2.0;
            double xm =xl+xd;
            if (f.apply(xl)*f.apply(xm)<0.0)
                xr = xm;
            else
                xl = xm;
        }
        System.out.println("Заданная точность: " + eps);
        System.out.println("Начальное приближение x1: " + a);
        System.out.println("Начальное приближение x2: " + b);
        System.out.println("Количество шагов n: " + n);
        System.out.println("Приближенное решение: " + xl);
        System.out.println("Величина невязки: " + Math.abs(f.apply(xl)));
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
    static double[] printTable(int m, double[] points, Function<Double, Double> f){
        double[] values = new double[points.length];
        System.out.println(" i|___xi__|__y___|");
        for (int i = 0; i < m+1; i++) {
            System.out.println(String.format("%2d|%5.5f|%5.5f",i,points[i],f.apply(points[i])));
            values[i] = f.apply(points[i]);
        }
        return values;
    }
    /*
    * Сортировка по дальности удаления от точки
    * */
    public static double[] sort(double[] points, double[] values, double polpoint){
        Set<Double> sorted = new TreeSet<>((o1, o2) -> Math.abs(polpoint-o1)<Math.abs(polpoint-o2)?-1:1);
        Arrays.stream(points).forEach(value -> sorted.add(value));
        double[] resultpoints = sorted.stream().mapToDouble(Double::doubleValue).toArray();
        double[] resultvalues = new double[values.length];
        for (int i = 0; i < resultpoints.length; i++) {
            for (int j = 0; j < points.length; j++) {
                if(resultpoints[i]==points[j])
                {
                    resultvalues[i] = values[j];
                }
            }
        }
        for (int i = 0; i < resultvalues.length; i++) {
                values[i]=resultvalues[i];
        }
        return resultpoints;
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
}

