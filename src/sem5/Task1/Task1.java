package sem5.Task1;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

/**
 * Created by andrey
 */
public class Task1 {

    public static void main(String[] args) {
        Function<Double,Double> f = new f();
        Function<Double,Double> df = new df();
        double a = -5.0;
        double b = 5.0;
        double eps = 0.00000001;
        double h = (b-a) / 10000.0;
        List<Point> points = createPoints(f,a,b,h);
        System.out.println("--------Бисекций-------------");
        for (Point point : points) {
            bisection(f, point.x, point.y, eps);
        }
        System.out.println("--------Ньютон-------------");
        for (Point point : points) {
            newton(f,df, point.x, eps);
        }
        System.out.println("---------МодНьютон------------");

        for (Point point : points) {
            newtonMod(f, df, point.x, point.x, eps);
        }
        System.out.println("----------Хорд-----------");

        for (Point point : points) {
            chord(f, point.x, point.y, eps);
        }
    }
    /*
    * Отделяем корни с помощию табулирования.
    *
    * */
    public static List<Point> createPoints(Function<Double, Double> f, double a, double b, double tab){
        ArrayList<Point> points = new ArrayList<>();
        for (double i = 1.0; i < 10000.0 && a+i*tab<b; i++) {
            if (f.apply(a + (i - 1.0)*tab) * f.apply(a + i*tab) < 0.0){
                points.add(new Point(a + (i - 1)*tab, a + i*tab));
            }
        }
        return points;
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
    * x Заданная точность.
    * Пусть вычисленно xi вычислим xi+1
    * Проведем касательную f(x) в точке xi и найдем точку пересечения этой касательной с точкой абсцисс
    * xi+1 положим найденной точке и повторим
    * Учитываем погрешность xi+1-xi
    * */
    public static void newton(Function<Double,Double> f, Function<Double,Double> df, double x,  double eps){
        double xi = x;
        double c = xi-(f.apply(xi)/df.apply(xi));
        int n = 0;
        while (Math.abs(c-xi)>eps){
            n++;
            xi = c;
            c = xi-(f.apply(xi)/df.apply(xi));
        }
        System.out.println("Заданная точность: " + eps);
        System.out.println("Начальное приближение: " + x);
        System.out.println("Количество шагов n: " + n);
        System.out.println("Приближенное решение: " + c);
        System.out.println("Величина невязки: " + Math.abs(f.apply(c)));
    }

    /*
    * x Заданная точность.
    * Пусть вычисленно xi вычислим xi+1
    * Проведем касательную f(x) в точке xi и найдем точку пересечения этой касательной с точкой абсцисс
    * xi+1 положим найденной точке и повторим
    * Учитываем погрешность xi+1-xi
    * Тоже самое только что f'(x0) не меняем и вычисляем один раз
    * */
    public static void newtonMod(Function<Double,Double> f, Function<Double,Double> df, double x0, double x,  double eps){
        double start = df.apply(x0);
        double xi = x;
        double c = xi-(f.apply(xi)/start);
        int n = 0;
        while (Math.abs(c-xi)>eps){
            n++;
            xi = c;
            c = xi-(f.apply(xi)/start);
        }
        System.out.println("Заданная точность: " + eps);
        System.out.println("Начальное приближение: " + x);
        System.out.println("Точка дифференцирования: " + x0);
        System.out.println("Количество шагов n: " + n);
        System.out.println("Приближенное решение: " + c);
        System.out.println("Величина невязки: " + Math.abs(f.apply(c)));
    }
    /*
    * Формулка
    * xi+1 = xi- f(xi)*(xi-x0))
*                 f(xi)-f(x0)
    *
    * */
    public static void chord(Function<Double,Double> f, double a, double b, double eps){
        double x0 = a;
        double x1 = b;
        double x2 = 0.0;
        int n = 0;
        while (Math.abs(x1 - x0)>eps){
            n++;
            x2 = x1 - f.apply(x1)*(x1-x0)/(f.apply(x1)-f.apply(x0));
            x0 = x1;
            x1 = x2;
        }
        System.out.println("Заданная точность: " + eps);
        System.out.println("Начальное приближение1: " + a);
        System.out.println("Начальное приближение2: " + b);
        System.out.println("Количество шагов n: " + n);
        System.out.println("Приближенное решение: " + x2);
        System.out.println("Величина невязки: " + Math.abs(f.apply(x2)));
    }

     static class f implements Function<Double, Double>{
        @Override
        public Double apply(Double x) {
            return x*x*x - Math.cos(x + 0.5);
        }
    }

    static class df implements Function<Double, Double>{
        @Override
        public Double apply(Double x) {
            return 3.0*x*x + Math.sin(x + 0.5);
        }
    }

    public static class Point{
        double x;
        double y;

        public Point(double x, double y) {
            this.x = x;
            this.y = y;
        }

        public double getX() {
            return x;
        }

        public double getY() {
            return y;
        }
    }
}
