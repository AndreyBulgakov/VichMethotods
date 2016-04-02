package sem5.Task5;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.function.Function;

/**
 * Created by andrey on 01.12.15.
 */
public class Task5 {
    public static void main(String[] args) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
//        Function<Double, Double> f = (x) -> 2 * Math.sin(x) - x / 2;
        Function<Double, Double> f = (x) -> 2 * Math.sin(x) - x / 2;
        Function<Double,Double> j = (x) -> -(x*x)/4-2*Math.cos(x);
//        Function<Double, Double> f = (x) -> (- x * Math.log(x))*Math.sin(x);
//        Function<Double,Double> j = (x) -> (x*x*x);

        System.out.print("Введите a:");
        double a = Double.parseDouble(reader.readLine());
        System.out.print("Введите b:");
        double b = Double.parseDouble(reader.readLine());
        System.out.print("Введите m:");
        double m = Double.parseDouble(reader.readLine());
        double h = (b - a) / m;
        double jx = j.apply(b) - j.apply(a);
        System.out.println("J=" + jx);
        System.out.println("------Формула левых прямоугольников------");
        double lefttrRes = leftRectangle(f, a, b, h,m);
        System.out.println("J(h)=" + lefttrRes);
        System.out.println("|J-J(h)|=" + Math.abs(jx-lefttrRes));

        System.out.println("\n------Формула средних прямоугольников------");
        double midtrRes =  midRectangle(f, a, b, h,m);
        System.out.println("J(h)=" + midtrRes);
        System.out.println("|J-J(h)|=" + Math.abs(jx-midtrRes));

        System.out.println("\n------Формула правых прямоугольников------");
        double righttrRes = rightRectangle(f, a, b, h,m);
        System.out.println("J(h)=" + righttrRes);
        System.out.println("|J-J(h)|=" + Math.abs(jx-righttrRes));

        System.out.println("\n------Формула трапеций------");
        double trapRes = trap(f, a, b, h, m);
        System.out.println("J(h)=" + trapRes);
        System.out.println("|J-J(h)|=" + Math.abs(jx-trapRes));

        System.out.println("\n------Формула Симпсона------");
        double simpsonRes = simpson(f, a, b, h, m);
        System.out.println("J(h)=" + simpsonRes);
        System.out.println("|J-J(h)|=" + Math.abs(jx-simpsonRes));
    }


    public static double leftRectangle(Function<Double, Double> f, double a, double b, double h, double m){
        double res = 0;
        for (int k = 1; k <= m; k++) {
            res += f.apply(a+(k-1)*h);
        }
        return res*h;
    }

    public static double midRectangle(Function<Double, Double> f, double a, double b, double h, double m){
        double res = 0;
        for (int k = 1; k <= m; k++) {
            res += f.apply((a+h/2)+(k-1)*h);
        }
        return res*h;
    }

    public static double rightRectangle(Function<Double, Double> f, double a, double b, double h, double m){
        double res = 0;
        for (int k = 1; k <= m; k++) {
            res += f.apply(a+h+(k-1)*h);
        }
        return res*h;
    }

    public static double trap(Function<Double, Double> f, double a, double b, double h, double m){
        double res = 0;
        for (int k = 1; k <= m-1; k++) {
            res += f.apply(a+k*h);
        }
        res *= 2;
        res += f.apply(a) + f.apply(b);
        return res*(b-a)/(2*m);
    }

    public static double simpson(Function<Double, Double> f, double a, double b, double h, double m){
        double res = 0;
        double subres1 = 0;
        double subres2 = 0;
        h /= 2;

        for (int k = 1; k <= 2*m-1; k++) {
            if (k % 2 != 0)
            subres1 += f.apply(a+k*h);
        }

        for (int k = 2; k <= 2*m-2; k++) {
            if (k % 2 == 0)
            subres2 += f.apply(a+k*h);
        }
        double fx0 = f.apply(a);
        double fxn = f.apply(b);
//        res +=  (Double.isNaN(fx0)?f.apply(a+h):fx0) + 4 * subres1 + 2 * subres2 + (Double.isNaN(fxn)?f.apply(a+h):fxn);
        res +=  (Double.isNaN(fx0)?0:fx0) + 4 * subres1 + 2 * subres2 + (Double.isNaN(fxn)?0:fxn);
        return res*(b-a)/(6*m);
    }
}
