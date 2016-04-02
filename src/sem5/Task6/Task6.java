package sem5.Task6;

import sem5.Task1.Task1;
import sem5.Task5.Task5;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

/**
 * Created by andrey on 07.12.15.
 */
public class Task6 {


    public static void main(String[] args) throws IOException {

        Function<Double, Double> f = Math::sin;
        Function<Double, Double> w =  x -> - x * Math.log(x);

        System.out.print("Введите a:");
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        double a = Double.parseDouble(reader.readLine());
        System.out.print("Введите b:");
        double b = Double.parseDouble(reader.readLine());
        System.out.print("Введите m:");
        double m = Double.parseDouble(reader.readLine());
        double h = (b - a) / m;

        //Вычисляем моменты весовой функции
        double m0 = Task5.simpson(w,a,b,h,m);
        double m1 = Task5.simpson((y -> w.apply(y)*y),a,b,h,m);
        double m2 = Task5.simpson((y -> w.apply(y)*y*y),a,b,h,m);
        double m3 = Task5.simpson((y -> w.apply(y)*y*y*y),a,b,h,m);
        //Коэффициенты
        double a1 = (m0*m3-m2*m1)/(m1*m1-m2*m0);
        double a2 = (m2*m2-m3*m1)/(m1*m1-m2*m0);

    //        IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(2,h,h,(int)m,(int)m+1);
    //        System.out.println(integrator.integrate(1000000,x-> Math.sin(x)*(-x*Math.log(x)

        //Многочлен
        Function<Double, Double> task = x -> x*x+a1*x+a2;
        System.out.println("Момент весовой функции: "+m0);
        System.out.println("Момент весовой функции: "+m1);
        System.out.println("Момент весовой функции: "+m2);
        System.out.println("Момент весовой функции: "+m3);

        System.out.println("Вид ортогонального многочлена: x^2+"+a1+"x+"+a2+"=0");

        //Ищем корни (узлы)
        List<Task1.Point> roots = Task1.createPoints(task,a,b,h);
        ArrayList<Double> doubles = new ArrayList<>();
        BisectionSolver solver = new BisectionSolver(h);
        for (Task1.Point root : roots) {
            doubles.add(solver.solve(10000,x->x*x+a1*x+a2,root.getX(),root.getY()));
            //doubles.add(bisection(task,root.getX(), root.getY(),0.000000001));
        }
        double A1 = (m1-doubles.get(1)*m0)/(doubles.get(0)-doubles.get(1));
        double A2 = (m1-doubles.get(0)*m0)/(doubles.get(1)-doubles.get(0));

        System.out.println("Узел КФ НАСТ:" + doubles.get(0));
        System.out.println("Узел КФ НАСТ:" + doubles.get(1));
        System.out.println("Коэффициент КФ НАСТ:"+A1);
        System.out.println("Коэффициент КФ НАСТ:"+A2);
        System.out.println((A1+A2)+":::"+m0);

        double res = A1 * f.apply(doubles.get(0))+A2*f.apply(doubles.get(1));
        System.out.println("Результат:"+res);
//        IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(2,h,h,(int)m,(int)m+1);
//        System.out.println(integrator.integrate(1000000,x-> Math.sin(x)*(-x*Math.log(x)),a,b));
//        System.out.println(Task5.simpson(x->Math.sin(x)*(-x*Math.log(x)),a,b,h,m));
    }

    public static double bisection(Function<Double,Double> f, double a, double b, double eps){
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
        return xl;
    }

}
