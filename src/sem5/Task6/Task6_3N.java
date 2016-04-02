package sem5.Task6;

import sem5.Task1.Task1;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;
import org.apache.commons.math3.linear.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

/**
 * Created by andrey on 07.12.15.
 */
public class Task6_3N {


    public static void main(String[] args) throws IOException {

        Function<Double, Double> f = Math::sin;
        Function<Double, Double> w =  x -> - x * Math.log(x);

        System.out.print("Введите a:");
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        double a = Double.parseDouble(reader.readLine());
        System.out.print("Введите b:");
        double b = Double.parseDouble(reader.readLine());
//        System.out.print("Введите m:");
//        double m = Double.parseDouble(reader.readLine());
        double h = (b - a) / 100;

        //Вычисляем моменты весовой функции
        //Вычисленно руками (а точнее в вольфраме)
        //Это взят интеграл пока без промежутка
        Function<Double,Double> v0 = x -> Math.pow(x,2) * (1.0-2.0*Math.log(x))/4.0;
        Function<Double,Double> v1 = x -> Math.pow(x,3) * (1.0-3.0*Math.log(x))/9.0;
        Function<Double,Double> v2 = x -> Math.pow(x,4) * (1.0-4.0*Math.log(x))/16.0;
        Function<Double,Double> v3 = x -> Math.pow(x,5) * (1.0-5.0*Math.log(x))/25.0;
        Function<Double,Double> v4 = x -> Math.pow(x,6) * (1.0-6.0*Math.log(x))/36.0;
        Function<Double,Double> v5 = x -> Math.pow(x,7) * (1.0-7.0*Math.log(x))/49.0;
        //Вычисляем по промежутку
        double m0 = v0.apply(b)-v0.apply(a+0.000000000001);
        double m1 = v1.apply(b)-v1.apply(a+0.000000000001);
        double m2 = v2.apply(b)-v2.apply(a+0.000000000001);
        double m3 = v3.apply(b)-v3.apply(a+0.000000000001);
        double m4 = v4.apply(b)-v4.apply(a+0.000000000001);
        double m5 = v5.apply(b)-v5.apply(a+0.000000000001);

        //Считаем коэфициенты ai
        //Их разрешалось вычислять используя сторонние либы
        //Кусок кода ниже ищет аi на основе моментов
        //Формулы смотреть в тетради
        //a1m2 + a2m1 + a3m0 = - m3
        //a1m3 + a2m2 + a3m1 = - m4
        //a1m4 + a2m3 + a3m2 = - m5
        RealMatrix coefficients = new Array2DRowRealMatrix(
                new double[][] {
                        { m2, m1, m0 },
                        { m3, m2, m1 },
                        { m4, m3, m2 } },
                        false);
        DecompositionSolver matrixSolver = new LUDecomposition(coefficients).getSolver();
        RealVector constants = new ArrayRealVector(new double[] { -m3, -m4, -m5 }, false);
        RealVector solution = matrixSolver.solve(constants);
       // System.out.println(solution.getEntry(0)+" "+ solution.getEntry(1)+" "+ solution.getEntry(2));

        //Коэфициенты ai
        //Просто вытаскиваем посчитанные
        double a1 = solution.getEntry(0);
        double a2 = solution.getEntry(1);
        double a3 = solution.getEntry(2);

        //Многочлен
        //Строим многочлен x^3 + a1*x^2 + a2*x + a3
        Function<Double, Double> task = x -> x*x*x+a1*x*x+a2*x+a3;
        System.out.println("Момент весовой функции 0: "+m0);
        System.out.println("Момент весовой функции 1: "+m1);
        System.out.println("Момент весовой функции 2: "+m2);
        System.out.println("Момент весовой функции 3: "+m3);
        System.out.println("Момент весовой функции 4: "+m4);
        System.out.println("Момент весовой функции 5: "+m5);

        System.out.println("Вид ортогонального многочлена: x^3+"+a1+"x^2+"+a2+"x"+a3+"=0");

        //Ищем корни (узлы)
        // Для поиска этого многочлена мы ранее вводили m
        List<Task1.Point> roots = Task1.createPoints(task,a,b,h);
        ArrayList<Double> doubles = new ArrayList<>();
        BisectionSolver bisectionSolver = new BisectionSolver(h);
        for (Task1.Point root : roots) {
            //doubles.add(bisectionSolver.solve(10000, x -> x*x*x+a1*x*x+a2*x+a3, root.getX(),root.getY()));
            doubles.add(bisection(task,root.getX(), root.getY(),0.000000001));
        }

        //Получаем узлы. По которым строится КФНАСТ
        double x1 = doubles.get(0);
        double x2 = doubles.get(1);
        double x3 = doubles.get(2);

        //Коэффициенты Ак посчитаны руками в тетраде.
        //Док-во того что эту задачу я не скатал.
        double A1 = (m2 - x3 * m1 - x2 * m1 + x2 * x3 * m0)/(x1 * x1 - x1 * x2 - x1 * x3 + x2 * x3);
        double A2 = (m2 - x3 * m1 - x1 * m1 + x1 * x3 * m0)/(x2 * x2 - x1 * x2 - x2 * x3 + x1 * x3);
        double A3 = (m2 - x1 * m1 - x2 * m1 + x1 * x2 * m0)/(x3 * x3 - x1 * x3 - x2 * x3 + x1 * x2);

        //Финальный вывод наших вычислений
        System.out.println("Узел КФ НАСТ:" + x1);
        System.out.println("Узел КФ НАСТ:" + x2);
        System.out.println("Узел КФ НАСТ:" + x3);
        System.out.println("Коэффициент КФ НАСТ A1:"+A1);
        System.out.println("Коэффициент КФ НАСТ A2:"+A2);
        System.out.println("Коэффициент КФ НАСТ A3:"+A3);

        //Проверка того что Ак вычисленны правильно и равны нулевому моменту
        // Смотри свойства коффициентов КФНАСТ
        System.out.println((A1+A2+A3)+"="+m0);

        //Ну и КФНАСТ в Действии
        double res = A1 * f.apply(x1) + A2 * f.apply(x2) + A3 * f.apply(x3);
        System.out.println("Результат:"+res);
//        IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(2,h,h,(int)m,(int)m+1);
//        System.out.println(integrator.integrate(1000000,x-> Math.sin(x)*(-x*Math.log(x)),a,b));
//        System.out.println(Task5.simpson(x->Math.sin(x)*(-x*Math.log(x)),a,b,h,m));
    }

    //Метод бисекций использованный когда-то давно
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
