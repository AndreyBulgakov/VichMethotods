package sem6.task6;

import java.util.function.DoubleFunction;
import java.util.function.Function;
import static sem6.task6.Task6.derivate;
/**
 * Created by andrey on 5/5/16.
 */
public class Ly_f {
    Function<Double, Double> p;
    Function<Double, Double> r;
    Function<Double, Double> f;
    double a1;
    double a2;
    double b1;
    double b2;
    TaskClass taskClass;

    public Ly_f(Function<Double, Double> p, Function<Double, Double> r, Function<Double, Double> f, double a1, double a2, double b1, double b2) {
        this.p = p;
        this.r = r;
        this.f = f;
        this.a1 = a1;
        this.a2 = a2;
        this.b1 = b1;
        this.b2 = b2;

        if (Math.abs(a1) + Math.abs(a2) == 0 || Math.abs(b1) + Math.abs(b2) ==0)
            throw new IllegalArgumentException("a1 = a2 = 0 or b1 = b2 = 0");
        if (a2 == 0)
            taskClass = TaskClass.I;
        else if (a1 == 0)
            taskClass = TaskClass.II;
        else
            taskClass = TaskClass.III;

    }

    public Function<Double, Double> Ly(Function<Double, Double> y,Function<Double, Double> dy,Function<Double, Double> ddy){
        double h = 1e-5;
        Function<Double, Double> res = x ->
                -(derivate(p,h).apply(x) * dy.apply(x) + (p.apply(x)) * ddy.apply(x))
                        + r.apply(x) * y.apply(x);
        return res;
    }
}
