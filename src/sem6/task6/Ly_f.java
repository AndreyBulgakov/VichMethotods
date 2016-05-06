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

    public Function<Double, Double> Ly(Function<Double, Double> y){
        double h = 1e-16;
//        Function<Double, Double> dy = Task6.derivate(y, h);
//        Function<Double, Double> dpdy = Task6.derivate(x -> p.apply(x)*dy.apply(x), h);
//        Function<Double, Double> res = x -> -dpdy.apply(x) + r.apply(x)*y.apply(x);
        Function<Double, Double> res = x ->
                - derivate(j -> p.apply(j)*derivate(y, h).apply(j), h).apply(x) +r.apply(x)*y.apply(x);
        return res;
    }
}
