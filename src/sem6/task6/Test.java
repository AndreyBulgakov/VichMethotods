package sem6.task6;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialsUtils;

import java.util.function.Function;

/**
 * Created by andrey on 5/6/16.
 */
public class Test {

    public static void main(String[] args) {
        int n = 5;
        //Заполняем wj
        Function<Double, Double>[] w = new Function[n];
        Function<Double, Double>[] dw = new Function[n];
        Function<Double, Double>[] derivatedw = new Function[n];
        Function<Double, Double>[] apcahedw = new Function[n];
        w[0] = x -> 1 + x;
        dw[0] = x -> 1.0;
        derivatedw[0] = x -> -2 * (0) * PolynomialsUtils.createJacobiPolynomial( 0, 0, 0).value(x);
        apcahedw[0] = x -> 1.0;
        for (int k = 1; k < n; k++) {
            // Лямбды в джаве...
            int finalK = k;
            w[k] = x ->
                    (1 - x*x) * PolynomialsUtils.createJacobiPolynomial(finalK - 1, 1, 1).value(x);
            // -2 ???
            dw[k] = x ->
                    -2 * (finalK - 1 + 1)*Math.pow(1-x*x, 1-1) * PolynomialsUtils.createJacobiPolynomial(finalK - 1 + 1, 0, 0).value(x);
            derivatedw[k] = derivate(w[k], 0.000001);
            apcahedw[k] = x -> PolynomialsUtils.createJacobiPolynomial(finalK - 1, 1, 1)
                    .multiply(new PolynomialFunction(new double[]{1.0, 0.0, -1.0})).polynomialDerivative().value(x);
        }

        for (int i = 0; i < n; i++) {
            for (double x = 0.0; x <= 1.0; x += 0.001){
                System.out.printf("%10.5f|%10.5f|%10.5f\n", dw[i].apply(x), derivatedw[i].apply(x), apcahedw[i].apply(x));
            }
        }
    }
    /**
     * Производная по формуле производных
     * */
    public static Function<Double, Double> derivate(Function<Double, Double> f, double h) {
//        Function<Double,Double> res =  x -> (f.apply(x + h) - f.apply(x - h)) / (2*h);
        Function<Double,Double> res =  x -> (3 * f.apply(x) - 4 * f.apply(x - h) + f.apply(x - 2*h)) / (2*h);
        return res;
    }
}
