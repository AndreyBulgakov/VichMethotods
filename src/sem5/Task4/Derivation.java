package sem5.Task4;

/**
 * Created by andrey on 25.11.15.
 */
public class Derivation {
    private double[] y;
    private double h;
    private int n;

    public Derivation(double[] y, double h, int n)
    {
        this.y = y;
        this.h = h;
        this.n = n;
    }

    public double derivate1(int i)
    {
        if (i != n && i != 0) return  (y[i+1]-y[i-1])/(2*h);
        if (i == n) return (y[i]-y[i-1])/h;
        else return (y[i+1]-y[i])/h;

    }

    public double derivate2(int i)
    {
        if (i == 0 || i == n) return Double.NaN;
        return (y[i+1] - 2 * y[i] + y[i-1]) / (h * h);
    }
}
