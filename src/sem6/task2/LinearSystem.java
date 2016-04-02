package sem6.task2;

/**
 * Created by andrey on 31.03.16.
 */
public class LinearSystem {


    double[][] A;
    int n ;
    public LinearSystem(double[][] a) {
        int m = a.length;
        int k = a[0].length;
        if (m != k)
            throw new IllegalArgumentException("Wrong size");
        A = a;
        n = m;
    }

    public void minus(int fromLine, int secondLine){
        for (int i = 0; i < n + 1; i++) {
            A[fromLine][i] = A[fromLine][i] - A[secondLine][i];
        }
    }
    public void minusWithCoefficient(int fromLine, int secondString){
        for (int i = 0; i < n + 1; i++) {
            double coefficient = A[fromLine][secondString]/A[secondString][secondString];
            A[fromLine][i] = A[fromLine][i] - A[secondString][i]*coefficient;
        }
    }
}
