package sem6.task4;

import sem6.task2.Task2;

import java.util.Arrays;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Created by andrey on 4/24/16.
 */
public class Task4 {

    public static void main(String[] args) {

    }


    public static void stepennoyMethod(double[][] A){
        double[] Y0 = new double[A.length];
        double eps = 10e-7;
        //Set all values to 1
        Arrays.fill(Y0, 1.0);
        boolean flag = true;
        int k = 0;

        double[] Yk_1 = Y0;
        double[] Yk = Task2.multiplyMatrixAndVector(A, Yk_1);
        double[] Lk = IntStream.range(0, Yk.length).mapToDouble(i -> Yk[i] / Yk_1[i]).toArray();



        // Uncomment Comment and there will be error
        while(flag){
//            Yk_1 = Y0;
//            Yk = Task2.multiplyMatrixAndVector(A, Yk_1);
//            Lk = IntStream.range(0, Yk.length).mapToDouble(i -> Yk[i] / Yk_1[i]).toArray();
        }


    }
}
