package sem6.task1;

import java.io.IOException;
import java.util.function.Function;
import java.util.stream.IntStream;

/**
 * Created by andrey on 29.03.16.
 */
public class Task1 {

    public static void main(String[] args) throws IOException {

//        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
//        System.out.print("n = ");
//        int n = Integer.valueOf(reader.readLine());
        int n = 100;
//
//        Function<Double,Double> p = x -> 1.0;
//        Function<Double,Double> q = x -> 1/(x+2);
//        Function<Double,Double> r = x -> -Math.sin((Math.PI*x)/2);
//        Function<Double,Double> f = x -> x;
//        double alfa0 = 0.8;
//        double alfa1 = -1;
//        double betta0 = 0.6;
//        double betta1 = 1;
//        double x0 = 0;
//        double xn = 1;
//        double A = 0.0;
//        double B = 0.0;

        //Новый
        double alfa = 0.1;
        Function<Double,Double> p = x -> 1.0;
        Function<Double,Double> q = x -> -(x*x+alfa);
        Function<Double,Double> r = x -> -2.0*x;
        Function<Double,Double> f = x -> 2*(3*x*x-alfa)/Math.pow((x*x+alfa),3);
        double alfa0 = 1.0;
        double alfa1 = -2.0;
        double betta0 = 1.0;
        double betta1 = 0.0;
        double x0 = 0.0;
        double xn = 1.0;
        double A = 1/alfa;
        double B = 1/(1+alfa);
        double h = (xn-x0)/n;


//        firstVariant(p,q,r,f,n,x0,xn,alfa0,alfa1,betta0,betta1,A,B);

        System.out.println("\n Второй способ \n");

        secondVariant(p,q,r,f,n,x0,xn,alfa0,alfa1,betta0,betta1,A,B);
        //Второй вариант
    }


    static void firstVariant(Function<Double,Double> p, Function<Double,Double> q, Function<Double,Double> r, Function<Double,Double> f,
            int n, double x0, double xn, double alfa0, double alfa1, double betta0, double betta1, double A, double B){

        //Создаем диагонали матрицы по формулам
        //Для ti i=0..n
        // Ai = p(ti)-0.5*q(ti)*h
        // Bi = -2.0 * p(ti) + r(ti) * h * h
        // Ci = p(ti) + 0.5 * q(ti) * h
        // Di = f(ti) * h * h

        double h = ( xn - x0 ) / n;

        double[] a = IntStream.range(0, n+1).mapToDouble(i -> i).map(i -> {
            if (i == 0) return 0.0;
            else if (i == n) return -1.0 * betta1;
            else return p.apply(x0 + i * h) - 0.5 * q.apply(x0 + i * h) * h;
        })
        .toArray();

        double[] b = IntStream.range(0, n+1).mapToDouble(i -> i).map(i -> {
            if (i == 0) return h * alfa0 - alfa1;
            else if (i == n) return h * betta0 + betta1;
            else return -2.0 * p.apply(x0 + i * h) + r.apply(x0 + i * h) * h * h;
        })
        .toArray();

        double[] c = IntStream.range(0, n+1).mapToDouble(i -> i).map(i -> {
            if (i == 0) return alfa1;
            else if (i == n) return 0.0;
            else return p.apply(x0 + i * h) + 0.5 * q.apply(x0 + i * h) * h;
        })
        .toArray();

        double[] d = IntStream.range(0, n+1).mapToDouble(i -> i).map(i -> {
            if (i == 0) return h * A;
            else if (i == n) return h * B;
            else return f.apply(x0 + i * h) * h * h;
        })
        .toArray();

        //Массивы s и t
        double[] m = new double[n+1];
        double[] k = new double[n+1];
        m[0] = -c[0]/b[0];
        k[0] = d[0]/b[0];

        //Создаем mi и ki
        for (int i = 0; i < n ; i++) {
            m[i + 1] = (-c[i]) / (a[i] * m[i] + b[i]);
            k[i + 1] = (d[i] - a[i] * k[i]) / (a[i] * m[i] + b[i]);
        }
        //Находим y. Обратной прогонкой
        double[] y = new double[n+1];
        y[n]=( d[n] - a[n] * k[n] ) / (a[n] * m[n] +b[n]);
        for (int i = n-1; i >= 0 ; i--) {
            y[i]=m[i+1]*y[i+1]+k[i+1];
        }

        //Выводим матрицу A|d
        System.out.println("Матрица A|d");
        double[][] matrix = getMatrixA(a,b,c,d,n);
        printMatrix(matrix);

        //Выводим прогоночные коэффициенты
        System.out.println("Прогоночные коэффициенты");
        System.out.println("|    mi    ||    ki    |");
        for (int i = 0; i < n + 1; i++) {
            System.out.printf("|%12.9f||%12.9f|\n",m[i],k[i]);
        }

        Function<Double,Double> resh = x -> 1/(x*x+0.1);

        //Вывод ветора неизвестных
        System.out.println("Вектор неизвестных");
        System.out.println("|    xi    ||  точное  |");
        for (int i = 0; i < n + 1; i++) {
            System.out.printf("|%12.9f||%12.9f|\n",y[i],resh.apply(x0+i*h));
        }
        System.out.println();

    }

    public static void secondVariant(Function<Double,Double> p, Function<Double,Double> q, Function<Double,Double> r, Function<Double,Double> f,
                              int n, double x0, double xn, double alfa0, double alfa1, double betta0, double betta1, double A, double B){
        //Создаем диагонали матрицы по формулам
        //Для ti i=0..n
        // Ai = p(ti)-0.5*q(ti)*h
        // Bi = -2.0 * p(ti) + r(ti) * h * h
        // Ci = p(ti) + 0.5 * q(ti) * h
        // Di = f(ti) * h * h

        double h = ( xn - x0 ) / n;


        double[] a = IntStream.range(0, n+1).mapToDouble(i -> i).map(i -> {
            return p.apply(x0 + i * h) - 0.5 * q.apply(x0 + i * h) * h;
        })
                .toArray();
        double[] b = IntStream.range(0, n+1).mapToDouble(i -> i).map(i -> {
            return -2.0 * p.apply(x0 + i * h) + r.apply(x0 + i * h) * h * h;
        })
                .toArray();
        double[] c = IntStream.range(0, n+1).mapToDouble(i -> i).map(i -> {
            return p.apply(x0 + i * h) + 0.5 * q.apply(x0 + i * h) * h;
        })
                .toArray();
        double[] d = IntStream.range(0, n+1).mapToDouble(i -> i).map(i -> {
            return f.apply(x0 + i * h) * h * h;
        })
                .toArray();

        //Расставляем коэффициенты
        a[0] = 0.0;
        b[0] = 2 * h * alfa0 + alfa1 * ((a[1] / c[1]) - 3);
        c[0] = alfa1 * ((b[1] / c[1]) + 4);
        d[0] = 2 * h * A + alfa1 * (d[1] / c[1]);

        a[n] = -betta1 * ((b[n - 1] / a[n - 1]) + 4);
        b[n] = 2 * h * betta0 + betta1 * (3 - (c[n - 1] / a[n - 1]));
        c[n] = 0.0;
        d[n] = 2 * h * B - betta1 * (d[n - 1] / a[n - 1]);

        //Массивы s и t
        double[] m = new double[n+1];
        double[] k = new double[n+1];
        m[0] = -c[0]/b[0];
        k[0] = d[0]/b[0];

        //Создаем mi и ki
        for (int i = 0; i < n ; i++) {
            m[i + 1] = (-c[i]) / (a[i] * m[i] + b[i]);
            k[i + 1] = (d[i] - a[i] * k[i]) / (a[i] * m[i] + b[i]);
        }
        //Находим y. Обратной прогонкой
        double[] y = new double[n+1];
        y[n]=( d[n] - a[n] * k[n] ) / (a[n] * m[n] +b[n]);
        for (int i = n-1; i >= 0 ; i--) {
            y[i]=m[i+1]*y[i+1]+k[i+1];
        }

//        //Выводим матрицу A|d
//        System.out.println("Матрица A|d");
//        double[][] matrix = getMatrixA(a,b,c,d,n);
//        printMatrix(matrix);

        //Выводим прогоночные коэффициенты
        System.out.println("Прогоночные коэффициенты");
        System.out.println("|    mi    ||    ki    |");
        for (int i = 0; i < n + 1; i++) {
            System.out.printf("|%12.9f||%12.9f|\n",m[i],k[i]);
        }

        Function<Double,Double> resh = x -> 1/(x*x+0.1);

        //Вывод ветора неизвестных
        System.out.println("Вектор неизвестных");
        System.out.println("|     yi     ||   точное   ||  невязка  |");
        double nev[] = nevyazka(a,b,c,d,y,n);
        for (int i = 0; i < n + 1; i++) {
            System.out.printf("|%12.9f||%12.9f||%19.18f|\n",y[i],resh.apply(x0+i*h),nev[i]);
        }
        System.out.println();


    }
    public static double[][] getMatrixA(double[] a, double[] b, double[] c, double[] d,int n){
        double[][] resultMatrix = new double[n+1][n+2];
        //Заполним диагонали
        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < n + 1; j++) {
                if (j == i + 1)
                    resultMatrix[i][j] = c[j];
                else
                if (j == i)
                    resultMatrix[i][j] = b[j];
                else
                if (j == i - 1)
                    resultMatrix[i][j] = a[j];
                else
                    resultMatrix[i][j] = 0.0;
            }
            resultMatrix[i][n+1] = d[i];
        }
        return resultMatrix;
    }

    public static void printMatrix(double[][] matrix){
        for (int i = 0; i < matrix.length; i++){
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.printf("|%12.9f|",matrix[i][j]);
            }
            System.out.println();
        }
    }

    public static double[] nevyazka(double[] a, double[] b, double[] c, double[] d, double[] y, int n){
        double[] result =  IntStream.range(0, n+1).mapToDouble(i -> i).map(k -> {
                int i = (int)k;
                if (i > 0 && i < n)
                    return (a[i] * y[i - 1] + b[i] * y[i] + c[i] * y[i + 1]) - d[i];
                else
                if (i == 0)
                    return b[i] * y[i] + c[i] * y[i + 1] - d[i];
                else
                    return a[i] * y[i - 1] + b[i] * y[i] - d[i];
        }).toArray();
    return result;
    }
}
