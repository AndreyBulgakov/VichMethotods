/**
 * Created by alexander on 29.02.16.
 */


fun main(args : Array<String>) {

    fun p(x: Double) = 1.0
    fun q(x: Double) = Math.pow(Math.E, -x)
    fun r(x: Double) = 1.0 / (1 + x) - x
    fun f(x: Double) = x + Math.pow(x, 3.0)

    val alfa0 = 0.25
    val alfa1 = -1
    val betta0 = 1.4
    val betta1 = 1

    println("write n")
    val n = readLine()!!.toInt() + 1
    println("write a")
    val x0: Double = 0.0//readLine()!!.toDouble()
    println("write b")
    val xn: Double = 1.0//readLine()!!.toDouble()
    val h: Double = (xn - x0) / (n - 1)

    var A = 0
    var B = 0

    val a = DoubleArray(n, { i ->
        if (i == 0) 0.0
        else if (i == n - 1) -1.0*betta1
        else p(x0 + i * h) - 0.5 * q(x0 + i * h) * h
    })

    val b = DoubleArray(n, { i ->
        if (i == 0) h * alfa0 - alfa1
        else if (i == n- 1) h * betta0 + betta1
        else -2.0 * p(x0 + i * h) + r(x0 + i * h) * h * h
    })

    val d = DoubleArray(n, { i ->
        if (i == 0) h * A
        else if (i == n - 1) h * B
        else f(x0 + i * h) * h * h
    })

    val c = DoubleArray(n, { i ->
        if (i == 0) alfa1.toDouble()
        else if (i == n - 1) 0.0
        else p(x0 + i * h) + 0.5 * q(x0 + i * h) * h
    })

    println("------------------")
    println(a.size)
    println(b.size)
    println(c.size)
    println(d.size)

    fun printingMatrix(): Unit {
        println("Matrix A | d")
        for (i in 0..n - 1) {
            for (k in 0..n - 1) {
                if (k > i) {
                    if (k != i + 1)
                        print ("| %12.9f | ".format(0.0))
                    else {
                        print("| %12.9f | ".format(c[i]))
                        //print("   ")
                    }
                } else {
                    if (k == i) {
                        print("| %12.9f | ".format(b[k]))
                    } else {
                        if (k == i - 1) {
                            print("| %12.9f | ".format(a[i]))
                        } else
                            print ("| %12.9f | ".format(0.0))
                    }
                }
            }
            println("| %12.9f".format(d[i]))
        }
    }

    val m = DoubleArray(n) //can't use lambda
    val k = DoubleArray(n) //can't use lambda
    val x = DoubleArray(n)

    fun mkGenerator() {
        println("vectors m | k")
        m[0] = -c[0] / b[0]
        k[0] =  d[0] / b[0]
        print ("| %12.9f | ".format(m[0]))
        println ("| %12.9f  ".format(k[0]))
        for (i in 0..n - 2) {
            m[i + 1] = (-c[i]) / (a[i] * m[i] + b[i]);
            k[i + 1] = (d[i] - a[i] * k[i]) / (a[i] * m[i] + b[i]);
            //m[i] = -c[i - 1] / (a[i - 1] * m[i - 1] + b[i - 1])
            //k[i] = (d[i - 1] - a[i - 1] * k[i - 1]) / (a[i - 1] * m[i - 1] + b[i - 1])
            print ("| %12.9f | ".format(m[i]))
            println ("| %12.9f  ".format(k[i]))
        }
    }

    fun xGenerator() {
        mkGenerator()
        x[n - 1] = (d[n - 1] - a[n - 1] * k[n - 1]) / (a[n - 1] * m[n - 1] + b[n - 1]);
        for (i in n-2 downTo 0)
            x[i] = m[i + 1] * x[i + 1] + k[i + 1];
        //for (i in 0..n - 1)
        //  x[i] = (d[i] - a[i] * k[i]) / (a[i] * m[i] + b[i])
    }

    fun nevXGenerator(){
        xGenerator()
        val nevyazka = DoubleArray(n, { i ->
            if (i > 0 && i < n - 1)
                (a[i] * x[i - 1] + b[i] * x[i] + c[i] * x[i + 1]) - d[i]
            else
                if (i == 0)
                    b[i] * x[i] + c[i] * x[i + 1] - d[i]
                else
                    a[i] * x[i - 1] + b[i] * x[i] - d[i]
        })
        println("vectors x | nevyazka")
        for (i in 0..n - 1) {
            print ("| %12.9f | ".format(x[i]))
            println ("| %20.18f  ".format(nevyazka[i]))
        }
    }

    printingMatrix()
    //xGenerator()
    nevXGenerator()


    a[0] = 0.0;
    b[0] = 2 * h * alfa0 + alfa1 * ((a[1] / c[1]) - 3);
    c[0] = alfa1 * ((b[1] / c[1]) + 4);
    d[0] = 2 * h * A + alfa1 * (d[1] / c[1]);
    a[n - 1] = -betta1 * ((b[n - 2] / a[n - 2]) + 4);
    b[n - 1] = 2 * h * betta0 + betta1 * (3 - (c[n - 2] / a[n - 2]));
    c[n - 1] = 0.0;
    d[n - 1] = 2 * h * B - betta1 * (d[n - 2] / a[n - 2]);


    println("Second Variant")
    printingMatrix()
    //xGenerator()
    nevXGenerator()
}