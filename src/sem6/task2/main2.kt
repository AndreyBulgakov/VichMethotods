package sem6.task2

import org.jetbrains.annotations.NotNull
import java.util.*

fun main(args : Array<String>) {
    //    val matrix = arrayOf(doubleArrayOf(2.0, -1.0, 3.0, 2.0), doubleArrayOf(3.0,3.0,3.0,2.0),doubleArrayOf(3.0,-1.0,-1.0,2.0),doubleArrayOf(3.0,-1.0,3.0,-1.0))
    //    var b = doubleArrayOf(4.0,6.0,6.0,6.0)

        val matrix = arrayOf(doubleArrayOf(6.9880E-06, -7.1240E-03, 4.37640), doubleArrayOf(6.3880E-03, -0.71240, 1.60640), doubleArrayOf(0.91880, 0.88769, 1.05640) )
        var b = doubleArrayOf(3.99895, 1.77052, 2.12361)
//    val matrix = arrayOf(
//            doubleArrayOf(3.0, 2.0, -5.0),
//            doubleArrayOf(2.0, -1.0, 3.0),
//            doubleArrayOf(1.0, 2.0, -1.0))
//    var b = doubleArrayOf(-1.0, 13.0, 9.0)

    var eps = 0.01
    val E = Array<DoubleArray>(matrix.size, {i -> DoubleArray(matrix.size, {j -> if (i == j) 1.0 else 0.0}) })

    fun opNorm(@NotNull x: Array<Array<Double>>): Double {
        var result = 0.0;
        for (i in 0..x.size - 1) {
            var res = 0.0;
            for (j in 0..x[0].size - 1)
                res += Math.abs(x[i][j])
            if (res > result)
                result = res
        }

        return result
    }


    fun simpleGauss(@NotNull m: Array<DoubleArray>, @NotNull res: DoubleArray): DoubleArray {
        //smth like deep copy
        val matrix = Array<DoubleArray>(m.size, {i -> DoubleArray(m.size, {j -> m[i][j]}) })
        val b = DoubleArray(res.size, {i -> res[i]})

        for (k in 0..matrix.size - 1) {
            if (Math.abs(matrix[k][k]) <= eps)
                println("May be mistakes, because element a[${k}][${k}] < epsilon")
            var tmp = matrix[k][k]
            for (j in k..matrix.size - 1)
                matrix[k][j] /= tmp
            b[k] /= tmp
            for(i in k + 1 .. matrix.size - 1){
                tmp = matrix[i][k]
                for (j in k + 1 .. matrix.size - 1)
                    matrix[i][j] -= matrix[k][j]*tmp

                b[i] -= b[k]*tmp
            }
        }

        var x = DoubleArray(matrix[0].size, { i -> b[i] })
        for (i in matrix[0].size - 1 downTo 0) {
            if (i != matrix[0].size - 1) {
                for (j in i + 1..matrix[i].size - 1)
                    x[i] += -x[j] * matrix[i][j]
            }
        }
        return x
    }

    fun gaussWithCol(@NotNull m:Array<DoubleArray>, @NotNull res:DoubleArray): DoubleArray {
        //smth like deep copy
        val matrix = Array<DoubleArray>(m.size, {i -> DoubleArray(m.size, {j -> m[i][j]}) })
        val b = DoubleArray(res.size, {i -> res[i]})
        val pairs =  ArrayList<Pair<Int, Int>>()

        fun swapColumns(matrix:Array<DoubleArray>, i:Int, j:Int){
            for (k in 0 .. matrix[0].size - 1) {
                var a = matrix[k][i]
                matrix[k][i] = matrix[k][j]
                matrix[k][j] = a
            }
        }

        fun reverseResults(x:DoubleArray, pairs: ArrayList<Pair<Int, Int>>): DoubleArray {
            for (i in pairs.size - 1 downTo 0){
                var pair = pairs[i]
                var a = x[pair.first]
                x[pair.first] = x[pair.second]
                x[pair.second] = a
            }
            return x
        }

        fun findNext(matrix:Array<DoubleArray>, i:Int, pairs:ArrayList<Pair<Int, Int>>){
            for (k in i+1 .. matrix[0].size - 1)
                if (Math.abs(matrix[k][i]) > eps) {
                    swapColumns(matrix, i, k)
                    pairs.add(Pair(i, k))
                    break
                }
        }

        for (k in 0..matrix.size - 1) {
            if (Math.abs(matrix[k][k]) <= eps)
                findNext(matrix, k, pairs)
            var tmp = matrix[k][k]
            for (j in k..matrix.size - 1)
                matrix[k][j] /= tmp
            b[k] /= tmp
            for(i in k + 1 .. matrix.size - 1){
                tmp = matrix[i][k]
                for (j in k + 1 .. matrix.size - 1)
                    matrix[i][j] -= matrix[k][j]*tmp
                b[i] -= b[k]*tmp
            }
        }

        var x = DoubleArray(matrix[0].size, { i -> b[i] })
        for (i in matrix[0].size - 2 downTo 0) {
            for (j in i + 1..matrix[i].size - 1)
                x[i] += -x[j] * matrix[i][j]
        }
        return reverseResults(x, pairs)
    }

    fun determinant(@NotNull m:Array<DoubleArray>): Double {
        //smth like deep copy
        val matrix = Array<DoubleArray>(m.size, {i -> DoubleArray(m.size, {j -> m[i][j]}) })
        val elements =  ArrayList<Double>()
        for (k in 0..matrix.size - 1) {
            var tmp = matrix[k][k]
            elements.add(tmp)
            for (j in k..matrix.size - 1)
                matrix[k][j] /= tmp
            for(i in k + 1 .. matrix.size - 1){
                tmp = matrix[i][k]
                for (j in k + 1 .. matrix.size - 1)
                    matrix[i][j] -= matrix[k][j]*tmp
            }
        }

        var res = 1.0
        for (i in elements)
            res *= i
        return res
    }

    fun multiMatrixVector(@NotNull matrix:Array<DoubleArray>, @NotNull b:DoubleArray):DoubleArray {
        val result = DoubleArray(matrix.size, {i ->
            var count = 0.0
            for(j in 0..matrix.size - 1)
                count += matrix[i][j]*b[j]
            count})
        return result
    }

    fun reverseMatrix(@NotNull matrix:Array<DoubleArray>, @NotNull E:Array<DoubleArray>): Array<DoubleArray> {
        //smth like deep copy
        val matrixT = Array<DoubleArray>(matrix.size, {i -> DoubleArray(matrix.size, {j -> matrix[j][i]}) })
        val pairs = ArrayList<Pair<Int, Int>>()
        val resultT = Array<DoubleArray>(matrix.size, { i -> gaussWithCol(matrix, E[i])})
        val result = Array<DoubleArray>(resultT.size, {i -> DoubleArray(resultT.size, {j -> resultT[j][i]}) })
        return result
    }

    fun vecNorm(@NotNull x: DoubleArray): Double {
        var result = 0.0
        for (i in 0..x.size - 1)
            result += Math.abs(x[i])
        return result
    }
    //pretty-printer
    println("matrix")
    for (i in 0..matrix.size - 1) {
        for (j in 0..matrix.size - 2)
            print("${matrix[i][j]}   ")
        println(matrix[i][matrix.size - 1])
    }
    println()
    println("result via simple Gauss's method")
    var a = simpleGauss(matrix, b)
    for (i in 0 .. a.size  - 1) {
        println(a[i])
    }
    //это x
    a = multiMatrixVector(matrix, a)
    for (i in 0.. a.size - 1)
        a[i] -= b[i]
    println("nevyazka is")
    println(vecNorm(a))
    println()
    println("result via Gauss's method with selection")
    var c = gaussWithCol(matrix, b)
    for (i in 0 .. c.size  - 1) {
        println(c[i])
    }
    c = multiMatrixVector(matrix, c)
    for (i in 0.. a.size - 1)
        c[i] -= b[i]
    println("nevyazka is")
    println(vecNorm(c))
    println()
    println("this is determinant matrix")
    println(determinant(matrix))
    //println("this is joke, I hadn't it, ho-ho-ho!")
    println()
    println("reverse matrix")
    var matrix_1 = reverseMatrix(matrix, E)
    for (i in 0..matrix_1.size - 1) {
        for (j in 0..matrix_1.size - 2)
        //need to cerrect format, but len'
            print("${matrix_1[i][j]}   ")
        println(matrix_1[i][matrix_1.size - 1])
    }
    println()
    println("result of x via A-1 matrix")
    var xnew = multiMatrixVector(matrix_1, b)
    for (i in 0..xnew.size - 1)
        print("${xnew[i]}   ")
    println()
    xnew = multiMatrixVector(matrix, xnew)
    for (i in 0.. a.size - 1)
        xnew[i] -= b[i]
    println("nevyazka is")
    println(vecNorm(xnew))

}