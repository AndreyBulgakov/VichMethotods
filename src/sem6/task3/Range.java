package sem6.task3;

/**
 * Created by andrey on 02.04.16.
 */
public class Range {
    double begin;
    double end;

    public Range(double begin, double end) {
        this.begin = begin;
        this.end = end;
    }

    @Override
    public String toString() {
        return begin + " < Лi < " +end;
    }
}
