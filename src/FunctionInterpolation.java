import java.util.Arrays;
import java.util.function.Function;

import static java.lang.Math.*;

class FunctionInterpolation {

    public static void main(String[] args){
        Function<Double, Double> task =  x -> tan(x) - cos(x) + 0.1;
        Function<Double, Double> taskModule = x -> task.apply(x)*abs(x);
        double[] myPoints = myPoints(3, 8, 10);
        double[] myMorePoints = myPoints(5,8,10);
        double[] recomendPoints = recomendPoints(2,8,10);
        double[] recomendMorePoints = recomendPoints(4,8,10);
        System.out.println("По равноотсаящим точкам:");
        double[] myPointsAnswer = createLagranzh(myPoints, task);
        double[] myMorePointsAnswer = createLagranzh(myMorePoints, task);
        double[] moduleAnswer = createLagranzh(myPoints, taskModule);
        double[] moduleMorePointsAnswer = createLagranzh(myMorePoints, taskModule);
        System.out.println(PolynomialToString(myPointsAnswer));
        System.out.println(PolynomialToString(myMorePointsAnswer));
        System.out.println(PolynomialToString(moduleAnswer));
        System.out.println(PolynomialToString(moduleMorePointsAnswer));
        System.out.println("По формуле (3.2):");
        System.out.println(PolynomialToString(createLagranzh(recomendPoints, task)));
        System.out.println(PolynomialToString(createLagranzh(recomendMorePoints, task)));
        System.out.println(PolynomialToString(createLagranzh(recomendPoints, taskModule)));
        System.out.println(PolynomialToString(createLagranzh(recomendMorePoints, taskModule)));
    }

    private static double[] createPolynomial(double[] array, int point){
        double[] result = new double[array.length-1];
        double[] firstPart = Arrays.copyOf(array, point);
        double[] secondPart = Arrays.copyOfRange(array, point + 1, array.length);
        System.arraycopy(firstPart, 0, result, 0, firstPart.length);
        System.arraycopy(secondPart, 0, result, firstPart.length, secondPart.length);
        return createPolynomial(result);
    }

    private static double[] createPolynomial(double[] points){
        int n = points.length;
        double[] result = new double[n+1];
        result[n] = 1;
        for (double point : points)
            for (int j = 0; j < n + 1; j++) {
                result[j] *= -point;
                if ( j != n ) result[j] += result[j + 1];
            }
        return result;
    }

    private static double[] createLagranzh(double[] points, Function<Double, Double> function){
        int n = points.length; double divider = 1;
        double[] result = new double[n];
        double[] current;
         for ( int i = 0; i < n; i++){
             for ( int j = 0; j < n; j++)
                 if( j != i) divider *= (points[i] - points[j]);
             divider = function.apply(points[i]) / divider;
             current = createPolynomial(points, i);
             for (int j = 0; j < n; j++)
                 result[j] += current[j]*divider;
             divider = 1;
         }
        return result;
    }

    private static double[] recomendPoints(int n, double start, double end){
        double[] answer = new double[n+1];
        for ( int i = 0; i <= n; i++)
            answer[i] = 0.5*(start + end + (end-start)*cos(PI*(2*i+1)/(2*(n+1))));
        return answer;
    }

    private static double[] myPoints(int n, double start, double end){
        double[] answer = new double[n];
        for ( int i = 0; i < n; i++)
            answer[i] = start + (end - start)/(n-1)*i;
        return answer;
    }

    static String PolynomialToString(double[] coefficients){
        StringBuilder sb = new StringBuilder();
        for ( int i = 0; i < coefficients.length; i++) {
            if(coefficients[i] > 0) sb.append("+");
            if (coefficients[i] != 0.0)
                sb.append(coefficients[i]).append("x^").append(coefficients.length - 1 - i);
        }
        return sb.toString();
    }
}
