import java.util.function.Function;
import static java.lang.Math.*;

public class FunctionApproximation {

    private final static double startPoint = -1;
    private final static double endPoint = 1;
    private static Function<Double, Double> Task = x -> x*x*x*sin(x);

    public static void main(String[] args) {
        double[] points = {-1, -0.5, 0, 0.5, 1};
        int degree = 4;
        System.out.println("LeastSquareMethod:");
        System.out.println(FunctionInterpolation.PolynomialToString(LeastSquareMethod(points, degree, Task)));
        System.out.println();
        System.out.println("LezhandrPolynom:");
        System.out.println(FunctionInterpolation.PolynomialToString(MethodLezhandrPolynom(degree-1)));
    }

    private static double[] LeastSquareMethod(double[] points, int degree, Function<Double, Double> function){
        int n = points.length;
        double[] system; double[] answer = new double[degree];
        double[][] y = new double[n][1];
        for ( int i = 0; i < n; i++)
            y[i][0] = function.apply(points[i]);
        double[][] Q = createQ(n, degree, points);
        double[][] transposeQ = new double[degree][n];
        QuadraticFunctionMinimization.transpose(Q,transposeQ);
        double[][] H = new double[degree][degree];
        QuadraticFunctionMinimization.multiplication(transposeQ, Q, H);
        double[][] b = new double[degree][1];
        QuadraticFunctionMinimization.multiplication(transposeQ, y, b);
        system = Gauss(H, b);
        for (int i = 0; i < degree; i++)
            answer[i] = system[degree -i-1];
        double sum = 0;
        for (double point : points) {
            for (int j = 0; j < degree; j++) {
                sum += answer[degree - 1 - j] * phiFunction(j, point);
            }
            System.out.println(abs(function.apply(point) - sum));
            sum = 0;
        }
        return answer;
    }

    private static double phiFunction(int k, double x){
        double result = 1;
        for ( int i = 0; i < k; i++)
            result *= x;
        return result;
    }

    private static double[][] createQ(int n, int m, double[] points){
        double[][] Q = new double[n][m];
        for (int i = 0; i < n; i++)
            for( int j = 0; j < m; j++)
                Q[i][j]  = phiFunction(j, points[i]);
        return Q;
    }

    private static double[] Gauss(double[][] A, double[][] B){
        double[] b = new double[B.length];
        for ( int i = 0; i < B.length; i++)
            b[i] = B[i][0];
        int n = b.length;
        double[] result = new double[n];
        for (int k = 0; k < n; k++) {
            for (int i = k + 1; i < n; i++)
                A[k][i] = A[k][i] / A[k][k];
            b[k] = b[k] / A[k][k];
            A[k][k] = 1;
            for (int i = k + 1; i < n; i++) {
                b[i] = b[i] - b[k] * A[i][k];
                for (int j = k; j < n; j++)
                    A[i][j] = A[i][j] - A[k][j] * A[i][k];
            }
        }
        double sum = 0;
        //Обратный ход
        for (int i = n - 1; i > -1; i--) {
            for (int j = n - 1; j > i; j--)
                sum += A[i][j] * result[j];
            result[i] = b[i] - sum;
            sum = 0;
        }
        return result;
    }

    private static Function<Double, Double> Lezhandr(int k){
        Function<Double, Double> lezhandr = x -> x;
        switch (k){
            case 0: lezhandr = x -> 1.0; break;
            case 1: lezhandr = x -> x; break;
            case 2: lezhandr = x -> 1.5*x*x - 0.5*x; break;
            default: break;
        }
        return lezhandr;
    }

    private static double[] LezhandrPoly(int k){
        double[] answer = new double[2];
        switch (k){
            case 0: answer = new double[]{0, 0, 1}; break;
            case 1: answer = new double[]{0, 1, 0}; break;
            case 2: answer = new double[]{1.5, 0, -0.5}; break;
            default: break;
        }
        return answer;
    }

    private static double integrate(Function<Double, Double> function1, Function<Double, Double> function2){
        double answer = 0;
        for (double i = startPoint; i < endPoint; i+= 0.001)
            answer += function1.apply(i)*function2.apply(i);
        return answer;
    }

    private static double[] MethodLezhandrPolynom(int degree){
        double[] answer = new double[degree];
        double[] coefficients = new double[degree];
        for (int i = 0; i < degree; i++)
            coefficients[i] = integrate(Task, Lezhandr(i))/integrate(Lezhandr(i),Lezhandr(i));
        for (int i = 0; i < degree; i++)
            for (int j = 0; j < degree; j++)
                answer[j] += coefficients[i]*LezhandrPoly(i)[j];
        return answer;
    }
}
