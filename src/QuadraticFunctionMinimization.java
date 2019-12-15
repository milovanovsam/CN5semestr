import java.util.Scanner;

public class QuadraticFunctionMinimization {

    private static int N = 12;
    private static double eps = 0.0000001;

    public static void main(String[] args) {
        Scanner in = new Scanner(System.in);
        int number = 1;
        while (number != 0) {
            double[][] A = {{4, 1, 1}, {1, 2 * (3 + 0.1 * N), -1}, {1, -1, 2 * (4 + 0.1 * N)}};
            double[][] b = {{1}, {-2}, {3}};
            int dimension = 3;
            System.out.println("МИНИМИЗАЦИЯ КВАДРАТИЧНОЙ ФУНКЦИИ\n" + "Выберите метод для тестирования:\n"
                    + "1 - МНГС\n" + "2 - МНПС");
            number = in.nextInt();
            switch (number){
                case 1: MNGS(A, b, dimension); break;
                case 2: MNPS(A, b, dimension); break;
                default: break;
            }
            System.out.println("Для окончания нажмите 0, иначе любую клавишу");
            number = in.nextInt();
        }
    }

    private static void MNGS (double[][] A, double[][] b, int dimension){
        double different;
        double[][] nyu = {{0}};
        double[][] q = new double[dimension][1];
        double[][] x = new double[dimension][1];
        x[0][0] = 1;
        different = function(x) + 1;
        System.out.println(function(x));
        while (Math.abs(function(x) - different) >= eps) {
            different = function(x);
            multiplication(A, x, q);
            sum(q, b);
            double[][] result = new double[dimension][1];
            multiplication(A, q, result);
            double[][] qnew = new double[1][dimension];
            transpose(q, qnew);
            multiplication(qnew, q, nyu);
            double k = -nyu [0][0];
            multiplication(qnew, result, nyu);
            k /= nyu[0][0];
            scalar(q, k);
            sum(x, q);
            System.out.println(function(x));
        }
    }

    private static void MNPS (double[][] A, double[][] b, int dimension){
        double different;
        double[][] nyu = {{0}};
        double[][] q = new double[dimension][1];
        double[][] x = new double[dimension][1];
        double[][] e = new double[dimension][1];
        double[][] Ae = new double[dimension][1];
        int coordination = 0;
        different = function(x) + 2;
        System.out.println(function(x));
        while (Math.abs(function(x) - different) >= eps) {
            e[coordination % dimension][0] = 1;
            e[(coordination - 1 + dimension) % dimension][0] = 0;
            different = function(x);
            multiplication(A, x, q);
            sum(q, b);
            double[][] enew = new double[1][dimension];
            transpose(e, enew);
            multiplication(A,e,Ae);
            multiplication(enew,Ae,nyu);
            double k = (-1)/nyu[0][0];
            multiplication(enew, q, nyu);
            k *= nyu[0][0];
            scalar(e,k);
            sum(x,e);
            System.out.println(function(x));
            coordination++;
        }
    }

    private static double function(double[][] vector) {
        double x = vector[0][0];
        double y = vector[1][0];
        double z = vector[2][0];
        return 2*x*x + y*y*(3+0.1*N) + (4+0.1*N)*z*z + x*y - y*z + x*z + x - 2*y + 3*z + N;
    }

    static void multiplication(double[][] first, double[][] second, double[][] result){
        double sum;
        for ( int i = 0; i < first.length; i++)
            for ( int j = 0; j < second[0].length; j++) {
                sum = 0;
                for (int k = 0; k < second.length; k++)
                    sum += first[i][k] * second[k][j];
                result[i][j] = sum;
            }
    }

    private static void sum(double[][] first, double[][] second){
        for ( int i = 0; i < first.length; i++)
            for ( int j = 0; j < first[0].length; j++)
                first[i][j] += second[i][j];
    }

    static void transpose (double[][] vector, double[][] result){
        for (int i = 0; i < result.length; i ++)
            for (int j = 0; j < result[0].length; j++)
                result[i][j] = vector[j][i];
    }

    private static void scalar (double[][] vector, double scalar){
        for ( int i = 0; i < vector.length; i++)
            for ( int j = 0; j < vector[0].length; j++)
                vector[i][j]*= scalar;
    }
}
