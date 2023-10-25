using System;
using MathNet.Numerics.LinearAlgebra;

namespace math_lab2
{
    internal class Program
    {
        public static void Print2DArray<T>(T[,] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    Console.Write(matrix[i, j] + "\t");
                }
                Console.WriteLine();
            }
        }

        public static void PrintArray<T>(T[] array)
        {
            for (int i = 0; i < array.Length; i++)
            {
                Console.Write(array[i] + "\t");
                Console.WriteLine();
            }
        }
        static void FillMatrix(Matrix<double> matrix, int N)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    matrix[i, j] = 1.0/((i+1) + (j+1) - 1);
                }
            }
        }

        static double ComputeConditionNumber(Matrix<double> matrix)
        {
            Matrix<double> inverseMatrix = matrix.Inverse();

            double norm = matrix.L2Norm();
            double normInverse = inverseMatrix.L2Norm();

            // Расчет числа обусловленности матрицы
            double conditionNumber = norm * normInverse;

            return conditionNumber;
        }

        static void FillF(double[] F, int N)
        {
            for (int i = 0; i < N; i++)
            {
                F[i] = ((double)N)/Math.Pow(i+1,2) ;
            }
        }

        static double[] SolveSystem(Matrix<double> matrix, double[] f, int n)
        {
            // Прямой ход метода Гаусса
            for (int k = 0; k < n - 1; k++)
            {
                for (int i = k + 1; i < n; i++)
                {
                    double factor = matrix[i, k] / matrix[k, k];
                    f[i] = f[i] - factor * f[k];
                    for (int j = k; j < n; j++)
                    {
                        matrix[i, j] = matrix[i, j] - factor * matrix[k, j];
                    }
                }

            }

            // Обратный ход метода Гаусса
            double[] x = new double[n];
            x[n - 1] = f[n - 1] / matrix[n - 1, n - 1];
            for (int i = n - 2; i >= 0; i--)
            {
                double sum = f[i];
                for (int j = i + 1; j < n; j++)
                {
                    sum = sum - matrix[i, j] * x[j];
                }
                x[i] = sum / matrix[i, i];
            }

            return x;
        }

        static double[] SolveSystemWithMainEl(Matrix<double> matrix, double[] f, int n)
        {
            double[] solution = new double[n];

            for (int i = 0; i < n; i++)
            {
                // Выбор главного элемента
                int maxRow = i;
                double maxValue = Math.Abs(matrix[i, i]);
                for (int k = i + 1; k < n; k++)
                {
                    if (Math.Abs(matrix[k, i]) > maxValue)
                    {
                        maxRow = k;
                        maxValue = Math.Abs(matrix[k, i]);

                    }
                }

                // Обмен строк
                for (int k = i; k < n; k++)
                {
                    double temp = matrix[maxRow, k];
                    matrix[maxRow, k] = matrix[i, k];
                    matrix[i, k] = temp;
                }

                double tempF = f[maxRow];
                f[maxRow] = f[i];
                f[i] = tempF;

                // Приведение матрицы к треугольному виду
                for (int j = i + 1; j < n; j++)
                {
                    double div = matrix[j, i] / matrix[i, i];
                    for (int k = i; k < n; k++)
                    {
                        matrix[j, k] -= div * matrix[i, k];
                    }
                    f[j] -= div * f[i];
                }
            }

            // Обратный ход
            for (int i = n - 1; i >= 0; i--)
            {
                solution[i] = f[i];
                for (int j = i + 1; j < n; j++)
                {
                    solution[i] -= matrix[i, j] * solution[j];
                }
                solution[i] /= matrix[i, i];
            }

            return solution;
        }
        static void Main(string[] args)
        {
            int N = 10;

            Matrix<double> matrix = Matrix<double>.Build.DenseOfArray(new double[N, N]);
            double[] f = new double[N];
            FillMatrix(matrix, N);

            FillF(f, N);
                
            Console.WriteLine("Matrix");
            Print2DArray(matrix.ToArray());

            Console.WriteLine("\n\nF");
            PrintArray(f);

            Console.WriteLine("\n\nЧисло обусловленности матрицы");
            double cond = ComputeConditionNumber(matrix);
            Console.WriteLine(cond);

            var result = SolveSystem(matrix,f,N);

            Console.WriteLine("\n\nРезульт решения методом Гаусса без выбора главного элемента:");
            PrintArray(result);

            var resultWithMainEl = SolveSystemWithMainEl(matrix, f, N);

            Console.WriteLine("\n\nРезульт решения методом Гаусса с выбором главного элемента:");
            PrintArray(resultWithMainEl);


            Console.WriteLine("\n\nЧужой алгоритм:");
            PrintArray(matrix.Solve(Vector<double>.Build.DenseOfArray(f)).ToArray());
        }
    }
}
