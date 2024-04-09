using System.ComponentModel;
using MathNet.Numerics;

namespace Algorithms
{
    public class DigitalFiltering
    {
        //FUNCTIONS OF SIGNALS
        //NameOfFunction(sizeOfArray)
        public static double[] AddNoise(double[] f)
        {
            double[] noisyFunc = (double[])f.Clone();
            var random = new Random();
            for (int i = 0; i < f.Length; i++)
            {
                double noise = (random.NextDouble() / 5.0) - 0.1;
                noisyFunc[i] = f[i] + noise;
            }
            return noisyFunc;
        }
        public static double[] CreateSinus(double[] f)
        {
            var func = (double[])f.Clone();
            for (int i = 0; i < f.Length; i++)
            {
                func[i] = Math.Sin(f[i]);
            }
            return func;
        }
        public static double[] CreateDustySignal(double[] f)
        {
            var func = (double[])f.Clone();
            for (int i = 0; i < f.Length; i++)
            {
                func[i] = Math.Atan(Math.Tan(f[i]));
            }
            return func;
        }
        public static double[] CreateTriangularSignal(double[] f)
        {
            var func = (double[])f.Clone();
            for (int i = 0; i < f.Length; i++)
            {
                func[i] = Math.Tan(Math.Cos(f[i]));
            }
            return func;
        }
        public static double[] CreateRectangularSignal(double[] f)
        {
            double amplitude = 1.0;
            double step = 1.0;
            var func = (double[])f.Clone();
            for (int i = 0; i < f.Length; i++)
            {
                double temp = Math.Cos(f[i] / step);
                if (temp >= 0)
                {
                    func[i] = amplitude;
                }
                else
                {
                    func[i] = 0;
                }
            }
            return func;
        }


        //WINDOWS
        //NameOfWindow(windowRadius)

        public static double[] DanielWindow(int m)
        {
            double[] f = new double[2 * m + 1];
            for (int i = 0; i < 2 * m + 1; i++)
            {
                f[i] = 1.0 / (2.0 * m + 1.0);
            }
            return f;
        }
        public static double[] BartlettWindow(int m)     //m >= 2
        {
            double[] f = new double[2 * m + 1];
            for (int j = -m; j <= m; j++)
            {
                f[j + m] = 1.0 / m - (double)Math.Abs(j) / (m * m);
            }
            return f;
        }
        public static double[] HemingWindow(int m)
        {
            double[] f = new double[2 * m + 1];
            for (int j = -m; j <= m; j++)
            {
                f[j + m] = (0.54 + 0.46 * Math.Cos(Math.PI * Math.Abs(j) / m)) / (0.08 + 1.08 * m);
            }
            return f;
        }
        public static double[] KaiserWindow(int m, double alpha)
        {
            if (alpha == 0)
            {
                return DanielWindow(m);
            }
            int M = 2 * m + 1;
            double[] f = new double[M];
            double beta = Math.PI * alpha;
            double sum = 0.0;
            for (int i = 0; i < M; i++)
            {
                double temp = beta * Math.Sqrt(1 - Math.Pow((2.0 * i / (M - 1)) - 1.0, 2));
                f[i] = SpecialFunctions.BesselI0(temp) / SpecialFunctions.BesselI0(beta);
                sum += f[i];
            }
            for (int i = 0; i < M; i++)
            {
                f[i] /= sum;
            }
            return f;
        }

        //ALGORITHMS
        //NameOfAlgorithm(noisyFunction, arrayOfCoefs, numOfIter, windowRadius)

        public static double[] Seq(double[] b, double[] f, int k, int m)
        {
            double[] a = (double[])b.Clone();
            double[] temp;
            for (int j = 0; j < k; j++)
            {
                temp = (double[])a.Clone();
                for (int i = 0; i < a.Length; i++)
                {
                    a[i] = 0;
                    for (int l = i - m; l <= i + m; l++)
                    {
                        if (l >= 0 && l < a.Length)
                        {
                            a[i] += temp[l] * f[l - i + m];
                        }
                    }
                }
            }
            return a;
        }
        public static double[] Par(double[] b, double[] f, int k, int m)
        {
            double[] a = (double[])b.Clone();
            double[] temp;
            for (int j = 0; j < k; j++)
            {
                temp = (double[])a.Clone();
                Parallel.For(0, a.Length, i =>
                {
                    a[i] = 0;
                    for (int l = i - m; l <= i + m; l++)
                    {
                        if (l >= 0 && l < a.Length)
                        {
                            a[i] += temp[l] * f[l - i + m];
                        }
                    }
                });
            }
            return a;
        }
        public static double[] ParBranch(double[] b, double[] f, int k, int m)
        {
            double[] a = (double[])b.Clone();
            double[] c = (double[])b.Clone();
            Parallel.For(0, a.Length, t =>
            {
                double[] temp = a;
                for (int j = 0; j < k; j++)
                {
                    int startIndex = Math.Max(0, (j - k + 1) * m + t);
                    int endIndex = Math.Min(a.Length - 1, (k - 1 - j) * m + t);
                    for (int i = startIndex; i <= endIndex; i++)
                    {
                        double sum = 0;
                        for (int l = i - m; l <= i + m; l++)
                        {
                            if (l >= 0 && l < a.Length)
                            {
                                sum += temp[l] * f[l - i + m];
                            }
                        }
                        temp[i] = sum;

                        if (j == k - 1)
                        {
                            c[t] = temp[i];
                        }
                    }
                }
            });
            return c;
        }
        public static double[] ParBranchLim(double[] b, double[] f, int k, int m)
        {
            double[] c = (double[])b.Clone();
            //int p = 3;
            int p = Environment.ProcessorCount;
            Parallel.For(0, p, t =>
            {
            //    int index1 = Math.Max(0, (t * b.Length / p)) - m * k;
            //    int index2 = Math.Min(b.Length - 1, ((t + 1) * b.Length / p) - 1) + m * k;
            //    double[] branch = new double[index2 - index1];
            //    for (int i = index1; i < index2; i++)
            //    {
            //        if (i >= 0 && i < b.Length)
            //        {
            //            branch[i] = b[i];
            //        }
            //        else
            //        {
            //            branch[i] = 0;
            //        }
            //    }
            //    double[] tempBranch = branch;

                double[] temp = c;
                double[] temp2;
                for (int j = 0; j < k; j++)
                {
                    temp2 = temp;
                    int startIndex = Math.Max(0, ((j - k + 1) * m) + (t * b.Length / p));
                    int endIndex = Math.Min(b.Length - 1, ((k - 1 - j) * m) + ((t + 1) * b.Length / p) - 1);
                    for (int i = startIndex; i <= endIndex; i++)
                    {
                        double sum = 0;
                        for (int l = i - m; l <= i + m; l++)
                        {
                            if (l >= 0 && l < b.Length)
                            {
                                sum += temp2[l] * f[l - i + m];
                            }
                        }
                        temp[i] = sum;

                        if (j == k - 1)
                        {
                            c[i] = temp[i];
                        }
                    }
                }
            });
            return c;
        }

        //new algs
        public static double[] SeqBranch(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }
            for (int t = 0; t < n; t++)
            {
                for (int j = 1; j <= k; j++)
                {
                    for (int i = Math.Max(0, (j - k) * m + t); i <= Math.Min(n - 1, (k - j) * m + t); i++)
                    {
                        x[j, i] = 0.0;
                        for (int s = i - m; s <= i + m; s++)
                        {
                            if (s >= 0 && s < n)
                            {
                                x[j, i] += x[j - 1, s] * f[s - i + m];
                            }
                        }
                    }
                }
            }
            double[] c = new double[n];
            for (int i = 0; i < n; i++)
            {
                c[i] = x[k, i];
            }
            return c;
        }
        public static double[] ParSeqBranch(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }
            Parallel.For(0, n, t =>
            {
                for (int j = 1; j <= k; j++)
                {
                    for (int i = Math.Max(0, (j - k) * m + t); i <= Math.Min(n - 1, (k - j) * m + t); i++)
                    {
                        double p = 0.0;
                        //x[j, i] = 0.0;
                        for (int s = i - m; s <= i + m; s++)
                        {
                            if (s >= 0 && s < n)
                            {
                                p += x[j - 1, s] * f[s - i + m];
                                //x[j, i] += x[j - 1, s] * f[s - i + m];
                            }
                        }
                        x[j, i] = p;
                    }
                }
            });
            double[] c = new double[n];
            for (int i = 0; i < n; i++)
            {
                c[i] = x[k, i];
            }
            return c;
        }
        public static double[] ParBranch2(double[] b, double[] f, int k, int m)
        {
            double[] c = (double[])b.Clone();
            Parallel.For(0, b.Length, t =>
            {
                double[] branch = new double[2 * m * k + 1];
                for (int i = 0; i < branch.Length; i++)
                {
                    int index = t - k * m + i;
                    if (index >= 0 && index < b.Length)
                    {
                        branch[i] = b[index];
                    }
                    else
                    {
                        branch[i] = 0;
                    }
                }
                double[] tempBranch = branch;

                for (int j = 0; j < k; j++)
                {
                    tempBranch = branch;
                    int startIndex = (j - k + 1) * m + m * k;
                    int endIndex = (k - 1 - j) * m + m * k;
                    for (int i = startIndex; i <= endIndex; i++)
                    {
                        if (i + t >= m * k && i + t <= b.Length - 1 + m * k)
                        {
                            double sum = 0.0;
                            for (int s = i - m; s <= i + m; s++)
                            {
                                sum += tempBranch[s] * f[s - i + m];
                            }
                            branch[i] = sum;
                        }
                    }
                    if (j == k - 1)
                    {
                        c[t] = branch[m * k];
                    }
                }
            });
            return c;
        }
        public static double[] ParLim(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }
            int p = Environment.ProcessorCount;
            Parallel.For(0, p, t =>
            {
                for (int j = 1; j <= k; j++)
                {
                    for (int i = Math.Max(0, (j - k) * m + (t * n / p)); i <= Math.Min(n - 1, (k - j) * m + ((t + 1) * n / p) - 1); i++)
                    {
                        double p2 = 0.0;
                        for (int s = i - m; s <= i + m; s++)
                        {
                            if (s >= 0 && s < n)
                            {
                                p2 += x[j - 1, s] * f[s - i + m];
                            }
                        }
                        x[j, i] = p2;
                    }
                }
            });
            double[] c = new double[n];
            for (int i = 0; i < n; i++)
            {
                c[i] = x[k, i];
            }
            return c;
        }
        public static double[] ParBranchLim2(double[] b, double[] f, int k, int m)
        {
            double[] c = (double[])b.Clone();
            //int p = 3;
            int p = Environment.ProcessorCount;
            Parallel.For(0, p, t =>
            {
                int index1 = (t * b.Length / p) - m * k;
                int index2 = ((t + 1) * b.Length / p) - 1 + m * k;
                double[] branch = new double[index2 - index1 + 1];
                for (int i = index1; i <= index2; i++)
                {
                    if (i >= 0 && i < b.Length)
                    {
                        branch[i - index1] = b[i];
                    }
                    else
                    {
                        branch[i - index1] = 0;
                    }
                }
                double[] tempBranch = (double[])branch.Clone();

                for (int j = 0; j < k; j++)
                {
                    tempBranch = branch;
                    int startIndex = m * (j + 1);
                    int endIndex = (branch.Length - 1) - (m * (j + 1));
                    for (int i = startIndex; i <= endIndex; i++)
                    {
                        if (t == 0 && i < k * m)
                        {
                            continue;
                        }
                        else if (t == p - 1 && i > branch.Length - 1 - k * m)
                        {
                            continue;
                        }
                        else
                        {
                            double sum = 0;
                            for (int l = i - m; l <= i + m; l++)
                            {
                                if (l >= 0 && l < b.Length)
                                {
                                    sum += tempBranch[l] * f[l - i + m];
                                }
                            }
                            branch[i] = sum;

                            if (j == k - 1)
                            {
                                c[i + index1] = branch[i];
                            }
                        }
                    }
                }
            });
            return c;
        }

        //my alg
        public static double[] MySeq(double[] b, int k, int m)
        {
            double[] c = (double[])b.Clone();
            Parallel.For(0, b.Length, t =>
            {
                double[] branch = new double[2 * m * k + 1];
                for (int i = 0; i < branch.Length; i++)
                {
                    int index = t - k * m + i;
                    if (index >= 0 && index < b.Length)
                    {
                        branch[i] = b[index];
                    }
                    else
                    {
                        branch[i] = 0;
                    }
                }
                double[] tempBranch = branch;

                for (int j = 0; j < k; j++)
                {
                    tempBranch = branch;
                    int startIndex = m * (j + 1);
                    int endIndex = tempBranch.Length - 1 - m * (j + 1);
                    for (int i = startIndex; i <= endIndex; i++)
                    {
                        double[] f = new double[2 * m + 1];
                        double avs = 0.0;
                        for (int s = i - m; s <= i + m; s++)
                        {
                            avs += tempBranch[s];
                        }
                        avs /= f.Length;

                        double sumf = 0.0;
                        for (int s = i - m; s <= i + m; s++)
                        {
                            f[s - i + m] = 1.0 / Math.Abs(tempBranch[s] - avs);
                            sumf += f[s - i + m];
                        }

                        double sum = 0.0;
                        for (int s = i - m; s <= i + m; s++)
                        {
                            sum += tempBranch[s] * (f[s - i + m] / sumf);
                        }
                        branch[i] = sum;
                    }
                    if (j == k - 1)
                    {
                        c[t] = branch[m * k];
                    }
                }
            });
            return c;
        }
    }
}
