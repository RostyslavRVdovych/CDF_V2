using System;
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


        //--------------------------------OLD ALGS


        //Послідовний алгоритм зі значеннями з попереднього кроку
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


        //Паралельний алгоритм (синхронна схема)
        public static double[] Par(double[] b, double[] f, int k, int m)
        {
            double[] a = (double[])b.Clone();
            double[] temp = new double[a.Length];
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
        public static double[] ParThreadPool(double[] b, double[] f, int k, int m)
        {
            double[] a = (double[])b.Clone();
            double[] temp = new double[a.Length];
            int workerThreads = Environment.ProcessorCount;
            ManualResetEvent[] doneEvents = new ManualResetEvent[workerThreads];

            for (int j = 0; j < k; j++)
            {
                temp = (double[])a.Clone();

                int itemsPerThread = (int)Math.Ceiling((double)a.Length / workerThreads);

                for (int t = 0; t < workerThreads; t++)
                {
                    doneEvents[t] = new ManualResetEvent(false);
                    int start = t * itemsPerThread;
                    int end = Math.Min(start + itemsPerThread, a.Length);

                    ThreadPool.QueueUserWorkItem((state) =>
                    {
                        int s = (int)state;
                        for (int i = start; i < end; i++)
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
                        doneEvents[s].Set();
                    }, t);
                }

                // Wait for all threads to complete
                WaitHandle.WaitAll(doneEvents);
            }

            return a;
        }
        public static double[] ParParallelFor(double[] b, double[] f, int k, int m)
        {
            double[] a = (double[])b.Clone();
            double[] temp = new double[a.Length];
            int workerThreads = Environment.ProcessorCount;

            for (int j = 0; j < k; j++)
            {
                temp = (double[])a.Clone();

                Parallel.For(0, workerThreads, t =>
                {
                    int itemsPerThread = (int)Math.Ceiling((double)a.Length / workerThreads);
                    int start = t * itemsPerThread;
                    int end = Math.Min(start + itemsPerThread, a.Length);

                    for (int i = start; i < end; i++)
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
                });
            }

            return a;
        }


        //Паралельний алгоритм з автономними гілками (не правильно працює)
        public static double[] OldParBranch(double[] b, double[] f, int k, int m)
        {
            double[] a = (double[])b.Clone();
            double[] c = (double[])b.Clone();
            Parallel.For(0, a.Length, t =>
            {
                double[] temp = (double[])a.Clone();
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
        //Паралельний алгоритм з обмеженим паралелізмом (не правильно працює)
        public static double[] OldParBranchLim(double[] b, double[] f, int k, int m)
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


        //--------------------------------NEW ALGS


        //Послідовний алгоритм зі значеннями з попереднього кроку з 2-вимірним масивом
        public static double[] Seq2(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }
            for (int j = 1; j <= k; j++)
            {
                for (int i = 0; i < n; i++)
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
            double[] c = new double[n];
            for (int i = 0; i < n; i++)
            {
                c[i] = x[k, i];
            }
            return c;
        }
        //Паралельний алгоритм (синхронна схема) з 2-вимірним масивом
        public static double[] Par2(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }
            for (int j = 1; j <= k; j++)
            {
                Parallel.For(0, n, i =>
                {
                    x[j, i] = 0.0;
                    for (int s = i - m; s <= i + m; s++)
                    {
                        if (s >= 0 && s < n)
                        {
                            x[j, i] += x[j - 1, s] * f[s - i + m];
                        }
                    }
                });
            }
            double[] c = new double[n];
            for (int i = 0; i < n; i++)
            {
                c[i] = x[k, i];
            }
            return c;
        }



        //Паралельний алгоритм з автономними гілками з масивами під кожну гілку
        public static double[] ParBranch(double[] b, double[] f, int k, int m)
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

                for (int j = 0; j < k; j++)
                {
                    double[] tempBranch = (double[])branch.Clone();
                    int startIndex = m * (j + 1);
                    int endIndex = (branch.Length - 1) - (m * (j + 1));
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
        public static double[] ParBranchThreadPool(double[] b, double[] f, int k, int m)
        {
            double[] c = (double[])b.Clone();

            CountdownEvent countdownEvent = new(b.Length);

            for (int i = 0; i < b.Length; i++)
            {
                ThreadPool.QueueUserWorkItem(state =>
                {
                    int t = (int)state;
                    double[] branch = new double[2 * m * k + 1];
                    for (int j = 0; j < branch.Length; j++)
                    {
                        int index = t - k * m + j;
                        if (index >= 0 && index < b.Length)
                        {
                            branch[j] = b[index];
                        }
                        else
                        {
                            branch[j] = 0;
                        }
                    }

                    for (int j = 0; j < k; j++)
                    {
                        double[] tempBranch = (double[])branch.Clone();
                        int startIndex = m * (j + 1);
                        int endIndex = (branch.Length - 1) - (m * (j + 1));
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

                    countdownEvent.Signal();
                }, i);
            }

            countdownEvent.Wait();

            return c;
        }




        //Послідовна реалізація паралельного алгоритму з автономними гілками з 2-вимірним масивом
        public static double[] SeqBranch2(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }


            //bool[,] isCalculated = new bool[k, n];
            //for (int j = 0; j < k; j++)
            //{
            //    for (int i = 0; i < n; i++)
            //    {
            //        isCalculated[j, i] = false;
            //    }
            //}

            for (int t = 0; t < n; t++)
            {
                for (int j = 1; j <= k; j++)
                {
                    for (int i = Math.Max(0, (j - k) * m + t); i <= Math.Min(n - 1, (k - j) * m + t); i++)
                    {

                        //if (!isCalculated[j - 1, i])
                        //{
                        //    isCalculated[j - 1, i] = true;

                        x[j, i] = 0.0;
                        for (int s = i - m; s <= i + m; s++)
                        {
                            if (s >= 0 && s < n)
                            {
                                x[j, i] += x[j - 1, s] * f[s - i + m];
                            }
                        }
                        //}
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
        //Паралельний алгоритм з автономними гілками з 2-вимірним масивом
        public static double[] ParBranch2(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }
            //ParallelOptions po = new()
            //{
            //    MaxDegreeOfParallelism = 2
            //};
            Parallel.For(0, n,/* po,*/ t =>
            {
                for (int j = 1; j <= k; j++)
                {
                    for (int i = Math.Max(0, (j - k) * m + t); i <= Math.Min(n - 1, (k - j) * m + t); i++)
                    {
                        double p = 0.0;
                        for (int s = i - m; s <= i + m; s++)
                        {
                            if (s >= 0 && s < n)
                            {
                                p += x[j - 1, s] * f[s - i + m];
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
        //Паралельний алгоритм з автономними гілками з 2-вимірним масивом з булевим масивом для уникнення дублювань
        public static double[] ParBranch2Bool(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }

            bool[,] isCalculated = new bool[k, n];
            for (int j = 0; j < k; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    isCalculated[j, i] = false;
                }
            }

            //ParallelOptions po = new()
            //{
            //    MaxDegreeOfParallelism = 2
            //};

            Parallel.For(0, n,/* po,*/ t =>
            {
                for (int j = 1; j <= k; j++)
                {
                    for (int i = Math.Max(0, (j - k) * m + t); i <= Math.Min(n - 1, (k - j) * m + t); i++)
                    {
                        // Умова для уникнення дублювань
                        if (!isCalculated[j - 1, i])
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
                            isCalculated[j - 1, i] = true;
                        }
                        else
                        {
                            continue;
                        }
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
        //Паралельний алгоритм з автономними гілками з 2-вимірним масивом, використовуючи Threads
        public static double[] ParBranch2Threads(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }

            Thread[] threads = new Thread[n];

            for (int t = 0; t < n; t++)
            {
                int threadIndex = t;

                threads[t] = new Thread(() =>
                {
                    for (int j = 1; j <= k; j++)
                    {
                        for (int i = Math.Max(0, (j - k) * m + threadIndex); i <= Math.Min(n - 1, (k - j) * m + threadIndex); i++)
                        {
                            double p = 0.0;
                            for (int s = i - m; s <= i + m; s++)
                            {
                                if (s >= 0 && s < n)
                                {
                                    p += x[j - 1, s] * f[s - i + m];
                                }
                            }
                            x[j, i] = p;
                        }
                    }
                });

                threads[t].Start();
            }

            foreach (Thread thread in threads)
            {
                thread.Join();
            }

            double[] c = new double[n];
            for (int i = 0; i < n; i++)
            {
                c[i] = x[k, i];
            }
            return c;
        }
        //Паралельний алгоритм з автономними гілками з 2-вимірним масивом, використовуючи Tasks
        public static double[] ParBranch2Tasks(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }

            Task[] tasks = new Task[n];

            for (int t = 0; t < n; t++)
            {
                int threadIndex = t;

                tasks[t] = Task.Run(() =>
                {
                    for (int j = 1; j <= k; j++)
                    {
                        for (int i = Math.Max(0, (j - k) * m + threadIndex); i <= Math.Min(n - 1, (k - j) * m + threadIndex); i++)
                        {
                            double p = 0.0;
                            for (int s = i - m; s <= i + m; s++)
                            {
                                if (s >= 0 && s < n)
                                {
                                    p += x[j - 1, s] * f[s - i + m];
                                }
                            }
                            x[j, i] = p;
                        }
                    }
                });
            }

            Task.WaitAll(tasks);

            double[] c = new double[n];
            for (int i = 0; i < n; i++)
            {
                c[i] = x[k, i];
            }
            return c;
        }
        //Паралельний алгоритм з автономними гілками з 2-вимірним масивом, використовуючи ThreadPool
        public static double[] ParBranch2ThreadPool(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }

            CountdownEvent countdownEvent = new(n);

            for (int t = 0; t < n; t++)
            {
                ThreadPool.QueueUserWorkItem(state =>
                {
                    int threadIndex = (int)state;

                    for (int j = 1; j <= k; j++)
                    {
                        for (int i = Math.Max(0, (j - k) * m + threadIndex); i <= Math.Min(n - 1, (k - j) * m + threadIndex); i++)
                        {
                            double p = 0.0;
                            for (int s = i - m; s <= i + m; s++)
                            {
                                if (s >= 0 && s < n)
                                {
                                    p += x[j - 1, s] * f[s - i + m];
                                }
                            }
                            x[j, i] = p;
                        }
                    }

                    countdownEvent.Signal();
                }, t);
            }

            countdownEvent.Wait();

            double[] c = new double[n];
            for (int i = 0; i < n; i++)
            {
                c[i] = x[k, i];
            }
            return c;
        }
        //Паралельний алгоритм з автономними гілками з 2-вимірним масивом з булевим масивом для уникнення дублювань, використовуючи Tasks
        public static double[] ParBranch2BoolTasks(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }

            bool[,] isCalculated = new bool[k, n];
            for (int j = 0; j < k; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    isCalculated[j, i] = false;
                }
            }

            List<Task> tasks = new List<Task>();

            for (int t = 0; t < n; t++)
            {
                int index = t;

                Task task = Task.Run(() =>
                {
                    for (int j = 1; j <= k; j++)
                    {
                        for (int i = Math.Max(0, (j - k) * m + index); i <= Math.Min(n - 1, (k - j) * m + index); i++)
                        {
                            // Умова для уникнення дублювань
                            if (!isCalculated[j - 1, i])
                            {
                                double p = 0.0;
                                for (int s = i - m; s <= i + m; s++)
                                {
                                    if (s >= 0 && s < n)
                                    {
                                        p += x[j - 1, s] * f[s - i + m];
                                    }
                                }
                                x[j, i] = p;
                                isCalculated[j - 1, i] = true;
                            }
                            else
                            {
                                continue;
                            }
                        }
                    }
                });

                tasks.Add(task);
            }

            Task.WaitAll(tasks.ToArray());

            double[] c = new double[n];
            for (int i = 0; i < n; i++)
            {
                c[i] = x[k, i];
            }
            return c;
        }
        //Паралельний алгоритм з автономними гілками з 2-вимірним масивом з булевим масивом для уникнення дублювань, використовуючи ThreadPool
        public static double[] ParBranch2BoolThreadPool(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }

            bool[,] isCalculated = new bool[k, n];
            for (int j = 0; j < k; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    isCalculated[j, i] = false;
                }
            }

            ManualResetEvent resetEvent = new ManualResetEvent(false);
            int tasksRemaining = n;

            for (int t = 0; t < n; t++)
            {
                int index = t;

                ThreadPool.QueueUserWorkItem((state) =>
                {
                    int taskIndex = (int)state;
                    try
                    {
                        for (int j = 1; j <= k; j++)
                        {
                            for (int i = Math.Max(0, (j - k) * m + taskIndex); i <= Math.Min(n - 1, (k - j) * m + taskIndex); i++)
                            {
                                // Умова для уникнення дублювань
                                if (!isCalculated[j - 1, i])
                                {
                                    double p = 0.0;
                                    for (int s = i - m; s <= i + m; s++)
                                    {
                                        if (s >= 0 && s < n)
                                        {
                                            p += x[j - 1, s] * f[s - i + m];
                                        }
                                    }
                                    x[j, i] = p;
                                    isCalculated[j - 1, i] = true;
                                }
                                else
                                {
                                    continue;
                                }
                            }
                        }
                    }
                    finally
                    {
                        if (Interlocked.Decrement(ref tasksRemaining) == 0)
                        {
                            resetEvent.Set();
                        }
                    }
                }, t);
            }

            resetEvent.WaitOne();

            double[] c = new double[n];
            for (int i = 0; i < n; i++)
            {
                c[i] = x[k, i];
            }
            return c;
        }



        //Паралельний алгоритм з обмеженим паралелізмом з масивами під кожну гілку
        public static double[] ParLim(double[] b, double[] f, int k, int m)
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
                    if (i >= 0 && i <= b.Length - 1)
                    {
                        branch[i - index1] = b[i];
                    }
                    else
                    {
                        branch[i - index1] = 0;
                    }
                }

                for (int j = 0; j < k; j++)
                {
                    double[] tempBranch = (double[])branch.Clone();
                    int startIndex = m * (j + 1);
                    int endIndex = (branch.Length - 1) - (m * (j + 1));
                    for (int i = startIndex; i <= endIndex; i++)
                    {
                        if (t == 0 && i < k * m)
                        {
                            continue;
                        }
                        else if (t == p - 1 && i > branch.Length - k * m - 1)
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
        //Паралельний алгоритм з обмеженим паралелізмом з масивами під кожну гілку, використовуючи ThreadPool
        public static double[] ParLimThreadPool(double[] b, double[] f, int k, int m)
        {
            double[] c = (double[])b.Clone();
            int p = Environment.ProcessorCount;
            // int p = 3;
            ManualResetEvent[] waitHandles = new ManualResetEvent[p];

            for (int t = 0; t < p; t++)
            {
                waitHandles[t] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(state =>
                {
                    int threadIndex = (int)state;
                    int index1 = (threadIndex * b.Length / p) - m * k;
                    int index2 = ((threadIndex + 1) * b.Length / p) - 1 + m * k;
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

                    for (int j = 0; j < k; j++)
                    {
                        double[] tempBranch = (double[])branch.Clone();
                        int startIndex = m * (j + 1);
                        int endIndex = (branch.Length - 1) - (m * (j + 1));
                        for (int i = startIndex; i <= endIndex; i++)
                        {
                            if (threadIndex == 0 && i < k * m)
                            {
                                continue;
                            }
                            else if (threadIndex == p - 1 && i > branch.Length - 1 - k * m)
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

                    waitHandles[threadIndex].Set(); // Сигналізуємо про завершення роботи потоку
                }, t);
            }

            WaitHandle.WaitAll(waitHandles); // Чекаємо завершення всіх потоків

            return c;
        }
        //Паралельний алгоритм з обмеженим паралелізмом з 2-вимірним масивом
        public static double[] ParLim2(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }
            int p = Environment.ProcessorCount;
            //int p = 12;
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
        //Паралельний алгоритм з обмеженим паралелізмом з 2-вимірним масивом, використовуючи Threads
        public static double[] ParLim2Threads(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }

            int p = Environment.ProcessorCount;
            //int p = 12;
            Thread[] threads = new Thread[p];
            int chunkSize = n / p;

            for (int t = 0; t < p; t++)
            {
                int start = t * chunkSize;
                int end = (t == p - 1) ? n - 1 : (t + 1) * chunkSize - 1;

                threads[t] = new Thread(() =>
                {
                    for (int j = 1; j <= k; j++)
                    {
                        for (int i = Math.Max(0, (j - k) * m + start); i <= Math.Min(n - 1, (k - j) * m + end); i++)
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
                threads[t].Start();
            }

            foreach (var thread in threads)
            {
                thread.Join();
            }

            double[] c = new double[n];
            for (int i = 0; i < n; i++)
            {
                c[i] = x[k, i];
            }
            return c;
        }
        //Паралельний алгоритм з обмеженим паралелізмом з 2-вимірним масивом, використовуючи ThreadPool
        public static double[] ParLim2ThreadPool(double[] b, double[] f, int k, int m)
        {
            int n = b.Length;
            double[,] x = new double[k + 1, n];
            for (int i = 0; i < n; i++)
            {
                x[0, i] = b[i];
            }

            int p = Environment.ProcessorCount;
            //int p = 12;
            int chunkSize = n / p;

            ManualResetEvent[] waitHandles = new ManualResetEvent[p];
            for (int t = 0; t < p; t++)
            {
                waitHandles[t] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(state =>
                {
                    int threadIndex = (int)state;
                    int start = threadIndex * chunkSize;
                    int end = (threadIndex == p - 1) ? n - 1 : (threadIndex + 1) * chunkSize - 1;

                    for (int j = 1; j <= k; j++)
                    {
                        for (int i = Math.Max(0, (j - k) * m + start); i <= Math.Min(n - 1, (k - j) * m + end); i++)
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
                    waitHandles[threadIndex].Set();
                }, t);
            }

            WaitHandle.WaitAll(waitHandles);

            double[] c = new double[n];
            for (int i = 0; i < n; i++)
            {
                c[i] = x[k, i];
            }
            return c;
        }


        //----------------my alg
        public static double[] MySeq(double[] b, int k, int m)
        {
            int n = b.Length;
            double[] a = (double[])b.Clone();
            for (int j = 0; j < k; j++)
            {
                double[] temp = (double[])a.Clone();
                for (int i = 0; i < n; i++)
                {
                    a[i] = 0;
                    double[] branch = new double[2 * m + 1];
                    for (int l = 0; l < 2 * m + 1; l++)
                    {
                        int index = i - m + l;
                        if (index < 0)
                        {
                            branch[l] = temp[0];
                        }
                        else if (index >= n)
                        {
                            branch[l] = temp[n - 1];
                        }
                        else
                        {
                            branch[l] = temp[index];
                        }
                    }
                    double mean = 0.0;
                    for (int l = 0; l < branch.Length; l++)
                    {
                        mean += branch[l];
                    }
                    mean /= branch.Length;
                    if (Math.Abs(mean - branch[m]) >= 0.02)
                    {
                        double[] coefs = GetCoefs(branch, mean);
                        for (int l = 0; l < 2 * m + 1; l++)
                        {
                            a[i] += branch[l] * coefs[l];
                        }
                    }
                    else
                    {
                        a[i] = branch[m];
                    }
                }
            }
            return a;
        }
        private static double[] GetCoefs(double[] branch, double mean)
        {
            double[] output = new double[branch.Length];
            double[] br = (double[])branch.Clone();
            double sum = 0.0;
            for (int i = 0; i < br.Length; i++)
            {
                if (mean - br[i] == 0)
                {
                    sum += br[i];
                }
                else
                {
                    br[i] = 1.0 / Math.Abs(mean - br[i]);
                    sum += br[i];
                }
            }
            for (int i = 0; i < br.Length; i++)
            {
                if (sum == 0.0)
                {
                    output[i] = br[i];
                }
                else
                {
                    output[i] = br[i] / sum;
                }
            }
            return output;
        }
    }
}
