using System.Diagnostics;
using System.Numerics;
using Algorithms;

class Program
{
    static void Main()
    {
        int n = 1200000; // кількість елементів
        int k = 5; // кількість переобчислень
        int m = 7; // рухоме вікно

        double xmin = 0;
        double xmax = 2 * Math.PI;
        double iteration_lenght = (xmax - xmin) / n;
        double[] x = new double[n + 1];
        x[0] = xmin;
        for (int i = 1; i <= n; i++)
        {
            x[i] = x[i - 1] + iteration_lenght;
        }

        var sw = new Stopwatch();
        sw.Start();
        sw.Stop();
        TimeSpan sum = sw.Elapsed;
        TimeSpan minTime = default, maxTime = default;

        double[] func = DigitalFiltering.CreateSinus(x);
        //double[] func = DigitalFiltering.CreateTriangularSignal(x);
        //double[] func = DigitalFiltering.CreateRectangularSignal(x);
        //double[] func = DigitalFiltering.CreateDustySignal(x);
        double[] noisyFunc = DigitalFiltering.AddNoise(func);
        //double[] noisyFunc =
        //    { 4,6,3,2,5 };

        double[] coefs = DigitalFiltering.DanielWindow(m);
        //double[] coefs = DigitalFiltering.BartlettWindow(m);
        //double[] coefs = DigitalFiltering.HemingWindow(m);
        //double alpha = 5.0;
        //double[] coefs = DigitalFiltering.KaiserWindow(m, alpha);
        double[] smoothedFunc;
        double[] smoothedFunc2;
        

        int iter = 100;
        for (int i = 0; i < iter; i++)
        {

            var stopwatch = new Stopwatch();
            stopwatch.Start();

            //smoothedFunc = DigitalFiltering.Seq(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParThreadPool(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParParallelFor(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParLimThreadPool(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParLim(noisyFunc, coefs, k, m);

            //smoothedFunc2 = DigitalFiltering.ParLim2(noisyFunc, coefs, k, m);
            smoothedFunc2 = DigitalFiltering.ParLimNew2(noisyFunc, coefs, k, m);


            //smoothedFunc = DigitalFiltering.Seq(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.Seq2(noisyFunc, coefs, k, m);

            //smoothedFunc2 = DigitalFiltering.Par(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParThreadPool(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParParallelFor(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.Par2(noisyFunc, coefs, k, m);

            //smoothedFunc2 = DigitalFiltering.ParBranch(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParBranchThreadPool(noisyFunc, coefs, k, m);

            //smoothedFunc2 = DigitalFiltering.SeqBranch2(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParBranch2(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParBranch2Threads(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParBranch2Tasks(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParBranch2ThreadPool(noisyFunc, coefs, k, m);

            //smoothedFunc2 = DigitalFiltering.ParBranch2Bool(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParBranch2BoolTasks(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParBranch2BoolThreadPool(noisyFunc, coefs, k, m);

            //smoothedFunc2 = DigitalFiltering.ParLim(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParLimThreadPool(noisyFunc, coefs, k, m);

            //smoothedFunc2 = DigitalFiltering.ParLim2Threads(noisyFunc, coefs, k, m);
            //smoothedFunc2 = DigitalFiltering.ParLim2ThreadPool(noisyFunc, coefs, k, m);

            //smoothedFunc = DigitalFiltering.MySeq(noisyFunc, k, m);

            stopwatch.Stop();
            TimeSpan elapsedTime = stopwatch.Elapsed;
            sum += elapsedTime;
            if (minTime == default)
            {
                minTime = elapsedTime;
            }
            if (elapsedTime < minTime)
            {
                minTime = elapsedTime;
            }
            if (maxTime == default)
            {
                maxTime = elapsedTime;
            }
            if (elapsedTime > maxTime)
            {
                maxTime = elapsedTime;
            }

            Console.WriteLine($"Час розпаралелювання: {elapsedTime}");

            //for (int j = 0; j < 5; j++)
            //{
            //    Console.WriteLine(smoothedFunc[j]);
            //}

            //double[] check = new double[func.Length];
            //for (int j = 0; j < func.Length; j++)
            //{
            //    check[j] = smoothedFunc[j] - smoothedFunc2[j];
            //    Console.WriteLine(check[j]);
            //}

            //double deviation;
            //double sumOfDevioations = 0.0;
            //for (int j = 0; j < func.Length; j++)
            //{
            //    sumOfDevioations += Math.Abs(func[j] - smoothedFunc[j]);
            //}
            //deviation = sumOfDevioations / func.Length;

            //Console.WriteLine("Середнє відхилення: " + deviation);

            //double meanSquareDev;
            //double sumOfMeanSquareDev = 0.0;
            //for (int j = 0; j < func.Length; j++)
            //{
            //    sumOfMeanSquareDev = Math.Pow(func[j] - smoothedFunc[j], 2);
            //}
            //meanSquareDev = Math.Sqrt(sumOfMeanSquareDev / func.Length);
            //Console.WriteLine("Середнє квадратичне відхилення: " + meanSquareDev);
        }
        sum /= iter;
        Console.WriteLine($"Середнє: {sum}");
        Console.WriteLine($"Мін: {minTime}");
        Console.WriteLine($"Макс: {maxTime}");
        Console.ReadKey();
    }
}