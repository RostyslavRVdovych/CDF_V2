using System.Diagnostics;
using Algorithms;

class Program
{
    static void Main()
    {
        int n = 1000000; // кількість елементів
        int k = 1; // кількість переобчислень
        int m = 1; // рухоме вікно

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

        double[] func = DigitalFiltering.CreateSinus(x);
        //double[] func = DigitalFiltering.CreateTriangularSignal(x);
        //double[] func = DigitalFiltering.CreateRectangularSignal(x);
        //double[] func = DigitalFiltering.CreateDustySignal(x);
        double[] noisyFunc = DigitalFiltering.AddNoise(func);

        double[] coefs = DigitalFiltering.DanielWindow(m);
        //double[] coefs = DigitalFiltering.BartlettWindow(m);
        //double[] coefs = DigitalFiltering.HemingWindow(m);
        //double alpha = 5.0;
        //double[] coefs = DigitalFiltering.KaiserWindow(m, alpha);
        double[] smoothedFunc;

        int iter = 1000;
        for (int i = 0; i < iter; i++)
        {

            var stopwatch = new Stopwatch();
            stopwatch.Start();

            //smoothedFunc = DigitalFiltering.Seq(noisyFunc, coefs, k, m);
            //smoothedFunc = DigitalFiltering.Par(noisyFunc, coefs, k, m);
            //smoothedFunc = DigitalFiltering.ParBranch(noisyFunc, coefs, k, m);
            //smoothedFunc = DigitalFiltering.ParLim(noisyFunc, coefs, k, m);
            //smoothedFunc = DigitalFiltering.Seq2(noisyFunc, coefs, k, m);
            smoothedFunc = DigitalFiltering.Par2(noisyFunc, coefs, k, m);
            //smoothedFunc = DigitalFiltering.SeqBranch2(noisyFunc, coefs, k, m);
            //smoothedFunc = DigitalFiltering.ParBranch2(noisyFunc, coefs, k, m);
            //smoothedFunc = DigitalFiltering.ParLim2(noisyFunc, coefs, k, m);
            //smoothedFunc = DigitalFiltering.MySeq(noisyFunc, k, m);

            stopwatch.Stop();
            TimeSpan elapsedTime = stopwatch.Elapsed;
            sum += elapsedTime;

            Console.WriteLine($"Час розпаралелювання: {elapsedTime}");

            double deviation;
            double sumOfDevioations = 0.0;
            for (int j = 0; j < func.Length; j++)
            {
                sumOfDevioations += Math.Abs(func[j] - smoothedFunc[j]);
            }
            deviation = sumOfDevioations / func.Length;

            Console.WriteLine("Середнє відхилення: " + deviation);

            double meanSquareDev;
            double sumOfMeanSquareDev = 0.0;
            for (int j = 0; j < func.Length; j++)
            {
                sumOfMeanSquareDev = Math.Pow(func[j] - smoothedFunc[j], 2);
            }
            meanSquareDev = Math.Sqrt(sumOfMeanSquareDev / func.Length);
            Console.WriteLine("Середнє квадратичне відхилення: " + meanSquareDev);
        }
        sum /= iter;
        Console.WriteLine($"Середнє: {sum}");
    }
}