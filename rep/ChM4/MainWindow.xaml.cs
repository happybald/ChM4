using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using LiveCharts;
using LiveCharts.Wpf;
using MathNet.Numerics;

namespace ChM4
{
    public partial class MainWindow : System.Windows.Window
    {
        public MainWindow()
        {
            InitializeComponent();

            Func();

        }

        public SeriesCollection SeriesCollection { get; set; }
        public string[] Labels { get; set; }
        public Func<double, string> YFormatter { get; set; }


        class MyTable
        {
            public MyTable(double x, double y1, double y2)
            {
                this.X = x.ToString("F3");
                this.Ay = y1.ToString("F10");
                this.Ym = y2.ToString("F10");
            }
            public string X { get; set; }
            public string Ay { get; set; }
            public string Ym { get; set; }
        }

        private void Func()
        {
            DataContext = null;
            int N = Convert.ToInt32(Ntext.Text);
            double ax = -1, bx = 1;
            double h = (bx - ax) / (N);
            double[] x = new double[N + 1];
            List<double> y = new List<double>();
            double[,] ACB = new double[N, N];
            double[] F = new double[N];
            double a = Convert.ToDouble(Atext.Text);
            double YAnalytic(double x) => (Trig.Cos(Math.PI * x) + x + ((x * SpecialFunctions.Erf(x / Math.Sqrt(2 * a)) + Math.Sqrt((2 * a) / Math.PI) * Math.Exp(-(x * x) / (2 * a))) / (SpecialFunctions.Erf(1 / Math.Sqrt(2 * a)) + Math.Sqrt(2 * a / Math.PI) * Math.Exp(-1 / (2 * a)))));


            for (int i = 0; i <= N; i++)
            {
                x[i] = ax + i * h;
            }

            for (int i = 0; i <= N; i++)
            {
                y.Add(YAnalytic(x[i]));
            }



            var ACBmatrix = MathNet.Numerics.LinearAlgebra.Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {a*Math.PI*Math.PI+1, Math.PI},
                {a*Math.PI*Math.PI+Math.PI/4-1, (-8*a*Math.PI*Math.PI)/Math.Sqrt(2)-1}
            });
            var Fvector = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.DenseOfArray(new double[]
            {
                -Math.PI/2,
              //  -Math.PI/2+((a*Math.PI*Math.PI/8)*Math.Sin(Math.PI/8))+(Math.PI/4)*Math.Cos(Math.PI/8)-1+2*Math.Sin(Math.PI/8),
               -(1+a*Math.PI*Math.PI)-Math.PI/4
                // -(1+a*Math.PI*Math.PI)-Math.PI/4+Math.Sqrt(2)*(((a*Math.PI*Math.PI)/8)*Math.Sin(5*Math.PI/16)-(Math.PI/8)*Math.Cos(5*Math.PI/16)-1+2*Math.Sin(5*Math.PI/16))
            });     
            var Cvector = ACBmatrix.Solve(Fvector);

            double c1 = Cvector[0];
            double c2 = Cvector[1];


            List<double> yres = new List<double>();

            for (int i = 0; i <= N; i++)
            {
                yres.Add(x[i] + c1 * (1 - Math.Pow(x[i], 2)) + c2 * (1 + x[i] - Math.Pow(x[i], 2) - Math.Pow(x[i], 3)));
            }

            SeriesCollection = new SeriesCollection
            {
                new LineSeries
                {
                    Title = "Analytic Y",
                    Values = new ChartValues<double>(y.ToList<double>())
                },
                new LineSeries
                {
                    Title = "Gilerkin",
                    Values = new ChartValues<double>(yres)
                }
            };
            Labels = x.Select(x => x.ToString("F2")).ToArray();
            YFormatter = value => value.ToString("F2");
            DataContext = this;
            List<MyTable> result = new List<MyTable>();
            for (int i = 0; i <= N; i++)
            {
                result.Add(new MyTable(x[i], y[i], yres[i]));
            }

            grid.ItemsSource = result;
        }

        private void Btn1_Click(object sender, RoutedEventArgs e)
        {
            Func();
        }

    }
}