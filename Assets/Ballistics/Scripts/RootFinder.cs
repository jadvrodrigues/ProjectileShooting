using System;
using System.Linq;
using System.Numerics;

namespace EquationSolver
{
    /// <summary>
    /// Provides static functions to calculate the roots of equations following the pattern ax^n + bx^(n-1) + cx^(n-2) +... = 0. 
    /// </summary>
    public static class RootFinder
    {
        const double Delta = 0.00000000000001;

        /// <summary>
        /// Finds the unique roots of the equation ax^4+bx^3+cx^2+dx+e=0.
        /// </summary>
        /// <param name="a">Coefficient of x^4.</param>
        /// <param name="b">Coefficient of x^3.</param>
        /// <param name="c">Coefficient of x^2.</param>
        /// <param name="d">Coefficient of x.</param>
        /// <param name="e">The constant term.</param>
        /// <returns>The unique roots.</returns>
        public static Complex[] FindRoots(double a, double b, double c, double d, double e)
        {
            Complex[] roots;
            if (a != 0.0)
            {
                roots = new Complex[4];
                (roots[0], roots[1], roots[2], roots[3]) = Quartic(b / a, c / a, d / a, e / a);
            }
            else if (b != 0.0)
            {
                roots = new Complex[3];
                (roots[0], roots[1], roots[2]) = Cubic(c / b, d / b, e / b);
            }
            else if (c != 0.0)
            {
                roots = new Complex[2];
                (roots[0], roots[1]) = Quadratic(d / c, e / c);
            }
            else if (d != 0.0)
            {
                roots = new Complex[1] { Linear(e / d) };
            }
            else
            {
                roots = Array.Empty<Complex>();
            }

            return roots.Distinct().ToArray();
        }

        /// <summary>
        /// Finds the unique real roots of the equation ax^4+bx^3+cx^2+dx+e=0.
        /// </summary>
        /// <param name="a">Coefficient of x^4.</param>
        /// <param name="b">Coefficient of x^3.</param>
        /// <param name="c">Coefficient of x^2.</param>
        /// <param name="d">Coefficient of x.</param>
        /// <param name="e">The constant term.</param>
        /// <returns>The unique real roots.</returns>
        public static double[] FindRealRoots(double a, double b, double c, double d, double e)
        {
            return FindRoots(a, b, c, d, e).Where(root => Math.Abs(root.Imaginary) <= Delta).Select(root => root.Real).ToArray();
        }

        /// <summary>
        /// Finds the real root of the equation x+a=0.
        /// </summary>
        /// <param name="a">The constant term.</param>
        /// <returns>The real root.</returns>
        public static double Linear(double a)
        {
            return -a;
        }

        /// <summary>
        /// Finds the roots of the equation x^2+ax+b=0.
        /// </summary>
        /// <param name="a">Coefficient of x.</param>
        /// <param name="b">The constant term.</param>
        /// <returns>The two roots.</returns>
        public static (Complex, Complex) Quadratic(double a, double b)
        {
            Complex c = Complex.Sqrt(a * a - 4.0 * b);
            return ((c - a) / 2.0, (a + c) / -2.0);
        }


        // Source: https://en.wikipedia.org/wiki/Cubic_equation#General_cubic_formula
        /// <summary>
        /// Finds the roots of the equation x^3+ax^2+bx+c=0.
        /// </summary>
        /// <param name="a">Coefficient of x^2.</param>
        /// <param name="b">Coefficient of x.</param>
        /// <param name="c">The constant term.</param>
        /// <returns>The three roots.</returns>
        public static (Complex, Complex, Complex) Cubic(double a, double b, double c)
        {
            const double SQRT2_DIV_3 = 0.86602540378443864676;

            double d = c;
            c = b;
            b = a;
            a = 1.0;

            double delta_0 = b * b - 3.0 * a * c;
            double delta_1 = 2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d;

            if (Math.Abs(delta_0) <= Delta && Math.Abs(delta_1) <= Delta)
            {
                Complex root = new(-b / (3.0 * a), 0.0);
                return (root, root, root);
            }

            Complex extra0 = Complex.Pow(delta_1 * delta_1 - 4.0 * delta_0 * delta_0 * delta_0, 1.0 / 2.0);

            Complex C = Complex.Pow((delta_1 + extra0) / 2.0, 1.0 / 3.0);
            if (Complex.Abs(C) <= Delta) C = Complex.Pow((delta_1 - extra0) / 2.0, 1.0 / 3.0); // delta_0 == 0

            Complex xi = new(-0.5, SQRT2_DIV_3);

            return (x_k(0), x_k(1), x_k(2));

            Complex x_k(int k)
            {
                Complex xiRaisedToK = xi;
                for (int i = 0; i < k; i++) xiRaisedToK *= xi;
                return -1.0 / (3.0 * a) * (b + xiRaisedToK * C + delta_0 / (xiRaisedToK * C));
            }
        }

        // Source: https://www.had2know.org/academics/quartic-equation-solver-calculator.html
        /// <summary>
        /// Finds the roots of the equation x^4+ax^3+bx^2+cx+d=0.
        /// </summary>
        /// <param name="a">Coefficient of x^3.</param>
        /// <param name="b">Coefficient of x^2.</param>
        /// <param name="c">Coefficient of x.</param>
        /// <param name="d">The constant term.</param>
        /// <returns>The four roots.</returns>
        public static (Complex, Complex, Complex, Complex) Quartic(double a, double b, double c, double d)
        {
            double e = d;
            d = c;
            c = b;
            b = a;
            a = 1.0;

            Complex root1 = new();
            Complex root2 = new();
            Complex root3 = new();
            Complex root4 = new();

            if (b != 0.0 || d != 0.0)
            {
                if ((c == b * b / 4.0 + 2.0 * Math.Sqrt(e) && d == b * Math.Sqrt(e)) || (c == b * b / 4.0 - 2.0 * Math.Sqrt(e) && d == -b * Math.Sqrt(e)))
                {
                    double aa = b / 2.0;
                    double kk = d / (2 * aa);
                    if (aa * aa - 4 * kk < 0)
                    {
                        double rad = Math.Sqrt(4.0 * kk - aa * aa) / 2.0;
                        double cad = -aa / 2.0;
                        return (new(cad, rad), new(cad, rad), new(cad, -rad), new(cad, -rad));
                    }
                    else if (aa * aa - 4.0 * kk >= 0.0)
                    {
                        double rad = Math.Sqrt(aa * aa - 4.0 * kk) / 2.0;
                        double cad = -aa / 2.0;
                        return (new(cad - rad, 0.0), new(cad - rad, 0.0), new(cad + rad, 0.0), new(cad + rad, 0.0));
                    }
                }
                else
                {
                    double y = Rescubic(b, c, d, e);
                    double A = b * b / 4.0 - c + 2.0 * y;
                    double r = (b * y - d) / (-2.0 * A);
                    
                    double qc1a;
                    double qc1b;
                    double qc2a;
                    double qc2b;
                    if (A > 0.0)
                    {
                        qc1a = b / 2.0 - Math.Sqrt(A);
                        qc1b = y + Math.Sqrt(A) * r;
                        qc2a = b / 2.0 + Math.Sqrt(A);
                        qc2b = y - Math.Sqrt(A) * r;
                    }
                    else
                    {
                        qc1a = b / 2.0;
                        qc1b = y + Math.Sqrt(y * y - e);
                        qc2a = b / 2.0;
                        qc2b = y - Math.Sqrt(y * y - e);
                    }

                    double jim = -qc1a / 2.0;
                    double bob = qc1a * qc1a - 4.0 * qc1b;
                    if (bob < 0.0)
                    {
                        double rad = Math.Sqrt(4.0 * qc1b - qc1a * qc1a) / 2.0;
                        root1 = new(jim, rad);
                        root2 = new(jim, -rad);
                    }
                    if (bob >= 0)
                    {
                        double rad = jim + Math.Sqrt(bob) / 2.0;
                        double rod = jim - Math.Sqrt(bob) / 2.0;
                        root1 = new(rad, 0.0);
                        root2 = new(rod, 0.0);
                    }

                    double jen = -0.5 * qc2a;
                    double bab = qc2a * qc2a - 4.0 * qc2b;
                    if (bab < 0.0)
                    {
                        double royd = 0.5 * Math.Sqrt(4.0 * qc2b - qc2a * qc2a);
                        root3 = new(jen, royd);
                        root4 = new(jen, -royd);
                    }
                    if (bab >= 0)
                    {
                        double royd = jen + 0.5 * Math.Sqrt(bab);
                        double rood = jen - 0.5 * Math.Sqrt(bab);
                        root3 = new(royd, 0.0);
                        root4 = new(rood, 0.0);
                    }
                }

            }
            else
            {
                if (c * c - 4.0 * e >= 0.0)
                {
                    if (-c + Math.Sqrt(c * c - 4.0 * e) >= 0.0)
                    {
                        root1 = new(Math.Sqrt(-c / 2 + 0.5 * Math.Sqrt(c * c - 4 * e)), 0.0);
                        root2 = new(-Math.Sqrt(-c / 2 + 0.5 * Math.Sqrt(c * c - 4 * e)), 0.0);
                    }
                    else
                    {
                        root1 = new(0.0, Math.Sqrt(c / 2 - 0.5 * Math.Sqrt(c * c - 4 * e)));
                        root2 = new(0.0, -Math.Sqrt(c / 2 - 0.5 * Math.Sqrt(c * c - 4 * e)));
                    }

                    if (-c - Math.Sqrt(c * c - 4.0 * e) >= 0.0)
                    {
                        root3 = new(Math.Sqrt(-c / 2 - 0.5 * Math.Sqrt(c * c - 4 * e)), 0.0);
                        root4 = new(-Math.Sqrt(-c / 2 - 0.5 * Math.Sqrt(c * c - 4 * e)), 0.0);
                    }
                    else
                    {
                        root3 = new(0.0, Math.Sqrt(c / 2 + 0.5 * Math.Sqrt(c * c - 4 * e)));
                        root4 = new(0.0, -Math.Sqrt(c / 2 + 0.5 * Math.Sqrt(c * c - 4 * e)));
                    }
                }

                if (c * c - 4.0 * e < 0.0)
                {
                    double Az = -c / 2.0;
                    double Bz = 0.5 * Math.Sqrt(4.0 * e - c * c);
                    double y0 = Math.Sqrt(0.5 * (Math.Sqrt(Az * Az + Bz * Bz) - Az));

                    root1 = new(0.5 * Bz / y0, y0);
                    root2 = new(0.5 * Bz / y0, -y0);
                    root3 = new(-0.5 * Bz / y0, y0);
                    root4 = new(-0.5 * Bz / y0, -y0);
                }
            }

            return (root1, root2, root3, root4);

            static double Rescubic(double r, double s, double t, double u)
            {
                double xx = 0.0;

                double bb = -s / 2.0;
                double cc = r * t / 4.0 - u;
                double dd = s * u / 2.0 - t * t / 8.0 - r * r * u / 8.0;

                double disc = (18.0 * bb * cc * dd) - (4.0 * bb * bb * bb * dd) + (bb * bb * cc * cc) - (4.0 * cc * cc * cc) - (27.0 * dd * dd);
                double pp = cc - (bb * bb / 3.0); 
                double qq = ((2.0 / 27.0) * (bb * bb * bb)) - (bb * cc / 3.0) + dd; 
                double ff = (27.0 / 2.0) * qq;

                if (disc > 0.0)
                {
                    double x1 = -bb / 3.0 + 2.0 * Math.Sqrt(-pp / 3.0) * Math.Cos((1.0 / 3.0) * Math.Acos((1.5 * qq / pp) * Math.Sqrt(-3.0 / pp)));
                    double x2 = -bb / 3.0 + 2.0 * Math.Sqrt(-pp / 3.0) * Math.Cos((1.0 / 3.0) * Math.Acos((1.5 * qq / pp) * Math.Sqrt(-3.0 / pp)) + (2.0 / 3.0) * Math.PI);
                    double x3 = -bb / 3.0 + 2.0 * Math.Sqrt(-pp / 3.0) * Math.Cos((1.0 / 3.0) * Math.Acos((1.5 * qq / pp) * Math.Sqrt(-3.0 / pp)) - (2.0 / 3.0) * Math.PI);
                    double g1 = Math.Max(x1, Math.Max(x2, x3));
                    double g3 = Math.Min(x1, Math.Min(x2, x3));
                    double g2 = x1 + x2 + x3 - g1 - g3;
                    if (r * r / 4.0 - s + 2.0 * g1 > 0.0) {
                        xx = g1;
                    }
                    else
                    {
                        if (r * r / 4 - s + 2 * g2 > 0) xx = g2;
                        else xx = g3;
                    }
                }

                if (disc == 0)
                {
                    double x1 = -bb / 3.0 - (2.0 / 3.0) * Math.Cbrt(ff);
                    double x2 = -bb / 3.0 + (1.0 / 3.0) * Math.Cbrt(ff);
                    double g1 = Math.Max(x1, x2);
                    double g2 = Math.Min(x1, x2);
                    if (r * r / 4.0 - s + 2.0 * g1 > 0.0) xx = g1;
                    else xx = g2;
                }

                if (disc < 0)
                {
                    double g1 = -bb / 3.0 - (1.0 / 3.0) * Math.Cbrt(ff + 0.5 * Math.Sqrt(-27.0 * disc)) - (1.0 / 3.0) * Math.Cbrt(ff - 0.5 * Math.Sqrt(-27.0 * disc));
                    xx = g1;
                }

                return xx;
            }
        }
    }
}