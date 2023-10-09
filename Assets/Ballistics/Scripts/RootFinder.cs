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

        // Source: https://www.1728.org/quartic2.htm
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

            double f = c - (3.0 * b * b / 8.0);
            double g = d + (b * b * b / 8.0) - (b * c / 2.0);
            double h = e - (3.0 * b * b * b * b / 256.0) + (b * b * c / 16.0) - (b * d / 4.0);

            Complex p;
            Complex q;

            var (y1, y2, y3) = Cubic(f / 2.0, (f * f - 4.0 * h) / 16.0, -g * g / 64.0);
            if (Math.Abs(y1.Imaginary) > Delta || Math.Abs(y2.Imaginary) > Delta || Math.Abs(y3.Imaginary) > Delta) // between y1, y2 and y3 there's one real and two complex numbers
            {
                p = Complex.Sqrt(Math.Abs(y1.Imaginary) > Delta ? y1 : y2);
                q = Complex.Sqrt(Math.Abs(y3.Imaginary) > Delta ? y3 : (Math.Abs(y2.Imaginary) > Delta ? y2 : y1));
            }
            else // y1, y2, y3 are real
            {
                (p, q, _) = SortByAbsoluteReal(y1, y2, y3);

                p = Complex.Sqrt(p);
                q = Complex.Sqrt(q);
            }

            Complex r = p * q != Complex.Zero ? -g / (8.0 * p * q) : Complex.Zero;
            double s = b / (4.0 * a);

            return (p + q + r - s, p - q - r - s, -p + q - r - s, -p - q + r - s);

            static (Complex first, Complex second, Complex third) SortByAbsoluteReal(Complex a, Complex b, Complex c)
            {
                (Complex first, Complex second, Complex third) result = (a, b, c);
                if (Math.Abs(result.first.Real) < Math.Abs(result.second.Real)) (result.first, result.second) = (result.second, result.first);
                if (Math.Abs(result.first.Real) < Math.Abs(result.third.Real)) (result.first, result.third) = (result.third, result.first);
                if (Math.Abs(result.second.Real) < Math.Abs(result.third.Real)) (result.second, result.third) = (result.third, result.second);

                return result;
            }
        }
    }
}