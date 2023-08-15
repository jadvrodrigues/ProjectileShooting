using System;

/// <summary>
/// Finds the real roots of equations following the pattern a * x^n + b * x^(n-1) + ... = 0. 
/// </summary>
public static class FindRoots
{
    /// <summary>
    /// Find all real roots of the equation ax+b=0.
    /// </summary>
    /// <param name="a">Coefficient of x.</param>
    /// <param name="b">The constant term.</param>
    /// <returns>The real roots (possibly with repeated values).</returns>
    public static double[] Linear(double a, double b)
    {
        if (a == 0d) return null;
        return new double[] { -b / a };
    }

    /// <summary>
    /// Find all real roots of the equation ax^2+bx+c=0.
    /// </summary>
    /// <param name="a">Coefficient of x^2.</param>
    /// <param name="b">Coefficient of x.</param>
    /// <param name="c">The constant term.</param>
    /// <returns>The real roots (possibly with repeated values).</returns>
    public static double[] Quadratic(double a, double b, double c)
    {
        if (a == 0d) return Linear(b, c);

        double disc = b * b - 4 * a * c;
        if (disc < 0) return null;
        disc = Math.Sqrt(disc);

        return new double[]
        {
            (-b + disc) / (2d * a),
            (-b - disc) / (2d * a)
        };
    }

    // C Source: http://www.realitypixels.com/turk/opensource/FindCubicRoots.c.txt
    /// <summary>
    /// Find all real roots of the equation ax^3+bx^2+cx+d=0.
    /// </summary>
    /// <param name="a">Coefficient of x^3.</param>
    /// <param name="b">Coefficient of x^2.</param>
    /// <param name="c">Coefficient of x.</param>
    /// <param name="d">The constant term.</param>
    /// <returns>The real roots (possibly with repeated values).</returns>
    public static double[] Cubic(double a, double b, double c, double d)
    {
        if (a == 0d) return Quadratic(b, c, d);

        b /= a;
        c /= a;
        d /= a;

        double Q = (b * b - 3 * c) / 9;
        double R = (2 * b * b * b - 9 * b * c + 27 * d) / 54;
        double Qcubed = Q * Q * Q;
        double D = Qcubed - R * R;

        /* Three real roots */
        if (D >= 0)
        {
            double theta = Math.Acos(R / Math.Sqrt(Qcubed));
            double sqrtQ = Math.Sqrt(Q);
            return new double[3]
            {
                -2 * sqrtQ * Math.Cos(theta / 3) - b / 3,
                -2 * sqrtQ * Math.Cos((theta + 2 * Math.PI) / 3) - b / 3,
                -2 * sqrtQ * Math.Cos((theta + 4 * Math.PI) / 3) - b / 3
            };
        }
        /* One real root */
        else
        {
            double e = Math.Pow(Math.Sqrt(-D) + Math.Abs(R), 1d / 3d);
            if (R > 0)
                e = -e;
            return new double[1] 
            { 
                (e + Q / e) - b / 3d 
            };
        }
    }

    // Java Source: https://github.com/fpsunflower/sunflow/blob/master/src/org/sunflow/math/Solvers.java
    /// <summary>
    /// Find all real roots of the equation ax^4+bx^3+cx^2+dx+e=0.
    /// </summary>
    /// <param name="a">Coefficient of x^4.</param>
    /// <param name="b">Coefficient of x^3.</param>
    /// <param name="c">Coefficient of x^2.</param>
    /// <param name="d">Coefficient of x.</param>
    /// <param name="e">The constant term.</param>
    /// <returns>The real roots (possibly with repeated values).</returns>
    public static double[] Quartic(double a, double b, double c, double d, double e)
    {
        if (a == 0d) return Cubic(b, c, d, e);

        double inva = 1 / a;
        double c1 = b * inva;
        double c2 = c * inva;
        double c3 = d * inva;
        double c4 = e * inva;
        // cubic resolvant
        double c12 = c1 * c1;
        double p = -0.375 * c12 + c2;
        double q = 0.125 * c12 * c1 - 0.5 * c1 * c2 + c3;
        double r = -0.01171875 * c12 * c12 + 0.0625 * c12 * c2 - 0.25 * c1 * c3 + c4;
        double z = CubicForQuartic(-0.5 * p, -r, 0.5 * r * p - 0.125 * q * q);
        double d1 = 2.0 * z - p;
        if (d1 < 0)
        {
            if (d1 > 1.0e-10)
                d1 = 0;
            else
                return null;
        }
        double d2;
        if (d1 < 1.0e-10)
        {
            d2 = z * z - r;
            if (d2 < 0)
                return null;
            d2 = Math.Sqrt(d2);
        }
        else
        {
            d1 = Math.Sqrt(d1);
            d2 = 0.5 * q / d1;
        }
        // setup useful values for the quadratic factors
        double q1 = d1 * d1;
        double q2 = -0.25 * c1;
        double pm = q1 - 4 * (z - d2);
        double pp = q1 - 4 * (z + d2);
        if (pm >= 0 && pp >= 0)
        {
            // 4 roots (!)
            pm = Math.Sqrt(pm);
            pp = Math.Sqrt(pp);
            double[] results = new double[4];
            results[0] = -0.5 * (d1 + pm) + q2;
            results[1] = -0.5 * (d1 - pm) + q2;
            results[2] = 0.5 * (d1 + pp) + q2;
            results[3] = 0.5 * (d1 - pp) + q2;
            // tiny insertion sort
            for (int i = 1; i < 4; i++)
            {
                for (int j = i; j > 0 && results[j - 1] > results[j]; j--)
                {
                    double t = results[j];
                    results[j] = results[j - 1];
                    results[j - 1] = t;
                }
            }
            return results;
        }
        else if (pm >= 0)
        {
            pm = Math.Sqrt(pm);
            double[] results = new double[2];
            results[0] = -0.5 * (d1 + pm) + q2;
            results[1] = -0.5 * (d1 - pm) + q2;
            return results;
        }
        else if (pp >= 0)
        {
            pp = Math.Sqrt(pp);
            double[] results = new double[2];
            results[0] = 0.5 * (d1 - pp) + q2;
            results[1] = 0.5 * (d1 + pp) + q2;
            return results;
        }
        return null;
    }

    // Java Source: https://github.com/fpsunflower/sunflow/blob/master/src/org/sunflow/math/Solvers.java
    /// <summary>
    /// Return only one root for the specified cubic equation. 
    /// This routine is only meant to be called by the quartic solver. 
    /// It assumes the cubic is of the form: x^3+px^2+qx+r.
    /// </summary>
    static double CubicForQuartic(double p, double q, double r)
    {
        double A2 = p * p;
        double Q = (A2 - 3.0 * q) / 9.0;
        double R = (p * (A2 - 4.5 * q) + 13.5 * r) / 27.0;
        double Q3 = Q * Q * Q;
        double R2 = R * R;
        double d = Q3 - R2;
        double an = p / 3.0;
        if (d >= 0)
        {
            d = R / Math.Sqrt(Q3);
            double theta = Math.Acos(d) / 3.0;
            double sQ = -2.0 * Math.Sqrt(Q);
            return sQ * Math.Cos(theta) - an;
        }
        else
        {
            double sQ = Math.Pow(Math.Sqrt(R2 - Q3) + Math.Abs
           (R), 1.0 / 3.0);
            if (R < 0)
                return (sQ + Q / sQ) - an;
            else
                return -(sQ + Q / sQ) - an;
        }
    }
}
