using UnityEngine;
using EquationSolver;
using System.Linq;

/// <summary>
/// A utility class for ballistic calculations related to projectile targeting and interception.
/// </summary>
public static class Ballistics
{
    /// <summary>
    /// Calculates the shoot velocity required to hit a moving target, given a certain projectile speed.
    /// </summary>
    /// <param name="targetStartPos">The starting position of the target.</param>
    /// <param name="targetVelocity">The velocity vector of the target.</param>
    /// <param name="targetAcceleration">The acceleration vector of the target.</param>
    /// <param name="projectileStartPos">The starting position of the projectile.</param>
    /// <param name="projectileAcceleration">The acceleration vector of the projectile.</param>
    /// <param name="projectileSpeed">The speed at which the projectile will be shot.</param>
    /// <param name="duration">The estimated time for the projectile to reach the target. Infinity if unreachable.</param>
    /// <returns>If reachable, the required velocity vector to hit the target; otherwise, the zero vector.</returns>
    public static Vector3 CalculateShootVelocity(Vector3 targetStartPos, Vector3 targetVelocity, Vector3 targetAcceleration,
        Vector3 projectileStartPos, Vector3 projectileAcceleration, float projectileSpeed, out float duration)
    {
        Vector3 a = targetAcceleration - projectileAcceleration;
        Vector3 relativeTargetStartPos = targetStartPos - projectileStartPos;

        double[] realRoots = RootFinder.FindRealRoots(
            Vector3.Dot(a, a) / 4.0f,
            Vector3.Dot(a, targetVelocity),
            Vector3.Dot(a, relativeTargetStartPos) + Vector3.SqrMagnitude(targetVelocity) - projectileSpeed * projectileSpeed,
            2.0f * Vector3.Dot(targetVelocity, relativeTargetStartPos),
            Vector3.SqrMagnitude(relativeTargetStartPos)
            );

        float? smallestPositiveRoot = realRoots.Where(root => root > 0.0).Select(root => (float?)root).DefaultIfEmpty(null).Min();

        if (smallestPositiveRoot.HasValue)
        {
            duration = smallestPositiveRoot.Value;
            return CalculateShootVelocity(targetStartPos, targetVelocity, targetAcceleration, projectileStartPos, projectileAcceleration, duration);
        }
        else
        {
            duration = Mathf.Infinity;
            return Vector3.zero;
        }
    }

    /// <summary>
    /// Calculates the shoot velocity required to hit a moving target, after a given duration.
    /// </summary>
    /// <param name="targetStartPos">The starting position of the target.</param>
    /// <param name="targetVelocity">The velocity vector of the target.</param>
    /// <param name="targetAcceleration">The acceleration vector of the target.</param>
    /// <param name="projectileStartPos">The starting position of the projectile.</param>
    /// <param name="projectileAcceleration">The acceleration vector of the projectile.</param>
    /// <param name="duration">The time needed for the projectile to reach the target.</param>
    /// <returns>If the duration is valid, the calculated velocity vector required to hit the target; otherwise, the zero vector.</returns>
    public static Vector3 CalculateShootVelocity(Vector3 targetStartPos, Vector3 targetVelocity, Vector3 targetAcceleration,
        Vector3 projectileStartPos, Vector3 projectileAcceleration, float duration) 
    {
        if (duration <= 0f) return Vector3.zero;

        Vector3 relativeTargetStartPos = targetStartPos - projectileStartPos;
        Vector3 travel = relativeTargetStartPos + targetVelocity * duration + duration * duration * (targetAcceleration / 2.0f);
        return travel / duration - (projectileAcceleration / 2.0f) * duration;
    }
}
