using UnityEngine;
using UnityEngine.UI;

public class ShootFromSpeed : MonoBehaviour
{
    [Header("Shooting")]
    [SerializeField] float shootSpeed;
    [SerializeField] KeyCode shootHotkey = KeyCode.Space;

    [Header("Projectile")]
    [SerializeField] Rigidbody projectile;
    [SerializeField] Vector3 projectileStartPos;
    Vector3 projectileAcceleration = Physics.gravity;

    [Header("Target")]
    [SerializeField] Rigidbody target;
    [SerializeField] Vector3 targetStartPos;
    [SerializeField] Vector3 targetStartVelocity;
    Vector3 targetAcceleration = Physics.gravity;

    [Header("UI")]
    [SerializeField] Text text;

    public void Shoot()
    {
        Vector3 projectileStartVelocity = Ballistics.CalculateShootVelocity(targetStartPos, targetStartVelocity,
            targetAcceleration, projectileStartPos, projectileAcceleration, shootSpeed, out _);

        ShootRigidbody(projectile, projectileStartPos, projectileStartVelocity);
        ShootRigidbody(target, targetStartPos, targetStartVelocity);

        static void ShootRigidbody(Rigidbody rigidbody, Vector3 position, Vector3 velocity)
        {
            rigidbody.position = position;
            rigidbody.velocity = velocity;
        }
    }

    void Update()
    {
        if (Input.GetKeyDown(shootHotkey)) Shoot();
    }

    void OnValidate()
    {
        if(text != null) text.text = shootSpeed + " m/s";
    }
}
