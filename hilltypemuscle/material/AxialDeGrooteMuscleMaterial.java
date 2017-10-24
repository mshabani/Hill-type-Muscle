package artisynth.models.hilltypemuscle.material;


/**
 * This class extends base class for Hill-type muscle material and implements
 * necessary curve characteristics. See
 * <a href="https://link.springer.com/article/10.1007/s10439-016-1591-9">
 * DeGroote et al. 2016</a> for details.
 * 
 * @author Mohammad
 */
public class AxialDeGrooteMuscleMaterial extends AxialHillTypeMuscleMaterial {

   // *********************************************************************/
   // *************** Characteristic curves constants *********************/
   // *********************************************************************/
   // Tendon Force-Length
   private static final double C1 = 0.200;
   private static final double C2 = 0.995;
   private static final double C3 = 0.250;
   private static final double KT = 35;

   // Active Muscle Force-Length
   private static final double B11 = 0.814483478343008;
   private static final double B21 = 1.05503342897057;
   private static final double B31 = 0.162384573599574;
   private static final double B41 = 0.0633034484654646;
   private static final double B12 = 0.433004984392647;
   private static final double B22 = 0.716775413397760;
   private static final double B32 = -0.0299471169706956;
   private static final double B42 = 0.200356847296188;
   private static final double B13 = 0.100;
   private static final double B23 = 1.000;
   private static final double B33 = 0.353553390593274;
   private static final double B43 = 0.000;

   // passive Muscle Force-Length
   private static final double KPE = 4.0;
   private static final double E0 = 0.6;
   private static final double E1 = -0.995172050006169;

   // Muscle Force-Velocity
   private static final double D1 = -0.318323436899127;
   private static final double D2 = -8.149156043475250;
   private static final double D3 = -0.374121508647863;
   private static final double D4 = 0.885644059915004;

   // ********************************************************************/
   // ************************* material constructors ********************/
   // ********************************************************************/

   /**
    * Default constructor, This constructor uses default values that are for the
    * experimental test case with tendon compliance
    */
   public AxialDeGrooteMuscleMaterial () {
      myOptLength = DEFAULT_OPT_LENGTH;
      myMaxForce = DEFAULT_MAX_FORCE;
      myTendonSlackLen = DEFAULT_TENDON_SLACK_LENGTH;
      myOptPenAngle = DEFAULT_OPTIMAL_PENNATION_ANGLE;
      myMaxFiberVel = DEFAULT_MAXIMUM_FIBER_Velocity;
      myDamping = DEFAULT_DAMPING;
      setMyDefaultActivation (DEFAULT_ACTIVATION);
   }

   /**
    * Main constructor of the material
    * 
    * @param optimalFiberLength
    * @param maximumForce
    * @param tendonSlackLength
    * @param optimalPennationAngle
    * @param maxFiberVelocity
    * @param damping
    * this is not muscle damping rather an extra damping parameter in force
    * equation to remove singularities from the equation See <a href=
    * "http://simtk-confluence.stanford.edu:8080/display/OpenSim/Millard+2012+Muscle+Models">
    * here </a> EQ 3 for more information about this dmaping parameter
    * @param defaultActivation
    */
   public AxialDeGrooteMuscleMaterial (double optimalFiberLength,
   double maximumForce, double tendonSlackLength, double optimalPennationAngle,
   double maxFiberVelocity, double damping, double defaultActivation) {

      myOptLength = optimalFiberLength;
      myMaxForce = maximumForce;
      myTendonSlackLen = tendonSlackLength;
      myOptPenAngle = optimalPennationAngle;
      myMaxFiberVel = maxFiberVelocity;
      myDamping = damping;
      setMyDefaultActivation (defaultActivation);
   }

   // ***********************************************************************/
   // * implementation of methods inherited from AxialHillTypeMuscleMaterial */
   // ***********************************************************************/

   /**
    * {@inheritDoc}
    */
   @Override
   public double computefactl (double normFiberLength) {

      double activeF_L =
         B11
         * Math.exp (
            ((-0.5 * Math.pow ((normFiberLength - B21), 2))
            / Math.pow ((normFiberLength * B41 + B31), 2)))
         + B12 * Math.exp (
            ((-0.5 * Math.pow ((normFiberLength - B22), 2))
            / Math.pow ((normFiberLength * B42 + B32), 2)))
         + B13 * Math.exp (
            ((-0.5 * Math.pow ((normFiberLength - B23), 2))
            / Math.pow ((normFiberLength * B43 + B33), 2)));

      return activeF_L;
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public double computefv (double normFiberVelocity) {
      return D1 * Math.log (
         (D2 * normFiberVelocity + D3)
         + Math.sqrt (Math.pow ((D2 * normFiberVelocity + D3), 2) + 1))
      + D4;
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public double computefpass (double normFiberLength) {
      return (Math.exp (KPE * (normFiberLength - 1) / E0) - 1 - E1)
      / (Math.exp (KPE) - 1);
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public double computeft (double normTendonLength) {
      return C1 * Math.exp (KT * (normTendonLength - C2)) - C3;
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public double computeDfactlmDnormlm (double normFiberLength) {
      double M, N, df;

      M = normFiberLength - B21;
      N = B31 + B41 * normFiberLength;
      df =
         B11 * ((M * M * B41 - M * N) / (N * N * N))
         * Math.exp (-0.5 * M * M / (N * N));

      M = normFiberLength - B22;
      N = B32 + B42 * normFiberLength;
      df +=
         B12 * ((M * M * B42 - M * N) / (N * N * N))
         * Math.exp (-0.5 * M * M / (N * N));

      M = normFiberLength - B23;
      N = B33 + B43 * normFiberLength;
      df +=
         B13 * ((M * M * B43 - M * N) / (N * N * N))
         * Math.exp (-0.5 * M * M / (N * N));

      return df;

   }

   /**
    * {@inheritDoc}
    */
   @Override
   public double computeDfvDnormvm (double normFiberVelocity) {
      return D1 * D2
      / (Math.sqrt (Math.pow (D2 * normFiberVelocity + D3, 2) + 1));
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public double computeDfpassDnormlm (double normFiberLength) {
      return KPE * Math.exp (KPE * (normFiberLength - 1) / E0)
      / (E0 * (Math.exp (KPE) - 1));
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public double computeInvfv (double forceVelocityMultiplier) {
      return (Math.sinh ((forceVelocityMultiplier - D4) / D1) - D3) / D2;
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public double computeDftDnormlt (double normTendonLength) {
      return C1 * KT * Math.exp (KT * (normTendonLength - C2));
   }

}
