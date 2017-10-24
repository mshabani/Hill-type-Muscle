package artisynth.models.hilltypemuscle.function;

import artisynth.models.hilltypemuscle.material.AxialHillTypeMuscleMaterial;

/**
 * This class implements fiber length dynamic function, which was used by
 * <a href="https://link.springer.com/article/10.1007/s10439-016-1591-9">
 * DeGroote et al. 2016</a>
 * 
 * @author Mohammad
 *
 */
public class DeGrooteFiberLengthDynamicfunction extends ODEFunction {

   // The muscle material that uses this activation dynamic function
   AxialHillTypeMuscleMaterial material;

   public DeGrooteFiberLengthDynamicfunction (AxialHillTypeMuscleMaterial mat) {
      super (1); // This function a first order ODE.
      material = mat;
   }

   /**
    * Function evaluation method.
    * 
    * @param in
    * A 2 element array consisting of {time, fiber length}
    * @return A double value of function evaluation
    */
   @Override
   public double eval (double[] in) {

      // define temp variables
      double H =
         material.getOptLength () * Math.sin (material.getMyOptPenAngle ());
      double fiberLength = in[1];
      double normFiberLength = fiberLength / material.getOptLength ();
      double muscleLength = material.getOwnerMuscle().getLength ();

      // calculate muscle force and fiber force
      double tendonLength =
         muscleLength - Math.sqrt (fiberLength * fiberLength - H * H);
      double tendonForce =
         material.computeft (tendonLength / material.getMyTendonSlackLen ());
      double cosAlpha = (muscleLength - tendonLength) / fiberLength;
      double muscleForce = tendonForce / cosAlpha;

      // calculate force-velocity multiplier
      double activation =
         (material.getCurrentStateInfo().getActivation () > 1) ? 1
            : material.getCurrentStateInfo().getActivation ();

      double forceVelMultiplier =
         (muscleForce - material.computefpass (normFiberLength))
         / (activation * material.computefactl (normFiberLength));

      double normFiberVelocity = material.computeInvfv (forceVelMultiplier);

      // return fiber velocity
      return normFiberVelocity * material.getMyMaxFiberVel ()
      * material.getOptLength ();
   }
}