package artisynth.models.hilltypemuscle.function;

import artisynth.models.hilltypemuscle.material.AxialHillTypeMuscleMaterial;

/**
 * This class implements activation dynamic function, which was used by
 * <a href="https://link.springer.com/article/10.1007/s10439-016-1591-9">
 * DeGroote et al. 2016</a>
 * 
 * @author Mohammad
 *
 */
public class DeGrooteActivationDynamicFunction extends ODEFunction {

   // The muscle material that uses this activation dynamic function
   protected AxialHillTypeMuscleMaterial material;

   public DeGrooteActivationDynamicFunction (AxialHillTypeMuscleMaterial mat) {
      super (1); // This function a first order ODE.
      material = mat;
   }

   /**
    * Function evaluation method.
    * 
    * @param in  A 2 element array consisting of {time, activation} 
    * @return A double value of function evaluation
    */
   @Override
   public double eval (double[] in) {

      double out = 0;
      double e = material.getOwnerMuscle().getNetExcitation ();
      double b = 0.1;
      double ta = 0.015; // activation time constant
      double td = 0.060; // deactivation time constant

      double f = 0.5 * Math.tanh (b * (e - in[1]));
      out += (f + 0.5) / (ta * (0.5 + 1.5 * in[1]));
      out += (0.5 + 1.5 * in[1]) * (-f + 0.5) / td;
      out *= (e - in[1]);
      return out;

   }
}