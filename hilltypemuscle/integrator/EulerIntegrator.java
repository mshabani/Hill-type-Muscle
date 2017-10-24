package artisynth.models.hilltypemuscle.integrator;

import artisynth.models.hilltypemuscle.function.ODEFunction;

/**
 * This class implements Euler integrator See
 * <a href="https://www.math.ubc.ca/~feldman/math/vble.pdf">here</a> for
 * mathematical details.
 * 
 * @author Mohammad
 */
public class EulerIntegrator extends ODEIntegrator {

   /* ************************** Constructors ****************************/
   /**
    * 
    * @see {@link artisynth.models.hilltypemuscle.integrator.ODEIntegrator #ODEIntegrator(ODEFunction, double, double)}
    */
   public EulerIntegrator (ODEFunction func, double stepSize, double tol) {
      super (func, stepSize, tol);
   }

   /**
    * @see {@link artisynth.models.hilltypemuscle.integrator.ODEIntegrator #ODEIntegrator(ODEFunction, double)}
    */
   public EulerIntegrator (ODEFunction func, double stepSize) {
      super (func, stepSize);
   }

   /**
    * @see {@link artisynth.models.hilltypemuscle.integrator.ODEIntegrator #ODEIntegrator(ODEFunction)}
    */
   public EulerIntegrator (ODEFunction func) {
      super (func);
   }

   /* ************************ Implementation ****************************/
   /**
    * {@inheritDoc}
    */
   @Override
   public double step (double[] in, double[] out, double h) {

      for (int i = 0; i < (order - 1); i++) {
         out[0] = in[i + 2] * h + in[i + 1];
      }
      out[order - 1] = function.eval (in) * h + in[order];

      // We don't know error, return 0
      return 0;
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public double[] adaptiveStep (double[] in, double[] out, double h) {
      double currentStep = h;
      double nextStep = h;
      if (tol == 0) {
         step (in, out, currentStep);
      }
      else {
         double err = 1;
         double fullApx, halfApx;

         double[] tmpIn = new double[in.length];

         do {
            currentStep = nextStep;
            tmpIn = in.clone ();
            step (tmpIn, out, currentStep);
            fullApx = out[0];
            step (tmpIn, out, currentStep / 2);
            halfApx = out[0];
            tmpIn[0] = tmpIn[0] + currentStep / 2;
            for (int i = 0; i < out.length; i++) {
               tmpIn[i + 1] = out[i];
            }
            step (tmpIn, out, currentStep / 2);
            halfApx = out[0];
            err = Math.abs (fullApx - halfApx) / currentStep;
            nextStep = safetyFactor * tol * currentStep / err;
         }
         while (err > tol);
         out[0] = 2 * halfApx - fullApx;
      }
      return new double[] { currentStep, nextStep };
   }
}
