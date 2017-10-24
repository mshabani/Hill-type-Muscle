package artisynth.models.hilltypemuscle.integrator;

import artisynth.models.hilltypemuscle.function.ODEFunction;

/**
 * This class implements Runge-Kutta-Merson integrator See
 * <a href="https://www.math.ubc.ca/~feldman/math/vble.pdf">here</a> for
 * mathematical details.
 * 
 * @author Mohammad
 */
public class RungeKuttaMersonIntegrator extends ODEIntegrator {

   /* ************************** Constructors ****************************/

   /**
    * @see {@link artisynth.models.hilltypemuscle.integrator.ODEIntegrator #ODEIntegrator(ODEFunction, double, double)}
    * @throws IllegalOrderException
    * Throws a checked exception if order is not one. That is because, this
    * class only integrates first order ODE.
    */
   public RungeKuttaMersonIntegrator (ODEFunction func, double stepSize,
   double tol) throws IllegalOrderException {
      super (func, stepSize, tol);
      if (func.getOrder () != 1) {
         String message = "Order should be one for RungeKuttaMersonIntegrator";
         throw new IllegalOrderException (message);
      }
   }

   /**
    * @see {@link artisynth.models.hilltypemuscle.integrator.ODEIntegrator #ODEIntegrator(ODEFunction, double)}
    * @throws IllegalOrderException
    * Throws a checked exception if order is not one. That is because, this
    * class only integrates first order ODE.
    */
   public RungeKuttaMersonIntegrator (ODEFunction func, double stepSize)
   throws IllegalOrderException {
      this (func, stepSize, 0);
   }

   /**
    * @see {@link artisynth.models.hilltypemuscle.integrator.ODEIntegrator #ODEIntegrator(ODEFunction)}
    * @throws IllegalOrderException
    * Throws a checked exception if order is not one. That is because, this
    * class only integrates first order ODE.
    */
   public RungeKuttaMersonIntegrator (ODEFunction func)
   throws IllegalOrderException {
      this (func, 0.01, 0);
   }

   /* ************************ Implementation ****************************/
   /**
    * {@inheritDoc}
    */
   @Override
   public double step (double[] in, double[] out, double h) {
      double[] tmpIn = new double[in.length];
      tmpIn = in.clone ();

      double k1 = function.eval (tmpIn);

      tmpIn[0] = in[0] + h / 3;
      tmpIn[1] = in[1] + k1 * h / 3;
      double k2 = function.eval (tmpIn);

      tmpIn[1] = in[1] + k1 * h / 6 + k2 * h / 6;
      double k3 = function.eval (tmpIn);

      tmpIn[0] = in[0] + h / 2;
      tmpIn[1] = in[1] + k1 * h / 8 + 3 * k3 * h / 8;
      double k4 = function.eval (tmpIn);

      tmpIn[0] = in[0] + h;
      tmpIn[1] = in[1] + k1 * h / 2 - 3 * k3 * h / 2 + 2 * h * k4;
      double k5 = function.eval (tmpIn);

      double A1 = in[1] + h * (k1 / 2 - 3 * k3 / 2 + 2 * k4);
      double A2 = in[1] + h * (k1 / 6 + 2 * k4 / 3 + k5 / 6);

      double E = (A1 - A2) / 5;
      double err = E / h;

      out[0] = A2 - E;

      return err;
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
         double[] tmpIn = new double[in.length];

         do {
            currentStep = nextStep;
            tmpIn = in.clone ();
            err = step (tmpIn, out, currentStep);
            nextStep = safetyFactor * currentStep * Math.pow (tol / err, 0.25);

         }
         while (err > tol);
      }
      nextStep = nextStep < maxStepSize ? nextStep : maxStepSize;
      return new double[] { currentStep, nextStep };
   }

   private class IllegalOrderException extends Exception {

      private static final long serialVersionUID = -8631010733740830267L;

      public IllegalOrderException (String message) {
         super (message);
      }
   }
}
