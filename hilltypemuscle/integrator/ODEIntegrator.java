package artisynth.models.hilltypemuscle.integrator;

import artisynth.models.hilltypemuscle.function.ODEFunction;
import maspack.util.DataBuffer;

/**
 * A base class for an ODE integrator.
 * 
 * @author Mohammad
 *
 */
public abstract class ODEIntegrator {

   private final static double DEFAULT_STEP_SIZE = 0.01;
   private final static double DEFAULT_TOLERANCE = 0.0;

   protected double maxStepSize;
   protected ODEFunction function;
   protected int order;
   protected double tol;
   protected double safetyFactor = 0.9;

   /************************* constructors ********************************/

   /**
    * Class main constructor
    * 
    * @param function
    * The ODE function to be integrated
    * @param stepSize
    * Integrator step size. If tolerance set to be zero, integration will be
    * performed by constant step size of <code> stepSize </code>. Otherwise,
    * this <code> 
    *                 stepSize </code> will be used as maximum stepSize and
    * integrator will adjust <code> stepSize </code> at each step to satisfy
    * tolerance.
    * @param tol
    * Integration error tolerance
    */
   public ODEIntegrator (ODEFunction function, double stepSize, double tol) {
      this.function = function;
      this.order = function.getOrder ();
      this.maxStepSize = stepSize;
      this.tol = tol;
   }

   /**
    * Alternative constructor
    * 
    * @param function
    * The ODE function to be integrated
    * @param stepSize
    * Integrator step size. Since tolerance is not specified, it will be set to
    * zero, and integration will be performed by constant step size of
    * <code> stepSize </code>.
    */
   public ODEIntegrator (ODEFunction function, double stepSize) {
      this (function, stepSize, DEFAULT_TOLERANCE);
   }

   /**
    * Alternative constructor. Since no <code> stepSize </code> is specified,
    * default step size of 0.01 will be used. And tolerance
    * 
    * @param function
    * The ODE function to be integrated
    */
   public ODEIntegrator (ODEFunction function) {
      this (function, DEFAULT_STEP_SIZE, DEFAULT_TOLERANCE);
   }

   /* *********************** setter and getters ****************************/

   public ODEFunction getFunction () {
      return function;
   }

   public void setFunction (ODEFunction function) {
      this.function = function;
   }

   public double getTolerance () {
      return tol;
   }

   public void setTolerance (double tol) {
      this.tol = tol;
   }

   public double getStepSize () {
      return maxStepSize;
   }

   public void setStepSize (double stepSize) {
      this.maxStepSize = stepSize;
   }

   public int getODEOrder () {
      return order;
   }

   public void setODEOrder (int order) {
      this.order = order;
   }

   /* ************************ implementation ****************************/

   /**
    * Integrates function from initial time to finalTime. Initial time is the
    * first element at <code> initCondition </code>.
    * 
    * @param initCondition
    * Integration initial condition with the following format:
    * {t,y,y_1,y_2,...,y_(n-1)}, in which <code> t </code> is the independent
    * variable(time) and y_K denotes Kth derivative of y (dependent variable)
    * @param finalTime
    * Integration final time
    * @return DataBuffer[] Array of DataBuffers for output. First buffer saves
    * time, second buffer is y, third buffer is y_1 (dy/dt) and so on.
    */
   public DataBuffer[] integrate (double[] initCondition, double finalTime) {

      // finalTime should be greater than initial time
      if (!(initCondition[0] < finalTime)) {
         return null;
      }

      // create output DataBuffer and input and output arrays for step method
      double[] funcIn = new double[initCondition.length];
      funcIn = initCondition.clone ();
      double[] funcOut = new double[order];
      double currentTime;
      currentTime = funcIn[0];
      DataBuffer[] outBuffer = new DataBuffer[order + 1];
      for (int i = 0; i < (order + 1); i++) {
         outBuffer[i] = new DataBuffer ();
         outBuffer[i].dput (funcIn[i]);
      }
      double[] stepData = new double[] { maxStepSize, maxStepSize };

      // integrate
      do {
         stepData = adaptiveStep (funcIn, funcOut, stepData[1]);

         if ((funcIn[0] + stepData[0]) > finalTime) {
            stepData[1] = (finalTime - funcIn[0]);
            stepData = adaptiveStep (funcIn, funcOut, stepData[1]);
         }

         funcIn[0] += stepData[0];
         currentTime = funcIn[0];

         outBuffer[0].dput (currentTime);
         for (int i = 0; i < order; i++) {
            outBuffer[i + 1].dput (funcOut[i]);
            funcIn[i + 1] = funcOut[i];
         }

      }
      while (currentTime < finalTime);

      return outBuffer;
   }

   /**
    * This method basically implements
    * {@link ODEIntegrator#integrate(double[], double)} with convenient input
    * and output format and a signature which is convenient to use. This method
    * advances the ODE from time t0 to t0+h (t0 is first element of in[] array),
    * but with an adaptive step. It is useful for components that implement
    * {@link artisynth.core.mechmodels.HasAuxState} and need to advance from
    * time t0 to t1 but with adaptive step. This method advances to the new time
    * regardless of the number of the steps it need to take. 
    * 
    * @param in
    * Input array with the following format: {t,y,y_1,y_2,...,y_(n-1)},
    * <code> t </code> is the independent variable(time) and y_K denotes Kth
    * derivative of y (dependent variable)
    * @param out
    * Output array with the following format {y,y_1,y_2,..y_(n-1)}
    * @param newTime
    * New time to advance to.
    * @param h
    * Step size to be used as initial step in
    * {@link ODEIntegrator#adaptiveStep(double[], double[], double)} In a set of
    * consecutive calls to this function, method can be invoked with some
    * initial h and after first step, return value of the method can be used as
    * next h.
    * @return Returns a stepSize for next step. This can be used as input
    * <code> h </code> for this method in the next step.
    */
   public double adaptiveStepToNewTime (
      double[] in, double[] out, double newTime, double h) {

      double[] funcIn = new double[in.length];
      funcIn = in.clone ();
      double[] stepData = new double[] { maxStepSize, h };
      if (h == -1) {
         stepData[1] = maxStepSize;
      }
      do {
         stepData = adaptiveStep (funcIn, out, stepData[1]);

         if ((funcIn[0] + stepData[0]) > newTime) {
            stepData[1] = (newTime - funcIn[0]);
            stepData = adaptiveStep (funcIn, out, stepData[1]);
         }
         funcIn[0] += stepData[0];

         for (int i = 0; i < order; i++) {
            funcIn[i + 1] = out[i];
         }

      }
      while (funcIn[0] < newTime);

      return stepData[1];
   }

   /**
    * Should advance ODE function from time t0 = in[0] to t0+h
    * 
    * @param in
    * Input array with the following format: {t,y,y_1,y_2,...,y_(n-1)},
    * <code> t </code> is the independent variable(time) and y_K denotes Kth
    * derivative of y (dependent variable)
    * @param out
    * Output array with the following format {y,y_1,y_2,..y_(n-1)}
    * @param h
    * Step size
    * @return If the method has a notion of stepping error, it should return the
    * error, otherwise it should return 0.
    */

   public abstract double step (double[] in, double[] out, double h);

   /**
    * Performs adaptive step with changing step size to satisfy tolerance. So,
    * the final step that is used to advance in time might be much smaller than
    * initial step (h)
    * 
    * @param in
    * Input array with the following format: {t,y,y_1,y_2,...,y_(n-1)},
    * <code> t </code> is the independent variable(time) and y_K denotes Kth
    * derivative of y (dependent variable)
    * @param out
    * Output array with the following format {y,y_1,y_2,..y_(n-1)}
    * @param h
    * Step size
    * @return Returns a 2 element array, with first element being the step that
    * is used to advance ODE function, and second element as a guess for next
    * step.
    */
   public abstract double[] adaptiveStep (double[] in, double[] out, double h);

}
