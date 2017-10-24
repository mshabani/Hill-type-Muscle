package artisynth.models.hilltypemuscle.function;

import maspack.function.MISOFunction;

/**
 * This is a base class for a general ordinary differential equation function in
 * form of (y_K denotes Kth derivative of y with respect to independent
 * variable): y_n = f(y_(n-1),y_(n-2),...,y,t)
 * 
 * for function y(t)
 * 
 * @author Mohammad
 *
 */
public abstract class ODEFunction implements MISOFunction {

   /*
    * order of the differential equation. This should be set at the
    * instantiation time.
    */
   private int order;

   /**
    * class constructor
    * 
    * @param order
    * order of the ordinary differential equation
    */
   public ODEFunction (int order) {
      this.order = order;
   }

   /**
    * @return Returns the input size of the ODE which should be equal to
    * order+1;
    */
   @Override
   public int getInputSize () {
      return order + 1;
   }

   public int getOrder () {
      return order;
   }

   /**
    * Abstract evaluation method.
    * 
    * @param in Array of doubles which is ODE input. Length if input array should be equal
    * to the order of the ODE, with the following format: in[] =
    * {t,y,y_1,y_2,...,y_(n-1)}, in which <code> t </code> is the independent
    * variable and y_K denotes Kth derivative of y (dependent variable)
    * @return A double value of function evaluation
    */
   @Override
   public abstract double eval (double[] in);

}
