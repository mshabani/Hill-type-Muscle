package artisynth.models.hilltypemuscle.integrator;

import artisynth.models.hilltypemuscle.function.ODEFunction;
import maspack.matrix.ImproperSizeException;
import maspack.util.DataBuffer;

/**
 * A test method for integrator. It integrates a function and calculates mean
 * absolute error of integration compared to analytical integration of the
 * function.
 * 
 * @author Mohammad
 *
 */
public class TestODEIntegrator {

   private class FunctionSecondOrder extends ODEFunction {
      // y = t*exp(-t)
      protected FunctionSecondOrder (int order) {
         super (order);
      }

      @Override
      public double eval (double[] in) throws ImproperSizeException {
         if (in.length != getInputSize ()) {
            String msg =
               "Input size for this function (" + this.getClass ().getName ()
               + ") is not compatible with order";
            throw new ImproperSizeException (msg);
         }
         return -Math.exp (-in[0]) - in[2];
      }
   }

   private class FunctionFirstOrder extends ODEFunction {
      // y = t*exp(-t)
      protected FunctionFirstOrder (int order) {
         super (order);
      }

      @Override
      public double eval (double[] in) throws ImproperSizeException {
         if (in.length != getInputSize ()) {
            String msg =
               "Input size for this function (" + this.getClass ().getName ()
               + ") is not compatible with order";
            throw new ImproperSizeException (msg);
         }
         return Math.exp (-in[0]) - in[1];
      }
   }

   public static void main (String args[]) {
      EulerIntegrator eulerIntegrator;
      RungeKuttaMersonIntegrator rungeKuttaMersonIntegrator;
      TestODEIntegrator instance = new TestODEIntegrator ();

      // test integration
      ODEFunction funcOrder2 = instance.new FunctionSecondOrder (2);
      ODEFunction funcOrder1 = instance.new FunctionFirstOrder (1);
      double[] initCondition = new double[] { 0, 0, 1 };
      double tol = 0.001;

      try {
         rungeKuttaMersonIntegrator =
            new RungeKuttaMersonIntegrator (funcOrder1, 0.01, tol);
         eulerIntegrator = new EulerIntegrator (funcOrder2, 0.01, tol);

         DataBuffer[] dataFirstOrder =
            rungeKuttaMersonIntegrator.integrate (new double[] { 0, 0 }, 1);
         DataBuffer[] dataSecondOrder =
            eulerIntegrator.integrate (initCondition, 1);

         int size1st = dataFirstOrder[0].dsize ();
         int size2nd = dataSecondOrder[0].dsize ();

         double error1st = 0;
         double error2nd = 0;

         for (int i = 0; i < size1st; i++) {
            double time = dataFirstOrder[0].dget ();
            double ny = dataFirstOrder[1].dget ();
            double na = time * Math.exp (-time);
            // System.out.println(time + " " + ny + " " + na);
            if (na != 0) {
               error1st += Math.abs (ny - na);
            }
         }
         double meanAbsoluteError1st = error1st / size1st;

         // System.out.println(meanAbsoluteError2nd);
         if (meanAbsoluteError1st < tol) {
            System.out.println (
               "Runge-Kutta-Merson integrator passed for sample function integration");
         }
         else {
            System.out.println (
               "Runge-Kutta-Merson integrator failed for sample function integration");
         }

         for (int i = 0; i < size2nd; i++) {
            double time = dataSecondOrder[0].dget ();
            double ny = dataSecondOrder[1].dget ();
            double na = time * Math.exp (-time);
            // System.out.println(time + " " + ny + " " + na);
            if (na != 0) {
               error2nd += Math.abs (ny - na);
            }
         }
         double meanAbsoluteError2nd = error2nd / size2nd;

         // System.out.println(meanAbsoluteError2nd);
         if (meanAbsoluteError2nd < tol) {
            System.out.println (
               "Euler integrator passed for sample function integration");
         }
         else {
            System.out.println (
               "Euler integrator failed for sample function integration");
         }

      }
      catch (Exception e) {
         System.out.println (e.getMessage ());
      }
   }

}
