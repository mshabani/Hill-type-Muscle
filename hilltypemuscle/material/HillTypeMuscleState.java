package artisynth.models.hilltypemuscle.material;

/**
 * An state object to hold muscle state information. It saves time , activation
 * at that time and fiber length at that time
 * 
 * @author Mohammad
 *
 */
public class HillTypeMuscleState {

   private double time;
   private double fiberLength;
   private double activation;

   public HillTypeMuscleState (double t, double a, double l) {
      setTime (t);
      setActivation (a);
      setFiberLength (l);
   }

   /**
    * a utility method to copy muscle state object
    * @param in an instance of the state
    */
   public void copy (HillTypeMuscleState in) {

      this.time = in.time;
      this.activation = in.activation;
      this.fiberLength = in.fiberLength;

   }

   public double getActivation () {
      return activation;
   }

   public void setActivation (double a) {
      activation = a;
   }

   public double getTime () {
      return time;
   }

   public void setTime (double t) {
      time = t;
   }

   public double getFiberLength () {
      return fiberLength;
   }

   public void setFiberLength (double length) {
      fiberLength = length;
   }

}
