package artisynth.models.hilltypemuscle.material;

import artisynth.core.materials.AxialMuscleMaterial;
import artisynth.core.mechmodels.HasAuxState;
import artisynth.models.hilltypemuscle.HillTypeMuscle;
import artisynth.models.hilltypemuscle.function.DeGrooteActivationDynamicFunction;
import artisynth.models.hilltypemuscle.function.DeGrooteFiberLengthDynamicfunction;
import artisynth.models.hilltypemuscle.function.ODEFunction;
import artisynth.models.hilltypemuscle.integrator.RungeKuttaMersonIntegrator;
import maspack.properties.PropertyList;
import maspack.properties.PropertyMode;
import maspack.properties.PropertyUtils;
import maspack.util.DataBuffer;

/**
 * This abstract class serves as a base for an equilibrium Hill-Type muscle
 * model. Equilibrium means that tendon force is equal to the fiber force along
 * tendon. Any Hill type muscle model should extend this class and implement
 * model characteristic curves: active force length(fa), force velocity(fv),
 * passive force length (fpass), and tendon force (ft) and their derivative.
 * 
 * @see {@link artisynth.core.materials.AxialMuscleMaterial}
 * @see {@link artisynth.models.hilltypemuscle.RBDeGrooteMuscle} for example usage
 * @author Mohammad
 */
public abstract class AxialHillTypeMuscleMaterial extends AxialMuscleMaterial
implements HasAuxState {

   // *******************************************************************/
   // ************************ default values ***************************/
   // *******************************************************************/

   protected static double DEFAULT_MAX_FORCE = 1.2758948506;
   protected static double DEFAULT_OPT_LENGTH = 0.0171;
   protected static double DEFAULT_TENDON_SLACK_LENGTH = 0.0171;
   protected static double DEFAULT_OPTIMAL_PENNATION_ANGLE = 0.10471975511966; // Unit:
                                                                               // Radian
   protected static double DEFAULT_MAXIMUM_FIBER_Velocity = 5.3156293513;
   protected static double DEFAULT_ACTIVATION = 0.01;

   // *******************************************************************/
   // ************************* properties ******************************/
   // *******************************************************************/

   protected double myTendonSlackLen;
   protected double myOptPenAngle;
   protected double myMaxFiberVel;
   protected double myDefaultActivation;

   protected PropertyMode myTendonSlackLenMode = PropertyMode.Inherited;
   protected PropertyMode myOptPenAngleMode = PropertyMode.Inherited;
   protected PropertyMode myMaxFiberVelMode = PropertyMode.Inherited;
   protected PropertyMode myDefaultActivationMode = PropertyMode.Inherited;

   public static PropertyList myProps =
      new PropertyList (
         AxialHillTypeMuscleMaterial.class, AxialMuscleMaterial.class);

   static {
      myProps.addInheritable (
         "myTendonSlackLen *", "Tendon slack length",
         DEFAULT_TENDON_SLACK_LENGTH);
      myProps.addInheritable (
         "myOptPenAngle * *", "Optimal pennation angle",
         DEFAULT_OPTIMAL_PENNATION_ANGLE);
      myProps.addInheritable (
         "myMaxFiberVel * *", "Maximum fiber evlocity",
         DEFAULT_MAXIMUM_FIBER_Velocity);
      myProps.addInheritable (
         "myDefaultActivation * *", "Default muscle activation",
         DEFAULT_ACTIVATION);
   }

   // *******************************************************************/
   // ************************* variables *******************************/
   // *******************************************************************/

   private HillTypeMuscleInfo muscleInfo;
   private HillTypeMuscleState nextStateInfo;
   private HillTypeMuscleState currentStateInfo;
   private ODEFunction activationDynamicFunction;
   private ODEFunction fiberLengthDynamicFunction;
   private RungeKuttaMersonIntegrator activationIntegrator;
   private RungeKuttaMersonIntegrator fiberLengthIntegrator;
   private double adaptiveActStep;
   private double adaptiveFiberStep;

   private boolean ignoreActDyn = false;
   private boolean ignoreTendonComp = false;
   private boolean useSolver = false;
   private HillTypeMuscle ownerMuscle;

   private double activationIntegratorTimeStep = 0.0001; // default values are
                                                         // good to use
   private double activationIntegratorTolerance = 0.00001;
   private double fiberLengthIntegratorTimeStep = 0.000001;
   private double fiberLengthIntegratorTolerance = 0.0000000001;

   // **************************************************************//
   // *************** Setter and Getter (properties) ***************//
   // **************************************************************//

   public PropertyList getAllPropertyInfo () {
      return myProps;
   }

   public double getMyTendonSlackLen () {
      return myTendonSlackLen;
   }

   public synchronized void setMyTendonSlackLen (double slackL) {
      myTendonSlackLen = slackL;
      myTendonSlackLenMode =
         PropertyUtils.propagateValue (
            this, "myTendonSlackLen", myTendonSlackLen, myTendonSlackLenMode);
      notifyHostOfPropertyChange ();
   }

   public PropertyMode getMyTendonSlackLenMode () {
      return myTendonSlackLenMode;
   }

   public void setMyTendonSlackLenMode (PropertyMode mode) {
      myTendonSlackLenMode =
         PropertyUtils.setModeAndUpdate (
            this, "myTendonSlackLenMode", myTendonSlackLenMode, mode);
   }

   public double getMyOptPenAngle () {
      return myOptPenAngle;
   }

   public synchronized void setMyOptPenAngle (double angle) {
      myOptPenAngle = angle;
      myOptPenAngleMode =
         PropertyUtils.propagateValue (
            this, "myOptPenAngle", myOptPenAngle, myOptPenAngleMode);
      notifyHostOfPropertyChange ();
   }

   public PropertyMode getMyOptPenAngleMode () {
      return myOptPenAngleMode;
   }

   public void setMyOptPenAngleMode (PropertyMode mode) {
      myOptPenAngleMode =
         PropertyUtils.setModeAndUpdate (
            this, "myOptPenAngleMode", myOptPenAngleMode, mode);
   }

   public double getMyDefaultActivation () {
      return myDefaultActivation;
   }

   public synchronized void setMyDefaultActivation (double vel) {
      myDefaultActivation = vel;
      myDefaultActivationMode =
         PropertyUtils.propagateValue (
            this, "myDefaultActivation", myDefaultActivation,
            myDefaultActivationMode);
      notifyHostOfPropertyChange ();
   }

   public PropertyMode getMyDefaultActivationMode () {
      return myDefaultActivationMode;
   }

   public void setMyDefaultActivationMode (PropertyMode mode) {
      myDefaultActivationMode =
         PropertyUtils.setModeAndUpdate (
            this, "myDefaultActivationMode", myDefaultActivationMode, mode);
   }

   public double getMyMaxFiberVel () {
      return myMaxFiberVel;
   }

   public synchronized void setMyMaxFiberVel (double vel) {
      myMaxFiberVel = vel;
      myMaxFiberVelMode =
         PropertyUtils.propagateValue (
            this, "maxFiberVel", myMaxFiberVel, myMaxFiberVelMode);
      notifyHostOfPropertyChange ();
   }

   public PropertyMode getMyMaxFiberVelMode () {
      return myMaxFiberVelMode;
   }

   public void setMyMaxFiberVelMode (PropertyMode mode) {
      myMaxFiberVelMode =
         PropertyUtils.setModeAndUpdate (
            this, "maxFiberVelMode", myMaxFiberVelMode, mode);
   }

   // **************************************************************//
   // *********** Setter and Getter (internal variables) ***********//
   // **************************************************************//

   public boolean isUseSolver () {
      return useSolver;
   }

   public void setUseSolver (boolean useSolver) {
      this.useSolver = useSolver;
   }

   public double getTendonLength () {
      return muscleInfo.getTendonLength();
   }

   public void setOwner (HillTypeMuscle mus) {
      setOwnerMuscle (mus);
   }

   public boolean isIgnoreTendonComp () {
      return ignoreTendonComp;
   }

   public void setIgnoreTendonComp (boolean ignoreTendonComp) {
      this.ignoreTendonComp = ignoreTendonComp;
      if (ignoreTendonComp) {
         setUseSolver (false); // we have rigid tendon, so solver is not
                               // required
      }
   }

   public boolean isIgnoreActDyn () {
      return ignoreActDyn;
   }

   public void setIgnoreActDyn (boolean ignore) {
      this.ignoreActDyn = ignore;
   }

   public HillTypeMuscle getOwnerMuscle () {
      return ownerMuscle;
   }

   public void setOwnerMuscle (HillTypeMuscle ownerMuscle) {
      this.ownerMuscle = ownerMuscle;
   }

   public HillTypeMuscleState getCurrentStateInfo () {
      return currentStateInfo;
   }

   public void setCurrentStateInfo (HillTypeMuscleState currentStateInfo) {
      this.currentStateInfo = currentStateInfo;
   }

   public HillTypeMuscleState getNextStateInfo () {
      return nextStateInfo;
   }

   public void setNextStateInfo (HillTypeMuscleState nextStateInfo) {
      this.nextStateInfo = nextStateInfo;
   }

   public double getActivationIntegratorTimeStep () {
      return activationIntegratorTimeStep;
   }

   public void setActivationIntegratorTimeStep (
      double activationIntegratorTimeStep) {
      this.activationIntegratorTimeStep = activationIntegratorTimeStep;
   }

   public double getActivationIntegratorTolerance () {
      return activationIntegratorTolerance;
   }

   public void setActivationIntegratorTolerance (
      double activationIntegratorTolerance) {
      this.activationIntegratorTolerance = activationIntegratorTolerance;
   }

   public double getFiberLengthIntegratorTimeStep () {
      return fiberLengthIntegratorTimeStep;
   }

   public void setFiberLengthIntegratorTimeStep (
      double fiberLengthIntegratorTimeStep) {
      this.fiberLengthIntegratorTimeStep = fiberLengthIntegratorTimeStep;
   }

   public double getFiberLengthIntegratorTolerance () {
      return fiberLengthIntegratorTolerance;
   }

   public void setFiberLengthIntegratorTolerance (
      double fiberLengthIntegratorTolerance) {
      this.fiberLengthIntegratorTolerance = fiberLengthIntegratorTolerance;
   }

   // ************************************************************//
   // ****************** utility functions***********************//
   // ************************************************************//

   /**
    * A utility function to initialize integrators and initial conditions This
    * method should be called prior to any call to computeF method
    * @param t time at which initialization should be done
    */
   public void initMaterialFromProperties (double t) {

      muscleInfo = new HillTypeMuscleInfo (this);
      
      double initAct =
      isIgnoreActDyn () ? getOwnerMuscle ().getNetExcitation ()
         : getMyDefaultActivation ();
      double initFiberLength = getOptLength ();
      
      setCurrentStateInfo (
         new HillTypeMuscleState (t, initAct, initFiberLength));
      setNextStateInfo (new HillTypeMuscleState (t, initAct, initFiberLength));
      fiberLengthDynamicFunction =
         new DeGrooteFiberLengthDynamicfunction (this);
      activationDynamicFunction = new DeGrooteActivationDynamicFunction (this);
      try {
         activationIntegrator =
            new RungeKuttaMersonIntegrator (
               activationDynamicFunction, getActivationIntegratorTimeStep (),
               activationIntegratorTolerance);
         fiberLengthIntegrator =
            new RungeKuttaMersonIntegrator (
               fiberLengthDynamicFunction, fiberLengthIntegratorTimeStep,
               fiberLengthIntegratorTolerance);
         adaptiveActStep = getActivationIntegratorTimeStep();
         adaptiveFiberStep = this.getFiberLengthIntegratorTimeStep();
      }
      catch (Exception e) {
         System.out.println (e.getMessage ());
      }
      equilibriateMaterial ();
      getNextStateInfo ().copy (getCurrentStateInfo ());
   }

   /**
    * Equilibriates material at initial time t
    */
   private void equilibriateMaterial () {
      muscleInfo.equilibriate ();
   }

   // ************************************************************//
   // implementation of methods inherited from AxialMuscleMaterial//
   // ************************************************************//

   /**
    * {@inheritDoc}
    */
   @Override
   public double computeF (double l, double ldot, double l0, double ex) {

      if (isUseSolver ()) {
         equilibriateMaterial ();
      }
      muscleInfo.calcMuscleDynamicInfo (l, ldot);

      
      double timeIndex = getCurrentStateInfo ().getTime () * 100 + 0.0000001;
      if ((timeIndex - Math.floor (timeIndex)) < 0.00001) {
         System.out.print (
            getCurrentStateInfo ().getTime () + "  " + muscleInfo.getFiberForce()
            * muscleInfo.getCosAlpha());
         System.out
            .print ("  " + getMaxForce () * computeft (muscleInfo.getNormTendonLength()));

         System.out.print ("   " + getCurrentStateInfo ().getActivation ());
         System.out.print ("    " + muscleInfo.getFiberVelocity());
         // System.out.println("vm " +muscleInfo.vm);
         System.out.print ("     " + muscleInfo.getTendonLength());
         System.out.print ("    " + muscleInfo.normFiberLength);
         System.out.println ("    " + l);
      }
      
      return muscleInfo.getTendonForce();
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public double computeDFdl (double l, double ldot, double l0, double ex) {

      if (isUseSolver ()) {
         equilibriateMaterial ();
      }
      muscleInfo.calcMuscleDynamicInfo (l, ldot);
      
      return muscleInfo.getMuscleStiffness();
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public double computeDFdldot (double l, double ldot, double l0, double ex) {

      if (isUseSolver ()) {
         equilibriateMaterial ();
      }
      muscleInfo.calcMuscleDynamicInfo (l, ldot);
      
      return muscleInfo.getMuscleDamping();
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public boolean isDFdldotZero () {
      return false;
   }

   // ************************************************************//
   // ************************abstract methods********************//
   // ************************************************************//

   /**
    * Computes active force length multiplier for normalized fiber length
    * 
    * @param normlm
    * normalized fiber length
    * @return active fiber force length multiplier
    */
   public abstract double computefactl (double normlm);

   /**
    * Computes force velocity multiplier for normalized fiber velocity
    * 
    * @param normvm
    * normalized fiber velocity
    * @return force velocity multiplier
    */
   public abstract double computefv (double normvm);

   /**
    * Computes passive fiber force
    * 
    * @param normlm
    * normalize fiber length
    * @return passive fiber force
    */
   public abstract double computefpass (double normlm);

   /**
    * Computes tendon force for normalized tendon length
    * 
    * @param normlt
    * normalized tendon length
    * @return normalized tendon force
    */
   public abstract double computeft (double normlt);

   /**
    * Computes active force length multiplier derivative for normalized fiber
    * length
    * 
    * @param normlm
    * normalized fiber length
    * @return derivative of active fiber force length multiplier
    */
   public abstract double computeDfactlmDnormlm (double normlm);

   /**
    * Computes force velocity multiplier derivative for normalized fiber
    * velocity
    * 
    * @param normvm
    * normalized fiber velocity
    * @return derivative of force velocity multiplier
    */
   public abstract double computeDfvDnormvm (double normvm);

   /**
    * Computes passive fiber force derivative for normalized fiber length
    * 
    * @param normlm
    * normalize fiber length
    * @return derivative of passive fiber force
    */
   public abstract double computeDfpassDnormlm (double normlm);

   /**
    * Computes tendon force derivative w.r.t normalized tendon length
    * 
    * @param normlm
    * normalize fiber length
    * @return derivative of passive fiber force
    */
   public abstract double computeDftDnormlt (double normlt);

   /**
    * Computes normalized fiber velocity
    * 
    * @param fv
    * fiber force velocity multiplier
    * @return normalized fiber velocity;
    */
   public abstract double computeInvfv (double fv);

   // ************************************************************//
   // *********************HasAuxState interface******************//
   // ************************************************************//

   /**
    * {@inheritDoc}
    */
   @Override
   public void advanceAuxState (double t0, double t1) {

      // nextstateInfo is calculated in previous advance, hence it becomes
      // current state info in this step, copy it. Solver will not use
      // nextStateStep
      getCurrentStateInfo ().copy (getNextStateInfo ());

      // if activation dynamic should be ignored set excitation directory to
      // state, otherwise integrate
      if (isIgnoreActDyn ()) {
         getNextStateInfo ()
            .setActivation (getOwnerMuscle ().getNetExcitation ());
      }
      else {
         double[] in = new double[2];
         in[0] = t0;
         in[1] = getNextStateInfo ().getActivation ();
         double[] out = new double[1];
         adaptiveActStep =
            activationIntegrator
               .adaptiveStepToNewTime (in, out, t1, adaptiveActStep);
         getNextStateInfo ().setActivation (out[0]);
      }

      // fiber length dynamic integration, if tendon is rigid, or solver should
      // be used instead of integrator this part should be skipped

      if (!isIgnoreTendonComp () & !isUseSolver ()) {
         double[] in = new double[2];
         in[0] = t0;
         in[1] = getNextStateInfo ().getFiberLength ();
         double[] out = new double[1];
         adaptiveFiberStep =
            fiberLengthIntegrator
               .adaptiveStepToNewTime (in, out, t1, adaptiveFiberStep);

         if (!Double.isNaN (out[0])) {
            getNextStateInfo ().setFiberLength (out[0]);
         }
         else {
            String muscleName = this.getOwnerMuscle ().getName ();
            equilibriateMaterial ();
            System.out.println (
               "Warning: Fiber length Integration for muscle material "
               + muscleName + " failed. Tring to Continue"
               + " with equilibriation");
         }

      }

      // set time
      getNextStateInfo ().setTime (t1);
   }

   /**
    * {@inheritDoc} <br/>
    * Advances the data buffer to skip over muscle state info. We are storing
    * state data in an object so here, we skip object buffer by 1
    */
   @Override
   public void skipAuxState (DataBuffer data) {
      data.oskip (1);
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public void getAuxState (DataBuffer data) {
      data.oput (getNextStateInfo ());
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public void getInitialAuxState (DataBuffer newData, DataBuffer oldData) {

      if (oldData != null) {
         newData.oput (oldData.oget ());
      }
      else {
         newData.oput (getNextStateInfo ());
      }
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public void setAuxState (DataBuffer data) {
      Object o = data.oget ();
      if (o instanceof HillTypeMuscleState) {
         setNextStateInfo ((HillTypeMuscleState)data.oget ());
      }
   }
}
