package artisynth.models.hilltypemuscle.material;

/**
 * This class stores Hill-type muscle information and do all the required
 * calculations prior to a force or stiffness or damping calculations. The idea
 * for creating this class is to hide internal calculation from the Hill-type
 * material, so material implementation is only focused on interacting and
 * providing force and stiffness for the rest of the system. </br>
 * Most of the methods in this class are directly translated from Millard muscle
 * model in OpenSim
 */
class HillTypeMuscleInfo {

   // reference to Hill-type material associated with this class
   private AxialHillTypeMuscleMaterial material;

   // Constants
   private double Epsilon;
   private double maxPenAngle;
   private double minFiberLength;
   private double minFiberLengthAT;
   private double H;
   private double minActivation;

   // length info
   private double tendonLength = Double.NaN;
   double fiberLength, normFiberLength, normTendonLength;

   // alpha is pennation angle
   private double alpha, cosAlpha, sinAlpha;

   // velocity info
   private double tendonVelocity = Double.NaN;
   private double fiberVelocity, normFiberVelocity, alphaVelocity;
   private double lMultiplier, vMultiplier;
   private double fiberForce;
   private double tendonForce;

   private double tendonStiffness, fiberStiffness, fiberStiffnessAT,
   muscleStiffness;
   private double muscleDamping;

   /**
    * Initializes min and max values for muscle parameters. These are used to
    * prevent muscle from getting too close to its singularity
    */
   HillTypeMuscleInfo (AxialHillTypeMuscleMaterial mat) {
      material = mat;
      Epsilon = Double.MIN_VALUE;
      maxPenAngle = Math.acos (0.1); // value is chosen based on Millard

      H = material.getOptLength () * Math.sin (material.getMyOptPenAngle ());
      minFiberLength =
         (maxPenAngle > Epsilon) ? (H / Math.sin (maxPenAngle))
            : (material.getOptLength () * 0.01);
      minFiberLengthAT = minFiberLength * Math.cos (maxPenAngle);
      minActivation = Epsilon;
   }

   /**
    * length related variables calculation
    * 
    * @param muscleLength
    */
   public void calcMuscleLengthInfo (double muscleLength) {
      if (material.isIgnoreTendonComp ()) {
         fiberLength =
            calcFiberLength (muscleLength, material.getMyTendonSlackLen ());
         normFiberLength = fiberLength / material.getOptLength ();
      }
      else {
         fiberLength = material.getCurrentStateInfo ().getFiberLength ();
         normFiberLength = fiberLength / material.getOptLength ();
      }

      alpha = calcPenAngle (fiberLength);
      sinAlpha = Math.sin (alpha);
      setCosAlpha (Math.cos (alpha));
      setTendonLength (muscleLength - fiberLength * getCosAlpha ());
      setNormTendonLength (
         getTendonLength () / material.getMyTendonSlackLen ());

   }

   /**
    * velocity related variables calculations
    * 
    * @param muscleVelocity
    */
   public void calcMuscleVelocityInfo (double muscleVelocity) {

      if (material.isIgnoreTendonComp ()) {

         setFiberVelocity (
            calcFiberVelocity (getCosAlpha (), muscleVelocity, 0));
         normFiberVelocity =
            getFiberVelocity ()
            / (material.getOptLength () * material.getMyMaxFiberVel ());

      }
      else {
         double activation =
            material.isIgnoreActDyn ()
               ? material.getOwnerMuscle ().getNetExcitation ()
               : material.getCurrentStateInfo ().getActivation ();
         double a = clampActivation (activation);

         if (material.getDamping () == 0) {
            // elastic tendon no damping

            if (getCosAlpha () < Epsilon) {
               throw new RuntimeException ("in calcMuscleLengthInfo, pennationa"
               + " angle is 90, causing sinularity in HillType muscle model");
            }

            double factl = material.computefactl (normFiberLength);
            double fpassl = material.computefpass (normFiberLength);
            double ft = material.computeft (getNormTendonLength ());
            if (factl < Epsilon) {
               throw new RuntimeException (
                  "in calcMuscleLengthInfo,active fiber"
                  + " force length is 0, causing sinularity in HillType muscle model");
            }
            double fv =
               Math.max (0, (ft / getCosAlpha () - fpassl) / (a * factl));
            normFiberVelocity = material.computeInvfv (fv);
            setFiberVelocity (
               normFiberVelocity * material.getOptLength ()
               * material.getMyMaxFiberVel ());
         }
         else {
            // elastic tendon with damping

            setFiberVelocity (-1);
            normFiberVelocity = -1;
            if (material.getDamping () <= -Epsilon) {
               throw new RuntimeException (
                  "in calcMuscleLengthInfo,Damping"
                  + " should be positive, causing sinularity in HillType muscle model");
            }
            double fiberVelocityV[] = calcDampedNormFiberVelocity (a);

            // If the Newton method converged, update the fiber velocity.
            if (fiberVelocityV[2] > 0.5) { // flag is set to 0.0 or 1.0
               normFiberVelocity = fiberVelocityV[0];
               setFiberVelocity (
                  normFiberVelocity * material.getOptLength ()
                  * material.getMyMaxFiberVel ());
            }
            else {
               // Throw an exception here because there is no point integrating
               // a muscle velocity that is invalid (it will end up producing
               // invalid fiber lengths and will ultimately cause numerical
               // problems). The idea is to produce an exception and catch this
               // early before it can cause more damage.
               throw new RuntimeException (
                  " Fiber velocity Newton method did not converge");
            }
         } // end of damping vs no damping
      } // end of rigid vs elastic

      // Compute the other velocity-related components.
      if (fiberLength < Epsilon) {
         throw new RuntimeException (
            "in calcMuscleLengthInfo, fiber length can not be zero");
      }
      alphaVelocity =
         (material.getMyOptPenAngle () > Epsilon)
            ? (-(getFiberVelocity () / fiberLength) * Math.tan (alpha)) : 0;

      setTendonVelocity (0);
      if (!material.isIgnoreTendonComp ()) {
         setTendonVelocity (
            muscleVelocity - getFiberVelocity () * getCosAlpha ()
            + fiberLength * sinAlpha * alphaVelocity);
      }
   }

   /**
    * dynamic related variables ( force, stiffness, damping)
    * 
    * @param muscleLength
    * @param muscleVelocity
    */
   public void calcMuscleDynamicInfo (

      double muscleLength, double muscleVelocity) {
      calcMuscleLengthInfo (muscleLength);
      calcMuscleVelocityInfo (muscleVelocity);

      double activation =
         material.isIgnoreActDyn ()
            ? material.getOwnerMuscle ().getNetExcitation ()
            : material.getCurrentStateInfo ().getActivation ();
      double a = clampActivation (activation);

      lMultiplier = material.computefactl (normFiberLength);
      vMultiplier = material.computefv (normFiberVelocity);
      setTendonForce (
         material.getMaxForce () * material.computeft (getNormTendonLength ()));
      double fa = material.getMaxForce () * (a * lMultiplier * vMultiplier);
      double fp1 =
         material.getMaxForce () * material.computefpass (normFiberLength);
      double fp2 =
         material.getMaxForce () * material.getDamping () * normFiberVelocity;
      setFiberForce (fa + (fp1 + fp2));

      fiberStiffness =
         calcFiberStiffness (
            material.getMaxForce (), a, material.computefv (normFiberVelocity),
            normFiberLength, material.getOptLength ());

      double dFmAT_dlm =
         calc_DFiberForceAT_DFiberLength (
            getFiberForce (), fiberStiffness, fiberLength, sinAlpha,
            getCosAlpha ());

      // dFmAT_dlmAT
      fiberStiffnessAT =
         calc_DFiberForceAT_DFiberLengthAT (
            dFmAT_dlm, sinAlpha, getCosAlpha (), fiberLength);

      if (!material.isIgnoreTendonComp ()) {
         tendonStiffness =
            material.computeDftDnormlt (getNormTendonLength ())
            * (material.getMaxForce () / material.getMyTendonSlackLen ());

         // Compute the stiffness of the whole musculotendon actuator.
         if (Math.abs (fiberStiffnessAT * tendonStiffness) > 0.0
         && Math.abs (fiberStiffnessAT + tendonStiffness) > Epsilon) {
            setMuscleStiffness (
               (fiberStiffnessAT * tendonStiffness)
               / (fiberStiffnessAT + tendonStiffness));
         }
      }
      else {
         tendonStiffness = Double.POSITIVE_INFINITY;
         setMuscleStiffness (fiberStiffnessAT);
      }

      double fact = material.computefactl (normFiberLength);
      double dfvdnormvm = material.computeDfvDnormvm (normFiberVelocity);
      setMuscleDamping (
         getCosAlpha () * getCosAlpha () * a * material.getMaxForce ()
         * (fact * dfvdnormvm + material.getDamping ())
         / (material.getMyMaxFiberVel () * material.getOptLength ()));

   }

   public void equilibriate () {
      if (material.isIgnoreTendonComp ()) { // rigid tendon
         return;
      }

      // Elastic tendon initialization routine.
      try {

         // Compute the fiber length where the fiber and tendon are in static
         // equilibrium. Fiber and tendon velocity are set to zero.
         int flag_status = -1;
         double fiberLength = Double.NaN;


         // tol is the desired tolerance in Newtons.
         double tol = 1e-3 * material.getMaxForce ();
         if (tol < Double.MIN_VALUE * 10) {
            tol = Double.MIN_VALUE * 10;
         }
         int maxIter = 400;
         double pathLength = material.getOwnerMuscle ().getLength ();
         double pathLengtheningSpeed =
            material.getOwnerMuscle ().getLengthDot ();
         double act =
            material.isIgnoreActDyn ()
               ? material.getOwnerMuscle ().getNetExcitation ()
               : material.getCurrentStateInfo ().getActivation ();
         double[] soln;
         soln =
            estimateMuscleFiberState (
               clampActivation (act), pathLength, pathLengtheningSpeed, tol,
               maxIter, false);
         flag_status = (int)soln[0];
         fiberLength = soln[3];
         fiberVelocity = soln[4];
         tendonForce = soln[5];

         switch (flag_status) {
            case 0: // converged
            {
               material.getCurrentStateInfo ().setFiberLength (fiberLength);
            }
               break;

            case 1: // lower bound on fiber length was reached
            {
               material.getCurrentStateInfo ().setFiberLength (fiberLength);
               String msg = "lower bound on fiber length was reached";
               System.out.println (msg);
            }
               break;

            case 2: // maximum number of iterations reached
            {
               material.getCurrentStateInfo ().setFiberLength (fiberLength);

               String msg = "maximum number of iterations reached";
               System.out.println (msg);
            }
               break;

            default:
               String msg =
                  "\n\nWARNING: invalid error flag returned from "
                  + "static solution for ."
                  + "Setting tendon force to 0.0 and fiber length to the "
                  + "optimal fiber length.";
               System.out.println (msg);
               material.getCurrentStateInfo ().setFiberLength (
                  material.getOptLength ());
         }

      }
      catch (Exception e) {
         // If the initialization routine fails in some unexpected way, tell the
         // user and continue with some valid initial conditions.

         material
            .getCurrentStateInfo ().setFiberLength (material.getOptLength ());
      }
   }

   private double[] estimateMuscleFiberState (
      double aActivation, double pathLength, double pathLengtheningSpeed,
      double aSolTolerance, int aMaxIterations, boolean staticSolution) {

      // Results vector format:
      // [0] flag: 0 = converged
      // 1 = diverged
      // 2 = no solution due to length singularity
      // 3 = no solution due to pennation angle singularity
      // [1] solution error (N)
      // [2] iterations
      // [3] fiber length (m)
      // [4] fiber velocity (N)
      // [5] tendon force (N)
      double[] results = new double[6];

      // Using short variable names to facilitate writing out long equations

      double ma = aActivation;
      double ml = pathLength;
      double dml = pathLengtheningSpeed;
      double tsl = material.getMyTendonSlackLen ();
      double ofl = material.getOptLength ();
      double fiso = material.getMaxForce ();
      double vmax = material.getMyMaxFiberVel ();

      double fse = 0.0; // normalized tendon (series element) force
      double fal = 0.0; // normalized active-force-length multiplier
      double fv = 0.0; // normalized force-velocity multiplier
      double fpe = 0.0; // normalized parallel element force

      // Position level
      double lce = 0.0;
      double tl =
         (staticSolution | Double.isNaN (getTendonLength ())) ? tsl * 1.01
            : getTendonLength (); // begin with small tendon


      lce = clampFiberLength (calcFiberLength (ml, tl));

      double phi = calcPenAngle (lce);
      double cosphi = Math.cos (phi);
      double sinphi = Math.sin (phi);
      double tlN = tl / tsl;
      double lceN = lce / ofl;

      double dtl = 0;// Double.isNaN (getTendonVelocity()) ? 0 :
                     // getTendonVelocity();

      double dlce =
         (staticSolution) ? 0.0 : calcFiberVelocity (cosphi, dml, dtl);
      double dlceN = (staticSolution) ? 0.0 : dlce / (vmax * ofl);

      // Internal variables for the loop

      double Fm = 0.0; // fiber force
      double FmAT = 0.0; // fiber force along tendon
      double Ft = 0.0; // tendon force
      double ferr = 10.0; // solution error
      double dFm_dlce = 0.0; // partial of muscle force w.r.t. lce
      double dFmAT_dlce = 0.0; // partial of muscle force along tl w.r.t. lce
      double dFmAT_dlceAT = 0.0; // partial of muscle force along tl w.r.t. lce
                                 // along the tendon
      double dFt_d_lce = 0.0; // partial of tendon force w.r.t. lce
      double dFt_d_tl = 0.0; // partial of tendon force w.r.t. tl
      double dferr_d_lce = 0.0; // partial of solution error w.r.t lce
      double delta_lce = 0.0; // change in lce
      double Ke = 0.0; // linearized local stiffness of the muscle

      // Initialize the loop
      int iter = 0;
      double[] fiberForceV = new double[4];
      // System.out.print(material.currentStateInfo.getTime()+ " " + tl + " " );
      while (Math.abs (ferr) > aSolTolerance && iter < aMaxIterations) {
         // Update the multipliers and their partial derivatives
         fal = material.computefactl (lceN);
         fpe = material.computefpass (lceN);
         fse = material.computeft (tlN);
         fv = (staticSolution) ? 1.0 : material.computefv (dlceN);

         // Compute the force error
         fiberForceV = calcFiberForce (fiso, ma, fal, fv, fpe, dlceN);
         Fm = fiberForceV[0];
         FmAT = Fm * cosphi;
         Ft = fse * fiso;
         ferr = FmAT - Ft;

         // Compute the partial derivative of the force error w.r.t. lce
         dFm_dlce = calcFiberStiffness (fiso, ma, fv, lceN, ofl);
         dFmAT_dlce =
            calc_DFiberForceAT_DFiberLength (Fm, dFm_dlce, lce, sinphi, cosphi);
         dFmAT_dlceAT =
            calc_DFiberForceAT_DFiberLengthAT (dFmAT_dlce, sinphi, cosphi, lce);
         dFt_d_tl = material.computeDftDnormlt (tlN) * fiso / tsl;
         dFt_d_lce =
            calc_DTendonForce_DFiberLength (dFt_d_tl, lce, sinphi, cosphi);

         // Error derivative
         dferr_d_lce = dFmAT_dlce - dFt_d_lce;

         if (Math.abs (ferr) > aSolTolerance) {
            if (Math.abs (dferr_d_lce) > Epsilon) {
               // Take a full Newton Step if the derivative is nonzero
               delta_lce = -ferr / dferr_d_lce;
               lce = lce + delta_lce;
            }
            else {
               // We've stagnated; perturb the current solution
               double perturbation = 2.0 * (Math.random ()) - 1.0;

               double lengthPerturbation =
                  0.5 * perturbation * material.getOptLength ();
               lce += lengthPerturbation;
            }

            // Update position level quantities only if they won't go singular
            phi = calcPenAngle (lce);
            sinphi = Math.sin (phi);
            cosphi = Math.cos (phi);
            tl = ml - lce * cosphi;// calcTendonLength(cosphi,lce,ml);
            lceN = lce / ofl;
            tlN = tl / tsl;

            /*
             * Update velocity-level quantities. Share the muscle velocity
             * between the tendon and the fiber according to their relative
             * stiffnesses:
             * 
             * Fm-Ft = 0 Equilibrium equation [1] d/dt Fm - d/dt Ft = 0 Time
             * derivative [2] lp = lm + lt Path definition [3] d/dt lp = d/dt lm
             * + d/dt lt Path derivative [4]
             * 
             * Computing a linearized model of [2]: Fm = Fm0 + Km*lceAT [5] Ft =
             * Ft0 Kt*xt [6]
             * 
             * Taking its time derivative: dFm_d_xm = Km*dlceAT + dKm_d_t*lceAT
             * (assume dKm_d_t = 0) [7] dFt_d_xt = Kt*dtl + dKt_d_t*dtl (assume
             * dKt_d_t = 0) [8]
             * 
             * Subtituting 7 and 8 into 2: Km dlceAT - Kt dtl = 0
             * 
             * Using Eqn 4, we have 2 equations in 2 unknowns. Can now solve for
             * tendon velocity, or the velocity of the fiber along the tendon.
             * 
             * This is a heuristic. The above assumptions are necessary since
             * computing the partial derivatives of Km or Kt requires
             * acceleration- level knowledge, which is not available in general.
             * 
             * Stiffness of the muscle is the stiffness of the tendon and the
             * fiber (along the tendon) in series.
             * 
             * The "if" statement here is to handle the special case where the
             * negative stiffness of the fiber (which happens in this model) is
             * equal to the positive stiffness of the tendon.
             */
            if (!staticSolution) {
               // Keep velocities set to zero if seeking a static solution.
               if (Math.abs (dFmAT_dlceAT + dFt_d_tl) > Epsilon && tlN > 1.0) {
                  Ke = (dFmAT_dlceAT * dFt_d_tl) / (dFmAT_dlceAT + dFt_d_tl);
                  dtl = (1 / dFt_d_tl) * Ke * dml;

               }
               else {
                  dtl = dml;
               }
               dlce = calcFiberVelocity (cosphi, dml, dtl);
               dlceN = dlce / (vmax * ofl);
            }
         }

         iter++;
      }

      // System.out.println(" ferr " +ferr);
      // Populate the results vector:
      // [0] flag: 0 = converged
      // 1 = diverged
      // 2 = no solution due to length singularity
      // 3 = no solution due to pennation angle singularity
      // [1] solution error (N)
      // [2] iterations
      // [3] fiber length (m)
      // [4] fiber velocity (N)
      // [5] tendon force (N)

      if (Math.abs (ferr) < aSolTolerance) {
         // The solution converged
         results[0] = 0;
         results[1] = ferr;
         results[2] = (double)iter;
         results[3] = lce;
         results[4] = dlce;
         results[5] = fse * fiso;

      }
      else {
         // The fiber length hit its lower bound
         if (iter < aMaxIterations) {
            lce = minFiberLength;
            phi = calcPenAngle (lce);
            cosphi = Math.cos (phi);
            tl = ml - lce * cosphi; // calcTendonLength(cosphi,lce,ml);
            lceN = lce / ofl;
            tlN = tl / tsl;
            fse = material.computeft (tlN);

            results[0] = 1.0;
            results[1] = ferr;
            results[2] = (double)iter;
            results[3] = lce;
            results[4] = 0;
            results[5] = fse * fiso;

         }
         else {
            // The solution diverged
            results[0] = 2.0;
            results[1] = ferr;
            results[2] = (double)iter;
            results[3] = lce;
            results[4] = Double.NaN;
            results[5] = Double.NaN;
         }
      }
      return results;
   }

   private double[] calcDampedNormFiberVelocity (double a) {

      double[] fiberForceV = new double[4];
      double[] result = new double[3];

      double fiso = material.getMaxForce ();
      double fal = material.computefactl (normFiberLength);
      double fpe = material.computefpass (normFiberLength);
      double fse = material.computeft (getNormTendonLength ());
      double beta = material.getDamping ();
      double cosPhi = getCosAlpha ();

      int maxIter = 20;
      double tol = 1.0e-10 * fiso;
      if (tol < Epsilon * 100) {
         tol = Epsilon * 100;
      }

      double perturbation = 0.0;
      double fiberForce = 0.0;
      double err = 1.0e10;
      double derr_d_dlceNdt = 0.0;
      double delta = 0.0;
      double iter = 0.0;

      // Get a really excellent starting position to reduce the number of
      // iterations. This reduces the simulation time by about 1%.
      double fv =
         (fse / Math.max (cosPhi, 0.01) - fpe)
         / (Math.max (a, 0.01) * Math.max (fal, 0.01));
      double dlceN_dt = material.computeInvfv (fv);

      // The approximation is poor beyond the maximum velocities.
      if (dlceN_dt > 1.0) {
         dlceN_dt = 1.0;
      }
      if (dlceN_dt < -1.0) {
         dlceN_dt = -1.0;
      }

      double df_d_dlceNdt = 0.0;

      while (Math.abs (err) > tol && iter < maxIter) {
         fv = material.computefv (dlceN_dt);
         fiberForceV = calcFiberForce (fiso, a, fal, fv, fpe, dlceN_dt);
         fiberForce = fiberForceV[0];

         err = fiberForce * cosPhi - fse * fiso;
         df_d_dlceNdt =
            fiso * (a * fal * material.computeDfvDnormvm (dlceN_dt) + beta);
         derr_d_dlceNdt = df_d_dlceNdt * cosPhi;

         if (Math.abs (err) > tol && Math.abs (derr_d_dlceNdt) > Epsilon) {
            delta = -err / derr_d_dlceNdt;
            dlceN_dt = dlceN_dt + delta;

         }
         else if (Math.abs (derr_d_dlceNdt) < Epsilon) {
            // Perturb the solution if we've lost rank. This should never happen
            // for this problem since dfv_d_dlceNdt > 0 and b > 0 (and so
            // derr_d_dlceNdt > 0).
            perturbation = 2.0 * (Math.random ()) - 1.0;
            dlceN_dt = dlceN_dt + perturbation * 0.05;
         }
         iter++;
      }

      double converged = 1.0;

      // If we failed to converge, it's probably because the fiber is at its
      // lower
      // bound. That decision is made further down the line, so if convergence
      // didn't happen, let the user know and return a NaN.
      if (Math.abs (err) > tol) {
         dlceN_dt = -1.0;
         converged = 0.0;
      }
      result[0] = dlceN_dt;
      result[1] = err;
      result[2] = converged;
      return result;
   }

   private double[] calcFiberForce (
      double fiso, double a, double fal, double fv, double fpe, double dlceN) {
      double beta = material.getDamping ();
      double fa = fiso * (a * fal * fv);
      double fp1 = fiso * fpe;
      double fp2 = fiso * beta * dlceN;
      double fm = fa + (fp1 + fp2);

      double[] fiberF = new double[4];
      fiberF[0] = fm;
      fiberF[1] = fa;
      fiberF[2] = fp1; // conservative passive force
      fiberF[3] = fp2;
      return fiberF;
   }

   private double clampActivation (double act) {
      return clamp (minActivation, act, 1);
   }

   private double clamp (double min, double val, double max) {
      val = (val < min) ? min : val;
      val = (val > max) ? max : val;
      return val;
   }

   private double calcFiberVelocity (
      double cosAlpha, double muscleVelocity, double tendonVelocity) {
      return (muscleVelocity - tendonVelocity) * cosAlpha;
   }

   private double calcPenAngle (double lm) {
      double alpha = 0;
      // This computation is only worth performing on pennated muscles.
      if (material.getMyOptPenAngle () > Epsilon) {
         if (lm > minFiberLength) {
            double sinAlpha = H / lm;
            alpha =
               (sinAlpha < Math.sin (maxPenAngle)) ? Math.asin (sinAlpha)
                  : maxPenAngle;
         }
         else {
            alpha = maxPenAngle;
         }
      }
      return alpha;
   }

   private double calcFiberLength (double muscleLength, double tendonLength) {
      double fiberLengthAT = muscleLength - tendonLength;
      double fiberLength = 0.0;

      if (fiberLengthAT >= minFiberLengthAT) {
         fiberLength = Math.sqrt (H * H + fiberLengthAT * fiberLengthAT);
      }
      else {
         fiberLength = minFiberLength;
      }
      return fiberLength;
   }

   private double clampFiberLength (double lce) {
      return Math.max (lce, minFiberLength);
   }

   private double calcFiberStiffness (
      double fiso, double a, double fv, double lceN, double optFibLen) {

      double DlceN_Dlce = 1.0 / optFibLen;
      double Dfal_Dlce = material.computeDfactlmDnormlm (lceN) * DlceN_Dlce;
      double Dfpe_Dlce = material.computeDfpassDnormlm (lceN) * DlceN_Dlce;

      // DFm_Dlce
      return fiso * (a * Dfal_Dlce * fv + Dfpe_Dlce);
   }

   private double calc_DFiberForceAT_DFiberLengthAT (
      double dFmAT_d_lce, double sinPhi, double cosPhi, double lce) {
      double dphi_d_lce = calc_DPennationAngle_DfiberLength (lce);

      double DlceAT_Dlce = cosPhi - lce * sinPhi * dphi_d_lce;

      return dFmAT_d_lce * (1.0 / DlceAT_Dlce);
   }

   private double calc_DFiberForceAT_DFiberLength (
      double fiberForce, double fiberStiffness, double lce, double sinPhi,
      double cosPhi) {
      double Dphi_Dlce = calc_DPennationAngle_DfiberLength (lce);
      double Dcosphi_Dlce = -sinPhi * Dphi_Dlce;

      return fiberStiffness * cosPhi + fiberForce * Dcosphi_Dlce;
   }

   private double calc_DPennationAngle_DfiberLength (double fiberLength) {
      if (fiberLength < H) {
         throw new RuntimeException (
            "in calc_DPennationAngle_DfiberLength,"
            + "Fiber length is below the lower bound for this muscle.");
      }

      double h_over_l = H / fiberLength;
      return (-h_over_l / fiberLength) / Math.sqrt (1.0 - h_over_l * h_over_l);
   }

   private double calc_DTendonForce_DFiberLength (
      double dFt_d_tl, double lce, double sinphi, double cosphi) {
      double dphi_d_lce = calc_DPennationAngle_DfiberLength (lce);
      double dtl_d_lce =
         calc_DTendonLength_DfiberLength (lce, sinphi, cosphi, dphi_d_lce);

      return dFt_d_tl * dtl_d_lce;
   }

   private double calc_DTendonLength_DfiberLength (
      double fiberLength, double sinPennationAngle, double cosPennationAngle,
      double DpennationAngle_DfiberLength) {
      if (fiberLength < H) {
         throw new RuntimeException (
            "in calc_DPennationAngle_DfiberLength,"
            + "Fiber length is below the lower bound for this muscle.");
      }

      return fiberLength * sinPennationAngle * DpennationAngle_DfiberLength
      - cosPennationAngle;
   }

   // **************************************************************//
   // *********** Setter and Getter (internal variables) ***********//
   // **************************************************************//
   public double getTendonLength () {
      return tendonLength;
   }

   public void setTendonLength (double tendonLength) {
      this.tendonLength = tendonLength;
   }

   public double getNormTendonLength () {
      return normTendonLength;
   }

   public void setNormTendonLength (double normTendonLength) {
      this.normTendonLength = normTendonLength;
   }

   public double getFiberForce () {
      return fiberForce;
   }

   public void setFiberForce (double fiberForce) {
      this.fiberForce = fiberForce;
   }

   public double getFiberVelocity () {
      return fiberVelocity;
   }

   public void setFiberVelocity (double fiberVelocity) {
      this.fiberVelocity = fiberVelocity;
   }

   public double getTendonVelocity () {
      return tendonVelocity;
   }

   public void setTendonVelocity (double tendonVelocity) {
      this.tendonVelocity = tendonVelocity;
   }

   public double getTendonForce () {
      return tendonForce;
   }

   public void setTendonForce (double tendonForce) {
      this.tendonForce = tendonForce;
   }

   public double getCosAlpha () {
      return cosAlpha;
   }

   public void setCosAlpha (double cosAlpha) {
      this.cosAlpha = cosAlpha;
   }

   public double getMuscleStiffness () {
      return muscleStiffness;
   }

   public void setMuscleStiffness (double muscleStiffness) {
      this.muscleStiffness = muscleStiffness;
   }

   public double getMuscleDamping () {
      return muscleDamping;
   }

   public void setMuscleDamping (double muscleDamping) {
      this.muscleDamping = muscleDamping;
   }

}