package artisynth.models.hilltypemuscle;

import java.awt.Color;
import maspack.render.*;
import artisynth.core.mechmodels.*;
import artisynth.core.workspace.RootModel;
import artisynth.models.hilltypemuscle.material.AxialDeGrooteMuscleMaterial;

/**
 * Simple demo of a Hill-type muscle. Left end is fixed, and the position and
 * velocity of the right end is controlled by a controller with data read from
 * file (to recreate the experiment).
 */
public class RBDeGrooteMuscle extends RootModel {

   public void build (String[] args) {

      // create MechModel and add to RootModel
      MechModel mech = new MechModel ("mech");
      addModel (mech);

      mech.setGravity (0, 0, 0);
      setMaxStepSize (0.0001);

      Particle p1 = new Particle ("p1", /* mass= */1, /* x,y,z= */0, 0, 0);
      p1.setDynamic (false);

      double muscleRestLength = 0.03382188;
      Particle p2 =
         new Particle ("p2", /* mass= */1, /* x,y,z= */muscleRestLength, 0, 0);
      p2.setDynamic (false);

      // Create muscle
      HillTypeMuscle muscle = new HillTypeMuscle ();
      muscle.setPoints (p1, p2);
      muscle.setName ("testMuscle");
      muscle.setExcitation (1);

      // Create material and set properties
      AxialDeGrooteMuscleMaterial mat;

      // OpenSim test uses different for rigid and non-rigid muscles
      // this if is intended to switch between parameters easily. If your muscle
      // is rigid set boolean variable rigid to be true
      boolean rigid = false;
      if (rigid) {
         mat =
            new AxialDeGrooteMuscleMaterial (
               0.0171 /* optimalFiberLength */, 1.2022943404 /* maximumForce */,
               0.0171 /* tendonSlackLength */,
               0.10471975511966 /* optimalPennationAngle */,
               9.8599297056 /* maxFiberVelocity */, 0.0 /* Damping */,
               0.01 /* default activation */);
      }
      else {
         mat =
            new AxialDeGrooteMuscleMaterial (
               0.0171 /* optimalFiberLength */, 1.2758948506 /* maximumForce */,
               0.0171 /* tendonSlackLength */,
               0/*0.10471975511966  optimalPennationAngle */,
               5.3156293513 /* maxFiberVelocity */, 0.0 /* Damping */,
               0.01 /* default activation */);
      }
      mat.setIgnoreTendonComp (rigid);
      mat.setIgnoreActDyn(true);

      // this is a solver based muscle implementation. It uses solver at each
      // step instead of integration.
      // Its performance is not good at current we will use integrator
      mat.setUseSolver (true);

      // You should set all the material properties before this call.
      // Material is copied in PointSpringBase class, hence any change to
      // material after this call will not affect the actual material of the
      // muscle
      muscle.setMaterial (mat);

      // add components to the mech model
      addController (new MyMuscleController (p2));
      mech.addParticle (p1);
      mech.addParticle (p2);
      mech.addAxialSpring (muscle);

      // set render properties for the components
      RenderProps.setSphericalPoints (p1, 0.001, Color.RED);
      RenderProps.setSphericalPoints (p2, 0.001, Color.RED);
      RenderProps.setCylindricalLines (muscle, 0.001, Color.BLUE);
   }
}
