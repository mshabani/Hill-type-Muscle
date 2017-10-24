package artisynth.models.hilltypemuscle;

import artisynth.core.mechmodels.Muscle;
import artisynth.models.hilltypemuscle.material.AxialHillTypeMuscleMaterial;
import maspack.matrix.Vector3d;
import maspack.render.Renderer;

public class HillTypeMuscle extends Muscle {

   /**
    * {@inheritDoc}
    */
   @Override
   public void initialize (double t) {

      // Set this muscle to be the owner of the material.
      // This way material would have access to muscle, length and velocity at
      // any place
      if (myAuxStateMat instanceof AxialHillTypeMuscleMaterial) {
         ((AxialHillTypeMuscleMaterial)myAuxStateMat).setOwner (this);
      }

      // This initialization is required to allocate internal variables.
      // Also, this calls material equilibriation, which chooses fiber length in
      // a way that fiber force along tendon and tendon force are equal, this
      // initial fiber length will be used as initial value for fiber length
      // integration
      AxialHillTypeMuscleMaterial mat =
         (AxialHillTypeMuscleMaterial)getEffectiveMaterial ();
      mat.initMaterialFromProperties (t);
   }

   /**
    * {@inheritDoc}
    */
   @Override
   public void render (Renderer renderer, int flags) {

      float[] midPoint = new float[3];
      getTendonPointForRenderer (midPoint);
      float[] myTendonColor = new float[] { 1, 1, 1 };

      renderer.drawLine (
         myRenderProps, myPnt0.myRenderCoords, midPoint, myTendonColor,
         /* capped= */false, isSelected ());
      renderer.drawLine (myRenderProps, midPoint, myPnt1.myRenderCoords,
         myRenderColor, /* capped= */false, isSelected ());
   }

   /**
    * This functions assumes pnt0 to be the start point of the tendon and
    * calculates the end point based on current tendon and muscle length and its
    * direction.
    * 
    * @param pointCoordinate
    * an array to store tendon end point coordinates
    */
   private void getTendonPointForRenderer (float[] pointCoordinate) {

      AxialHillTypeMuscleMaterial mat =
         (AxialHillTypeMuscleMaterial)getEffectiveMaterial ();
      float tendonLength = (float)mat.getTendonLength ();

      updateU();
      float[] p0 = myPnt0.myRenderCoords;
      Vector3d dir = new Vector3d (myU);

      pointCoordinate[0] = p0[0] + tendonLength * (float)dir.get (0);
      pointCoordinate[1] = p0[1] + tendonLength * (float)dir.get (1);
      pointCoordinate[2] = p0[2] + tendonLength * (float)dir.get (2);
   }
}
