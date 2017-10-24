package artisynth.models.hilltypemuscle;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;
import org.python.google.common.base.CharMatcher;
import artisynth.core.mechmodels.Particle;
import artisynth.core.modelbase.ControllerBase;
import artisynth.core.util.ArtisynthPath;
import maspack.matrix.Point3d;
import maspack.matrix.Vector3d;

public class MyMuscleController extends ControllerBase {

   private Particle myParticle;
   private Point3d myPos0;
   private Vector3d myVel0;
   private static final double dataTimeStep = 0.0001;

   private double[][] lengthData;
   private double[][] velocityData;
   private int[] lengthRowNum = new int[1];
   private int[] velocityRowNum = new int[1];

   public MyMuscleController (Particle p) {
      myParticle = p;
      myPos0 = new Point3d (p.getPosition());
      myVel0 = new Vector3d(p.getVelocity()); 
            
      lengthData = readData("Trial6FreeEndDataLength.txt",lengthRowNum);
      velocityData = readData("Trial6FreeEndDataVelocity.txt",velocityRowNum);      
   }

   @Override
   public void apply (double t0, double t1) {
      Point3d pos = new Point3d (myPos0);
      Vector3d vel = new Vector3d(myVel0);
      
      int idx = (int)((t1 * (1/dataTimeStep))+0.000001);
      if (idx >= lengthRowNum[0]) {
         pos.x = lengthData[lengthRowNum[0] - 1][1];
         vel.x = velocityData[velocityRowNum[0] - 1][1];
      }
      else {
         pos.x = lengthData[idx][1];
         vel.x = velocityData[idx][1];
      }

      myParticle.setPosition (pos);
      myParticle.setVelocity(vel);
   }
   

   /**
    * A private utility method to read data ( experimental length and velocity data of the muscle)
    * and return it in a 2D array.
    * @param fileName file name to read data from
    * @param rowNum number of rows (data)
    * @return returns 2D array of data (each row consist of the {time, data at time })
    */
   private double[][] readData(String fileName, int[] rowNum){
      
      
      double[][] storage = null;
      try {
         String fullFileName =
            ArtisynthPath.getSrcRelativePath (
               RBDeGrooteMuscle.class, fileName );
         Scanner scanner = new Scanner (new File (fullFileName));
         scanner.useDelimiter ("\n");
         ArrayList<String> lines = new ArrayList<String> ();

         while (scanner.hasNext ()) {
            lines.add (scanner.next ());
         }
         scanner.close ();

         rowNum[0] = lines.size ();
         int ncol = CharMatcher.is (',').countIn (lines.get (0)) + 1;
         storage = new double[rowNum[0]][ncol];

         for (int i = 0; i < rowNum[0]; i++) {
            String st = lines.get (i);
            int ind = st.indexOf (',');
            int j = 0;
            while (ind > 0) {
               storage[i][j++] = Double.parseDouble (st.substring (0, ind));
               st = st.substring (ind + 1);
               ind = st.indexOf (',');
            }
            storage[i][j] = Double.parseDouble (st);
         }
      }
      catch (FileNotFoundException e) {
         e.printStackTrace ();
      }
      return storage;
   }

}
