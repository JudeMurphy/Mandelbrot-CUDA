import java.io.*;
import java.util.Scanner;
import edu.rit.color.HSB;
import edu.rit.image.PJGColorImage;
import edu.rit.image.PJGImage;

class imgview
{
  static int height, width;
  static int matrix[][];
  static PJGColorImage image;
  static File outputFileName;
	
  public static void main(String args[]) throws Exception
	{
    
    height = Integer.parseInt(args[0]);
    width = Integer.parseInt(args[1]);
    
    matrix = new int[height][width];
    image = new PJGColorImage(height, width, matrix);
    
    outputFileName = new File(args[3]);
    Scanner in = new Scanner(new FileReader(args[2]));
    
    for (int i=0; i<height; i++)
      for (int j=0; j<width; j++)
      {
        float r = in.nextFloat();
        float g = in.nextFloat();
        float b = in.nextFloat();
        matrix[i][j] = HSB.pack(r, g, b);
      }
    
    PJGImage.Writer writer = image.prepareToWrite(new BufferedOutputStream(new FileOutputStream(outputFileName)));
    writer.write();
    writer.close();
	}
}