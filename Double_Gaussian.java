import ij.*;
import ij.process.*;
import ij.gui.*;
import static ij.measure.Measurements.MEAN;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import java.util.ArrayList;
//import java.util.Arrays;

public class Double_Gaussian implements PlugIn 
{
    static Point loc = null;
    ImagePlus im;
    ImagePlus imp;
    ImagePlus filter=null;
    ImageProcessor ip;
    ImagePlus result;
    ImageStack is;
    ImageStack is2;
    ImagePlus imp2;
    FHT fht;
    ImageProcessor PS;
    double sigma1=0.4;
    double sigma2=0.2;
    double sig1=0.4;
    double sig2=0.2;
    double weight=0.3;
    boolean show_filter=false;
	boolean withFFT=false;
	boolean bnorm = false;
    String imtitle;
    
    private boolean padded;
    private int originalWidth;
    private int originalHeight;
    private int stackSize = 1;
    private int slice = 1;
    private boolean doFFT;
    private int maxN=2;

    float[][] kernel;
    float[] kernel1d;

    public void run(String arg) 
    {
        im= IJ.getImage();
        pad(im.getProcessor());
        imtitle=im.getTitle();
        int frNb;
        if(im.getNFrames()!=0)
            frNb=im.getNFrames();
        else
            frNb=1;
        int slNb;
        if(im.getNSlices()!=0)
            slNb=im.getNSlices();
        else
            slNb=1;
        is = new ImageStack(im.getWidth(),im.getHeight(),slNb*frNb);
        
        result=new ImagePlus();
        result.setCalibration(im.getCalibration());
        result.setFileInfo(im.getFileInfo());
        imp2= new ImagePlus();
        result.setTitle("Double-Gaussian-Filter-of-"+imtitle);
        while(doDialog())
        {
            imp=WindowManager.getImage(imtitle);
            get_kernel(sig1,sig2,weight);
            int currentSlice = imp.getCurrentSlice();
            for (int i=1; i<=slNb*frNb; i++)
            {
                do_Process(i);
                if (i == currentSlice && show_filter)
                {
                    byte[] PSPixels = (byte[])PS.getPixels();
                    float[] resPS = new float[PSPixels.length];
                    for (int k=0; k<PSPixels.length; k++)
                    {   
                        short PSPixel = PSPixels[k];
                        if (PSPixel < 0)
                            PSPixel += 256;
                        if (!withFFT)
                        	PSPixel=1;
                        resPS[k] = PSPixel*kernel1d[k];
                    }
                    ImageStack fSt;
                    filter = WindowManager.getImage("DG Filter");
                    if (filter==null)
                    {
                        filter = new ImagePlus();
                        fSt = new ImageStack(maxN, maxN, 1);
                    }
                    else
                    {
                        fSt = filter.getImageStack();
                    }
                    fSt.setPixels(resPS, 1);
                    filter.setStack(fSt);
                    filter.setTitle("DG Filter");
                    if (WindowManager.getImage("DG Filter") == null)
                        filter.show();
                    filter.updateAndRepaintWindow();
                    
                    
                }
            }
            
            result.setStack(is);
            result.show();
			/* result.setStack(is);
			if (WindowManager.getImage("Double Gaussian Filter of "+imtitle)==null)
			{
            result.show();
			}
			else
				result.updateAndRepaintWindow();
			} */
        }
    }
    
    void do_Process(int i)
    {
        
        imp.setSliceWithoutUpdate(i);
        ip= imp.getProcessor();
        
        fht=newFHT(ip);
        is2= new ImageStack(maxN,maxN,1);
        apply_filter(kernel, i);
        //is=is2.crop(0, 0,0, originalWidth-1, originalHeight-1,1);
        is.setPixels(crop((float[])is2.getPixels(1),maxN,originalWidth,originalHeight), i); 
    }
    
    public boolean doDialog()
    {
        
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog(getClass().getSimpleName());
        gd.setSmartRecording(true);
        if(loc != null)
        {
            gd.centerDialog(false);
            gd.setLocation(loc.x+10,loc.y+10); //compensate the drift
        }
        loc = gd.getLocation();
        
        gd.addNumericField("sigma1 (large)", sigma1, 2);
        gd.addNumericField("sigma2 (small)", sigma2, 2);
        gd.addNumericField("weight", weight, 2);
        gd.addCheckbox("show filter", show_filter);
        gd.addCheckbox("add FFT?", withFFT);
	gd.addCheckbox("keep mean?", bnorm);
        gd.showDialog();
        
        sigma1=gd.getNextNumber();
        sigma2=gd.getNextNumber();
        weight=gd.getNextNumber();
        show_filter=gd.getNextBoolean();
        withFFT=gd.getNextBoolean();
	bnorm=gd.getNextBoolean();
        sig1=sigma1*sigma1;
        sig2=sigma2*sigma2;
        if(gd.wasCanceled())
        {
            return false;
        }
        return true;
        
    }
    public void get_kernel(double sigma1, double sigma2, double weight)
    {
        int size=maxN;
        kernel=new float[size][size];
        kernel1d=new float[size*size];
	float sum = 0f;
	float val;
        for(int row=0;row<size;row++)
        {
            for(int col=0;col<size;col++)
            {
                double dist=Math.sqrt((col-size/2)*(col-size/2)+(row-size/2)*(row-size/2))/(size/2);
                if (sigma1==0)
		{
			val = (float)(1-(1.0-weight)*Math.exp(-0.5*dist*dist/(sigma2*sigma2)));
		}
                else
                {
			val = (float)(Math.exp(-0.5*dist*dist/(sigma1*sigma1))-(1.0-weight)*Math.exp(-0.5*dist*dist/(sigma2*sigma2)));
                }
		kernel[row][col] = val;
                kernel1d[row*size+col] = val;
                sum+=val;
            }
            
        }
	for(int row=0;row<size;row++)
        {
            for(int col=0;col<size;col++)
            {
		kernel[row][col]/=sum;
		kernel1d[row*size+col]/=sum;   
            }
            
        }
	
	    
	if (bnorm)
	{
        	kernel[size/2][size/2]=1;
        	kernel1d[size*size/2]=1;    
	}
    }
    
    public void apply_filter(float[][] kernel, int slice_number)
    {
        
        fht.transform();
        PS = fht.getPowerSpectrum();
        fht.swapQuadrants();
        fht=fht.multiply(new FHT(new FloatProcessor(kernel)));      
        fht.swapQuadrants();
        fht.inverseTransform();
        is2.setPixels((float[])fht.getPixels(), 1);
        //is.setPixels((float[])fht.getPixels(), slice_number);       

    }
         
    
    ImageProcessor pad(ImageProcessor ip) {
        originalWidth = ip.getWidth();
        originalHeight = ip.getHeight();
        maxN = Math.max(originalWidth, originalHeight);
        int i = 2;
        while(i<maxN) i *= 2;
        if (i==maxN && originalWidth==originalHeight) {
            padded = false;
            return ip;
        }
        maxN = i;
        //showStatus("Padding to "+ maxN + "x" + maxN);
        ImageStatistics stats = ImageStatistics.getStatistics(ip, MEAN, null);
        ImageProcessor ip2 = ip.createProcessor(maxN, maxN);
        ip2.setValue(stats.mean);
        ip2.fill();
        ip2.insert(ip, 0, 0);
        padded = true;
        Undo.reset();
        //new ImagePlus("padded", ip2.duplicate()).show();
        return ip2;
    }
    
 
    
    FHT newFHT(ImageProcessor ip) 
    {
        FHT fht;
        if (ip instanceof ColorProcessor) 
        {
            //showStatus("Extracting brightness");
            ImageProcessor ip2 = ((ColorProcessor)ip).getBrightness();
            fht = new FHT(pad(ip2));
            fht.rgb = (ColorProcessor)ip.duplicate(); // save so we can later update the brightness
        } 
        else
            fht = new FHT(pad(ip));
        if (padded) {
            fht.originalWidth = originalWidth;
            fht.originalHeight = originalHeight;
        }
        int bitDepth = imp.getBitDepth();
        fht.originalBitDepth = bitDepth;
        if (bitDepth!=24)
        	fht.originalColorModel = ip.getColorModel();
        return fht;
    }
    
    float[] crop(float[] ar1, int w1, int w2, int h2)
    {
        float[] ar2= new float[w2*h2];
        //int h1=ar1.length/w1;
        int c=0;
        for(int i=0;i<ar1.length;++i)
        {
            int x=(i%w1);
            int y=(int)(i/w1);
            if(x>=w2 || y>=h2)
            {
                continue;
            }
            ar2[c]=ar1[i];
            c++;
        }
        
        return ar2;
    }
    
    void remove_NaNs(float[] imarray)
    {
    	float sum=0;
    	int counter=0;
    	int index=0;
    	ArrayList<Integer> nanlist=new ArrayList<Integer>();
    	for(float el : imarray)
    	{
    		if (Float.isNaN(el))
    			nanlist.add(index);
    		else
    		{
    			sum+=el;
    			counter+=1;
    		}
    		index+=1;
    	}
    	float mean=sum/counter;
    for (Integer ind : nanlist)
    {
    	imarray[ind]=mean;
    }
    
    	   
    }
  
    

}
