package plugins.mitiv.myEzPlug;

public class MyMetadata {
    //Just a public object with all psf values inside
    public double dxy     = 64.5;
    public double dz      = 160;
    public double nxy     = 256;
    public double nz      = 128;
    public double na      = 1.4;
    public double lambda  = 542;
    public double ni      = 1.518;

    public String toString(){
        return new String("dxy: "+dxy+" dz: "+dz+" nxy: "+nxy+" nz: "+nz+" na "+na+" lambda "+lambda+" ni "+ni);
    }
}
