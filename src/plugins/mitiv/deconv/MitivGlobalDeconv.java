package plugins.mitiv.deconv;

import icy.gui.main.GlobalSequenceListener;
import icy.image.IcyBufferedImage;
import icy.main.Icy;
import icy.sequence.Sequence;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;

import commands.TotalVariationDeconvolution;
import loci.formats.ome.OMEXMLMetadata;
import loci.formats.ome.OMEXMLMetadataImpl;
import mitiv.array.Double1D;
import mitiv.array.Double2D;
import mitiv.array.Double3D;
import mitiv.array.DoubleArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.invpb.ReconstructionJob;
import mitiv.invpb.ReconstructionViewer;
import mitiv.linalg.WeightGenerator;
import mitiv.linalg.shaped.DoubleShapedVector;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import mitiv.microscopy.MicroscopyModelPSF1D;
import mitiv.microscopy.PSF_Estimation;
import mitiv.utils.FFTUtils;
import mitiv.utils.MathUtils;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.mitiv.io.IcyBufferedImageUtils;

public class MitivGlobalDeconv extends EzPlug implements GlobalSequenceListener, EzStoppable {
    /***************************************************/
    /**                  MyMetaData                   **/
    /***************************************************/
    public class tvViewer implements ReconstructionViewer{
        @Override
        public void display(ReconstructionJob job) {
            setResult(job);
        }
    }

    /***************************************************/
    /**                  MyMetaData                   **/
    /***************************************************/
    public class myMetaData {
        //Just a public object with all psf values inside
        public double dxy     = 64.5;
        public double dz      = 160;
        public double nxy     = 256;
        public double nz      = 128;
        public double na      = 1.4;
        public double lambda  = 542;
        public double ni      = 1.518;
    }

    /***************************************************/
    /**                  MyJcomboBox                  **/
    /***************************************************/
    public class myComboBox extends JPanel{
        private static final long serialVersionUID = 1L;
        private JLabel filler;
        private JComboBox<String> jcb;

        myComboBox(String name, String[] inputs){
            filler = new JLabel(name);
            setLayout(new FlowLayout());
            add(filler);
            jcb = new JComboBox<String>();
            for (int j = 0; j < inputs.length; j++) {
                jcb.addItem(inputs[j]);
            }
            add(jcb);
        }

        public void updateData(String[] newInputs){
            for (int i = 0; i < newInputs.length; i++) {
                jcb.removeAllItems();
                for (int j = 0; j < newInputs.length; j++) {
                    jcb.addItem(newInputs[j]);
                }
            }
        }

        public String getValue(){
            if (jcb.getSelectedItem() == null) {
                return "None";
            }
            return jcb.getSelectedItem().toString();
        }

        public void setValue(String value){
            int size = jcb.getItemCount();
            for (int i = 0; i < size; i++) {
                String tmp = jcb.getItemAt(i);
                if (tmp.compareTo(value) == 0) {
                    jcb.setSelectedIndex(i);
                    return;
                }
            }
        }

        public void addActionListener(ActionListener l){
            jcb.addActionListener(l);
        }
    }

    /***************************************************/
    /**                  MyDouble                     **/
    /***************************************************/
    public class myDouble extends JPanel{
        private static final long serialVersionUID = 1L;
        private JLabel filler;
        private JTextField field;
        private double mult;
        public myDouble(String name, double input) {
            this(name, input, 1.0);
        }

        //Here we give the the multiplication factor, the result will be multiply by this factor
        public myDouble(String name, double input, double valueMult) {
            mult = valueMult;
            filler = new JLabel(name);
            field = new JTextField();
            field.setText(String.valueOf(input));
            field.setColumns(10);
            field.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent e) {
                    try {
                        Double.valueOf(field.getText());
                    } catch (Exception e2) {
                        field.setText("0.0");
                    }
                }
            });
            setLayout(new FlowLayout());
            add(filler);
            add(field);
        }



        public double getValue(){
            return mult*Double.valueOf(field.getText());
        }

        public void setValue(double value){
            field.setText(String.valueOf(value));
        }
    }

    /***************************************************/
    /**                  MyBoolean                    **/
    /***************************************************/
    public class myBoolean extends JPanel{
        private static final long serialVersionUID = 1L;
        private JLabel filler;
        private JCheckBox box;

        public myBoolean(String name, boolean defaultValue) {
            filler = new JLabel(name);
            box = new JCheckBox();
            box.setSelected(defaultValue);
            setLayout(new FlowLayout());
            add(filler);
            add(box);
        }

        public boolean getValue(){
            return box.isSelected();
        }

        public void addListener(ActionListener l){
            box.addActionListener(l);
        }
    }

    /***************************************************/
    /**                  All variables                **/
    /***************************************************/
    private ArrayList<myComboBox> listChoiceList = new ArrayList<myComboBox>();
    private myDouble dxy, dz, nxy, nz, na, lambda, ni, nbAlphaCoef, nbBetaCoef;    //PSF
    private MicroscopyModelPSF1D pupil;
    private int psfInitFlag = 0;
    private myDouble mu, epsilon, grtol, nbIteration, zeroPadding;          //Deconvolution
    private myDouble gain,noise;                          //VARIANCE
    private myDouble grtolPhase, grtolModulus, grtolDefocus, bDecTotalIteration;          //BDec
    private myComboBox image, canalImage, psf, weights, deadPixel, deconvOptions;
    private myBoolean deadPixGiven, restart;
    private String[] seqList;
    private final String[] deconvStringOptions = new String[]{"Total Variation","Tichonov"}; 
    private String[] canalImageOptions = new String[]{"None"}; 

    private JButton psfShow, showModulus, showPhase;

    private JPanel psfGlob, imageGlob, varianceGlob, deconvGlob, bdecGlob; 

    private JTabbedPane tabbedPane;

    private Shape shape;
    boolean run = true;

    /*********************************/
    /**            DEBUG            **/
    /*********************************/
    private boolean debug = true; //Show psf steps 
    private boolean verbose = false;    //show some values, need debug to true
    
    //Global variables for the algorithms
    TotalVariationDeconvolution tvDec;
    PSF_Estimation PSFEstimation;
    //Global variable for the deconvolution
    Sequence sequence; //The reference to the sequence we use to plot 
    int width, height, sizeZ;
    //Update all sequence with the new sequence remove
    private void updateAllList(){
        for (int i = 0; i < listChoiceList.size(); i++) {
            myComboBox box = listChoiceList.get(i);
            String current = box.getValue();
            box.updateData(seqList);
            box.setValue(current);
        }
    }

    //Mini factory for ComboBox => EzVarSequences
    private myComboBox createChoiceList(String name, String[] inputs){
        myComboBox jcb = new myComboBox(name, inputs);
        return jcb;
    }

    //Mini factory for double panel => EzVarDouble
    private myDouble createDouble(String name, double input){
        return new myDouble(name, input);
    }

    private myDouble createDouble(String name, double input, double mult){
        return new myDouble(name, input, mult);
    }

    //Mini Jlabel factory
    private JLabel createLabel(String input){
        JLabel tmp = new JLabel(input);
        tmp.setAlignmentX(Component.CENTER_ALIGNMENT);
        return tmp;
    }

    //Get from Icy all the sequences  presents
    private String[] getSequencesName(){
        ArrayList<Sequence> list = getSequences();
        String[] listString = new String[list.size()+1];
        listString[0] = "None";
        for (int i = 0; i < list.size(); i++) {
            listString[i+1] = list.get(i).getName();
        }
        return listString;
    }

    //Will convert the input String to the corresponding sequence
    private Sequence getSequence(myComboBox box){
        String seqName = box.getValue();
        ArrayList<Sequence> list = getSequences();
        for (Sequence sequence : list) {
            if (sequence.getName().compareTo(seqName) == 0) {
                return sequence;
            }
        }
        return null;
    }

    private void throwError(String s){
        throw new IllegalArgumentException(s);
    }
    @Override
    public void clean() {

    }

    @Override
    protected void initialize() {

        Icy.getMainInterface().addGlobalSequenceListener(this);
        seqList = getSequencesName();
        tabbedPane = new JTabbedPane();

        /****************************************************/
        /**                    IMAGE TAB                   **/
        /****************************************************/
        //Creation of the inside of IMAGE TAB
        imageGlob = new JPanel(new BorderLayout()); //Border layout to be sure that the images are stacked to the up
        JPanel imagePan = new JPanel(false);
        imagePan.setLayout(new BoxLayout(imagePan, BoxLayout.Y_AXIS));
        imagePan.add((image = createChoiceList(     "<html><pre>Sequence: </pre></html>", seqList)));
        imagePan.add((canalImage = createChoiceList("<html><pre>Canal:    </pre></html>", canalImageOptions)));
        image.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Sequence seq = getSequence(image);
                if (seq == null) {
                    canalImageOptions = new String[]{"None"};
                } else {
                    //Update channel availables in the image
                    int nbChan = seq.getSizeC();
                    canalImageOptions = new String[nbChan];
                    for (int i = 0; i < nbChan; i++) {
                        canalImageOptions[i] = seq.getChannelName(i);
                    }
                    canalImage.updateData(canalImageOptions);
                    //Update PSF metadata
                    myMetaData meta = getMetaData(seq); //FIXME Do it only one time
                    dxy.setValue(    meta.dxy);
                    dz.setValue(     meta.dz);
                    nxy.setValue(    meta.nxy);
                    nz.setValue(     meta.nz);
                    na.setValue(     meta.na);
                    lambda.setValue( meta.lambda);
                    ni.setValue(     meta.ni);
                }
            }
        });

        //Little hack to have a minimal size at startup
        //JLabel sizeTab = new JLabel("                                                                    ");
        //imagePan.add(sizeTab);

        //Creation of IMAGE TAB
        imageGlob.add(imagePan, BorderLayout.NORTH);
        tabbedPane.addTab("FILE", null, imageGlob, "Selecting input data");

        /****************************************************/
        /**                    PSF TAB                     **/
        /****************************************************/
        //Creation of the inside of PSF TAB
        psfGlob = new JPanel(new BorderLayout()); //Border layout to be sure that the images are stacked to the up
        JPanel psfPannel = new JPanel(false);
        psfPannel.setLayout(new BoxLayout(psfPannel, BoxLayout.Y_AXIS));
        psfPannel.add((psf = createChoiceList("<html><pre>Load PSF:       </pre></html>", seqList)));
        psfPannel.add((na = createDouble(     "<html><pre>NA:        </pre></html>", 1.4)));
        psfPannel.add((ni = createDouble(     "<html><pre>ni:        </pre></html>", 1.518)));
        psfPannel.add((lambda = createDouble( "<html><pre>\u03BB(nm):     </pre></html>", 542, 1E-9)));     //Here we give the the multiplication factor, the result will be multiply by this factor
        psfPannel.add((nxy = createDouble(    "<html><pre>Nxy:       </pre></html>", 256)));
        psfPannel.add((nz = createDouble(     "<html><pre>Nz:        </pre></html>", 128)));
        psfPannel.add((dxy = createDouble(    "<html><pre>dxy(nm):   </pre></html>", 64.5, 1E-9)));
        psfPannel.add((dz = createDouble(     "<html><pre>dz(nm):    </pre></html>", 160, 1E-9)));
        psfPannel.add((nbAlphaCoef = createDouble(     "<html><pre>N\u03B2:        </pre></html>", 76)));
        psfPannel.add((nbBetaCoef = createDouble(     "<html><pre>N\u03B1:        </pre></html>", 22)));
        psfPannel.add((psfShow = new JButton("Show PSF")));
        psfPannel.add((showPhase = new JButton("Show phase of the pupil")));
        psfPannel.add((showModulus = new JButton("Show modulus of the pupil")));

        psfShow.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // Show the initial PSF
                PSF0Clicked();
                System.out.println("First PSF compute");
            }
        });

        showPhase.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                phaseClicked();
                System.out.println("Show phase");
            }
        });

        showModulus.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                modulusClicked();
                System.out.println("Show modulus");
            }
        });


        //Creation of IMAGE TAB
        psfGlob.add(psfPannel, BorderLayout.NORTH);
        tabbedPane.addTab("PSF", null, psfGlob, "Choice of the PSF, visuakization of theoritical PSF");

        /****************************************************/
        /**                 VARIANCE TAB                   **/
        /****************************************************/
        //Creation of the inside of VARIANCE TAB
        varianceGlob = new JPanel(new BorderLayout()); //Border layout to be sure that the images are stacked to the up
        JPanel varianceTab = new JPanel(false);
        varianceTab.setLayout(new BoxLayout(varianceTab, BoxLayout.Y_AXIS));
        varianceTab.add((weights = createChoiceList(  "<html><pre>Variance MAP:     </pre></html>", seqList)));
        varianceTab.add(createLabel(                  "Variance Map:"));
        varianceTab.add((gain = createDouble(         "<html><pre>Gain:             </pre></html>", 0.0)));
        varianceTab.add((noise = createDouble(        "<html><pre>Readout Noise:    </pre></html>", 0.0)));
        varianceTab.add((deadPixGiven = new myBoolean("<html><pre>Dead Pixel Map ?  </pre></html>", false)));
        deadPixel = createChoiceList(                 "<html><pre>Dead pixel Map:   </pre></html>", seqList);
        deadPixGiven.addListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                deadPixel.setVisible(deadPixGiven.getValue());
            }
        });
        deadPixel.setVisible(false);
        varianceTab.add(deadPixel);
        //Creation of VARIANCE TAB
        varianceGlob.add(varianceTab, BorderLayout.NORTH);
        tabbedPane.addTab("Variance", null, varianceGlob, "Selecting the weights, the variance and the dead pixels");

        /****************************************************/
        /**                    DECONV TAB                  **/
        /****************************************************/
        //Creation of the inside of DECONVOLUTION TAB
        deconvGlob = new JPanel(new BorderLayout()); //Border layout to be sure that the images are stacked to the up
        JPanel deconvTab = new JPanel(false);
        deconvTab.setLayout(new BoxLayout(deconvTab, BoxLayout.Y_AXIS));
        deconvTab.add((deconvOptions = createChoiceList("<html><pre>Method:                  </pre></html>", deconvStringOptions)));
        deconvTab.add((mu = new myDouble(               "<html><pre>Mu:                      </pre></html>", 5E-4)));
        deconvTab.add((epsilon = new myDouble(          "<html><pre>Epsilon:                 </pre></html>", 1E-2)));
        deconvTab.add((grtol = new myDouble(            "<html><pre>Grtol:                   </pre></html>", 1E-2)));
        deconvTab.add((zeroPadding = new myDouble(      "<html><pre>Padding multiplication:  </pre></html>", 1.0)));
        deconvTab.add((nbIteration = new myDouble(      "<html><pre>Number of iterations:    </pre></html>", 50)));
        deconvTab.add((restart = new myBoolean(         "<html><pre>Start from last result:  </pre></html>", false)));

        //Creation of DECONVOLUTION TAB
        deconvGlob.add(deconvTab, BorderLayout.NORTH);
        tabbedPane.addTab("Deconvolution", null, deconvGlob, "Methods and options of the deconvolution");


        /****************************************************/
        /**                      BDEC TAB                  **/
        /****************************************************/
        //Creation of the inside of BDec TAB
        bdecGlob = new JPanel(new BorderLayout()); //Border layout to be sure that the images are stacked to the up
        JPanel bdecTab = new JPanel(false);
        bdecTab.setLayout(new BoxLayout(bdecTab, BoxLayout.Y_AXIS));
        bdecTab.add((grtolDefocus = new myDouble("<html><pre>Grtol defocus:           </pre></html>", 0.1)));
        bdecTab.add((grtolPhase = new myDouble("<html><pre>Grtol phase:           </pre></html>", 0.1)));
        bdecTab.add((grtolModulus = new myDouble("<html><pre>Grtol modulus:               </pre></html>", 0.1)));
        bdecTab.add((bDecTotalIteration = new myDouble("<html><pre>Number of total iterations:        </pre></html>", 2)));
        //Creation of BDec TAB
        bdecGlob.add(bdecTab, BorderLayout.NORTH);
        tabbedPane.addTab("BDec", null, bdecGlob,    "All the options for the blind deconvolution");

        /****************************************************/
        /**                      ToolTips                  **/
        /****************************************************/
        image.setToolTipText(ToolTipText.sequenceImage);
        psf.setToolTipText(ToolTipText.sequencePSF);
        weights.setToolTipText(ToolTipText.sequenceWeigth);
        deadPixel.setToolTipText(ToolTipText.sequencePixel);

        canalImage.setToolTipText(ToolTipText.textCanal);
        deconvOptions.setToolTipText(ToolTipText.textMethod);

        dxy.setToolTipText(ToolTipText.doubleDxy);
        dz.setToolTipText(ToolTipText.doubleDz);
        nxy.setToolTipText(ToolTipText.doubleNxy);
        nz.setToolTipText(ToolTipText.doubleNz);
        na.setToolTipText(ToolTipText.doubleNa);
        lambda.setToolTipText(ToolTipText.doubleLambda);
        ni.setToolTipText(ToolTipText.doubleNi);
        gain.setToolTipText(ToolTipText.doubleGain);
        noise.setToolTipText(ToolTipText.doubleNoise);
        mu.setToolTipText(ToolTipText.doubleMu);
        epsilon.setToolTipText(ToolTipText.doubleEpsilon);
        grtol.setToolTipText(ToolTipText.doubleGrtoll);
        nbIteration.setToolTipText(ToolTipText.doubleMaxIter);
        zeroPadding.setToolTipText(ToolTipText.doublePadding);

        grtolPhase.setToolTipText(ToolTipText.doubleGrtolPhase);
        grtolModulus.setToolTipText(ToolTipText.doubleGrtolModulus);
        grtolDefocus.setToolTipText(ToolTipText.doubleGrtolDefocus);
        bDecTotalIteration.setToolTipText(ToolTipText.doubleBDecTotalIteration);

        restart.setToolTipText(ToolTipText.booleanRestart);
        // Adding image, canalImage, psf, weights, deadPixel to auto refresh when sequence is added/removed
        listChoiceList.add(image);
        listChoiceList.add(psf);
        listChoiceList.add(weights);
        listChoiceList.add(deadPixel);
        //Add the tabbed pane to this panel.
        addComponent(tabbedPane);
    }

    public void launchDeconvolution(DoubleArray imgArray, DoubleArray psfArray, double[] weight){
        if (deconvOptions.getValue() == deconvStringOptions[0]) { //Total variation
            if (tvDec != null &&  restart.getValue()) {
                tvDec.setResult(tvDec.getResult());
            } else {
                tvDec = new TotalVariationDeconvolution();
                tvDec.setResult(null);
                tvDec.setAbsoluteTolerance(0.0);
                tvDec.setWeight(weight);
                tvDec.setData(imgArray);
                tvDec.setPsf(psfArray);
                tvDec.setViewer(new tvViewer());
            }
            tvDec.setRegularizationWeight(mu.getValue()); //TODO : faire comme PSFEstimationInit
            tvDec.setRegularizationThreshold(epsilon.getValue());
            tvDec.setRelativeTolerance(grtol.getValue());
            tvDec.setMaximumIterations((int)nbIteration.getValue());
            tvDec.deconvolve(shape);
            setResult(tvDec);
        } else if (deconvOptions.getValue() == deconvStringOptions[1]) { //tichonov
            System.out.println("Ticho ON");
        } else {
            throw new IllegalArgumentException("Unknow deconvolution option");
        }
    }

    @Override
    protected void execute() {
        if (debug && verbose) {
            System.out.println("-------------IMAGE-------------------");
            System.out.println("File: "+image.getValue());              //Used
            System.out.println("Canal: "+canalImage.getValue());
            System.out.println("--------------PSF------------------");
            System.out.println("PSF: "+psf.getValue());                 //Used
            System.out.println("dxy: "+dxy.getValue());
            System.out.println("dz: "+dz.getValue());
            System.out.println("nxy: "+nxy.getValue());
            System.out.println("nz: "+nz.getValue());
            System.out.println("NA: "+na.getValue());
            System.out.println("\u03BB: "+lambda.getValue());
            System.out.println("ni: "+ni.getValue());
            System.out.println("--------------Variance------------------");
            System.out.println("Weights: "+weights.getValue());
            System.out.println("Gain: "+gain.getValue());
            System.out.println("Noise: "+noise.getValue());
            System.out.println("deadPix: "+deadPixel.getValue());
            System.out.println("--------------DECONV------------------");
            System.out.println("Methode: "+deconvOptions.getValue());   //Used
            System.out.println("zeroPad: "+zeroPadding.getValue());
            System.out.println("nbIter: "+nbIteration.getValue());
            System.out.println("Restart: "+restart.getValue());
            System.out.println("--------------BDEC------------------");
            System.out.println("nbIter: "+nbIteration.getValue());
            System.out.println("zeroPad: "+zeroPadding.getValue());
            System.out.println("nbIterZern: "+grtolPhase.getValue());
            System.out.println("module: "+grtolModulus.getValue());
            System.out.println("defoc: "+grtolDefocus.getValue());
            System.out.println("Number of total iterations: "+bDecTotalIteration.getValue());
            System.out.println("");
        }
        run = true;

        // Preparing parameters and testing input
        Sequence imgSeq = getSequence(image);
        Sequence psfSeq = getSequence(psf);
        if (imgSeq == null)
        {
            throwError("An image/sequence of images should be given");
        }
        ArrayList<IcyBufferedImage> listImg = imgSeq.getAllImage();

        BufferedImage img = imgSeq.getFirstNonNullImage();
        // Set the informations about the input
        width = img.getWidth();
        height = img.getHeight();
        sizeZ = imgSeq.getSizeZ();
        if (sizeZ == 1) { //2D
            shape = Shape.make(width, height);
        } else {    //3D
            shape = Shape.make(width, height, sizeZ);
        }
        DoubleArray imgArray, psfArray;
        double coef = zeroPadding.getValue();
        boolean runBdec = true; 
        if (tabbedPane.getSelectedComponent() == deconvGlob) {  //If the deconv panel is selected we run only the deconvolution
            runBdec = false;
        }
        // If no PSF is loaded -> creation of a PSF
        if ( psfSeq == null) {
            if (runBdec) {
                
                //if(psfInitFlag == 0) {
                
                    psf0Init();
                    pupil.computePSF();
                    psfInitFlag = 1;
                //}
                if (shape.rank() == 2) {
                    psfArray =  Double2D.wrap(MathUtils.uint16(MathUtils.fftShift1D(pupil.getPSF(), width, height)) , shape);
                } else {
                    psfArray =  Double3D.wrap(MathUtils.uint16(MathUtils.fftShift3D(pupil.getPSF(), width, height, sizeZ)) , shape);
                }
            } else {
                throw new IllegalArgumentException("You should give a PSF or run from BDEC tab");
            }
        } else {
            ArrayList<IcyBufferedImage> listPSf = psfSeq.getAllImage();
            if (shape.rank() == 2) {
                double[] psfTmp = IcyBufferedImageUtils.icyImage3DToArray1D(listPSf, psfSeq.getWidth(), psfSeq.getHeight(), sizeZ, false);
                psfArray =  Double2D.wrap(psfTmp, Shape.make(psfSeq.getWidth(), psfSeq.getHeight()));
            } else {
                double[] psfTmp = IcyBufferedImageUtils.icyImage3DToArray1D(listPSf, width, height, sizeZ, false);
                psfArray =  Double3D.wrap(psfTmp, Shape.make(psfSeq.getWidth(), psfSeq.getHeight(), psfSeq.getSizeZ()));
            }
        }

        double[] image = IcyBufferedImageUtils.icyImage3DToArray1D(listImg, width, height, sizeZ, false);
        double[] weight = createWeight(image);
        if (shape.rank() == 2) {
            imgArray =  Double2D.wrap(image, shape);
        } else {
            imgArray =  Double3D.wrap(image, shape);
        }
        //BEWARE here we change the value to match the new padded image size
        width = FFTUtils.bestDimension((int)(width*coef));
        height = FFTUtils.bestDimension((int)(height*coef));
        sizeZ = FFTUtils.bestDimension((int)(sizeZ*coef));
        if (shape.rank() == 2) {
            shape = Shape.make(width, height);
        } else {
            shape = Shape.make(width, height, sizeZ);
        }
        /*---------------------------------------*/
        /*            OPTIMISATION               */
        /*---------------------------------------*/
        if (runBdec) {
            double[] alpha = new double[(int)nbAlphaCoef.getValue()];
            double[] beta = new double[(int)nbBetaCoef.getValue()];
            beta[0] = 1;
            double[] defocus = {ni.getValue()/lambda.getValue(), 0., 0.};
            DoubleShapedVectorSpace defocuSpace = new DoubleShapedVectorSpace(new int[]{defocus.length});
            DoubleShapedVector defocusVector = defocuSpace.wrap(defocus);
            DoubleShapedVectorSpace alphaSpace = new DoubleShapedVectorSpace(new int[]{alpha.length});
            DoubleShapedVector alphaVector = alphaSpace.create();
            DoubleShapedVectorSpace betaSpace = new DoubleShapedVectorSpace(new int[]{beta.length});
            DoubleShapedVector betaVector = betaSpace.wrap(beta);

            PSFEstimation = new PSF_Estimation();
            PSFEstimationInit(); //TODO : peut-être foutre les autres plus bas dedans
            PSFEstimation.setWeight(weight); //TODO : peut-être le mettre avant la boucle, inutile de réinitialiser les poids et data

            PSFEstimation.setData(imgArray);

            for(int i = 0; i < bDecTotalIteration.getValue(); i++) {
                /* OBJET ESTIMATION (by the current PSF) */
                launchDeconvolution(imgArray, psfArray, weight);

                /* PSF ESTIMATION (by the current objet) */


                PSFEstimation.setPupil(pupil);
                PSFEstimation.setPsf(tvDec.getData());

                /* Defocus estimation */
                if (debug && verbose) {
                    System.out.println("------------------");
                    System.out.println("Defocus estimation");
                    System.out.println("------------------");
                }
                PSFEstimation.setRelativeTolerance(0.1); //TODO : grtolDefocus..
                PSFEstimation.fitPSF(defocusVector, PSF_Estimation.DEFOCUS);

                /* Phase estimation */
                if (debug && verbose) {
                    System.out.println("Phase estimation");
                    System.out.println("------------------");
                }
                PSFEstimation.setResult(null);
                PSFEstimation.setRelativeTolerance(0.1);
                PSFEstimation.fitPSF(alphaVector, PSF_Estimation.ALPHA);

                /* Modulus estimation */
                if (debug && verbose) {
                    System.out.println("Modulus estimation");
                    System.out.println("------------------");
                }
                PSFEstimation.setResult(null);
                PSFEstimation.setRelativeTolerance(0.1);
                PSFEstimation.fitPSF(betaVector, PSF_Estimation.BETA);
                MathUtils.normalise(betaVector.getData());
                if (debug) {
                    showResult(i);
                }
                //If we want a emergency stop
                if (!run) {
                    return;
                }
            }
        } else {
            launchDeconvolution(imgArray, psfArray, weight);
        }
        setMetaData(imgSeq, sequence);
        sequence = null; //In any cases the next image will be in a new sequence
    }

    private double[] createWeight(double[] data){
        WeightGenerator weightGen = new WeightGenerator();
        ShapedArray array;
        ShapedArray deadPixMap = null;
        Sequence tmp;
        if (deadPixGiven.getValue() && (tmp = getSequence(deadPixel)) != null) {
            deadPixMap = IcyBufferedImageUtils.imageToArray(tmp.getAllImage());
        }
        if((tmp = getSequence(weights)) != null) {//Variance map
            //If the user give a varianceMap map we give it to weightGen
            array = IcyBufferedImageUtils.imageToArray(tmp.getAllImage());
            weightGen.setWeightMap(array);//Gain + readOut

        } else if(gain.getValue() != 0.0 &&  noise.getValue() != 0.0) {
            //In the case of computed variance: we give gain, readout noise, and the image
            array = Double1D.wrap(data, data.length);
            weightGen.setComputedVariance(array, gain.getValue(), noise.getValue());

        } else {
            //Last case, we create an array of 1, that correspond to the image
            double[] weight = new double[width*height*sizeZ];
            for (int i = 0; i < weight.length; i++) {
                weight[i] = 1;
            }
            weightGen.setWeightMap(Double1D.wrap(weight, weight.length));
        }
        weightGen.setPixelMap(deadPixMap);
        return weightGen.getWeightMap().toDouble().flatten();
    }

    //Get the results from a reconstruction and plot the intermediate result
    private void setResult(ReconstructionJob tvDec){
        if (sequence == null || (sequence != null && sequence.isEmpty())) {
            sequence = new Sequence();
            addSequence(sequence);
        }
        sequence.beginUpdate();
        double[] in = tvDec.getResult().toDouble().flatten();
        for (int j = 0; j < sizeZ; j++) {
            double[] temp = new double[width*height];
            for (int i = 0; i < width*height; i++) {
                temp[i] = in[i+j*width*height];
            }
            sequence.setImage(0,j, new IcyBufferedImage(width, height, temp));
        }
        sequence.setName("TV mu:"+mu.getValue()+" Iteration:"+tvDec.getIterations()+" grToll: "+tvDec.getRelativeTolerance());
        sequence.endUpdate();
    }

    private void showResult(int num)
    {
        Sequence psf3DSequence = new Sequence();
        double[] PSF_shift = MathUtils.fftShift3D(pupil.getPSF(), (int)nxy.getValue(), (int)nxy.getValue(), (int)nz.getValue());
        psf3DSequence.setName("PSF Estimated - " + num);
        for (int k = 0; k < (int)nz.getValue(); k++)
        {
            psf3DSequence.setImage(0, k, new IcyBufferedImage((int)nxy.getValue(), (int)nxy.getValue(),
                    MathUtils.getArray(PSF_shift, (int)nxy.getValue(), (int)nxy.getValue(), k)));
        }
        addSequence(psf3DSequence);
    }

    private void psf0Init()
    {

        double ns = 0;
        double zdepth = 0;
        int use_depth_scaling = 0;
        pupil = new MicroscopyModelPSF1D(na.getValue(), lambda.getValue(), ni.getValue(), ns, zdepth, dxy.getValue(),
                dz.getValue(), (int)nxy.getValue(), (int)nxy.getValue(), (int)nz.getValue(), use_depth_scaling);
    }

    private void PSF0Clicked()
    {
        /* PSF0 initialisation */
        if(psfInitFlag == 0)
        {
            psf0Init();
            pupil.computePSF();
            psfInitFlag = 1;
        }

        /* PSF0 Sequence */
        Sequence PSF0Sequence = new Sequence();
        PSF0Sequence.setName("Initial PSF");
        double[] PSF_shift = MathUtils.fftShift3D(pupil.getPSF(), (int)nxy.getValue(), (int)nxy.getValue(), (int)nz.getValue());
        for (int k = 0; k < (int)nz.getValue(); k++)
        {
            PSF0Sequence.setImage(0, k, new IcyBufferedImage((int)nxy.getValue(), (int)nxy.getValue(),
                    MathUtils.getArray(PSF_shift, (int)nxy.getValue(), (int)nxy.getValue(), k)));
        }
        addSequence(PSF0Sequence);
        psfInitFlag = 1;
    }

    private void phaseClicked()
    {
        /* PSF0 initialisation */
        if(psfInitFlag == 0)
        {
            psf0Init();
            pupil.computePSF();
            psfInitFlag = 1;
        }
        /* Phase Sequence */
        Sequence phaseSequence = new Sequence();
        phaseSequence.setName("Phase of the pupil");
        double[] phase_shift = MathUtils.fftShift1D(pupil.getPhi(), (int)nxy.getValue(), (int)nxy.getValue());
        phaseSequence.addImage(new IcyBufferedImage((int)nxy.getValue(), (int)nxy.getValue(), phase_shift));
        addSequence(phaseSequence);
    }

    private void modulusClicked()
    {
        /* PSF0 initialisation */
        if(psfInitFlag == 0)
        {
            psf0Init();
            pupil.computePSF();
            psfInitFlag = 1;
        }
        /* Modulus Sequence */
        Sequence modulusSequence = new Sequence();
        modulusSequence.setName("Modulus of the pupil");
        double[] modulus_shift = MathUtils.fftShift1D(pupil.getRho(), (int)nxy.getValue(), (int)nxy.getValue());
        modulusSequence.addImage(new IcyBufferedImage((int)nxy.getValue(), (int)nxy.getValue(), modulus_shift));
        addSequence(modulusSequence);
    }

    private void PSFEstimationInit()
    {
        PSFEstimation.setRegularizationWeight(0.1);
        PSFEstimation.setRegularizationThreshold(0.01);
        //pupilEstim.setRelativeTolerance(1.);
        PSFEstimation.setAbsoluteTolerance(0.);
        PSFEstimation.setMaximumIterations(10);
    }

    private myMetaData getMetaData(Sequence seq){
        OMEXMLMetadata metDat = seq.getMetadata();
        myMetaData data = new myMetaData();
        //If no instrument found, at least we have the right image size
        data.nxy     = seq.getSizeX(); //We suppose X and Y equal
        data.nz      = seq.getSizeZ();
        if (metDat.getInstrumentCount() > 0) {
            data.dxy     = metDat.getPixelsSizeX(0).getValue().doubleValue()*1E9;   //FIXME In doubt I suppose it will give the size in meters
            data.dz      = metDat.getPixelsSizeZ(0).getValue().doubleValue()*1E9;   //FIXME In doubt I suppose it will give the size in meters
            data.na      = metDat.getObjectiveLensNA(0, 0);
            data.lambda  = metDat.getChannelEmissionWavelength(0, 0).getValue().doubleValue()*1E9;  //FIXME In doubt I suppose it will give the size in meters
            //data.ni      = metDat.getObjectiveImmersion(0, 0).getValue().doubleValue(); //STRANGE why a string
        } else {
            if (debug && verbose) {
                System.out.println("INFO: Metadata: No instrument so no metadata.");
            }
        }
        return data;
    }

    private void setMetaData(Sequence seqOld, Sequence seqNew) {
        OMEXMLMetadataImpl metDat = seqOld.getMetadata();
        seqNew.setMetaData(metDat);
    }

    private void update(){
        seqList = getSequencesName();
        updateAllList();
    }

    @Override
    public void sequenceOpened(Sequence sequence) {
        update();
    }

    @Override
    public void sequenceClosed(Sequence sequence) {
        update();
    }

    @Override
    public void stopExecution() {
        if (tvDec != null) {
            tvDec.stop();
        }
        if (PSFEstimation != null) {
            PSFEstimation.stop();
        }
        run = false;
    }
}
