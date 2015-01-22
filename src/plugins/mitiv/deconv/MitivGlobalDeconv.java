package plugins.mitiv.deconv;

import icy.gui.frame.progress.AnnounceFrame;
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
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import loci.formats.ome.OMEXMLMetadata;
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
import mitiv.utils.reconstruction.ReconstructionThread;
import mitiv.utils.reconstruction.ReconstructionThreadToken;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.mitiv.io.IcyBufferedImageUtils;
import plugins.mitiv.reconstruction.TotalVariationJobForIcy;

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

        public String toString(){
            return new String("dxy: "+dxy+" dz: "+dz+" nxy: "+nxy+" nz: "+nz+" na "+na+" lambda "+lambda+" ni "+ni);
        }
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


        public double getValue(boolean withMult){
            if (withMult) {
                return mult*Double.valueOf(field.getText());
            } else {
                return Double.valueOf(field.getText());
            }
        }

        public double getValue(){
            return getValue(true);
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
    private myDouble dxy, dz, nxy, nz, na, lambda, ni;    //PSF
    private MicroscopyModelPSF1D pupil;
    private boolean psfInitFlag = false;
    private myDouble mu, epsilon, grtol, nbIteration, zeroPadding;          //Deconvolution
    private myDouble gain,noise;                          //VARIANCE
    private myDouble grtolPhase, grtolModulus, grtolDefocus, bDecTotalIteration;          //BDec
    private myComboBox image, canalImage, psf, weightsMethod, weights, deadPixel, deconvOptions, nbAlphaCoef, nbBetaCoef;
    private myBoolean deadPixGiven, restart, positivity;
    private String[] seqList;           //Global list given to all ComboBox that should show the actual image
    private final String[] deconvStringOptions = new String[]{"Total Variation","Tichonov"};
    private final String[] weightOptions = new String[]{"None","Personnalized map","Variance map","Computed variance"}; 
    private final String[] nAlphaOptions = new String[]{"-2","1","8","19","34","53","76","103","134","169"}; 
    private final String[] nBetaOptions = new String[]{"1","4","11","22","37","56","79","106","137","172"}; 
    private String[] canalImageOptions = new String[]{"None"}; 
    private myMetaData meta = null;     //The image metadata that we will move from one image to another
    private JButton psfShow, psfShow2, showModulus, showPhase;

    private JPanel psfGlob, imageGlob, varianceGlob, deconvGlob, bdecGlob; 
    private boolean canRunBdec = true;      //In the case where a psf is given we will not allow to run bdec
    private JTabbedPane tabbedPane;

    private Shape shape;
    boolean run = true;

    /*********************************/
    /**            Job              **/
    /*********************************/
    private ReconstructionThreadToken token;
    ReconstructionThread thread;

    /*********************************/
    /**            DEBUG            **/
    /*********************************/
    private boolean debug = false; //Show psf steps 
    private boolean verbose = true;    //show some values, need debug to true

    //Global variables for the algorithms
    TotalVariationJobForIcy tvDec;
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
    @SuppressWarnings("unused")
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
        if (token != null) {
            token.stop();
            token.exit();
        }
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
                    meta = getMetaData(seq);
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
        psfPannel.add((psfShow = new JButton("Show PSF")));

        psf.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                canRunBdec = psf.getValue().equals("None"); // If no psf given, we can run bdec tab, else we disable the run button in bdec tab
                na.setVisible(canRunBdec);                  // AND if we can run BDEC we can show the options of bdec
                ni.setVisible(canRunBdec);
                lambda.setVisible(canRunBdec);
                nxy.setVisible(canRunBdec);
                nz.setVisible(canRunBdec);
                dxy.setVisible(canRunBdec);
                dz.setVisible(canRunBdec);
                psfShow.setVisible(canRunBdec);
            }
        });

        psfShow.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // Show the initial PSF
                PSF0Clicked();
                if (debug) {
                    System.out.println("First PSF compute");
                }
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
        varianceTab.add((weightsMethod = createChoiceList(  "<html><pre> Weighting:     </pre></html>", weightOptions)));
        varianceTab.add((weights = createChoiceList(        "<html><pre> Map:     </pre></html>", seqList)));
        varianceTab.add((gain = createDouble(               "<html><pre>Gain:             </pre></html>", 0.0)));
        varianceTab.add((noise = createDouble(              "<html><pre>Readout Noise:    </pre></html>", 0.0)));
        varianceTab.add((deadPixGiven = new myBoolean(      "<html><pre>Dead Pixel Map ?  </pre></html>", false)));
        varianceTab.add((deadPixel = createChoiceList(      "<html><pre>Map: </pre></html>", seqList)));

        weightsMethod.addActionListener(new ActionListener() {
            //weightOptions =  new String[]{"None","Personnalized map","Variance map","Computed variance"};
            @Override
            public void actionPerformed(ActionEvent e) {
                if (weightsMethod.getValue() == weightOptions[0]) { //None
                    weights.setVisible(false);
                    gain.setVisible(false);
                    noise.setVisible(false);
                } else if (weightsMethod.getValue() == weightOptions[1] || weightsMethod.getValue() == weightOptions[2]) {  //Personnalized map ou Variance map
                    weights.setVisible(true);
                    gain.setVisible(false);
                    noise.setVisible(false);
                } else if (weightsMethod.getValue() == weightOptions[3]) {  //Computed variance
                    weights.setVisible(false);
                    gain.setVisible(true);
                    noise.setVisible(true);
                } else {
                    throw new IllegalArgumentException("Invalid argument passed to weight method");
                }
            }
        });

        deadPixGiven.addListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                deadPixel.setVisible(deadPixGiven.getValue());
            }
        });
        weights.setVisible(false);
        gain.setVisible(false);
        noise.setVisible(false);
        deadPixel.setVisible(false);
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
        deconvTab.add((deconvOptions = createChoiceList("<html><pre>Method:                           </pre></html>", deconvStringOptions)));
        deconvTab.add((mu = new myDouble(               "<html><pre>Mu:                               </pre></html>", 5E-4)));
        deconvTab.add((epsilon = new myDouble(          "<html><pre>Epsilon:                          </pre></html>", 1E-2)));
        deconvTab.add((grtol = new myDouble(            "<html><pre>Grtol:                            </pre></html>", 1E-2)));
        deconvTab.add((zeroPadding = new myDouble(      "<html><pre>Number of lines to add (padding): </pre></html>", 0)));
        deconvTab.add((nbIteration = new myDouble(      "<html><pre>Number of iterations:             </pre></html>", 50)));
        deconvTab.add((positivity = new myBoolean(      "<html><pre>Enable positivity:                </pre></html>", false)));
        deconvTab.add((restart = new myBoolean(         "<html><pre>Start from last result:           </pre></html>", false)));

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
        bdecTab.add((nbAlphaCoef = createChoiceList(    "<html><pre>N\u03B1:                          </pre></html>", nAlphaOptions)));
        bdecTab.add((nbBetaCoef = createChoiceList(     "<html><pre>N\u03B2:                          </pre></html>", nBetaOptions)));
        bdecTab.add((grtolDefocus = new myDouble(   "<html><pre>Grtol defocus:               </pre></html>", 0.1)));
        bdecTab.add((grtolPhase = new myDouble(     "<html><pre>Grtol phase:                 </pre></html>", 0.1)));
        bdecTab.add((grtolModulus = new myDouble(   "<html><pre>Grtol modulus:               </pre></html>", 0.1)));
        bdecTab.add((bDecTotalIteration = new myDouble("<html><pre>Number of total iterations:  </pre></html>", 2)));
        bdecTab.add((psfShow2 = new JButton(        "Show PSF"))); //Already created in psf tab
        bdecTab.add((showPhase = new JButton(       "Show phase of the pupil")));
        bdecTab.add((showModulus = new JButton(     "Show modulus of the pupil")));

        psfShow2.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // Show the initial PSF
                PSF0Clicked();
                if (debug) {
                    System.out.println("First PSF compute");
                }
            }
        });

        showPhase.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                phaseClicked();
                if (debug) {
                    System.out.println("Show phase");
                }
            }
        });

        showModulus.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                modulusClicked();
                if (debug) {
                    System.out.println("Show modulus");
                }
            }
        });
        //Creation of BDec TAB
        bdecGlob.add(bdecTab, BorderLayout.NORTH);
        tabbedPane.addTab("BDec", null, bdecGlob,    "All the options for the blind deconvolution");

        /****************************************************/
        /**                      ToolTips                  **/
        /****************************************************/
        image.setToolTipText(ToolTipText.sequenceImage);
        psf.setToolTipText(ToolTipText.sequencePSF);
        weights.setToolTipText(ToolTipText.sequenceWeigth);
        weightsMethod.setToolTipText(ToolTipText.sequenceWeigth);
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

        nbAlphaCoef.setToolTipText(ToolTipText.doubleNalpha);
        nbBetaCoef.setToolTipText(ToolTipText.doubleNbeta);
        grtolPhase.setToolTipText(ToolTipText.doubleGrtolPhase);
        grtolModulus.setToolTipText(ToolTipText.doubleGrtolModulus);
        grtolDefocus.setToolTipText(ToolTipText.doubleGrtolDefocus);
        bDecTotalIteration.setToolTipText(ToolTipText.doubleBDecTotalIteration);

        restart.setToolTipText(ToolTipText.booleanRestart);
        positivity.setToolTipText(ToolTipText.booleanPositivity);
        // Adding image, canalImage, psf, weights, deadPixel to auto refresh when sequence is added/removed
        listChoiceList.add(image);
        listChoiceList.add(psf);
        listChoiceList.add(weights);
        listChoiceList.add(deadPixel);
        //Add the tabbed pane to this panel.
        tabbedPane.addChangeListener(new ChangeListener() {
            //This is the part where we disable the run button if we are not in the good tab (bdec or deconv)
            //Also if there is a psf we canno't run the BDEC
            @Override
            public void stateChanged(ChangeEvent e) {
                if (tabbedPane.getSelectedComponent() == deconvGlob || (tabbedPane.getSelectedComponent() == bdecGlob && canRunBdec )) {
                    getUI().setRunButtonEnabled(true);
                }else{
                    getUI().setRunButtonEnabled(false);
                }
            }
        });
        getUI().setRunButtonEnabled(false); //Disable start button if we are not on bdec or deconv tab
        addComponent(tabbedPane);

        token = new ReconstructionThreadToken(new double[]{mu.getValue(),epsilon.getValue(),0.0,grtol.getValue()});
        thread = new ReconstructionThread(token);
        thread.start();
    }

    public void launchDeconvolution(DoubleArray imgArray, DoubleArray psfArray, DoubleArray weight){
        if (deconvOptions.getValue() == deconvStringOptions[0]) { //Total variation
            if (tvDec != null &&  restart.getValue()) {
                tvDec.setResult(tvDec.getResult());
            } else {
                tvDec = new TotalVariationJobForIcy(token);
                tvDec.setResult(null);
                tvDec.setAbsoluteTolerance(0.0);
                tvDec.setWeight(weight);
                tvDec.setData(imgArray);
                tvDec.setPsf(psfArray);
                tvDec.setViewer(new tvViewer());
                thread.setJob(tvDec);
            }
            tvDec.setPositivity(positivity.getValue());
            tvDec.setRegularizationWeight(mu.getValue());
            tvDec.setRegularizationThreshold(epsilon.getValue());
            tvDec.setRelativeTolerance(grtol.getValue());
            tvDec.setMaximumIterations((int)nbIteration.getValue());
            tvDec.setOutputShape(shape);
            token.start();
            setResult(tvDec);
        } else if (deconvOptions.getValue() == deconvStringOptions[1]) { //tichonov
            System.out.println("Ticho ON");
        } else {
            throw new IllegalArgumentException("Unknow deconvolution option");
        }
    }

    @Override
    protected void execute() {
        try {
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
                System.out.println("Positivity: "+positivity.getValue());
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
            //ArrayList<IcyBufferedImage> listImg = imgSeq.getAllImage();

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

            // Guessing the canal
            String canalToUse = canalImage.getValue();
            int numCanal = -1;
            for (int ii = 0; ii < imgSeq.getSizeC(); ii++) {
                if (canalToUse.equals(imgSeq.getChannelName(ii))) {
                    numCanal = ii;
                }
            }

            DoubleArray imgArray, psfArray;
            if (zeroPadding.getValue() < 0.0) {
                throw new IllegalArgumentException("Padding value canno't be inferior to the image size");
            }
            double coef = (width + zeroPadding.getValue())/width;
            boolean runBdec = (tabbedPane.getSelectedComponent() == bdecGlob); //If the BDEC panel is selected we the blind deconvolution
            // If no PSF is loaded -> creation of a PSF
            if (psfSeq == null) {
                psf0Init();
                pupil.computePSF();
                psfInitFlag = true;
                if (shape.rank() == 2) {
                    psfArray =  Double2D.wrap(MathUtils.uint16(MathUtils.fftShift1D(pupil.getPSF(), width, height)) , shape);
                } else {
                    psfArray =  Double3D.wrap(MathUtils.uint16(MathUtils.fftShift3D(pupil.getPSF(), width, height, sizeZ)) , shape);
                }
            } else {
                if (shape.rank() == 2) {
                    psfArray = (DoubleArray) IcyBufferedImageUtils.imageToArray(psfSeq, Shape.make(psfSeq.getWidth(), psfSeq.getHeight()), numCanal);
                } else {
                    psfArray = (DoubleArray) IcyBufferedImageUtils.imageToArray(psfSeq, Shape.make(psfSeq.getWidth(), psfSeq.getHeight(), psfSeq.getSizeZ()), numCanal);
                }
            }

            imgArray = (DoubleArray) IcyBufferedImageUtils.imageToArray(imgSeq, shape, numCanal);

            DoubleArray weight = createWeight(imgArray).toDouble();
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
                double[] alpha = new double[Integer.parseInt(nbAlphaCoef.getValue())];
                double[] beta = new double[Integer.parseInt(nbBetaCoef.getValue())];
                beta[0] = 1;
                double[] defocus = {ni.getValue()/lambda.getValue(), 0., 0.};
                DoubleShapedVectorSpace defocuSpace = new DoubleShapedVectorSpace(new int[]{defocus.length});
                DoubleShapedVector defocusVector = defocuSpace.wrap(defocus);
                DoubleShapedVectorSpace alphaSpace = new DoubleShapedVectorSpace(new int[]{alpha.length});
                DoubleShapedVector alphaVector = alphaSpace.create();
                DoubleShapedVectorSpace betaSpace = new DoubleShapedVectorSpace(new int[]{beta.length});
                DoubleShapedVector betaVector = betaSpace.wrap(beta);

                PSFEstimation = new PSF_Estimation();
                PSFEstimationInit();
                PSFEstimation.setWeight(weight.flatten());

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
                    PSFEstimation.setRelativeTolerance(0.1);
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
        } catch (IllegalArgumentException e) {
            new AnnounceFrame("Oops, Error: "+ e.getMessage());
        }
    }

    private ShapedArray createWeight(ShapedArray data){
        WeightGenerator weightGen = new WeightGenerator();
        ShapedArray array;
        ShapedArray deadPixMap = null;
        Sequence tmp;
        //If a dead pixel map was given
        if (deadPixGiven.getValue() && (tmp = getSequence(deadPixel)) != null) {
            deadPixMap = IcyBufferedImageUtils.imageToArray(tmp.getAllImage());
        }
        //We check the values given
        if (weightsMethod.getValue() == weightOptions[0]) { //None
            //we create an array of 1, that correspond to the image
            double[] weight = new double[width*height*sizeZ];
            for (int i = 0; i < weight.length; i++) {
                weight[i] = 1;
            }
            weightGen.setWeightMap(Double1D.wrap(weight, weight.length));
        } else if (weightsMethod.getValue() == weightOptions[1]) {  //Personnalized map
            if((tmp = getSequence(weights)) != null) {      
                array = IcyBufferedImageUtils.imageToArray(tmp.getAllImage());
                weightGen.setWeightMap(array);
            }
        } else if (weightsMethod.getValue() == weightOptions[2]) {  //Variance map
            if((tmp = getSequence(weights)) != null) {//Variance map
                array = IcyBufferedImageUtils.imageToArray(tmp.getAllImage());
                weightGen.setVarianceMap(array);
            }
        }else if (weightsMethod.getValue() == weightOptions[3]) {   //Computed variance
            weightGen.setComputedVariance(data, gain.getValue(), noise.getValue());
        }
        weightGen.setPixelMap(deadPixMap);
        return weightGen.getWeightMap(data.getShape()).toDouble();
    }

    //Get the results from a reconstruction and plot the intermediate result
    private void setResult(ReconstructionJob tvDec){
        try{
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
            sequence.endUpdate();
            sequence.setName("TV mu:"+mu.getValue()+" Iteration:"+tvDec.getIterations()+" grToll: "+tvDec.getRelativeTolerance());
        } catch (NullPointerException e) {
            //Here in case of brutal stop the sequence can become null but it's not important as it's an emergency stop
            //So we do nothing
            System.out.println("INFO: Emergency stop detected in setResult");
        }
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
        if(!psfInitFlag)
        {
            psf0Init();
            pupil.computePSF();
            psfInitFlag = true;
        }

        /* PSF0 Sequence */
        Sequence PSF0Sequence = new Sequence();
        PSF0Sequence.setName("PSF");
        double[] PSF_shift = MathUtils.fftShift3D(pupil.getPSF(), (int)nxy.getValue(), (int)nxy.getValue(), (int)nz.getValue());
        for (int k = 0; k < (int)nz.getValue(); k++)
        {
            PSF0Sequence.setImage(0, k, new IcyBufferedImage((int)nxy.getValue(), (int)nxy.getValue(),
                    MathUtils.getArray(PSF_shift, (int)nxy.getValue(), (int)nxy.getValue(), k)));
        }
        addSequence(PSF0Sequence);
        psfInitFlag = true;
    }

    private void phaseClicked()
    {
        /* PSF0 initialisation */
        if(!psfInitFlag)
        {
            psf0Init();
            pupil.computePSF();
            psfInitFlag = true;
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
        if(!psfInitFlag)
        {
            psf0Init();
            pupil.computePSF();
            psfInitFlag = true;
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
        if (meta == null) {
            meta = new myMetaData();
            meta.dxy     = metDat.getPixelsPhysicalSizeX(0).getValue().doubleValue()*1E3;  //FIXME In doubt I suppose it will give the size in umeters
            meta.dz      = metDat.getPixelsPhysicalSizeZ(0).getValue().doubleValue()*1E3;  //FIXME In doubt I suppose it will give the size in umeters
            if (metDat.getInstrumentCount() > 0) {
                meta.na      = metDat.getObjectiveLensNA(0, 0);
                meta.lambda  = metDat.getChannelEmissionWavelength(0, 0).getValue().doubleValue()*1E9;  //FIXME In doubt I suppose it will give the size in meters
                //data.ni      = metDat.getObjectiveImmersion(0, 0).getValue().doubleValue(); //STRANGE why a string
                if (debug && verbose) {
                    System.out.println("INFO: Metadata: No instrument so no metadata.");
                }
            } 
        }
        //If no instrument found, at least we have the right image size
        meta.nxy     = seq.getSizeX(); //We suppose X and Y equal
        meta.nz      = seq.getSizeZ();
        //We keep the existing values for n(xyz) and d(xyz)
        meta.dxy     = dxy.getValue()*1E9;  
        meta.dz      = dz.getValue()*1E9;
        meta.na      = na.getValue();
        meta.lambda  = lambda.getValue(false);
        meta.ni      = ni.getValue();
        return meta;
    }

    //TODO overwrite all data that we want to keep
    private void setMetaData(Sequence seqOld, Sequence seqNew) {
        //OMEXMLMetadataImpl metDat = seqOld.getMetadata();
        //FIXME AT this point if we modify the name we will modify the input and output, we need to duplicate/clone metDat 
        //seqNew.setMetaData(metDat);
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
        if (sequence == this.sequence) {
            this.sequence = null;
        }
    }

    @Override
    public void stopExecution() {
        if (token != null) {
            token.stop();
        }
        if (PSFEstimation != null) {
            PSFEstimation.stop();
        }
        run = false;
    }
}
