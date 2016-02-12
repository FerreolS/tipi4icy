package plugins.mitiv.deconv;

import icy.gui.frame.progress.AnnounceFrame;
import icy.gui.main.GlobalSequenceListener;
import icy.image.IcyBufferedImage;
import icy.main.Icy;
import icy.sequence.MetaDataUtil;
import icy.sequence.Sequence;
import icy.util.OMEUtil;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import loci.common.services.ServiceException;
import loci.formats.ome.OMEXMLMetadata;
import loci.formats.ome.OMEXMLMetadataImpl;
import mitiv.array.ArrayUtils;
import mitiv.array.Double1D;
import mitiv.array.Double3D;
import mitiv.array.DoubleArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.invpb.ReconstructionJob;
import mitiv.invpb.ReconstructionViewer;
import mitiv.linalg.WeightGenerator;
import mitiv.linalg.shaped.DoubleShapedVector;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import mitiv.microscopy.WideFieldModel;
import mitiv.microscopy.PSF_Estimation;
import mitiv.utils.FFTUtils;
import mitiv.utils.MathUtils;
import mitiv.utils.reconstruction.ReconstructionThread;
import mitiv.utils.reconstruction.ReconstructionThreadToken;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.mitiv.io.IcyBufferedImageUtils;
import plugins.mitiv.myEzPlug.MyBoolean;
import plugins.mitiv.myEzPlug.MyComboBox;
import plugins.mitiv.myEzPlug.MyDouble;
import plugins.mitiv.myEzPlug.MyMetadata;
import plugins.mitiv.reconstruction.TotalVariationJobForIcy;

/**
 * MiTivGlobalDeconv is a blind deconvolution tool built on the same basis than
 * MiTivTotalVariation. The blind deconvolution process is trying to guess the PSF 
 * and then the it is a standard deconvolution. By iterating on the result we can affine 
 * the PSF until the result is good enough for the user.
 * 
 * @author light
 *
 */
public class MitivBlindDeconvolution extends EzPlug implements GlobalSequenceListener, EzStoppable {
    /***************************************************/
    /**               Viewer Update result            **/
    /***************************************************/
    public class TvViewer implements ReconstructionViewer{
        @Override
        public void display(ReconstructionJob job) {
            setResult(job);
        }
    }

    /***************************************************/
    /**                  All variables                **/
    /***************************************************/
    private ArrayList<MyComboBox> listChoiceList = new ArrayList<MyComboBox>(); //Contain all the list that will be updated
    private MyDouble dxy, dz, nxy, nz, na, lambda, ni;    //PSF
    private MyDouble ns;
    private MyDouble zdepth;
    private boolean use_depth_scaling = false;
    private WideFieldModel pupil=null;
    // private boolean psfInitFlag = false;
    private MyDouble mu, epsilon, nbIteration, zeroPaddingxy, zeroPaddingz;          //Deconvolution
    private MyDouble gain,noise;                          //VARIANCE
    private MyDouble bDecTotalIteration,DefocusMaxIter,PhaseMaxIter,ModulusMaxIter;          //BDec
    private MyComboBox image, canalImage, psf, restart, weightsMethod, weights, deadPixel, nbAlphaCoef, nbBetaCoef;
    private MyBoolean deadPixGiven, positivity, crop, resetPSF;
    private String[] seqList;           //Global list given to all ComboBox that should show the actual image
    private final String[] weightOptions = new String[]{"None","Inverse covariance map","Variance map","Computed variance"}; 
    private final String[] nAlphaOptions = new String[]{"1","8","19","34","53","76","103","134","169"}; 
    private final String[] nBetaOptions = new String[]{"1","4","11","22","37","56","79","106","137","172"}; 
    private String[] canalImageOptions = new String[]{"None"}; 
    private MyMetadata meta = null;     //The image metadata that we will move from one image to another
    private JButton saveMetaData, showPSF, psfShow2, showWeight, showModulus, showPhase;
    private JLabel resultPad, resultCostPrior, resultDefocus, resultPhase, resultModulus;

    private JPanel psfGlob, imageGlob, varianceGlob, deconvGlob, bdecGlob, resultGlob; 
    private boolean canRunBdec = true;      //In the case where a psf is given we will not allow to run bdec
    private JTabbedPane tabbedPane;

    private double grtol = 0.0;
    private int nbAlpha=0, nbBeta=1;


    DoubleShapedVectorSpace defocuSpace = null, alphaSpace=null, betaSpace=null;
    DoubleShapedVector defocusVector = null, alphaVector = null, betaVector =null;


    private Shape shape, shapePad;
    boolean run = true;
    boolean runBdec;

    /*********************************/
    /**            Job              **/
    /*********************************/
    private ReconstructionThreadToken token;
    ReconstructionThread thread;

    /*********************************/
    /**            DEBUG            **/
    /*********************************/
    private boolean debug = true;      //Show psf steps 
    private boolean verbose = false;     //show some values, need debug to true

    //Global variables for the algorithms
    TotalVariationJobForIcy tvDec;
    PSF_Estimation PSFEstimation;
    //Global variable for the deconvolution
    Sequence sequence; //The reference to the sequence we use to plot 
    int width, height, sizeZ, widthPad, heightPad, sizeZPad;

    //Update all sequence with the new sequence remove
    private void updateAllList(){
        for (int i = 0; i < listChoiceList.size(); i++) {
            MyComboBox box = listChoiceList.get(i);
            String current = box.getValue();
            box.updateData(seqList);
            box.setValue(current);
        }
    }

    /*********************************/
    /**            FACTORY          **/
    /*********************************/
    //Mini factory for ComboBox => EzVarSequences
    private MyComboBox createChoiceList(String name, String[] inputs){
        MyComboBox jcb = new MyComboBox(name, inputs);
        return jcb;
    }

    //Mini factory for double panel => EzVarDouble
    private MyDouble createDouble(String name, double input){
        return new MyDouble(name, input);
    }

    private MyDouble createDouble(String name, double input, double mult){
        return new MyDouble(name, input, mult);
    }

    //Mini Jlabel factory
    @SuppressWarnings("unused")
    private JLabel createLabel(String input){
        JLabel tmp = new JLabel(input);
        tmp.setAlignmentX(Component.CENTER_ALIGNMENT);
        return tmp;
    }

    /*********************************/
    /**    Communication with Icy   **/
    /*********************************/
    //Get, from Icy, all the existings sequences
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
    private Sequence getSequence(MyComboBox box){
        String seqName = box.getValue();
        ArrayList<Sequence> list = getSequences();
        for (Sequence sequence : list) {
            if (sequence.getName().compareTo(seqName) == 0) {
                return sequence;
            }
        }
        return null;
    }

    // Guessing the canal
    private int getNumCanal(Sequence imgSeq) {
        String canalToUse = canalImage.getValue();
        for (int ii = 0; ii < imgSeq.getSizeC(); ii++) {
            if (canalToUse.equals(imgSeq.getChannelName(ii))) {
                return ii;
            }
        }
        return -1;
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

    protected void updateOutputPad() {
        int xy = FFTUtils.bestDimension((int)(nxy.getValue() + zeroPaddingxy.getValue()));
        int z = FFTUtils.bestDimension((int)(nz.getValue() + zeroPaddingz.getValue()));
        String text = "<html><pre>Output: <font color=\"red\">"+xy+"x"+xy+"x"+z+"</font></pre></html>";
        resultPad.setText(text);
    }

    /*********************************/
    /**      Initialization         **/
    /*********************************/
    @Override
    protected void initialize() {
        Icy.getMainInterface().addGlobalSequenceListener(this);
        getUI().setParametersIOVisible(false);
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
        nxy = createDouble(             "<html><pre>Nxy:        </pre></html>", 512);
        nz = createDouble(              "<html><pre>Nz:         </pre></html>",128);
        resultPad = new JLabel(         "<html><pre>Output: <font color=\"red\">"+(int)nxy.getValue()+"x"+(int)nxy.getValue()+"x"+(int)nz.getValue()+"</font></pre></html>");
        dxy = createDouble(             "<html><pre>dxy(nm):    </pre></html>", 64);
        dz = createDouble(              "<html><pre>dz(nm):     </pre></html>", 160);
        zeroPaddingxy = new MyDouble(   "<html><pre>padding xy: </pre></html>", 0);
        zeroPaddingz = new MyDouble(    "<html><pre>padding z : </pre></html>", 0);

        nxy.setEditable(false);
        nz.setEditable(false);

        imagePan.add(nxy);
        imagePan.add(nz);
        imagePan.add(zeroPaddingxy);
        imagePan.add(zeroPaddingz);
        imagePan.add(resultPad);
        imagePan.add(dxy);
        imagePan.add(dz);

        imagePan.add((saveMetaData = new JButton("Save metadata")));

        saveMetaData.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (debug) {
                    System.out.println("Saving metadata");
                    updateMetaData();
                }
            }
        });

        dxy.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Sequence seq = getSequence(image);
                if (seq != null)  {              
                    setMetaData(seq) ;
                }
            }
        });

        ActionListener zeroPadActionListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                updateOutputPad();
            }
        };
        zeroPaddingxy.addActionListener(zeroPadActionListener);
        zeroPaddingz.addActionListener(zeroPadActionListener);

        dz.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Sequence seq = getSequence(image);
                if (seq != null)  {              
                    setMetaData(seq) ;
                }
            }
        });

        image.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Sequence seq = getSequence(image);
                if (seq == null) {
                    canalImageOptions = new String[]{"None"};
                } else {
                    //Update channel available in the image
                    int nbChan = seq.getSizeC();
                    canalImageOptions = new String[nbChan];
                    for (int i = 0; i < nbChan; i++) {
                        canalImageOptions[i] = seq.getChannelName(i);
                    }
                    canalImage.updateData(canalImageOptions);
                    //Update PSF metadata
                    if(seq.getSizeX() != seq.getSizeY()){
                        throwError("Height and width of the image must be equal ");
                    }
                    meta = getMetaData(seq);
                    dxy.setValue(    meta.dxy);
                    dz.setValue(     meta.dz);
                    nxy.setValue(    meta.nxy);
                    nz.setValue(     meta.nz);
                    na.setValue(     meta.na);
                    lambda.setValue( meta.lambda);
                    ni.setValue(     meta.ni);

                    updateOutputPad();
                }
            }
        });
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
        psfPannel.add((showPSF = new JButton("Show PSF")));

        ns =createDouble(     "<html><pre>ns:        </pre></html>", 0.0);
        zdepth = createDouble(     "<html><pre>ns:        </pre></html>", 0.0);



        psf.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                canRunBdec = psf.getValue().equals("None"); // If no psf given, we can run bdec tab, else we disable the run button in bdec tab
                na.setVisible(canRunBdec);                  // AND if we can run BDEC we can show the options of bdec
                ni.setVisible(canRunBdec);
                lambda.setVisible(canRunBdec);
                showPSF.setVisible(canRunBdec);
            }
        });

        showPSF.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // Show the initial PSF
                psfClicked();
                if (debug) {
                    System.out.println("First PSF compute");
                }
            }
        });

        //Creation of IMAGE TAB
        psfGlob.add(psfPannel, BorderLayout.NORTH);
        tabbedPane.addTab("PSF", null, psfGlob, "Choice of the PSF, visualization of theoritical PSF");

        /****************************************************/
        /**                 VARIANCE TAB                   **/
        /****************************************************/
        //Creation of the inside of VARIANCE TAB
        varianceGlob = new JPanel(new BorderLayout()); //Border layout to be sure that the images are stacked to the up
        JPanel varianceTab = new JPanel(false);
        varianceTab.setLayout(new BoxLayout(varianceTab, BoxLayout.Y_AXIS));
        varianceTab.add((weightsMethod = createChoiceList(  "<html><pre>Weighting:      </pre></html>", weightOptions)));
        varianceTab.add((weights = createChoiceList(        "<html><pre>Map:      </pre></html>", seqList)));
        varianceTab.add((gain = createDouble(               "<html><pre>Gain:             </pre></html>", 1.0)));
        varianceTab.add((noise = createDouble(              "<html><pre>Readout Noise:    </pre></html>", 1.0)));
        varianceTab.add((deadPixGiven = new MyBoolean(      "<html><pre>Data Map?        </pre></html>", false)));
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
                    throwError("Invalid argument passed to weight method");
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
        varianceTab.add((showWeight = new JButton("Show weight map")));	

        showWeight.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // Show the initial PSF
                showWeightClicked();
                if (debug) {
                    System.out.println("Weight compute");
                }
            }
        });
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
        deconvTab.add((mu = new MyDouble(               "<html><pre>Regularization level:             </pre></html>", 5E-4)));
        deconvTab.add((epsilon = new MyDouble(          "<html><pre>Threshold level:                  </pre></html>", 1E-2)));
        //deconvTab.add((zeroPadding = new MyDouble(      "<html><pre>Number of lines to add (padding): </pre></html>", 0)));
        deconvTab.add((nbIteration = new MyDouble(      "<html><pre>Number of iterations:             </pre></html>", 50)));
        deconvTab.add((positivity = new MyBoolean(      "<html><pre>Enforce nonnegativity:            </pre></html>", true)));
        deconvTab.add((crop = new MyBoolean(            "<html><pre>Crop output to match input:       </pre></html>", false)));
        deconvTab.add((restart = createChoiceList(      "<html><pre>Start from last result:           </pre></html>", seqList)));

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
        bdecTab.add((nbAlphaCoef = createChoiceList(    "<html><pre>Number of phase coefs N\u03B1:                          </pre></html>", nAlphaOptions)));
        bdecTab.add((nbBetaCoef = createChoiceList(     "<html><pre>Number of modulus coefs N\u03B2:                          </pre></html>", nBetaOptions)));    
        bdecTab.add((resetPSF = new MyBoolean(      "<html><pre>Reset PSF:            </pre></html>", true)));
        bdecTab.add((DefocusMaxIter = new MyDouble(       "<html><pre>Max. nb. of iterations for defocus:               </pre></html>", 0)));
        bdecTab.add((PhaseMaxIter = new MyDouble(         "<html><pre>Max. nb. of iterations for  phase:                 </pre></html>", 0)));
        bdecTab.add((ModulusMaxIter = new MyDouble(       "<html><pre>Max. nb. of iterations for  modulus:               </pre></html>", 0)));
        /*        bdecTab.add((grtolDefocus = new MyDouble(       "<html><pre>Grtol defocus:               </pre></html>", 0.001)));
        bdecTab.add((grtolPhase = new MyDouble(         "<html><pre>Grtol phase:                 </pre></html>", 0.001)));
        bdecTab.add((grtolModulus = new MyDouble(       "<html><pre>Grtol modulus:               </pre></html>", 0.001)));*/
        bdecTab.add((bDecTotalIteration = new MyDouble( "<html><pre>Number of loops:             </pre></html>", 1)));
        bdecTab.add((psfShow2 = new JButton(        "Show PSF"))); //Already created in psf tab
        bdecTab.add((showPhase = new JButton(       "Show phase of the pupil")));
        bdecTab.add((showModulus = new JButton(     "Show modulus of the pupil")));

        psfShow2.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // Show the initial PSF
                psfClicked();
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
        tabbedPane.addTab("BDec", null, bdecGlob,    "Options for the blind deconvolution");

        /****************************************************/
        /**                    RESULT TAB                  **/
        /****************************************************/
        //Creation of the inside of BDec TAB
        String empty = "                    ";
        resultGlob = new JPanel(new BorderLayout()); //Border layout to be sure that the images are stacked to the up
        JPanel resultTab = new JPanel(false);
        resultTab.setLayout(new BoxLayout(resultTab, BoxLayout.Y_AXIS));
        resultTab.add((resultCostPrior = new JLabel(    "<html><pre>"+empty+"No results yet   </pre></html>")));
        resultTab.add((resultDefocus = new JLabel(       "<html><pre>"+empty+"No results yet   </pre></html>")));
        if (debug) {
            resultTab.add((resultModulus = new JLabel(      "<html><pre>"+empty+"No results yet   </pre></html>")));
            resultTab.add((resultPhase = new JLabel(        "<html><pre>"+empty+"No results yet   </pre></html>")));
        }

        resultGlob.add(resultTab, BorderLayout.NORTH);
        tabbedPane.addTab("Results", null, resultGlob,    "Results of the deconvolution");

        /****************************************************/
        /**                      ToolTips                  **/
        /****************************************************/
        image.setToolTipText(ToolTipText.sequenceImage);
        psf.setToolTipText(ToolTipText.sequencePSF);
        weights.setToolTipText(ToolTipText.sequenceWeigth);
        weightsMethod.setToolTipText(ToolTipText.sequenceWeigth);
        deadPixel.setToolTipText(ToolTipText.sequencePixel);

        canalImage.setToolTipText(ToolTipText.textCanal);

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
        nbIteration.setToolTipText(ToolTipText.doubleMaxIter);
        zeroPaddingxy.setToolTipText(ToolTipText.doublePadding);
        zeroPaddingz.setToolTipText(ToolTipText.doublePadding);

        nbAlphaCoef.setToolTipText(ToolTipText.doubleNalpha);
        nbBetaCoef.setToolTipText(ToolTipText.doubleNbeta);
        //grtolPhase.setToolTipText(ToolTipText.doubleGrtolPhase);
        //grtolModulus.setToolTipText(ToolTipText.doubleGrtolModulus);
        //grtolDefocus.setToolTipText(ToolTipText.doubleGrtolDefocus);
        bDecTotalIteration.setToolTipText(ToolTipText.doubleBDecTotalIteration);

        restart.setToolTipText(ToolTipText.booleanRestart);
        positivity.setToolTipText(ToolTipText.booleanPositivity);
        crop.setToolTipText(ToolTipText.booleanCrop);
        resultTab.setToolTipText(ToolTipText.textOutput);
        // Adding image, canalImage, psf, weights, deadPixel to auto refresh when sequence is added/removed
        listChoiceList.add(image);
        listChoiceList.add(psf);
        listChoiceList.add(weights);
        listChoiceList.add(deadPixel);
        listChoiceList.add(restart);
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

        token = new ReconstructionThreadToken(new double[]{mu.getValue(),epsilon.getValue(),0.0,grtol});
        thread = new ReconstructionThread(token);
        thread.start();
    }

    public void launchDeconvolution(DoubleArray imgArray, DoubleArray psfArray, DoubleArray weight){
        if (tvDec == null) {
            tvDec = new TotalVariationJobForIcy(token);
            tvDec.setResult(null);
        }
        tvDec.setAbsoluteTolerance(0.0);
        tvDec.setWeight(weight);
        tvDec.setData(imgArray);
        tvDec.setPsf(psfArray);
        tvDec.setViewer(new TvViewer());
        thread.setJob(tvDec);

        if (restart.getValue() != "None") {
            Sequence restartSeq = getSequence(restart);
            // We verify that the previous result is conform to our expectations: !Null and same dim as input
            if (restartSeq != null) {
                DoubleArray tmpDoubleArray = (DoubleArray) IcyBufferedImageUtils.imageToArray(restartSeq, 0);
                for (int i = 0; i < imgArray.getOrder(); i++) {
                    if (imgArray.getOrder() != tmpDoubleArray.getOrder() || imgArray.getDimension(i) != tmpDoubleArray.getDimension(i)) {
                        throwError("The previous result does not have the same dimensions as the input image");
                    }
                }
                tvDec.setResult(tmpDoubleArray);
            }
        } else {
            tvDec.setResult(null);
        }
        tvDec.setPositivity(positivity.getValue());
        tvDec.setRegularizationWeight(mu.getValue());
        tvDec.setRegularizationThreshold(epsilon.getValue());
        tvDec.setRelativeTolerance(grtol);
        tvDec.setMaximumIterations((int)nbIteration.getValue());
        tvDec.setOutputShape(shapePad);
        token.start();
        setResult(tvDec);
    }

    @Override
    protected void execute() {
        try {
            if (debug) {
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
                System.out.println("zeroPad xy: "+zeroPaddingxy.getValue());
                System.out.println("zeroPad z: "+zeroPaddingz.getValue());
                System.out.println("nbIter: "+nbIteration.getValue());
                System.out.println("Restart: "+restart.getValue());
                System.out.println("Positivity: "+positivity.getValue());
                System.out.println("--------------BDEC------------------");
                System.out.println("nbIter: "+nbIteration.getValue());
                System.out.println("zeroPad: "+zeroPaddingxy.getValue());
                /*System.out.println("nbIterZern: "+grtolPhase.getValue());
                System.out.println("module: "+grtolModulus.getValue());
                System.out.println("defoc: "+grtolDefocus.getValue());*/
                System.out.println("Number of total iterations: "+bDecTotalIteration.getValue());
                System.out.println("------------------------------------");
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

            widthPad  = FFTUtils.bestDimension((int)(width + zeroPaddingxy.getValue()));
            heightPad = FFTUtils.bestDimension((int)(height + zeroPaddingxy.getValue()));
            sizeZPad  = FFTUtils.bestDimension((int)(sizeZ + zeroPaddingz.getValue()));
            if (sizeZ == 1) {
                throwError("Input data must be 3D");
            } 
            if (zeroPaddingxy.getValue() < 0.0) {
                throwError("Padding value cannot be negative");
            }
            if (zeroPaddingz.getValue() < 0.0) {
                throwError("Padding value cannot be negative");
            }
            //double coef = (width + zeroPadding.getValue())/width;

            shapePad = Shape.make(widthPad, heightPad, sizeZPad);    
            shape    = Shape.make(width, height, sizeZ);

            int numCanal = getNumCanal(imgSeq);

            DoubleArray imgArray, psfArray;
            runBdec = (tabbedPane.getSelectedComponent() == bdecGlob); //If the BDEC panel is selected we the blind deconvolution
            // If no PSF is loaded -> creation of a PSF
            if (psfSeq == null) {
                if (pupil == null) {
                    buildpupil();
                }
                psfArray =  Double3D.wrap(MathUtils.fftShift3D(pupil.getPSF(), width, height, sizeZ) , shape);
            } else {
                psfArray = (DoubleArray) IcyBufferedImageUtils.imageToArray(psfSeq, Shape.make(psfSeq.getWidth(), psfSeq.getHeight(), psfSeq.getSizeZ()), numCanal);
            }

            imgArray = (DoubleArray) IcyBufferedImageUtils.imageToArray(imgSeq, shape, numCanal);

            DoubleArray weight = createWeight(imgArray).toDouble();

            psfArray = (DoubleArray) ArrayUtils.pad(psfArray, shapePad);
            imgArray = (DoubleArray) ArrayUtils.pad(imgArray, shapePad);
            weight   = (DoubleArray) ArrayUtils.pad(weight  , shapePad);

            //BEWARE here we change the value to match the new padded image size
            // FIXME we pad every input so for nown we use xxxPad

            /*---------------------------------------*/
            /*            OPTIMISATION               */
            /*---------------------------------------*/

            if (runBdec) {
                if (Integer.parseInt(nbAlphaCoef.getValue()) != nbAlpha){
                    nbAlpha = Integer.parseInt(nbAlphaCoef.getValue());
                    alphaSpace = new DoubleShapedVectorSpace(new int[]{nbAlpha});
                    alphaVector = alphaSpace.create();
                }
                if (Integer.parseInt(nbBetaCoef.getValue()) != nbBeta){
                    nbBeta = Integer.parseInt(nbBetaCoef.getValue());
                    double[] beta = new double[nbBeta];
                    beta[0] = 1;
                    betaSpace = new DoubleShapedVectorSpace(new int[]{beta.length});
                    betaVector = betaSpace.wrap(beta);
                }
                if(resetPSF.getValue())
                {
                    defocuSpace = null;
                    defocusVector= null;
                    alphaSpace = null;
                    alphaVector = null;
                    betaSpace = null;
                    betaVector = null;
                }
                if (defocuSpace==null){
                    double[] defocus = {ni.getValue()/lambda.getValue(), 0., 0.};
                    defocuSpace = new DoubleShapedVectorSpace(new int[]{defocus.length});
                    defocusVector = defocuSpace.wrap(defocus);
                }
                if (alphaSpace==null){
                    alphaSpace = new DoubleShapedVectorSpace(new int[]{nbAlpha});
                    alphaVector = alphaSpace.create();
                }
                if (betaSpace==null){
                    double[] beta = new double[Integer.parseInt(nbBetaCoef.getValue())];
                    beta[0] = 1;
                    betaSpace = new DoubleShapedVectorSpace(new int[]{beta.length});
                    betaVector = betaSpace.wrap(beta);
                }

                PSFEstimation = new PSF_Estimation();
                PSFEstimation.setWeight(weight);
                PSFEstimation.setData(imgArray);

                PSFEstimation.enablePositivity(false);
                PSFEstimation.setAbsoluteTolerance(0.0);

                for(int i = 0; i < bDecTotalIteration.getValue(); i++) {
                    /* OBJET ESTIMATION (by the current PSF) */
                    launchDeconvolution(imgArray, psfArray, weight);
                    /* PSF ESTIMATION (by the current objet) */
                    PSFEstimation.setPupil(pupil);
                    //   PSFEstimation.setPsf(tvDec.getData());
                    PSFEstimation.setObj(tvDec.getResult());

                    /* Defocus estimation */
                    if (DefocusMaxIter.getValue()!=0){
                        if (debug && verbose) {
                            System.out.println("------------------");
                            System.out.println("Defocus estimation");
                            System.out.println("------------------");
                        }
                        PSFEstimation.setRelativeTolerance(0.);
                        PSFEstimation.setMaximumIterations((int)DefocusMaxIter.getValue());
                        PSFEstimation.fitPSF(defocusVector, PSF_Estimation.DEFOCUS);
                    }

                    /* Phase estimation */
                    if(PhaseMaxIter.getValue()!=0){
                        if (debug && verbose) {
                            System.out.println("Phase estimation");
                            System.out.println("------------------");
                        }
                        PSFEstimation.setResult(null);                    
                        PSFEstimation.setMaximumIterations((int)PhaseMaxIter.getValue());
                        PSFEstimation.fitPSF(alphaVector, PSF_Estimation.ALPHA);
                    }


                    /* Modulus estimation */
                    if(ModulusMaxIter.getValue() !=0){
                        if (debug && verbose) {
                            System.out.println("Modulus estimation");
                            System.out.println("------------------");
                        }
                        PSFEstimation.setResult(null);           
                        PSFEstimation.setMaximumIterations((int)ModulusMaxIter.getValue());
                        PSFEstimation.fitPSF(betaVector, PSF_Estimation.BETA);
                        // MathUtils.normalise(betaVector.getData());
                    }
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
            sequence = null; //In any cases the next image will be in a new sequence
        } catch (IllegalArgumentException e) {
            new AnnounceFrame("Oops, Error: "+ e.getMessage());
            if (debug) {
                e.printStackTrace();
            }
        }
    }

    /**
     * The goal is to create an weight array, but it will be created depending
     * the user input so we will have to test each cases:
     *  	-None
     *  	-A given map
     *  	-A variance map
     *  	-A computed variance
     * Then we apply the dead pixel map
     * 
     * @param data
     * @return
     */
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

    /**
     * Get the results from a reconstruction and plot the intermediate result
     * 
     * @param tvDec
     */
    private void setResult(ReconstructionJob tvDec){
        try{
            //Here we will update the sequence
            if (sequence == null || (sequence != null && sequence.isEmpty())) {
                sequence = new Sequence();
                setMetaData(getSequence(image), sequence);
                addSequence(sequence);
            }
            sequence.beginUpdate();
            if (crop.getValue()) {

                ShapedArray croppedArray = ArrayUtils.crop(tvDec.getResult(), shape);
                double[] in = croppedArray.toDouble().flatten();
                for (int j = 0; j < sizeZ; j++) {
                    double[] temp = new double[width*height];
                    for (int i = 0; i < width*height; i++) {
                        temp[i] = in[i+j*width*height];
                    }
                    sequence.setImage(0,j, new IcyBufferedImage(width, height, temp));
                }
            } else {
                double[] in = tvDec.getResult().toDouble().flatten();
                for (int j = 0; j < sizeZPad; j++) {
                    double[] temp = new double[widthPad*heightPad];
                    for (int i = 0; i < widthPad*heightPad; i++) {
                        temp[i] = in[i+j*widthPad*heightPad];
                    }
                    sequence.setImage(0,j, new IcyBufferedImage(widthPad, heightPad, temp));
                }
            }

            sequence.endUpdate();
            sequence.setName("TV mu:"+mu.getValue()+" Iteration:"+tvDec.getIterations());
            update();
            //Then we will update the result tab pannel
            if (runBdec) {
                String empty = "      ";
                //               resultCostData.setText( "<html><pre>"+empty+"FCostData  "+tvDec.getCost()                       +"</pre></html>");
                resultCostPrior.setText("<html><pre>"+empty+"Cost "+tvDec.getCost()                       +"</pre></html>");
                resultDefocus.setText(   "<html><pre>"+empty+"Defocus   "+Arrays.toString(pupil.getDefocusMultiplyByLambda())   +"</pre></html>");
                if (debug) {
                    resultModulus.setText(  "<html><pre>"+empty+"Modulus    "+pupil.getRho()[0]                     +"</pre></html>");
                    resultPhase.setText(    "<html><pre>"+empty+"Phase      "+pupil.getPhi()[0]                     +"</pre></html>");
                }
            }
        } catch (NullPointerException e) {
            //Here in case of brutal stop the sequence can become null but it's not important as it's an emergency stop
            //So we do nothing
            System.out.println("INFO: Emergency stop detected in setResult");
            e.printStackTrace();
        }
    }

    /**
     * A debug function linked to the PSF
     * 
     * @param num Number of PSF already generated
     */
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

    /*****************************************/
    /** All the PSF buttons call are here   **/
    /*****************************************/

    private void buildpupil()
    {
        pupil = new WideFieldModel(na.getValue(), lambda.getValue(), ni.getValue(), ns.getValue(), zdepth.getValue(), dxy.getValue()*1E-9,
                dz.getValue()*1E-9, (int)nxy.getValue(), (int)nxy.getValue(), (int)nz.getValue(), use_depth_scaling);
    }	

    private void psfClicked()
    {
        /* PSF0 initialisation */
        if(pupil==null)
        {
            buildpupil();
        }

        /* PSF0 Sequence */
        Sequence PSF0Sequence = new Sequence();
        double[] PSF_shift = MathUtils.fftShift3D(pupil.getPSF(), (int)nxy.getValue(), (int)nxy.getValue(), (int)nz.getValue());
        for (int k = 0; k < (int)nz.getValue(); k++)
        {
            PSF0Sequence.setImage(0, k, new IcyBufferedImage((int)nxy.getValue(), (int)nxy.getValue(),
                    MathUtils.getArray(PSF_shift, (int)nxy.getValue(), (int)nxy.getValue(), k)));
        }
        setMetaData(PSF0Sequence) ;
        PSF0Sequence.setName("PSF");
        addSequence(PSF0Sequence);
    }

    private void phaseClicked()
    {
        /* PSF0 initialisation */ 
        if(pupil==null)
        {
            buildpupil();
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
        if(pupil==null)
        {
            buildpupil();
        }
        /* Modulus Sequence */
        Sequence modulusSequence = new Sequence();
        modulusSequence.setName("Modulus of the pupil");
        double[] modulus_shift = MathUtils.fftShift1D(pupil.getRho(), (int)nxy.getValue(), (int)nxy.getValue());
        modulusSequence.addImage(new IcyBufferedImage((int)nxy.getValue(), (int)nxy.getValue(), modulus_shift));
        addSequence(modulusSequence);
    }

    private void showWeightClicked()
    {

        //    DoubleArray weight = createWeight(...).toDouble();
        /* PSF0 Sequence */

        Sequence img = getSequence(image);
        if (img == null) {
            new AnnounceFrame("No image was chosen");
            return;
        }
        Shape myShape; 
        if (img.getSizeZ() == 1) { //2D
            myShape = Shape.make(img.getSizeX(), img.getSizeY());
        } else {    //3D
            myShape = Shape.make(img.getSizeX(), img.getSizeY(), img.getSizeZ());
        }
        width = img.getSizeX();
        height = img.getSizeY();
        sizeZ = img.getSizeZ();

        int numCanal = getNumCanal(img);
        DoubleArray input = (DoubleArray) IcyBufferedImageUtils.imageToArray(img, myShape, numCanal);
        Sequence WeightSequence = new Sequence();
        WeightSequence.setName("Weight");
        double[] wght = createWeight(input).toDouble().flatten();
        for (int k = 0; k < (int)nz.getValue(); k++)
        {
            WeightSequence.setImage(0, k, new IcyBufferedImage((int)nxy.getValue(), (int)nxy.getValue(),
                    MathUtils.getArray(wght, (int)nxy.getValue(), (int)nxy.getValue(), k)));
        }
        addSequence(WeightSequence);
    }

    /**
     * Here we get the informations given by the users but not all.
     * In fact we trust only a few data that we know that are given by Icy.
     * Else we are trying to keep them for the next run.
     * 
     * Remember: if users may lie, they will !
     * 
     * @param seq
     * @return
     */
    private MyMetadata getMetaData(Sequence seq){
        OMEXMLMetadata metDat = seq.getMetadata();
        if (meta == null) {
            meta = new MyMetadata();
            if (metDat.getInstrumentCount() > 0) {
                try {
                    meta.na      = metDat.getObjectiveLensNA(0, 0);
                    //meta.lambda  = metDat.getChannelEmissionWavelength(0, 0).getValue().doubleValue()*1E6;  //I suppose the value I will get is in um
                } catch(Exception e){
                    System.out.println("Failed to get some metadatas, will use default values for na, lambda");
                }
            } else {
                if (debug && verbose) {
                    System.out.println("INFO: Metadata: No instrument so no metadata.");
                }
            }
        }
        //If no instrument found, at least we have the right image size
        meta.nxy     = seq.getSizeX(); //We suppose X and Y equal
        meta.nz      = seq.getSizeZ();
        meta.dxy     = seq.getPixelSizeX()*1E3;
        meta.dz      = seq.getPixelSizeZ()*1E3;
        meta.na      = na.getValue();
        meta.lambda  = lambda.getValue(false);
        meta.ni      = ni.getValue();
        return meta;
    }

    //Copy the input metadata to the output. We may want to change some with our values
    //So it should be done here
    private void setMetaData(Sequence seqOld, Sequence seqNew) {
        OMEXMLMetadataImpl newMetdat = OMEUtil.createOMEMetadata(seqOld.getMetadata());
        //newMetdat.setImageDescription("MyDescription", 0);
        newMetdat.setPixelsPhysicalSizeX(OMEUtil.getLength(dxy.getValue()*1E-3), 0);
        newMetdat.setPixelsPhysicalSizeY(OMEUtil.getLength(dxy.getValue()*1E-3), 0);
        newMetdat.setPixelsPhysicalSizeZ(OMEUtil.getLength(dz.getValue()*1E-3), 0);
        seqNew.setMetaData(newMetdat);
    }


    //Copy the input metadata to the output. We may want to change some with our values
    //So it should be done here
    private void setMetaData(Sequence seqNew) {
        OMEXMLMetadataImpl newMetdat = OMEUtil.createOMEMetadata();
        //newMetdat.setImageDescription("MyDescription", 0);
        newMetdat.setPixelsPhysicalSizeX(OMEUtil.getLength(dxy.getValue()*1E-3), 0);
        newMetdat.setPixelsPhysicalSizeY(OMEUtil.getLength(dxy.getValue()*1E-3), 0);
        newMetdat.setPixelsPhysicalSizeZ(OMEUtil.getLength(dz.getValue()*1E-3), 0);
        seqNew.setMetaData(newMetdat);
    }

    private void updateMetaData() {
        Sequence seq = getSequence(image);
        try {
            OMEXMLMetadata newMetdat = MetaDataUtil.generateMetaData(seq, false);
            newMetdat.setPixelsPhysicalSizeX(OMEUtil.getLength(dxy.getValue()*1E-3), 0);
            newMetdat.setPixelsPhysicalSizeY(OMEUtil.getLength(dxy.getValue()*1E-3), 0);
            newMetdat.setPixelsPhysicalSizeZ(OMEUtil.getLength(dz.getValue()*1E-3), 0);
            seq.setMetaData((OMEXMLMetadataImpl) newMetdat);
        } catch (ServiceException e) {
            e.printStackTrace();
        }

    }

    /**
     * This function update all the names of the available sequences contains 
     * in the MyComboBox AND that have been added to the update list
     */
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
