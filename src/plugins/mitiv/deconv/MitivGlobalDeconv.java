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
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;

import ome.xml.model.primitives.PositiveInteger;
import commands.TotalVariationDeconvolution;
import loci.formats.ome.OMEXMLMetadata;
import loci.formats.ome.OMEXMLMetadataImpl;
import mitiv.array.Double1D;
import mitiv.array.Double3D;
import mitiv.array.DoubleArray;
import mitiv.array.ShapedArray;
import mitiv.invpb.ReconstructionJob;
import mitiv.invpb.ReconstructionViewer;
import mitiv.io.IcyBufferedImageUtils;
import mitiv.linalg.WeightGenerator;
import mitiv.utils.CommonUtils;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;

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
        public double dxy     = 1;
        public double dz      = 2;
        public double nxy     = 3;
        public double nx      = 4;
        public double no      = 5;
        public double lambda  = 6;
        public double ni      = 7;
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

        public myDouble(String name, double input) {
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
            return Double.valueOf(field.getText());
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
    private myDouble dxy, dz, nxy, nx, no,lambda,ni;    //PSF
    private myDouble mu, epsilon, grtol, nbIteration, zeroPadding;          //Deconvolution
    private myDouble gain,noise;                          //VARIANCE
    private myDouble nbIterationZern, moduleZern, defocusZern, loopZern;          //BDec
    private myComboBox image, canalImage, psf, weights, deadPixel, deconvOptions;
    private myBoolean deadPixGiven, restart;
    private String[] seqList;
    private final String[] deconvStringOptions = new String[]{"Total Variation","Tichonov"}; 
    private String[] canalImageOptions = new String[]{"None"}; 

    //Global variables for the algorithms
    TotalVariationDeconvolution tvDec;
    
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
        JTabbedPane tabbedPane = new JTabbedPane();

        /****************************************************/
        /**                    IMAGE TAB                   **/
        /****************************************************/
        //Creation of the inside of IMAGE TAB
        JPanel imageGlob = new JPanel(new BorderLayout()); //Border layout to be sure that the images are stacked to the up
        JPanel imagePan = new JPanel(false);
        imagePan.setLayout(new BoxLayout(imagePan, BoxLayout.Y_AXIS));
        imagePan.add((image = createChoiceList( "Sequence: ", seqList)));
        imagePan.add((canalImage = createChoiceList("Canal:    ", canalImageOptions)));
        image.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Sequence seq = getSequence(image);
                if (seq == null) {
                    canalImageOptions = new String[]{"None"};
                } else {
                    int nbChan = seq.getSizeC();
                    canalImageOptions = new String[nbChan];
                    for (int i = 0; i < nbChan; i++) {
                        canalImageOptions[i] = seq.getChannelName(i);
                    }
                    canalImage.updateData(canalImageOptions);
                }
            }
        });
        //Little hack to have a minimal size at startup
        JLabel sizeTab = new JLabel("                                                                    ");
        imagePan.add(sizeTab);

        //Creation of IMAGE TAB
        imageGlob.add(imagePan, BorderLayout.NORTH);
        tabbedPane.addTab("FILE", null, imageGlob, "Selecting input data");

        /****************************************************/
        /**                    PSF TAB                     **/
        /****************************************************/
        //Creation of the inside of PSF TAB
        JPanel psfGlob = new JPanel(new BorderLayout()); //Border layout to be sure that the images are stacked to the up
        JPanel psfPannel = new JPanel(false);
        psfPannel.setLayout(new BoxLayout(psfPannel, BoxLayout.Y_AXIS));
        psfPannel.add((psf = createChoiceList("PSF:       ", seqList)));
        psfPannel.add((dxy = createDouble(    "dxy:       ", 0.0)));
        psfPannel.add((dz = createDouble(     "dz:        ", 1.0)));
        psfPannel.add((nxy = createDouble(    "Nxy:       ", 2.0)));
        psfPannel.add((nx = createDouble(     "Nx:        ", 3.0)));
        psfPannel.add((no = createDouble(     "Nteta:     ", 4.0)));
        psfPannel.add((lambda = createDouble( "lambda:    ", 5.0)));
        psfPannel.add((ni = createDouble(     "ni:        ", 6.0)));

        psf.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Sequence seq = getSequence(psf);
                if (seq != null) {
                    myMetaData meta = getMetaData(seq);
                    dxy.setValue(    meta.dxy);
                    dz.setValue(     meta.dz);
                    nxy.setValue(    meta.nxy);
                    nx.setValue(     meta.nx);
                    no.setValue(     meta.no);
                    lambda.setValue( meta.lambda);
                    ni.setValue(     meta.ni);
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
        JPanel varianceGlob = new JPanel(new BorderLayout()); //Border layout to be sure that the images are stacked to the up
        JPanel varianceTab = new JPanel(false);
        varianceTab.setLayout(new BoxLayout(varianceTab, BoxLayout.Y_AXIS));
        varianceTab.add((weights = createChoiceList(  "Variance MAP:     ", seqList)));
        varianceTab.add(createLabel(                  "Variance Map:     "));
        varianceTab.add((gain = createDouble(         "Gain:             ", 0.0)));
        varianceTab.add((noise = createDouble(        "Readout Noise:    ", 0.0)));
        varianceTab.add((deadPixGiven = new myBoolean("Dead Pixel Map ?  ", false)));
        deadPixel = createChoiceList(                 "Dead pixel Map:   ", seqList);
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
        JPanel deconvGlob = new JPanel(new BorderLayout()); //Border layout to be sure that the images are stacked to the up
        JPanel deconvTab = new JPanel(false);
        deconvTab.setLayout(new BoxLayout(deconvTab, BoxLayout.Y_AXIS));
        deconvTab.add((deconvOptions = createChoiceList("Method:                  ", deconvStringOptions)));
        deconvTab.add((mu = new myDouble(               "Mu:                      ", 1.0)));
        deconvTab.add((epsilon = new myDouble(          "Epsilon:                 ", 1.0)));
        deconvTab.add((grtol = new myDouble(            "Grtol:                   ", 1.0)));
        deconvTab.add((zeroPadding = new myDouble(      "Padding multiplication:  ", 1.0)));
        deconvTab.add((nbIteration = new myDouble(      "Number of iterations:    ", 50)));
        deconvTab.add((restart = new myBoolean(         "Start from last result:  ", false)));
        //Creation of DECONVOLUTION TAB
        deconvGlob.add(deconvTab, BorderLayout.NORTH);
        tabbedPane.addTab("Deconvolution", null, deconvGlob, "Methods and options of the deconvolution");

        /****************************************************/
        /**                      BDEC TAB                  **/
        /****************************************************/
        //Creation of the inside of BDec TAB
        JPanel bdecGlob = new JPanel(new BorderLayout()); //Border layout to be sure that the images are stacked to the up
        JPanel bdecTab = new JPanel(false);
        bdecTab.setLayout(new BoxLayout(bdecTab, BoxLayout.Y_AXIS));
        bdecTab.add((nbIterationZern = new myDouble("Zernike iterations:           ", 50)));
        bdecTab.add((moduleZern      = new myDouble("Modulus:                      ", 50)));
        bdecTab.add((defocusZern     = new myDouble("Defocus:                      ", 50)));
        bdecTab.add((loopZern        = new myDouble("Loop:                         ", 50)));
        //Creation of BDec TAB
        bdecGlob.add(bdecTab, BorderLayout.NORTH);
        tabbedPane.addTab("BDec", null, bdecGlob,    "All the options for the blind deconvolution");

        // Addingimage, canalImage, psf, weights, deadPixel to auto refresh when sequence is added/removed
        listChoiceList.add(image);
        listChoiceList.add(psf);
        listChoiceList.add(weights);
        listChoiceList.add(deadPixel);
        //Add the tabbed pane to this panel.
        addComponent(tabbedPane);
    }
    
    @Override
    protected void execute() {
        System.out.println("-------------IMAGE-------------------");
        System.out.println("File: "+image.getValue());              //Used
        System.out.println("Canal: "+canalImage.getValue());
        System.out.println("--------------PSF------------------");
        System.out.println("PSF: "+psf.getValue());                 //Used
        System.out.println("dxy: "+dxy.getValue());
        System.out.println("dz: "+dz.getValue());
        System.out.println("nxy: "+nxy.getValue());
        System.out.println("nx: "+nx.getValue());
        System.out.println("no: "+no.getValue());
        System.out.println("l: "+lambda.getValue());
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
        System.out.println("nbIterZern: "+nbIterationZern.getValue());
        System.out.println("module: "+moduleZern.getValue());
        System.out.println("defoc: "+defocusZern.getValue());
        System.out.println("loop: "+loopZern.getValue());
        System.out.println("");
        
        // Preparing parameters and testing input
        Sequence imgSeq = getSequence(image);
        Sequence psfSeq = getSequence(psf);
        if (imgSeq == null || psfSeq == null) {
            throwError("A PSF/Image should be given");
        }
        ArrayList<IcyBufferedImage> listImg = imgSeq.getAllImage();
        ArrayList<IcyBufferedImage> listPSf = psfSeq.getAllImage();

        BufferedImage img = imgSeq.getFirstNonNullImage();
        // Set the informations about the input
        width = img.getWidth();
        height = img.getHeight();
        sizeZ = imgSeq.getSizeZ();
        double coef = zeroPadding.getValue();
        //3D Only
        DoubleArray imgArray, psfArray;
        double[] image = CommonUtils.icyImage3DToArray1D(listImg, width, height, sizeZ, false);
        double[] psfTmp = CommonUtils.icyImage3DToArray1D(listPSf, width, height, sizeZ, false);
        double[] weight = createWeight(image);
        weight = CommonUtils.imagePad(weight, width, height, sizeZ, coef);
        image = CommonUtils.imagePad(image, width, height, sizeZ, coef);
        psfTmp = CommonUtils.imagePad(psfTmp, width, height, sizeZ, coef);
        int[] shape = new int[]{(int)(width*coef), (int)(height*coef), (int)(sizeZ*coef)};

        imgArray =  Double3D.wrap(image, shape);
        double[] psfShift = new double[shape[0]*shape[1]*shape[2]];
        CommonUtils.fftShift3D(psfTmp, psfShift , shape[0], shape[1], shape[2]);
        psfArray =  Double3D.wrap(psfShift, shape);

        //BEWARE here we change the value to match the new padded image size
        width = (int)(width*coef);
        height = (int)(height*coef);
        sizeZ = (int)(sizeZ*coef);

        // Launching the method
        if (deconvOptions.getValue() == deconvStringOptions[0]) { //Total variation
            if (tvDec != null &&  restart.getValue()) {
                tvDec.setResult(tvDec.getResult());
            } else {
                tvDec = new TotalVariationDeconvolution();
                tvDec.setAbsoluteTolerance(0.0);
                tvDec.setWeight(weight);
                tvDec.setData(imgArray);
                tvDec.setPsf(psfArray);
                tvDec.setViewer(new tvViewer());
            }
            tvDec.setRegularizationWeight(mu.getValue());
            tvDec.setRegularizationThreshold(epsilon.getValue());
            tvDec.setRelativeTolerance(grtol.getValue());
            tvDec.setMaximumIterations((int)nbIteration.getValue());
            //Computation HERE
            tvDec.deconvolve();
            //Showing the results
            setResult(tvDec);
            
        } else if (deconvOptions.getValue() == deconvStringOptions[1]) { //tichonov
            System.out.println("Ticho ON");
        } else {
            throw new IllegalArgumentException("Unknow deconvolution option");
        }
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
        return weightGen.getWeightMap().toDouble().flatten();//FIXME WRONG but Good for now
    }

    //Get the results from a reconstruction and plot the intermediate result
    private void setResult(ReconstructionJob tvDec){
        if (sequence == null) {
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
        setMetaData(sequence);
        sequence.setName("TV mu:"+mu.getValue()+" Iteration:"+tvDec.getIterations()+" grToll: "+tvDec.getRelativeTolerance());
        sequence.endUpdate();
    }

    private myMetaData getMetaData(Sequence seq){
        OMEXMLMetadata metDat = seq.getMetadata();
        myMetaData data = new myMetaData();
        //TODO Add dxy, dz, nxy, nx, no,lambda,ni
        //dxy lateral pixel size
        //dz pixel size axial
        //nxy nb Pixel W et H
        //nx nb Pixels
        //NA Numerical Aperture
        //lambda getChannelEmissionWavelength
        //ni getObjectiveImmersion
        data.dxy = metDat.getPixelsSizeC(0).getValue().doubleValue();
        return data;
    }
    
    private void setMetaData(Sequence seq){
        OMEXMLMetadataImpl metDat = seq.getMetadata();
        //TODO Add dxy, dz, nxy, nx, no,lambda,ni
        metDat.setPixelsSizeC(new PositiveInteger((int)dxy.getValue()),0);
        seq.setMetaData(metDat);
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
        
    }

}
