package plugins.mitiv.myEzPlug;

import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

public class MyDouble extends JPanel{
    private static final long serialVersionUID = 1L;
    private JLabel filler;
    private JTextField field;
    private double mult;
    public MyDouble(String name, double input) {
        this(name, input, 1.0);
    }

    //Here we give the the multiplication factor, the result will be multiply by this factor
    public MyDouble(String name, double input, double valueMult) {
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

    public void addActionListener(ActionListener l){
        field.addActionListener(l);
    }
}