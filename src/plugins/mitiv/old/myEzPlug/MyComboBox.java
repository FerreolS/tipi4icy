package plugins.mitiv.old.myEzPlug;

import java.awt.FlowLayout;
import java.awt.event.ActionListener;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;

public class MyComboBox extends JPanel{
    private static final long serialVersionUID = 1L;
    private JLabel filler;
    private JComboBox jcb;

    public MyComboBox(String name, String[] inputs){
        filler = new JLabel(name);
        setLayout(new FlowLayout());
        add(filler);
        jcb = new JComboBox();
        for (int j = 0; j < inputs.length; j++) {
            jcb.addItem(inputs[j]);
        }
        add(jcb);
    }

    public void updateData(String[] newInputs){             //FIXME Pushing pixel exception try using swing: SwingUtilities.invokeLater
        jcb.removeAllItems();
        for (int j = 0; j < newInputs.length; j++) {
            jcb.addItem(newInputs[j]);
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
            String tmp = (String) jcb.getItemAt(i);
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