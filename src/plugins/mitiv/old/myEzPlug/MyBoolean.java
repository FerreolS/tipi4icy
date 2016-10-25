package plugins.mitiv.old.myEzPlug;

import java.awt.FlowLayout;
import java.awt.event.ActionListener;

import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;

public class MyBoolean extends JPanel{
    private static final long serialVersionUID = 1L;
    private JLabel filler;
    private JCheckBox box;

    public MyBoolean(String name, boolean defaultValue) {
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