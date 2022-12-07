/*
 * TreeSpacePathViewer.java

    Copyright (C) 2013  Tom M. W. Nye

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact the author at:  <tom.nye@ncl.ac.uk>
                            <http://www.mas.ncl.ac.uk/~ntmwn/>
 */

package geodesicgraphics;

/**
 * Provide a view and a slider for viewing a TreeSpacePath.
 */


import treebase.TreeAsSplits;
import geodesics.Geodesic;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import treegraphics.*;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JMenuItem;


public class TreeSpacePathViewer extends JPanel implements javax.swing.event.ChangeListener, java.awt.event.ActionListener {

    /* Instance variables */
    private double[] lambdaVals;
//    private int sliderMiddle;
    protected int tickSpace;
    private TreeSpacePath thePath;
    private PathPlotter thePlotter;

    /* GUI components */
    protected TreeViewer graphicsPane;
    private JSlider theSlider;
    private JMenuItem printItem;

    /** Constructor.
     Inputs:
     - a TreeSpacePath
     - number of points on slider  */
    public TreeSpacePathViewer(TreeSpacePath p, int numFrames) {
        super();
        tickSpace = (int)Math.round(0.05*numFrames);
        thePath = p;

        // Sort out max and min values on the slider
        double[] b = thePath.getPathBounds();
        double minLambda = b[0];
        double maxLambda = b[1];

        double width = (maxLambda-minLambda)/numFrames;

        lambdaVals = new double[numFrames];
         for (int i=0; i<numFrames; i++) {
            lambdaVals[i] = i*width;
        }

        // Create the plotter: object responsible for drawing individual trees
        thePlotter = new PathPlotter(thePath, lambdaVals);

        // Deal with the GUI

        //Make the slider
        theSlider = new JSlider(JSlider.HORIZONTAL, 0, numFrames-1, 0);
        theSlider.setMajorTickSpacing(tickSpace);
        theSlider.setPaintTicks(true);
        theSlider.addChangeListener(this);

        // Create the graphics pane
        graphicsPane = new TreeViewer(thePlotter);
        printItem = new JMenuItem("Print...");
        printItem.addActionListener(this);
        graphicsPane.addToBackgroundMenu(printItem);

        // Finish off layout
        this.setLayout(new javax.swing.BoxLayout(this, javax.swing.BoxLayout.Y_AXIS));
        this.add(graphicsPane);
        this.add(theSlider);
        graphicsPane.refreshView();
        setSize(600,400);
        setVisible(true);
    }

    /** Deal with slider changes */
    public void stateChanged(javax.swing.event.ChangeEvent e) {
//        if (!theSlider.getValueIsAdjusting()) {
            int k = theSlider.getValue();
            thePlotter.changeTree(k);
            graphicsPane.refreshView();
 //       }
    }

    /** Get current tree */
    public TreeAsSplits getCurrentTree() {
        int k = theSlider.getValue();
        TreeAsSplits t = thePath.getTreeOnPath(lambdaVals[k]);
        return t;
    }

    /** Handle clicks on the menu items */
    public void actionPerformed(java.awt.event.ActionEvent e){
        graphicsPane.actionPerformed(e);
        JMenuItem source = (JMenuItem)(e.getSource());
        if (source == printItem) {
            graphicsPane.print();
        }
    }

    /** Cheeky menu stuff -- add an item.
     Useful for customizing the viewer. */
    public void addToBackgroundMenu(JMenuItem theItem) {
        graphicsPane.addToBackgroundMenu(theItem);
    }

    public static TreeSpacePathViewer viewGeodesic(Geodesic g) {
        TreeSpacePathViewer theViewer = new TreeSpacePathViewer(g, 50);

        //Create and set up window.
        javax.swing.JFrame frame = new javax.swing.JFrame("Geodesic Viewer");
        frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);

        //Set up the content pane.
//        theViewer.setOpaque(true); //content panes must be opaque???
        frame.setContentPane(theViewer);

        // Set position and size
        frame.setLocation(100,100);

        //Display the window.
        frame.pack();
        frame.setVisible(true);

        return theViewer;
    }

    
    /* Application
    At the command line, pass in the name of a file containing two Newick strings. 
    Construct the BHV geodesic between these and create a TreeSpacePathViewer using this. */

    public static void main(String[] args) throws treebase.AlgorithmException, java.io.IOException {

        /* Get newick strings */
        File theFile = new File(args[0]);
        BufferedReader br = new BufferedReader(new FileReader(theFile));

        String strA = br.readLine();
        strA = strA.trim();
        String strB = br.readLine();
        strB = strB.trim();
        br.close();
        
        treebase.Tree tA = new treebase.Tree(strA);
        treebase.Tree tB = new treebase.Tree(strB);
        TreeAsSplits treeA = new TreeAsSplits(tA);
        TreeAsSplits treeB = new TreeAsSplits(tB);
        System.out.println(treeA.toString());
        System.out.println(treeB.toString());
        System.out.println();

        Geodesic h = new Geodesic(treeA, treeB);
        System.out.println(h.summary());

        TreeSpacePathViewer theViewer = new TreeSpacePathViewer(h, 50); // 50 intermediate trees

        //Create and set up window.
        javax.swing.JFrame frame = new javax.swing.JFrame("Geodesic Viewer");
        frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);

 
        frame.setContentPane(theViewer);

        // Set position and size
        frame.setLocation(100,100);

        //Display the window.
        frame.pack();
        frame.setVisible(true);


    }



}
