/*
 * PathFromArrayOfTrees.java

    Copyright (C) 2016  Tom M. W. Nye

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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.TreeAsSplits;
import treedatasets.TreeAsSplitsDataSet;

/**
 * Enables viewing of any arbitrary array of trees as a path in tree-space
 */

public class PathFromArrayOfTrees implements TreeSpacePath {
    
    private TreeAsSplits[] theTrees;
    
    public PathFromArrayOfTrees(TreeAsSplits[] s) {
        theTrees = s;
    }
    
    public PathFromArrayOfTrees(ArrayList<TreeAsSplits> s) {
        theTrees = (TreeAsSplits[]) s.toArray();
    }
    
    public PathFromArrayOfTrees(TreeAsSplitsDataSet s) {
        theTrees = (TreeAsSplits[]) s.theTrees.toArray();
    }

    @Override
    public TreeAsSplits getTreeOnPath(double s) {
        
        final int n = theTrees.length;
        if (s>=n) {
            return theTrees[n-1];
        }
        
        int ind = ((int) Math.floor(s));
        return theTrees[ind];
    }

    @Override
    public double[] getPathBounds() {
        double[] lim = new double[2];
        lim[0]=0.0; lim[1]=theTrees.length;
        return lim;
    }

    
    @Override
    /** Get rotational ordering */
    public ArrayList<String> getRotationalOrder() {

        int numTaxa = theTrees[0].getNumTaxa();
        ArrayList<String> rotationalOrder = new ArrayList(), referenceOrder = new ArrayList();
        rotationalOrder.addAll(theTrees[0].getTaxa());
        Collections.sort(rotationalOrder);
        referenceOrder.addAll(rotationalOrder);
        //Collections.shuffle(rotationalOrder); // A way to avoid local minima?
        
        /* Get all the non-triv splits */
        HashSet<Split> allNonTrivSplits = new HashSet();
        for (int i=0; i<theTrees.length; i++) {
            allNonTrivSplits.addAll(theTrees[i].getNonTrivialSplits());
        }

        int score = scoreRotationalOrder(rotationalOrder, allNonTrivSplits);
        int oldScore = score+1;
        int testScore, optScore, optPos;
        String taxon;

        while (score<oldScore) {

            oldScore = score;

            for (int k=0; k<numTaxa; k++) {
                taxon = referenceOrder.get(k);
                optScore = score;
                optPos = rotationalOrder.indexOf(taxon);
                // Try re-inserting taxon k in each position
                for (int i=0; i<numTaxa; i++) {
                    rotationalOrder.remove(taxon);
                    rotationalOrder.add(i, taxon);
                    testScore = scoreRotationalOrder(rotationalOrder, allNonTrivSplits);
                    if (testScore<optScore) {
                        optScore = testScore;
                        optPos = i;
                    }
                }
                score = optScore;
                rotationalOrder.remove(taxon);
                rotationalOrder.add(optPos, taxon);
                // Completes search for position for taxon k
            }

        }
//        for (int i=0; i<rotationalOrder.size(); i++) System.out.println(rotationalOrder.get(i));
        return rotationalOrder;
    }
    
    /* allSplits should contain all non-triv splits in the collections of trees */
    private int scoreRotationalOrder(ArrayList<String> theOrder, HashSet<Split> allSplits) {

        String tA="", tB="";
        Iterator<Split> it = allSplits.iterator();
        int score=0;
        try {
            while (it.hasNext()) {
                Split s = it.next();
                for (int i=0; i<theOrder.size(); i++) {
                    tA = theOrder.get(i);
                    if ((i+1)<theOrder.size()) {
                        tB = theOrder.get(i+1);
                    }
                    else {
                        tB = theOrder.get(0);
                    }
                    if (s.separatesTaxa(tA, tB)) score++;
                }
            }
        }
        catch (AlgorithmException anErr) {
            System.out.println("Error finding a rotational order for a geodesic. "+anErr.getMessage());
        }
        return score;
    }


    
    @Override
    public Split getRootSplit() {
        Split root = null;
                
        /* Find any split contained in all the trees */
        TreeAsSplits common = theTrees[0];
        for (int i=1; i<theTrees.length; i++) {
            TreeAsSplits intersect = TreeAsSplits.intersection(common, theTrees[i]);
            common = intersect;
        }
        
        // Find common internal splits
        HashSet<Split> sharedInternalSplits = common.getNonTrivialSplits();
        if (sharedInternalSplits.size()>0) {
            /* Now find a "balanced" root c.f. calculation in Geodesic */
            double balance = 2.0, b;
            int numTaxa = common.getNumTaxa();
            Iterator<Split> it = sharedInternalSplits.iterator();
            for (int i=0; i<sharedInternalSplits.size(); i++) {
                Split s =  it.next();
                b = ((double)s.getTaxonSubset(null).size()) / ((double) numTaxa);
                if (Math.abs(b-0.5)<balance) {
                    balance = Math.abs(b-0.5);
                    root = s;
                }
            }
        }
        
        if (root==null) {
            // No non-trivial splits in common, so use a pendant split
            ArrayList<Split> leafSplits = new ArrayList();
            leafSplits.addAll(common.getSplits());
            Collections.sort(leafSplits);
            root = leafSplits.get(0);
        }
        
        return root;
    }
    
    
    /* Example: generate a path from a random walk.  */
    
/*    

    public static void main(String[] args) throws treebase.AlgorithmException {


        String strA = "((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):2);";
        String strB = "((B:1,C:1):2,(A:1,E:1):1,(D:1,F:1):1);";
        
        strA = "((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);";
        
        treebase.Tree tA = new treebase.Tree(strA);
        treebase.Tree tB = new treebase.Tree(strB);
        TreeAsSplits treeA = new TreeAsSplits(tA);
        TreeAsSplits treeB = new TreeAsSplits(tB);
        System.out.println(treeA.toString());
        System.out.println(treeB.toString());
        System.out.println();
        
        final int k=50;
        TreeAsSplits[] theTrees = new TreeAsSplits[k+1];
        
        // Create an array of trees via random walk
        Random.setEngine((int)System.currentTimeMillis());
        theTrees[0] = treeA;
        for (int i=1; i<=k; i++) {
            theTrees[i] = RandomWalkSimulator.sampleRandomWalk(theTrees[i-1], 0.04, 1, false);
        }     
        
        PathFromArrayOfTrees thePath = new PathFromArrayOfTrees(theTrees);
        TreeSpacePathViewer theViewer = new TreeSpacePathViewer(thePath,(k+1));

        // Create and set up window.
        javax.swing.JFrame frame = new javax.swing.JFrame("Tree Space Path Viewer");
        frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);

 
        frame.setContentPane(theViewer);

        // Set position and size
        frame.setLocation(100,100);

        //Display the window.
        frame.pack();
        frame.setVisible(true);


    }

*/
    
}
