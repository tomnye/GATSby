/*
 * Geodesic.java

    Copyright (C) 2012  Tom M. W. Nye

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

package geodesics;

/**
 * Geodesic path between two trees.
 * Calculation based on Owen's 2010 polynomial time algorithm.
 * Developed from version of 5/1/2011 but based on trees as sets of splits.
 */


import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.ListIterator;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;


public class Geodesic implements geodesicgraphics.TreeSpacePath {

    public static final boolean DEBUG_ON = true;

    protected HashSet<Split> sharedSplits; /* Set of splits in both end points */
    protected TreeAsSplits treeA, treeB; /* Map from splits to lengths for end points */
    /* The geodesic goes from A to B */

    protected ArrayList<HashSet<Split>> partitionA, partitionB; /* The sets A_i B_i of splits defining the orthonts along the geodesic */
    protected ArrayList<MinWeightCover> covers; /* Order list of the cover objects for each set of splits in the partition. Null if not yet computed */

    protected double length=-1.0; // What it says!
    protected double sharedLength=-1.0, unsharedLength; // Lengths corresponding to shared and disjoint splits respectively -- useful for debugging
    protected double internalLength=-1.0; // Length corresponding to internal splits only


    /** Construct geodesic between two trees using Owen's 2010 algorithm.
     Constructor from two sets of splits.
     */
    public Geodesic(TreeAsSplits tA, TreeAsSplits tB) throws AlgorithmError {
        build(tA, tB);
    }
    public void build(TreeAsSplits tA, TreeAsSplits tB) throws AlgorithmError {
        treeA = tA.efficientClone();
        treeB = tB.efficientClone();
 
        treeA.removeZeroLengthSplits();
        treeB.removeZeroLengthSplits();

        // Check taxon sets match
        if (DEBUG_ON) {
            HashSet<String> taxaA = tA.getTaxa();
            HashSet<String> taxaB = tB.getTaxa();
            if (!((taxaA.containsAll(taxaB))&&(taxaB.containsAll(taxaA)))) {
                throw new AlgorithmError("Cannot build geodesic between trees with different taxon sets. ");
            }
        }

        build();
    }


    /** Used by constructors to build the geodesic */
    protected void build() {

        /* Get intersection of sets of splits -- sharedSplits. */
        sharedSplits = new HashSet();
        HashSet<Split> splitsA = treeA.getSplits();
        HashSet<Split> splitsB = treeB.getSplits();
       
        sharedSplits = symmetricDifference(splitsA, splitsB, treeA.getNumTaxa());
        
        // Get the initial partition corresponding to cone path
        //splitsA.removeAll(sharedSplits);
        //splitsB.removeAll(sharedSplits);
        partitionA = new ArrayList();
        partitionB = new ArrayList();
        covers = new ArrayList();

        // Take care of one tree being subset of other
        if (splitsA.size()==0) {
            sharedSplits.addAll(splitsB);
            splitsB.clear();
        }
        else if (splitsB.size()==0) {
            sharedSplits.addAll(splitsA);
            splitsA.clear();
        }

        // If there are some non-shared splits iterate thru' extending partitions
        if ((splitsA.size()>0)&&(splitsB.size()>0)) {
            partitionA.add(splitsA);
            partitionB.add(splitsB);
            covers.add(null);
            boolean done = false;
            while (!done) {
                done = refinePathSpace();
            }
        }

        // Tidy up and compute total length
        //computeLength();
        // No! don't compute length till you actually need it!
        length = -1;
    }
    
    
    private HashSet<Split> symmetricDifference(HashSet<Split> a, HashSet<Split> b, int numTaxa) {
        boolean resolved = ((a.size()==(2*numTaxa-3))&&(b.size()==(2*numTaxa-3)));
        HashSet<Split> shared = new HashSet();
        Iterator<Split> it = a.iterator(); // Could be more efficient to loop through smaller set...?
        while (it.hasNext()) {
            Split p = it.next();
            if (p.isTerminal()!=null) {
                shared.add(p);
                b.remove(p);
            }
            else if (b.contains(p)) {
                shared.add(p);
                b.remove(p);
            }
        }
        a.removeAll(shared);
        
        if (resolved) return shared;
         
        /* Little bit of extra work to do for unresolved ends */
        HashSet<Split> sharedA = new HashSet();
        Iterator<Split> itA = a.iterator(); 
        while (itA.hasNext()) {
            Split pA = itA.next();
            boolean addA = true;
            Iterator<Split> itB = b.iterator(); 
            while (itB.hasNext()) {
                Split pB = itB.next();
                if (!pA.isCompatible(pB)) {
                    addA = false;
                    break;
                }
            }
            if (addA) sharedA.add(pA);
        }
        
        HashSet<Split> sharedB = new HashSet();
        Iterator<Split> itB = b.iterator(); 
        while (itB.hasNext()) {
            Split pB = itB.next();
            boolean addB = true;
            itA = a.iterator(); 
            while (itA.hasNext()) {
                Split pA = itA.next();
                if (!pA.isCompatible(pB)) {
                    addB = false;
                    break;
                }
            }
            if (addB) sharedB.add(pB);
        }
        
        a.removeAll(sharedA);
        b.removeAll(sharedB);
        shared.addAll(sharedA);
        shared.addAll(sharedB);
        return shared;
    }  


    /** Given a path space compute all min weight covers and check weights>= 1.0 */
    private boolean refinePathSpace() {

        if (DEBUG_ON) {
            for (int i=1; i<partitionA.size(); i++) {
                HashSet<Split> pA = partitionA.get(i);
                HashSet<Split> pB = partitionB.get(i);
                double na = edgeLenNorm(pA, treeA);
                double nb = edgeLenNorm(pB, treeB);
                double y = (na/nb);
                pA = partitionA.get(i-1);
                pB = partitionB.get(i-1);
                na = edgeLenNorm(pA, treeA);
                nb = edgeLenNorm(pB, treeB);
                double x = (na/nb);
                if (x>y) System.out.println("Error building a geodesic: non-proper path found");
            }
        }

        int numNewCovers = 0;
        // Loop thru partition: first pass to compute new min weight covers
        for (int i=0; i<partitionA.size(); i++) {
            HashSet<Split> pA = partitionA.get(i);
            HashSet<Split> pB = partitionB.get(i);
            if (( covers.get(i) == null ) && ( (pA.size()>1 )||( pB.size()>1 ) )) {
                MinWeightCover mwc = new MinWeightCover(partitionA.get(i), partitionB.get(i));
                if (mwc.getWeight() < 1.0) numNewCovers++;
                covers.set(i, mwc);
            }
        }

        // Loop thru partition: second pass to refine partition
        // For each new cover...
        for (int j=0; j<numNewCovers; j++) {
            // ... find the corresponding part of the partition
            for (int i=0; i<partitionA.size(); i++) {
                if (covers.get(i)!=null) {
                    MinWeightCover mwc = covers.get(i);
                    if (mwc.getWeight() < 1.0) {
                        // Replace A_i and B_i with subdivision corresponding to min weight vertex cover
                        partitionA.set(i, mwc.newA1());
                        partitionA.add(i+1, mwc.newA2());
                        partitionB.set(i, mwc.newB1());
                        partitionB.add(i+1, mwc.newB2());
                        covers.set(i, null);
                        covers.add(i+1, null);
                        break; // We've changed the size of the partitions, so loop thru' again
                    }
                }
            }
        } // End second pass extending covers

        return (numNewCovers == 0);
    }

    /** Given the correct path space, now compute the total length */
    protected void computeLength() {

        // Work out distance for shared splits.
        // This is a simple euclidean distance
        double total = 0.0;
        double u,v;
        Iterator<Split> it = sharedSplits.iterator();
        double leafLength = 0.0;
        for (int i=0; i<sharedSplits.size(); i++) {
            Split s = it.next();
            if (treeA.contains(s)) u = treeA.getSplitLength(s);
            else u = 0.0;
            if (treeB.contains(s)) v = treeB.getSplitLength(s);
            else v = 0.0;
            total += (u-v)*(u-v);
            if (s.isTerminal()!=null) {
                // Terminal leaf
                leafLength += (u-v)*(u-v);
            }
        }
        sharedLength = Math.sqrt(total);
        if (partitionA.size()==0) {
            length = sharedLength;
            internalLength = Math.sqrt(total-leafLength);
            return;
        }

        // Work out total distance of the "unshared" splits
//        double[] normsA = new double[partitionA.size()];
//        double[] normsB = new double[partitionB.size()];
        total = 0.0;
        HashSet<Split> ai, bi;
        for (int i=0; i<partitionA.size(); i++) {
            ai = partitionA.get(i);
            bi = partitionB.get(i);
            u = Math.sqrt(edgeLenNorm(ai, treeA));
            v = Math.sqrt(edgeLenNorm(bi, treeB));
//            normsA[i] = u;
//            normsB[i] = v;
            total += (u+v)*(u+v);
        }
        unsharedLength = Math.sqrt(total);

        // Now combine shared and unshared pieces: straightforward Pythagoras
        length = Math.sqrt(sharedLength*sharedLength + unsharedLength*unsharedLength);
        internalLength = Math.sqrt(sharedLength*sharedLength + unsharedLength*unsharedLength-leafLength);

    }

    /** Create a string description */
    public String summary() {
        if (length<0.0) computeLength();
        java.text.DecimalFormat df = new java.text.DecimalFormat("####.######");
        String s = "Geodesic total length = "+df.format(length) +"\n";
        s += "Length from non-consensus splits = "+df.format(unsharedLength) +"\n";
        s += "Length from consensus splits = "+df.format(sharedLength) +"\n\n";
        s += "Euclidean distance = " +df.format(this.calcEuclideanDistance())+"\n";
        s += "Cone path distance = " +df.format(this.calcConePathDistance())+"\n";
        s += "\nDetails:-\n";
        s += "Geodesic passes through "+(partitionA.size()+1)+" orthants.\n";

        boolean simple = true;
        for (int i=0; i<partitionA.size(); i++) {
            HashSet<Split> h = partitionA.get(i);
            if (h.size()>1) simple = false;
        }
        for (int i=0; i<partitionB.size(); i++) {
            HashSet<Split> h = partitionB.get(i);
            if (h.size()>1) simple = false;
        }
        if (simple) s += "Geodesic is simple.\n";
        else s += "Geodesic is non-simple.\n";

        return s;
    }

    /** Return the length */
    public double getLength() {
        if (length<0.0) computeLength();
        return length;
    }

    /** Return the length for internal splits*/
    public double getInternalLength() {
        if (length<0.0) computeLength();
        return internalLength;
    }

    public int numOrthants() {
        return partitionA.size()+1;
    }

    public boolean isSimple() {
        boolean simple = true;
        for (int i=0; i<partitionA.size(); i++) {
            HashSet<Split> h = partitionA.get(i);
            if (h.size()>1) simple = false;
        }
        for (int i=0; i<partitionB.size(); i++) {
            HashSet<Split> h = partitionB.get(i);
            if (h.size()>1) simple = false;
        }
        return simple;
    }
    
    public int[][] getPartitionSizes() {
        int[][] ps = new int[2][partitionA.size()];
        for (int i=0; i<partitionA.size(); i++) {
            ps[0][i] = partitionA.get(i).size();
            ps[1][i] = partitionB.get(i).size();
        }
        return ps;
    }
    
    /* -------------------------------------------------------------------- */


    /** Get a tree on the geodesic.
     Needs s between 0 (=A) and 1 (=B) */
    public TreeAsSplits getTree(double s) {

        TreeAsSplits theTree = null;

        // Deal with trivial cases
        if (s==0.0) {
            theTree = treeA.clone();
            return theTree;
        }
        if (s==1.0) {
            theTree = treeB.clone();
            return theTree;
        }

        theTree = new TreeAsSplits(treeA.getNumTaxa());
        
        /* Start by getting the lambda values which separate orthonts along the
           geodesic and the norm of each partition.
           Simultaneously find which orthont we should be in. */
        int k = partitionA.size();
        int orthontIndex = k;  // Unless we find s to be less than a lambda critical value, then the orthont must be the last one i.e. k
        double lowerLam = 0.0;
        double[] normsA = new double[k];
        double[] normsB = new double[k];
        double[] criticalLambda = new double[k];
        HashSet<Split> ai, bi;
        for (int i=0; i<k; i++) {
            ai = partitionA.get(i);
            bi = partitionB.get(i);
            normsA[i] = Math.sqrt(edgeLenNorm(ai, treeA));
            normsB[i] = Math.sqrt(edgeLenNorm(bi, treeB));
            criticalLambda[i] = normsA[i]/(normsA[i]+normsB[i]);
            if ((lowerLam<=s)&&(s<criticalLambda[i])) orthontIndex = i;
            lowerLam = criticalLambda[i];
        }

        /* OK, assemble topology (i.e. the set of splits):
           first "unshared" splits, which come from the partition A_i and B_i */
        Iterator<Split> it;
        Split p;
        double nai, nbi; // norms of Ai and Bi
        double e, u, v;
        /* Loop thru parts of the B partition which we need to add:
           add in B_1 up to B_(orthontIndex) */
        for (int i=0; i<orthontIndex; i++) {
            bi = partitionB.get(i);
            // Loop thru Bi: work out length and add to treeMap
            nai = normsA[i];
            nbi = normsB[i];
            it = bi.iterator();
            for (int j=0; j<bi.size(); j++) {
                p = it.next();
                e = treeB.getSplitLength(p);
                u = (s*nbi-(1.0-s)*nai)*e/nbi;
                theTree.add(p, u);
            }
        }
        /* Loop thru parts of the A partition which we need to add:
           add in A_(orthontIndex+1) up to A_(k) */
        for (int i=orthontIndex; i<k; i++) {
            ai = partitionA.get(i);
            // Loop thru Ai: work out length and add to treeMap
            nai = normsA[i];
            nbi = normsB[i];
            it = ai.iterator();
            for (int j=0; j<ai.size(); j++) {
                p = it.next();
                e = treeA.getSplitLength(p);
                u = ((1.0-s)*nai-s*nbi)*e/nai;
                theTree.add(p, u);
            }
        }

        /* Remains to do the shared splits -- easier! */
        it = sharedSplits.iterator();
        for (int i=0; i<sharedSplits.size(); i++) {
            p = (Split) it.next();
            if (treeA.contains(p)) u = treeA.getSplitLength(p);
            else u = 0.0;
            if (treeB.contains(p)) v = treeB.getSplitLength(p);
            else v = 0.0;
            e = (1.0-s)*u+s*v; // Linear scale on shared splits!
            theTree.add(p, e);
        }
        
        return theTree;
     }

    /** Get a collection of trees along the geodesic. Useful for debugging!
     Pass in an array of doubles -- values between 0 and 1 where you'd like a tree.
     Alternatively, pass in null for s. The method computes s values corresponding
     to one tree per quadrant which the geodesic traverses. */
    public Tree[] getCollectionOfTrees(double[] s, int numPts) {
        if (s==null) {
            if (numPts==0) {
                // No set of coord values have been passed in so compute crucial ones.
                if (partitionA.size()==0) { // This is if both end points are in the same quadrant.
                    numPts = 3;
                    s = new double[numPts];
                    s[0] = 0.0;
                    s[1] = 0.5;
                    s[2] =1.0;
                }
                else {
                    // Non-trivial path through several orthonts.
                    // Get one s value in each orthont.
                    // Start by computing the critical lambda values
                    int k = partitionA.size();
                    double[] normsA = new double[k];
                    double[] normsB = new double[k];
                    double[] criticalLambda = new double[k];
                    HashSet<Split> ai, bi;
                    for (int i=0; i<k; i++) {
                        ai = partitionA.get(i);
                        bi = partitionB.get(i);
                        normsA[i] = Math.sqrt(edgeLenNorm(ai, treeA));
                        normsB[i] = Math.sqrt(edgeLenNorm(bi, treeB));
                        criticalLambda[i] = normsA[i]/(normsA[i]+normsB[i]);
                    }

                    numPts = k+1;
                    // Work out s vals
                    s = new double[numPts];
                    s[0] = 0.5*criticalLambda[0];
                    s[k] = 0.5*(1.0+criticalLambda[(k-1)]);
                    for (int i=0; i<(k-1); i++) {
                        s[i+1] = 0.5*(criticalLambda[i]+criticalLambda[i+1]);
                    }
                }
            } // End no info passed in at all
            else {
                s = new double[numPts];
                double delta = 1.0/(numPts-1);
                for (int i=0; i<numPts; i++) {
                    s[i] = i*delta;
                }
            }
        }
        else {
            numPts = s.length;
        }
        
        // OK, got a valid list of s values
        Tree[] theTrees = new Tree[numPts];
        for (int i=0; i<numPts; i++) {
            TreeAsSplits myTreeMap = this.getTree(s[i]);
            try {
                Tree myTree = new Tree(myTreeMap.getMap());
                theTrees[i] = myTree;
            }
            catch (AlgorithmError anErr) {
                System.out.println("AlgorithmError computing trees along a geodesic. "+anErr.getMessage());
            }

        }
        return theTrees;
    }  

    /* -------------------------------------------------------------------- */

    /* Methods for shifting the end points */

    public void shiftSingleEnd(boolean whichEnd, TreeAsSplits newEnd) {
        if (whichEnd) {
            treeB = newEnd.efficientClone();
        }
        else {
            treeA = newEnd.efficientClone();
        }
        length = -1; // Reset length to be re-computed
        build();
    }

    /** Create a new geodesic with different end points, specified as hashmaps.
     This is useful, for example, when iteratively moving one end point (see PCA).
     NB: assumes X and Y based on <same sets of splits> as A, B */
    public Geodesic shiftEndPoints(TreeAsSplits treeX, TreeAsSplits treeY) throws AlgorithmError {

        Geodesic g = null;

        try {
            /* 1. Check end point topologies match originals */
            HashSet<Split> splitsX = treeX.getSplits();
            HashSet<Split> splitsY = treeY.getSplits();
            if (treeX.getNumSplits()!=treeA.getNumSplits()) throw new AlgorithmException("Topology mismatch");
            if (treeY.getNumSplits()!=treeB.getNumSplits()) throw new AlgorithmException("Topology mismatch");
            if (!treeA.containsAll(splitsX)) throw new AlgorithmException("Topology mismatch");
            if (!treeB.containsAll(splitsY)) throw new AlgorithmException("Topology mismatch");
 
            /* 2. Check branch length condition. (Owen, property P2) */
            int k = partitionA.size();
            double[] normsX = new double[k];
            double[] normsY = new double[k];
            HashSet<Split> ai, bi;
            for (int i=0; i<k; i++) {
                ai = partitionA.get(i);
                bi = partitionB.get(i);
                normsX[i] = Math.sqrt(edgeLenNorm(ai, treeX));
                normsY[i] = Math.sqrt(edgeLenNorm(bi, treeY));
            }
            boolean good = true;
            for (int i=1; i<k; i++) {
                good = (good)&&( (normsX[i-1]/normsY[i-1]) <= (normsX[i]/normsY[i]) );
            }
            if (!(good)) throw new AlgorithmException("Proper path space condition failed");

            // OK create new geodesic with existing path -- conditon P3
            g = new Geodesic(treeX, treeY, partitionA, partitionB);
            return g;
        }
        catch (AlgorithmException anEx) {
            // OK -- we've failed some criterion -- return a new geodesic
            g = new Geodesic(treeX, treeY); // 0 means no checking!
            return g;
        }

    }

    /** Constructor where the path must be provided. */
    private Geodesic(TreeAsSplits treeX, TreeAsSplits treeY, ArrayList<HashSet<Split>> pA, ArrayList<HashSet<Split>> pB) throws AlgorithmException {

        sharedSplits = new HashSet();
        sharedSplits.addAll(treeX.getSplits());
        sharedSplits.addAll(treeY.getSplits());
        for (int i=0; i<pA.size(); i++) {
            HashSet<Split> o = pA.get(i);
            sharedSplits.removeAll(o);
            o = pB.get(i);
            sharedSplits.removeAll(o);
        }
        
        partitionA = new ArrayList(pA);
        partitionB = new ArrayList(pB);
        int k = partitionA.size();

        /* 3. Check min weight covers for each partition */
        // First replace A by X and B by Y
        treeA = treeX; 
        treeB = treeY; 
        int numNewCovers = 0;
        for (int i=0; i<k; i++) {
            MinWeightCover mwc = new MinWeightCover(partitionA.get(i), partitionB.get(i));
            if (mwc.getWeight() < 1.0) numNewCovers++;
        }
        if (numNewCovers>0) throw new AlgorithmException("Partition condition failed");

        // Path space is good! Re-calculate length of geodesic
        //computeLength();
        length = -1;
    }

    /* -------------------------------------------------------------------- */

    /* Util methods */

    /** Edge length norm for a set of splits */
    protected static double edgeLenNorm(ArrayList<Split> x, TreeAsSplits tree) {
        double s = 0.0;
        double u;
        ListIterator<Split> it = x.listIterator();
        for (int i=0; i<x.size(); i++) {
            Split p = it.next();
//            if (!(tree.contains(p))) {
//                System.out.println("Warning in Geodesic: tree doesn't contain split!");
//            }
            u = tree.getSplitLength(p);
            s += u*u;
        }
        return s;
    }
    protected static double edgeLenNorm(HashSet<Split> x, TreeAsSplits tree) {
       double s = 0.0;
        double u;
        Iterator<Split> it = x.iterator();
        for (int i=0; i<x.size(); i++) {
            Split p = it.next();
//            if (!(tree.contains(p))) {
//                System.out.println("Warning in Geodesic: tree doesn't contain split!");
//            }
            u = tree.getSplitLength(p);
            s += u*u;
        }
        return s;
    }

    /** Euc distance */
    public double calcEuclideanDistance() {
        double total = sharedLength*sharedLength;
        total += treeA.sumSquaredLengths(sharedSplits);
        total += treeB.sumSquaredLengths(sharedSplits);
        return Math.sqrt(total);
    }

    /** Cone path distance */
    public double calcConePathDistance() {
        double abranch = Math.sqrt(treeA.sumSquaredLengths(sharedSplits));
        double bbranch = Math.sqrt(treeB.sumSquaredLengths(sharedSplits));

        double x = abranch + bbranch;
        return Math.sqrt(x*x+sharedLength*sharedLength);
    }


    /** Get split map for simple geodesic */
    public HashMap<Split, Split> getSimpleSplitMap() throws AlgorithmException {
        
        if (!this.isSimple()) {
            throw new AlgorithmException("Geodesic is not simple.");
        }
        
        HashMap<Split, Split> res = new HashMap();
        Iterator<Split> it = sharedSplits.iterator();
        while (it.hasNext()) {
            Split p = it.next();
            res.put(p, p);
        }
        
        for (int i=0; i<partitionA.size(); i++) {
            HashSet<Split> pA = partitionA.get(i);
            HashSet<Split> pB = partitionB.get(i);
            Iterator<Split> itA = pA.iterator();
            Iterator<Split> itB = pB.iterator();
            res.put(itA.next(), itB.next());
        }
        
        if (res.size()!=treeA.getNumSplits()) {
            throw new AlgorithmError("Error constructing split map for simple geodesic.");
        }
        
        return res;
    }
        
    /* -------------------------------------------------------------------- */

    /** Private inner class: MinWeightCover
     This represents a min weight cover on the incompatibility graph of two sets of splits. */
    private class MinWeightCover {

        private HashSet<Split> coverA, coverAbar, coverB, coverBbar;
        /* These are the splits in the min weight cover*/
        private double weight;

        public MinWeightCover(HashSet<Split> aSet, HashSet<Split> bSet) {

            ArrayList<Split> A = new ArrayList(aSet);
            ArrayList<Split> B = new ArrayList(bSet);

            /* Work out vertex weights */
            double[] weightsA;
            double[] weightsB;
            weightsA = new double[A.size()];
            weightsB = new double[B.size()];
            double normA = edgeLenNorm(A, treeA);
            double normB = edgeLenNorm(B, treeB);
            Split p;
            double u;
            Iterator<Split> itA = A.iterator();
            for (int i=0; i<A.size(); i++) {
                p = itA.next();
                u = treeA.getSplitLength(p);
                weightsA[i] = u*u/normA;
            }
            Iterator<Split> itB = B.iterator();
            for (int i=0; i<B.size(); i++) {
                p = itB.next();
                u = treeB.getSplitLength(p);
                weightsB[i] = u*u/normB;
            }

            /* Get the max flow */
            EKMaxFlow theFlow = null;
            theFlow = new EKMaxFlow(A, B, weightsA, weightsB);

            /* Calculate the weight */
            weight = theFlow.max_flow;
            /* Extract the result from theFlow */
            if (weight<1.0) {
                coverA = theFlow.getACover();
                coverB = theFlow.getBCover();
                coverAbar = new HashSet(aSet);
                coverAbar.removeAll(coverA);
                coverBbar = new HashSet(bSet);
                coverBbar.removeAll(coverB);
            }

        }

        /** Return the cover weight */
        public double getWeight() {
            return weight;
        }

        /** Provide access to cover partitions */
        public HashSet<Split> newA1() {return coverA;}
        public HashSet<Split> newA2() {return coverAbar;}
        public HashSet<Split> newB1() {return coverBbar;}
        public HashSet<Split> newB2() {return coverB;}

    }

    /* -------------------------------------------------------------------- */

    /* Graphics methods */

    /** Get "slider" limits: 0 and 1 */
    public double[] getPathBounds() {
        double[] lim = new double[2];
        lim[0]=0.0; lim[1]=1.0;
        return lim;
    }

    private int scoreRotationalOrder(ArrayList<String> theOrder) {
        // Get a list of all the nontrivial splits
        HashSet<Split> allSplits = new HashSet();
        allSplits.addAll(treeA.getNonTrivialSplits());
        allSplits.addAll(treeB.getNonTrivialSplits());

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

    /** Get rotational ordering */
    public ArrayList<String> getRotationalOrder() {

        int numTaxa = treeA.getNumTaxa();
        ArrayList<String> rotationalOrder = new ArrayList(), referenceOrder = new ArrayList();
        rotationalOrder.addAll(treeA.getTaxa());
        Collections.sort(rotationalOrder);
        referenceOrder.addAll(rotationalOrder);
        //Collections.shuffle(rotationalOrder); // A way to avoid local minima?

        int score = scoreRotationalOrder(rotationalOrder);
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
                    testScore = scoreRotationalOrder(rotationalOrder);
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


    /** Get a shared root split */
    public Split getRootSplit() {

        /* Find a suitable root split for graphics.
         Find non-terminal shared splits and select one which is "most balanced" */

        HashSet<Split> leafSplits = new HashSet();
        Split root = null, s;
        double balance = 2.0, b;
        int numTaxa = treeA.getNumTaxa();
        Iterator<Split> it = sharedSplits.iterator();
        for (int i=0; i<sharedSplits.size(); i++) {
            s =  it.next();
            if (!(s.isTerminal()==null)) {
                leafSplits.add(s);
            }
            else {
                b = ((double)s.getTaxonSubset(null).size()) / ((double) numTaxa);
                if (Math.abs(b-0.5)<balance) {
                    balance = Math.abs(b-0.5);
                    root = s;
                }
            }
        }
        if (root==null) {
            // No shared internal splits
            ArrayList<Split> orderedLeaves = new ArrayList();
            orderedLeaves.addAll(leafSplits);
            Collections.sort(orderedLeaves);
            root = orderedLeaves.get(0);
        }

        return root; 
    }

    public TreeAsSplits getTreeOnPath(double s) {
        return getTree(s);
    }

    /* -------------------------------------------------------------------- */

    /* Debug and testing area */

    /** Prepare a representative set of trees along the geodesic
        and output as newick strings. */
    public void outputTreesOnPath() {
        Tree[] theTrees = this.getCollectionOfTrees(null, 0);
        Tree t;
        String s;
        for (int i=0; i<theTrees.length; i++) {
            t = theTrees[i];
            s = t.toOrderedString();
            System.out.println(s);
        }
    }
   

}