/*
 * TreeAsSplits.java

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

package treebase;

/**
 * Representation of a tree as a weighted set of splits.
 */

import java.util.*;

public class TreeAsSplits implements java.io.Serializable {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 0L;


    /* Instance data: the map */
    private HashMap<Split,Double> theMap;
    private int numTaxa;

    /** Constructor from a tree */
    public TreeAsSplits(Tree theTree) {
        theMap = theTree.getSplitBranchLengths();
        numTaxa = theTree.numTaxa();
        checkNaNs();
    }

    /** Empty constructor */
    public TreeAsSplits(int n) {
        theMap = new HashMap();
        numTaxa = n;
    }

    /** Constructor from a HashMap */
    public TreeAsSplits(HashMap<Split,Double> t) {
        theMap = t;
        Split p = t.keySet().iterator().next();
        numTaxa = p.numTaxa();
        checkNaNs();
    }

    public TreeAsSplits clone() {
        HashMap<Split,Double> t = new HashMap();
        Iterator<Split> it = theMap.keySet().iterator();
        for (int i=0; i<theMap.size(); i++) {
            Split p = it.next();
            double l = this.getSplitLength(p);
            t.put(p, new Double(l));
        }
        return new TreeAsSplits(t);
    }

    /** Clone but retain the Doubles */
    public TreeAsSplits efficientClone() {
        HashMap<Split,Double> t = new HashMap();
        t.putAll(theMap);
        return new TreeAsSplits(t);
    }

    /* Util methods --------------------------------------------------------- */

    /** Iterator thru' the splits */
    public Iterator<Split> getSplitIterator() {
        return theMap.keySet().iterator();
    }

    public Iterator<Map.Entry<Split,Double>> getEntryIterator() {
        return theMap.entrySet().iterator();
    }

    /** Get split length */
    public double getSplitLength(Split p) {
        if (theMap.containsKey(p)) return theMap.get(p).doubleValue();
        else return 0.0;
    }

    /** Set split length */
    public void setSplitLength(Split p, double x) throws AlgorithmError {
        if (theMap.containsKey(p)) {
            theMap.remove(p);
            theMap.put(p, x);
        }
        else {
            try {
                this.addWithCheck(p, x);
            }
            catch (AlgorithmException anErr) {
                throw new AlgorithmError("SetSplitLength to bad split in TreeAsMap.");
            }
        }
    }

    /** Scale a split */
    public void scaleSplitLength(Split p, double x) throws AlgorithmError {
        if (!theMap.containsKey(p)) throw new AlgorithmError("Request to change length of a split that is not contained in a TreeAsMap");
        double l = this.getSplitLength(p);
        theMap.remove(p);
        theMap.put(p, new Double(x*l));
    }
    public void scaleAllLengths(double x) {
        HashSet<Split> splits = new HashSet();
        splits.addAll(theMap.keySet());
        Iterator<Split> it = splits.iterator();
        while (it.hasNext()) {
            Split p = it.next();
            double l = this.getSplitLength(p);
            theMap.remove(p);
            theMap.put(p, new Double(x*l));
        }
    }


    public void scaleViaMap(HashMap<Split,Double> transfac, boolean invert) {
        Iterator<Split> itS = transfac.keySet().iterator();
        Split s;
        Double f;
        try {
            while (itS.hasNext()) {
                s = itS.next();
                f = transfac.get(s);
                if (f!=null) {
                    if (this.contains(s)) {
                        if (invert) this.scaleSplitLength(s, 1.0/f.doubleValue());
                        else this.scaleSplitLength(s, f.doubleValue());
                    }
                }
            }
        }
        catch (AlgorithmError anErr) {
            System.out.println("Error scaling tree. "+anErr.getMessage());
        }
    }

    /** Add to edge length */
    public void addToEdgeLength(Split p, double x) throws AlgorithmException {
        double l = this.getSplitLength(p);
        if ((l+x)<0.0) throw new AlgorithmException("Negative branch length");
        else this.setSplitLength(p, l+x);
    }

    /** Remove split */
    public void remove(Split p) {
        theMap.remove(p);
    }
    public void removeAll(Collection<Split> c) {
        Iterator<Split> it = c.iterator();
        for (int i=0; i<c.size(); i++) {
            theMap.remove(it.next());
        }
    }

    /** Remove zero length INTERNAL splits */
    public void removeZeroLengthSplits() {
        Iterator<Map.Entry<Split,Double>> it = getEntryIterator();
        HashSet<Split> toRemove = new HashSet();
        Map.Entry<Split,Double> e;
        Split p;
        double x;
        for (int i=0; i<theMap.size(); i++) {
            e = it.next();
            p = e.getKey();
            x = e.getValue();
            if ((x<1.0E-20)&&(p.isTerminal()==null)) {
                toRemove.add(p);
            }
        }
        Iterator<Split> itR = toRemove.iterator();
        for (int i=0; i<toRemove.size(); i++) {
            p = itR.next();
            remove(p);
        }
    }

    /** Remove INTERNAL splits with length below some threshold */
    public void removeZeroLengthSplits(double threshold) {
        Iterator<Map.Entry<Split,Double>> it = getEntryIterator();
        HashSet<Split> toRemove = new HashSet();
        Map.Entry<Split,Double> e;
        Split p;
        double x;
        for (int i=0; i<theMap.size(); i++) {
            e = it.next();
            p = e.getKey();
            x = e.getValue();
            if ((x<threshold)&&(p.isTerminal()==null)) {
                toRemove.add(p);
            }
        }
        Iterator<Split> itR = toRemove.iterator();
        for (int i=0; i<toRemove.size(); i++) {
            p = itR.next();
            remove(p);
        }
    }

    /** Add a split: it's up to the programmer to ensure compatibility! */
    public void add(Split p, double x) {
        theMap.put(p, x);
    }

    /** Add a split carefully */
    public void addWithCheck(Split p, double x) throws AlgorithmException {
        if (!checkCompatibility(p)) throw new AlgorithmException("Added split not compatible with a TreeAsMap");
        theMap.put(p, x);
    }

    /** Access the map */
    public HashMap getMap() {
        return theMap;
    }

    /** Check for presence of a split */
    public boolean contains(Split p) {
        return theMap.containsKey(p);
    }
    public boolean containsAll(HashSet<Split> h) {
        Iterator<Split> it = h.iterator();
        for (int i=0; i<h.size(); i++) {
            if (!theMap.containsKey(it.next())) return false;
        }
        return true;
    }

    public boolean matchesTopology(TreeAsSplits t) {
        Iterator<Split> it = getSplitIterator();
        for (int i=0; i<theMap.size(); i++) {
            Split p = it.next();
            if (!t.contains(p)) return false;
        }
        it = t.getSplitIterator();
        for (int i=0; i<t.getNumSplits(); i++) {
            Split p = it.next();
            if (!this.contains(p)) return false;
        }
        return true;
    }

    /** Check whether a split is compatible with the current split set */
    public boolean checkCompatibility(Split p) {
       Iterator<Split> it = theMap.keySet().iterator();
        for (int i=0; i<getNumSplits(); i++) {
            Split q = it.next();
            if (!p.isCompatible(q)) return false;
        }
       return true;
    }

    public int getNumSplits() {
        return theMap.size();
    }

    public int getNumInternalSplits() {
        return theMap.size()-numTaxa;
    }

    public HashSet<Split> getSplits() {
        HashSet<Split> s = new HashSet();
        s.addAll(theMap.keySet());
        return s;
    }

    public HashSet<Split> getNonTrivialSplits() {
        HashSet<Split> s = new HashSet();
        Iterator<Split> it = theMap.keySet().iterator();
        for (int i=0; i<theMap.size(); i++) {
            Split p = it.next();
            if (p.isTerminal()==null) {
                s.add(p);
            }
        }
        return s;
    }

    /** Check splits are compatible */
    public void sanityCheck() {
        Iterator<Split> itA = theMap.keySet().iterator();
        for (int i=0; i<getNumSplits(); i++) {
            Split a = itA.next();
            Iterator<Split> itB = theMap.keySet().iterator();
            for (int j=0; j<getNumSplits(); j++) {
                Split b = itB.next();
                if (!a.isCompatible(b)) {
                    System.out.println("Error: compatibility check failed on a TreeAsMap");
                    System.out.println(getNumSplits());
                    return;
                }
            }
        }
    }

    /** Check for NaNs */
    public void checkNaNs() {
       Iterator<Split> it = getSplitIterator();
        Split p;
        for (int i=0; i<getNumSplits(); i++) {
            p = it.next();
            if (Double.isNaN(getSplitLength(p))) {
                System.out.println("NaN in TreeAsMap.");
            }
        }
    }


    /** Get the number of taxa */
    public int getNumTaxa() {
        return numTaxa;
    }

    public HashSet<String> getTaxa() {
        Split p = theMap.keySet().iterator().next();
        HashSet<String> t = p.allTaxa();
        return t;
    }

    /** Check whether fully resolved */
    public boolean fullyResolved() {
        return ((2*this.numTaxa-3)==getNumSplits());
    }

    /** Build a tree */
    public Tree getTree() throws AlgorithmError {
        return new Tree(theMap);
    }

    /** Sum squared lengths */
    public double sumSquaredLengths(HashSet<Split> ignore) {
        double s = 0.0, e;
        Iterator<Split> it = theMap.keySet().iterator();
        for (int i=0; i<theMap.size(); i++) {
            Split p = it.next();
            if (ignore==null) {
                e = getSplitLength(p);
                s += e*e;
            }
            else {
                if (!ignore.contains(p)) {
                    e = getSplitLength(p);
                    s += e*e;
                }
            }
        }
        return s;
    }

    public double sumSquaredLengths(boolean ignorePendants) {
        double s = 0.0, e;
        Iterator<Split> it = theMap.keySet().iterator();
        for (int i=0; i<theMap.size(); i++) {
            Split p = it.next();
            if (p.isTerminal()==null) {
                e = getSplitLength(p);
                s += e*e;
            }
        }
        return s;
    }

    
    /** Translate taxa */
    public void translateTaxa(HashMap<String, String> translationInfo) throws AlgorithmError {
        HashMap<Split,Double> toAdd = new HashMap();
        Iterator<Split> it = theMap.keySet().iterator();
        for (int i=0; i<theMap.size(); i++) {
            Split p = it.next();
            Double x = theMap.get(p);
            //theMap.remove(p);
            Split q = p.translate(translationInfo);
            toAdd.put(q, x);
        }
        theMap.clear();
        theMap.putAll(toAdd);
    }

    /** Tricky: write to string */
    public String toString() {
        return toString(true);
    }
    public String toString(boolean includeEdgeLengths) {
        ArrayList<Split> orderedSplits = new ArrayList();
        orderedSplits.addAll(theMap.keySet());
        Collections.sort(orderedSplits);

        Split p=null;
        for (int i=0; i<orderedSplits.size(); i++) {
            p = orderedSplits.get(i);
            if (p.isTerminal()==null) {
                break;
            }
        }

        String s;
        if (p==null) {
            // Star tree
            s = "(";
            for (int i=0; i<orderedSplits.size(); i++) {
                p = orderedSplits.get(i);
                s += p.isTerminal();
                if (includeEdgeLengths) {
                    s += ":";
                    s += String.format("%7.7f", this.getSplitLength(p));
                }
                if (i==(orderedSplits.size()-1)) {
                    s += ");";
                }
                else {
                    s += ",";
                }
            }
        } // End star tree
        else {
            // Recursive string output
            orderedSplits.remove(p);
            HashSet<String> setA = p.getTaxonSubset(null);
            HashSet<String> setB = p.allTaxa();
            setB.removeAll(setA);

            s = "(";
            s += recursiveToString(orderedSplits,setA,includeEdgeLengths);
            if (includeEdgeLengths) {
                s += ":";
                s += String.format("%7.7f", 0.5*this.getSplitLength(p));
            }
            s += ",";
            s += recursiveToString(orderedSplits,setB,includeEdgeLengths);
            if (includeEdgeLengths) {
                s += ":";
                s += String.format("%7.7f", 0.5*this.getSplitLength(p));
            }
            s += ");";
        }
        return s;
    }

    private String recursiveToString(ArrayList<Split> splits, HashSet<String> envelope, boolean includeEdgeLengths) {

        /* Find "maximal" splits */
        HashSet<Split> toAdd = new HashSet<Split>();

        Iterator<Split> itQ, itP;
        Split q, p;
        boolean maximal, test;

        // loop thru' splits
        itQ = splits.iterator();
        for (int i=0; i<splits.size(); i++) {
            q = itQ.next();
            if (q.getTaxonSubset(envelope)!=null) {
                // q sits below the root clade
                maximal = true;

                itP = splits.iterator();
                for (int j=0; j<splits.size(); j++) {
                    p = itP.next();
                    if (i!=j) {
                        // Compare the splits
                        try {
                            test = q.partialOrder(p, envelope);
                            if (!test) maximal = false;
                        }
                        catch (AlgorithmException incomparable) {
                            // splits aren't nested -- nothing to do.
                        }
                    }
                } // End comparison of q with all others

                if (maximal) {
                    toAdd.add(q);
                }
            }
        }

        // toAdd contains the top level of splits

        // Remove these from the pool
        splits.removeAll(toAdd);

        // Initialize
        Split s;
        String theLabel, res="(";
        double x;
        HashSet<String> newEnvelope;

        // Create edges and vertices corresponding to toAdd
        Iterator<Split> it = toAdd.iterator();
        for (int i=0; i<toAdd.size(); i++) {
            s = it.next();
            x = this.getSplitLength(s);
            theLabel = s.isTerminal();

            if (theLabel==null) {
                newEnvelope = s.getTaxonSubset(envelope);
                res += recursiveToString(splits, newEnvelope, includeEdgeLengths);
            }
            else {
                res += theLabel;
            }
            if (includeEdgeLengths) {
                res += ":";
                res += String.format("%7.7f", x);
            }
            if (i<(toAdd.size()-1)) res += ",";
        }
        res += ")";
        return res;

    }

    /* METRICS -------------------------------------------------------------- */

    public static double[] computeEuclideanAndConeMetrics(TreeAsSplits treeA, TreeAsSplits treeB, boolean includeLeaves) {

        HashSet<Split> sharedSplits = Split.intersection(treeA.getSplits(), treeB.getSplits());
        Iterator<Split> it = sharedSplits.iterator();
        double total = 0.0, u, v, sharedLength;
        Split s;
        while (it.hasNext()) {
            s = it.next();
            if ((includeLeaves)||(s.isTerminal()==null)) {
                u = treeA.getSplitLength(s);
                v = treeB.getSplitLength(s);
                total += (u-v)*(u-v);                
            }
        }
        sharedLength = Math.sqrt(total);

        double res[] = new double[2];

        double abranch = Math.sqrt(treeA.sumSquaredLengths(sharedSplits));
        double bbranch = Math.sqrt(treeB.sumSquaredLengths(sharedSplits));
        //double x = abranch + bbranch;
        double x = sharedLength*sharedLength + abranch*abranch + bbranch*bbranch;

        res[0] = Math.sqrt(x); // Euc
        res[1] = Math.sqrt(x+Math.abs(2.0*abranch*bbranch)); // Cone

        return res;

    }

    /** Compute Robinson Foulds distance */
    public static int rfDistance(TreeAsSplits treeA, TreeAsSplits treeB) {
        int res = 0;
        Iterator<Split> it = treeA.getSplitIterator();
        while (it.hasNext()) {
            if (!(treeB.contains(it.next()))) res++;
        }
        it = treeB.getSplitIterator();
        while (it.hasNext()) {
            if (!(treeA.contains(it.next()))) res++;
        }
        return res;
    }
    
    /* STATIC UTILS --------------------------------------------------------- */
    
    /** Find consensus / intersection.
     NB make sure both tA and tB have same taxon sets! */
    public static TreeAsSplits intersection(TreeAsSplits tA, TreeAsSplits tB) {
        TreeAsSplits res = new TreeAsSplits(tA.getNumTaxa());
        Iterator<Split> it = tA.getSplitIterator();
        while (it.hasNext()) {
            Split p = it.next();
            if (tB.contains(p)) {
                res.add(p, 0.5*(tA.getSplitLength(p)+tB.getSplitLength(p)));
            }
        }
        return res;
    }
    
    /** Get symmetric difference */
    public static HashSet<Split> symmetricDifference(TreeAsSplits tA, TreeAsSplits tB){
        HashSet<Split> res = new HashSet();
        Iterator<Split> it = tA.getSplitIterator();
        while (it.hasNext()) {
            Split p = it.next();
            if (!tB.contains(p)) {
                res.add(p);
            }
        }
        it = tB.getSplitIterator();
        while (it.hasNext()) {
            Split p = it.next();
            if (!tA.contains(p)) {
                res.add(p);
            }
        }
        return res;
    }

}
