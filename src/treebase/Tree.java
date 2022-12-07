/**
 * Tree 
 *
 * Copyright (C) 2011 Tom M. W. Nye
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact the author at:  <tom.nye@ncl.ac.uk>
 * <http://www.mas.ncl.ac.uk/~ntmwn/>
 */
package treebase;

/**
 * Unrooted tree described as an extension of the Graph class i.e. in terms of
 * edges and vertices.
 * 
 */
import java.util.*;

public class Tree extends Graph {

    /**
     * Version number for serialization - increment when structural changes are
     * made
     */
    private static final long serialVersionUID = 0L;

    protected HashSet<String> taxa; // List of strings
    protected HashMap<Edge, Split> splits; // Map from edges to splits

    /* -------------------------------------------------------------------------- */

    /* Constructors */
    protected Tree() {
        super();
        taxa = new HashSet<String>();
    }

    /**
     * Make a tree from a string
     */
    public Tree(String theString) throws AlgorithmException {

        super();
        taxa = new HashSet<String>();

        theString = NewickStringUtils.removeTrailingBrackets(theString);

        // Create the first vertex
        Vertex v = addNewVertex("");

        // Now parse the string
        generateConnectionsFromNewick(v, theString);
    }

    /**
     * Same as constructor from a string -- just read string from file first
     */
    public Tree(java.io.File theFile) throws AlgorithmException {
        super();
        taxa = new HashSet<String>();

        // Read in a string from a file
        String theString = new String();
        String s;
        try {
            java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(theFile));
            // Read a string
            while ((s = br.readLine()) != null) {
                theString += s;
            }
        } catch (java.io.IOException ioex) {
            throw new AlgorithmException("Error reading newick file.");
        }

        theString = NewickStringUtils.removeTrailingBrackets(theString);

        // Create the first vertex
        Vertex v = addNewVertex("");

        // Now parse the string
        generateConnectionsFromNewick(v, theString);

    }

    /**
     * Generate vertices and edges from a newick string
     */
    private void generateConnectionsFromNewick(Vertex v, String theString) throws AlgorithmException {
        /* Split the string up into comma separated blocks, each with a branch length.
         For each block generate a new vertex and edges, and recursively parse each block */
        HashMap<String, Double> blocks = NewickStringUtils.parseToSubstrings(theString);

        if (blocks.size() == 0) {
            // The string is already minimal -- i.e. just a taxon name
            v.label = new String(theString);
            taxa.add(v.label);
            return;
        }

        // Loop thru' blocks
        HashSet<String> vertexStrings = new HashSet(blocks.keySet());
        Iterator<String> it = vertexStrings.iterator();
        for (int i = 0; i < vertexStrings.size(); i++) {
            String s = it.next();
            double x = (blocks.get(s)).doubleValue();
            Vertex w = addNewVertex("");
            connect(v, w, x);
            generateConnectionsFromNewick(w, s);
        }
    }

    /**
     * Constructor from a hashmap from splits to lengths
     */
    public Tree(HashMap<Split, Double> splitsToLengths) throws AlgorithmError {
        super();
        taxa = new HashSet<String>();
        if (splitsToLengths.size() == 0) {
            Vertex v = addNewVertex("");
            return;
        }
        // Pick *any* split to be root split
        Iterator<Split> itS = splitsToLengths.keySet().iterator();
        Split rootSplit = itS.next();
        buildFromSplits(splitsToLengths, rootSplit);
    }

    /**
     * Clone a tree
     */
    public Tree clone() {
        return this.clone(null);
    }

    /**
     * Clone a tree but also return the correspondence from old to new vertices.
     * Pass in null if this information is not needed, otherwise an empty map
     */
    public Tree clone(HashMap<Vertex, Vertex> corresp) {
        Tree t = new Tree();
        t.buildFromTemplate(this, corresp);
        return t;
    }

    /**
     * Copy from a template tree Pass in corresp=empty map get get info on how
     * vertices match between original and clone, otherwise a null
     */
    protected void buildFromTemplate(Tree t, HashMap<Vertex, Vertex> corresp) {
        buildFromTemplate((Graph) t, corresp);
        taxa = t.getTaxa();
    }

    protected void buildFromTemplate(Tree t) {
        buildFromTemplate(t, null);
    }

    /**
     * Direct copy method: use with caution!!!
     */
    private void setVerticesAndEdges(HashSet<Vertex> v, HashSet<Edge> e, HashSet<String> t) {
        vertices = v;
        edges = e;
        taxa = t;
    }

    public void forceCopy(Tree t) {
        t.setVerticesAndEdges(vertices, edges, taxa);
    }

    public static Tree copy(Tree t) {
        Tree newTree = new Tree();
        newTree.buildFromTemplate(t);
        return newTree;
    }

    /* -------------------------------------------------------------------------- */

    /* Util Methods */
    public int numTaxa() {
        return taxa.size();
    }

    public HashSet<String> getTaxa() {
        HashSet<String> h = new HashSet<String>();
        h.addAll(taxa);
        return h;
    }

    public Iterator<String> getTaxaIterator() {
        return taxa.iterator();
    }

    /**
     * Output as a string: but pick "central" vertex by deterministic method.
     * Central vertex forms the "root" of the Newick string.
     */
    public String toOrderedString() {
        /* Pick a vertex as the "root" of the string.
         Needs to be "deterministic" in order to produce same result each time. */
        Vertex v = pickCentralVertex(true);
        return toString(v, null, new CompareVerticesByOrdering());
    }

    /**
     * Output as a string, ignoring branch lengths. This will produce the same
     * string every time, and so can be used to compare topologies.
     */
    public String toTopologyString() {
        Vertex v = pickCentralVertex(true);
        return toString(v, null, false, new CompareVerticesByOrdering());
    }

    /**
     * Output as a string
     */
    public String toString() {
        /* Pick a vertex as the "root" of the string.
         This might produce a different result each time because the vertex might change. */
        Iterator<Vertex> it = vertices.iterator();
        Vertex v = null;
        while (it.hasNext()) {
            Vertex w = it.next();
            if (w.degree() == 2) {
                v = w;
                break;
            }
            if (w.degree() > 2) {
                v = w;
            }
        }
        if (v == null) {
            // No degree 2 or >2 vertices
            it = vertices.iterator();
            v = it.next();
        }
        return toString(v, null, null);
    }

    /**
     * Output as a string -- from a specific vertex and edge. fromV=null for
     * full string rooted at v. Otherwise fromV is a vertex connected to v: the
     * corresponding edge serves as root.
     */
    protected String toString(Vertex v, Vertex fromV, boolean showEdgeLengths, Comparator vertexComp) {

        String theString = "";
        double x;

        // Deal with case that current vertex is a leaf
        if (fromV == null) {
            if (v.degree() == 0) {
                return v.label + ";"; // In case tree is a single vertex!!!
            }
            if (v.degree() == 1) {
                // Horrible case: introduce false node half-way along branch
                Iterator<? extends Vertex> itN = v.neighboursIterator();
                Vertex w = itN.next();
                x = w.getNeighbourDistance(v);
                if (showEdgeLengths) {
                    theString = "(" + v.label + ":" + String.format("%7.7f", 0.5 * x) + "," + toString(w, v, showEdgeLengths, vertexComp) + ":" + String.format("%7.7f", 0.5 * x) + ");";
                } else {
                    theString = "(" + v.label + "," + toString(w, v, showEdgeLengths, vertexComp) + ");";
                }
                return theString;
            }
        } else {
            if (v.degree() == 1) {
                return v.label;
            }
        }

        /* OK, in situation where v is not a leaf */
        // Sort if vertex comp!=null. You need an arraylist of vertices.
        ArrayList<Vertex> orderedVertices = new ArrayList<Vertex>();
        orderedVertices.addAll(v.getNeighbours());
        if (fromV != null) {
            orderedVertices.remove(fromV);
        }
        if (vertexComp != null) {
            Collections.sort(orderedVertices, vertexComp);
        }

        // Loop thru neighbours other than fromV which we've already dealt with
        theString = "(";
        for (int i = 0; i < orderedVertices.size(); i++) {
            Vertex w = orderedVertices.get(i);
            String s = this.toString(w, v, showEdgeLengths, vertexComp);
            theString += s;
            if (showEdgeLengths) {
                theString += ":";
                theString += String.format("%7.7f", w.getNeighbourDistance(v));
            }
            if (i < (orderedVertices.size() - 1)) {
                theString += ",";
            }
        }
        theString += ")";
        if (fromV == null) {
            theString += ";";
        }

        return theString;
    }

    /**
     * Version which automatically includes branch lengths
     */
    protected String toString(Vertex v, Vertex fromV, Comparator vertexComp) {
        return toString(v, fromV, true, vertexComp);
    }

    /* Get an ordered list of edges. The order corresponds to the left-to-right
     order of edges in a string produced by toOrderedString. This should produce 
     the same order for two trees with the same topology. */
    public ArrayList<Graph.Edge> getOrderedListOfEdges() {
        Vertex v = pickCentralVertex(true);
        ArrayList<Graph.Edge> orderedEdges = new ArrayList<Graph.Edge>();
        getOrderedListOfEdges(v, null, new CompareVerticesByOrdering(), orderedEdges);
        return orderedEdges;
    }

    /**
     * Get ordered set of edges and return them in the return object "theEdges".
     * fromV=null for full string rooted at v. Otherwise fromV is a vertex
     * connected to v: the corresponding edge serves as root.
     */
    private void getOrderedListOfEdges(Vertex v, Vertex fromV, Comparator vertexComp, ArrayList<Graph.Edge> theEdges) {

        // Deal with case that current vertex is a leaf
        if (fromV == null) {
            if (v.degree() == 0) {
                return; // In case tree is a single vertex!!!
            }
            if (v.degree() == 1) {
                // Only one edge
                Iterator<Vertex> itN = v.neighboursIterator();
                Vertex w = itN.next();
                try {
                    theEdges.add(v.getEdge(w));
                } catch (AlgorithmException ex) {
                    System.out.println("Error ordering edges. This should not be possible. " + ex.getMessage());
                }
                return;
            }
        } else {
            if (v.degree() == 1) {
                return;
            }
        }

        /* OK, in situation where v is not a leaf */
        // Sort. You need an arraylist of vertices.
        ArrayList<Vertex> orderedVertices = new ArrayList<Vertex>();
        orderedVertices.addAll(v.getNeighbours());
        if (fromV != null) {
            orderedVertices.remove(fromV);
        }
        Collections.sort(orderedVertices, vertexComp);

        // Loop thru neighbours other than fromV which we've already dealt with
        for (int i = 0; i < orderedVertices.size(); i++) {
            Vertex w = orderedVertices.get(i);
            this.getOrderedListOfEdges(w, v, vertexComp, theEdges);
            try {
                theEdges.add(w.getEdge(v));
            } catch (AlgorithmException ex) {
                System.out.println("Error ordering edges. This should not be possible. " + ex.getMessage());
            }
        }
    }

    /**
     * Find out whether an edge is terminal (ends in a leaf)
     */
    public static boolean isEdgeTerminal(Edge e) {
        Vertex[] w = e.getVertices();
        if (w[0].degree() == 1) {
            return true;
        }
        if (w[1].degree() == 1) {
            return true;
        }
        return false;
    }

    /* Return internal vertices in a HashSet */
    public HashSet<Vertex> getInternalVertices() {
        HashSet<Vertex> vtces = new HashSet<Vertex>();
        Iterator<Vertex> it = vertices.iterator();
        for (int i = 0; i < vertices.size(); i++) {
            Vertex v = it.next();
            if (v.degree() != 1) {
                vtces.add(v);
            }
        }
        return vtces;
    }
    /* Return internal edgess in a HashSet */

    public HashSet<Edge> getInternalEdges() {
        HashSet<Edge> h = new HashSet<Edge>();
        Iterator<Edge> it = edges.iterator();
        for (int i = 0; i < edges.size(); i++) {
            Edge e = it.next();
            if (!isEdgeTerminal(e)) {
                h.add(e);
            }
        }
        return h;
    }

    /**
     * Recursively add leaves to a hash set
     */
    protected void getLeaves(Vertex v, Vertex fromw, HashSet<String> leaves) {
        if (v.degree() == 1) {
            leaves.add(v.label);
        }

        Iterator<Vertex> it = v.neighboursIterator();
        for (int i = 0; i < v.degree(); i++) {
            Vertex w = it.next();
            if (w != fromw) {
                getLeaves(w, v, leaves);
            }
        }
    }

    /* -------------------------------------------------------------------------- */

    /* Methods relating to sets of splits.
     Return the root edge corresponding to the root split passed in. */
    protected Edge buildFromSplits(HashMap<Split, Double> splitsToLengths, Split rootSplit) throws AlgorithmError {

        // Pick *any* split to start
        HashSet<Split> splitSet = new HashSet<Split>();
        splitSet.addAll(splitsToLengths.keySet());

        splitSet.remove(rootSplit); // Remove split from consideration
        taxa = rootSplit.allTaxa();

        // Sanity check on number of trivial splits
        Iterator<Split> itCheck = splitsToLengths.keySet().iterator();
        Split trivSplit;
        int checkSum = 0;
        for (int i = 0; i < splitsToLengths.size(); i++) {
            trivSplit = itCheck.next();
            if (trivSplit.isTerminal() != null) {
                checkSum++;
            }
        }
        if (checkSum != taxa.size()) {
            throw new AlgorithmError("Attempt to construct tree from splits failed: wrong number of terminal splits.");
        }

        // Sanity check on pairwise compatibility
        ArrayList<Split> testList = new ArrayList(splitsToLengths.keySet());
        if (!Split.testPairwiseCompatibility(testList)) {
            throw new AlgorithmError("Attempt to construct tree from splits failed: incompatible splits");
        }

        Vertex v1, v2;
        double x;

        v1 = addNewVertex("");
        v2 = addNewVertex("");
        x = (splitsToLengths.get(rootSplit)).doubleValue();
        connect(v1, v2, x);
        Edge rootEdge = null;
        try {
            rootEdge = v1.getEdge(v2);
        } catch (AlgorithmException anEx) {
            throw new AlgorithmError("Missing edge between two vertices just joined.");
        }

        HashSet<String> left = rootSplit.getTaxonSubset(null);
        HashSet<String> right = this.getTaxa();
        right.removeAll(left);

        recursiveBuildFromSplits(v1, left, splitSet, splitsToLengths);
        recursiveBuildFromSplits(v2, right, splitSet, splitsToLengths);

        // Deal with vertex labels
        if (left.size() == 1) {
            Iterator<String> it = left.iterator();
            String newLabel = it.next();
            v1.label = newLabel;
        }
        if (right.size() == 1) {
            Iterator<String> it = right.iterator();
            String newLabel = it.next();
            v2.label = newLabel;
        }

        return rootEdge;

    }

    /**
     * Recursively build from splits
     */
    private void recursiveBuildFromSplits(Vertex v, HashSet<String> rootClade, HashSet<Split> splitSet, HashMap<Split, Double> splitsToLengths) throws AlgorithmError {
        /* Find "maximal" splits */
        HashSet<Split> toAdd = new HashSet<Split>();

        Iterator<Split> itQ, itP;
        Split q, p;
        boolean maximal, test;

        // loop thru' splits
        itQ = splitSet.iterator();
        for (int i = 0; i < splitSet.size(); i++) {
            q = itQ.next();
            if (q.getTaxonSubset(rootClade) != null) {
                // q sits below the root clade
                maximal = true;

                itP = splitSet.iterator();
                for (int j = 0; j < splitSet.size(); j++) {
                    p = itP.next();
                    if (!p.equals(q)) {
                        // Compare the splits
                        try {
                            test = q.partialOrder(p, rootClade);
                            if (!test) {
                                maximal = false;
                            }
                        } catch (AlgorithmException incomparable) {
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
        splitSet.removeAll(toAdd);

        // Initialize
        Split s;
        Vertex w;
        double x;
        HashSet<String> newRootClade;

        // Create edges and vertices corresponding to toAdd
        Iterator<Split> it = toAdd.iterator();
        for (int i = 0; i < toAdd.size(); i++) {
            s = it.next();
            x = (splitsToLengths.get(s)).doubleValue();

            String label = s.isTerminal();
            if (label == null) {
                w = addNewVertex("");
            } else {
                w = addNewVertex(label);
            }

            connect(v, w, x);

            // Recursively call next level
            newRootClade = s.getTaxonSubset(rootClade);
            recursiveBuildFromSplits(w, newRootClade, splitSet, splitsToLengths);
        }

    }

    /**
     * Compute the hashmap from edges to splits
     */
    public void generateSplits() {
        /* Pick any vertex as the "root" of the process. */
        Iterator<Vertex> it = vertices.iterator();
        Vertex v = it.next();
        splits = new HashMap<Edge, Split>();
        generateSplits(v, null);

        /* Sort out repeated splits */
        ArrayList<Split> sortedSplits = new ArrayList<Split>();
        sortedSplits.addAll(splits.values());
        Collections.sort(sortedSplits);

        int k, pos;
        Split s;
        Edge e;

        /* Loop thru' edges
         If corresponding split occurs more than once then replace it*/
        Iterator<Edge> itE = edges.iterator();
        for (int i = 0; i < edges.size(); i++) {
            e = itE.next();
            s = splits.get(e);
            k = Collections.frequency(sortedSplits, s);
            if (k > 1) {
                // Occurs more than once!
                pos = sortedSplits.indexOf(s); // Find first occurance
                splits.remove(e);
                splits.put(e, sortedSplits.get(pos));
            }
        }
        // Done replacing
    }

    /**
     * Recursively generate splits down the tree
     */
    private HashSet<String> generateSplits(Vertex v, Vertex w) {

        Vertex nbr;
        Edge en;
        HashSet<String> h, comp;

        // Deal with case of leaf
        if ((v.degree() == 1) && (w != null)) {
            h = new HashSet<String>();
            h.add(v.label);
            return h;
        }

        // Deal with general case
        HashSet<String> hangingTaxa = new HashSet<String>();
        Iterator<Vertex> it = v.neighboursIterator();
        try {
            for (int i = 0; i < v.degree(); i++) {
                nbr = it.next();
                if ((w == null) || (nbr != w)) {
                    en = v.getEdge(nbr);
                    h = generateSplits(nbr, v);
                    comp = getTaxa();
                    comp.removeAll(h);
                    Split s = new Split(h, comp);
                    splits.put(en, s);
                    hangingTaxa.addAll(h);
                }
            }
        } catch (AlgorithmException anErr) {
            javax.swing.JOptionPane.showMessageDialog(null, "Serious problem generating splits from a tree: " + anErr.getMessage(), "Algorithm Error", javax.swing.JOptionPane.ERROR_MESSAGE);
            System.exit(1);
        }

        return hangingTaxa;

    }

    /**
     * Return split corresponding to a particular edge
     */
    public Split getSplit(Edge e) throws AlgorithmError {
        if (!edges.contains(e)) {
            throw new AlgorithmError("Request for split corresponding to an invalid edge. ");
        }
        if (splits == null) {
            generateSplits();
        }
        return splits.get(e);
    }

    /**
     * Get split associated to an edge without computing any other splits
     */
    public Split getSingleSplit(Edge e) throws AlgorithmError {
        if (!edges.contains(e)) {
            throw new AlgorithmError("Request for split corresponding to an invalid edge. ");
        }
        Vertex[] x = e.getVertices();
        HashSet<String> a = new HashSet();
        HashSet<String> b = new HashSet();
        getLeaves(x[0], x[1], a);
        getLeaves(x[1], x[0], b);
        return new Split(a, b);
    }

    /**
     * Get hashset of all splits
     */
    public HashSet<Split> getSplits() {
        if (splits == null) {
            generateSplits();
        }
        HashSet<Split> h = new HashSet<Split>();
        h.addAll(splits.values());
        return h;
    }

    /**
     * Finds the edge associated with a split
     */
    public Edge findEdge(Split s) {
        if (splits == null) {
            generateSplits();
        }
        Edge e;
        Iterator<Edge> it = edges.iterator();
        while (it.hasNext()) {
            e = it.next();
            if (splits.get(e).equals(s)) {
                return e;
            }
        }
        return null;
    }

    /**
     * Get a hashmap from splits to branch lengths. Useful for generating input
     * to tree constructors from splits
     */
    public HashMap<Split, Double> getSplitBranchLengths() {
        HashMap<Split, Double> res = new HashMap<Split, Double>();

        if (splits == null) {
            generateSplits();
        }

        Edge e;
        Split s;
        double x;

        // Loop thru edges
        Iterator<Edge> itE = edges.iterator();
        for (int i = 0; i < edges.size(); i++) {
            e = itE.next();
            s = splits.get(e);
            if (res.containsKey(s)) {
                x = ((Double) res.get(s)).doubleValue();
                res.remove(s);
                res.put(s, new Double(x + e.getLength())); // Add branch length to existing sum
            } else {
                res.put(s, new Double(e.getLength()));
            }
        }

        return res;
    }

    /** Indicate whether splits already computed */
    public boolean gotSplits() {
        return (splits == null);
    }


    /** Reset splits -- eg. call after topological operation */
    public void resetSplits() {
        // Force recalculation of all splits. 
        splits = null;
    }


    /* -------------------------------------------------------------------------- */
    
    /* Methods relating to distance matrices */
    
    /** Return distances between taxa. Taxa in lexographical order if orderedTaxa is null */
    public double[][] getDistanceMatrix (ArrayList<String> orderedTaxa) {
        
        if (orderedTaxa==null) {
            orderedTaxa = new ArrayList();
            orderedTaxa.addAll(taxa);
            Collections.sort(orderedTaxa, null);
        }
        
        int n = taxa.size();
        
        double[][] dist = new double[n][n];
        try {
            for (int i=0; i<n; i++) {
                String tA = orderedTaxa.get(i);
                for (int j=0; j<i; j++) {
                    String tB = orderedTaxa.get(j);

                    /* Loop thru' splits */
                    for (Graph.Edge e : edges) {
                        Split s = getSplit(e);
                        if (s.separatesTaxa(tA, tB)) {
                            dist[i][j] += e.getLength();
                            dist[j][i] = dist[i][j];
                        }
                    }
                }
            }
        }
        catch (AlgorithmException anErr) {
            System.out.println("Error computing distance matrix. This should not be possible.");
        }
        return dist;
    }
    
    /** Return derivative of distance matrix. */
    public double[][] getDistanceMatrixDerivative (Edge e, ArrayList<String> orderedTaxa) throws AlgorithmException {
        
        if (!edges.contains(e)) throw new AlgorithmException("Request to differentiate distance matrix with respect to bad edge length.");
        
        if (orderedTaxa==null) {
            orderedTaxa = new ArrayList();
            orderedTaxa.addAll(taxa);
            Collections.sort(orderedTaxa, null);
        }
        int n = taxa.size();
        Split s = getSplit(e);
        
        double[][] d = new double[n][n];
        for (int i=0; i<n; i++){
            String tA = orderedTaxa.get(i);
            for (int j=0; j<i; j++) {
                String tB = orderedTaxa.get(j); 
                if (s.separatesTaxa(tA, tB)) {
                        d[i][j] = 1;
                        d[j][i] = 1;
                }
            }
        }

        return d;
    }
    

    /* -------------------------------------------------------------------------- */

    /* Methods relating to paths */
    /**
     * Get max path length (used for graphics). Put topology flag true to treat
     * all branches as length 1.
     */
    public double getMaxLeafDistance(Vertex v, boolean topologyFlag) {
        return recursiveGetMaxLeafDistance(v, null, topologyFlag);
    }

    private double recursiveGetMaxLeafDistance(Vertex v, Vertex fromV, boolean topologyFlag) {
        double pathLength = 0.0;
        // Loop thru' neighbours
        HashSet<Vertex> nh = v.getNeighbours();
        if (fromV != null) {
            nh.remove(fromV);
        }
        Iterator<Vertex> it = nh.iterator();
        for (int i = 0; i < nh.size(); i++) {
            Vertex w = it.next();
            // max path length = time to child + max path length for child
            double t = recursiveGetMaxLeafDistance(w, v, topologyFlag);
            if (topologyFlag) {
                t += 1;
            } else {
                t += v.getNeighbourDistance(w);
            }
            if (t > pathLength) {
                pathLength = t;
            }
        }
        return pathLength;
    }

    /**
     * Compute path length between taxa
     */
    public double pathLength(String t1, String t2) throws AlgorithmException {
        // Algorithm exception thrown if species not found.
        if (splits == null) {
            generateSplits();
        }
        Edge e;
        Split s;
        double x = 0.0;
        Iterator<Edge> it = edges.iterator();
        for (int i = 0; i < edges.size(); i++) {
            e = it.next();
            s = splits.get(e);
            if (s.separatesTaxa(t1, t2)) {
                x += e.getLength();
            }
        }
        return x;
    }

    /**
     * A test method
     */
    public double debugPathLength(String t1, String t2) throws AlgorithmException {
        LinkedList edgePath = new LinkedList();
        // Find vertices associated to t1, t2
        Vertex v1 = null;
        Vertex v2 = null;
        Iterator it = vertices.iterator();
        for (int i = 0; i < vertices.size(); i++) {
            Vertex v = (Vertex) it.next();
            if (v.label.equals(t1)) {
                v1 = v;
            }
            if (v.label.equals(t2)) {
                v2 = v;
            }
        }
        getPath(v1, v2, null, edgePath);
        double x = 0.0;
        ListIterator itE = edgePath.listIterator();
        for (int i = 0; i < edgePath.size(); i++) {
            Edge e = (Edge) itE.next();
            x += e.getLength();
        }
        return x;
    }

    /**
     * Compute path between two vertices. vertexPath and edgePath return a
     * linked list of objects from v to w. Pass in null if you don't want either
     * of these.
     */
    public void getPath(Vertex v, Vertex w, LinkedList<Vertex> vertexPath, LinkedList<Edge> edgePath) throws AlgorithmException {

        // Check v,w are in the tree!?
        if ((!(vertices.contains(v))) || (!(vertices.contains(w)))) {
            throw new AlgorithmException("Attempt to find path between two vertices which don't both lie in a tree.");
        }

        // Initialize
        if (vertexPath != null) {
            vertexPath.clear();
        }
        if (edgePath != null) {
            edgePath.clear();
        }

        Iterator<Vertex> it = v.neighboursIterator();
        Vertex u;
        for (int i = 0; i < v.degree(); i++) {
            u = it.next();
            recursiveGetPath(u, v, w, vertexPath, edgePath);
        }

        if (vertexPath != null) {
            vertexPath.addFirst(v);
        }
    }

    /**
     * Recursive call to find path: used by getPath
     */
    private boolean recursiveGetPath(Vertex v, Vertex from, Vertex target, LinkedList<Vertex> vertexPath, LinkedList<Edge> edgePath) {

        if (v.equals(target)) {
            // We're done: at end of the path
            // Add in vertex and edge, and return
            if (vertexPath != null) {
                vertexPath.addLast(v);
            }
            if (edgePath != null) {
                try {
                    edgePath.addLast(v.getEdge(from));
                } catch (AlgorithmException anErr) {
                    javax.swing.JOptionPane.showMessageDialog(null, "Serious problem finding path on a tree: " + anErr.getMessage(), "Algorithm Error", javax.swing.JOptionPane.ERROR_MESSAGE);
                    System.exit(1);
                }
            }
            return true;
        }

        // Otherwise seek thru' vertices attached to v
        Iterator<Vertex> it = v.neighboursIterator();
        Vertex w;
        for (int i = 0; i < v.degree(); i++) {
            w = it.next();
            if (w != from) {
                if (recursiveGetPath(w, v, target, vertexPath, edgePath)) {
                    // w lies on the path: add in edge and vertex, and return
                    if (vertexPath != null) {
                        vertexPath.addFirst(v);
                    }
                    if (edgePath != null) {
                        try {
                            edgePath.addFirst(v.getEdge(from));
                        } catch (AlgorithmException anErr) {
                            javax.swing.JOptionPane.showMessageDialog(null, "Serious problem finding path on a tree: " + anErr.getMessage(), "Algorithm Error", javax.swing.JOptionPane.ERROR_MESSAGE);
                            System.exit(1);
                        }
                    }
                    return true;
                }
            }
        }
        return false;
    }

    /* -------------------------------------------------------------------------- */

    /* Core methods to support all tree topological operation classes */
    /**
     * Add a vertex -- check for new taxon
     */
    @Override
    public Vertex addNewVertex(String theLabel) {
        Vertex v = new Vertex(theLabel);
        vertices.add(v);
        if (theLabel != null) {
            if (!theLabel.isEmpty()) {
                taxa.add(theLabel);
            }
        }
        return v;
    }

    /**
     * Remove a vertex -- delete taxon label if necessary
     */
    public void remove(Vertex v) throws AlgorithmError {
        super.remove(v);
        if (v.label != null) {
            if (!v.label.isEmpty()) {
                taxa.remove(v.label);
            }
        }
    }

    /**
     * Get a list of vertices either side of an edge. Use to obtain vertices for
     * NNI. The return array is a set of vertices. w[0] are vertices on one side
     * of the edge, w[1] on the other. w[0][0] is the vertex on the edge while
     * w[0][1], w[0][2] etc are the adjacent vertices. This assumes e is a fully
     * resolved internal edge so w will have dimension [2][3]
     */
    public Vertex[][] getVerticesAdjacentToEdge(Edge e) throws AlgorithmException {
        Vertex[] v = e.getVertices();
        int d1 = v[0].degree();
        int d2 = v[1].degree();
        if ((d1 != 3) || (d2 != 3)) {
            throw new AlgorithmException("Request to get adjacent vertices for non resolved edge.");
        }

        Vertex[][] w = new Vertex[2][3];
        w[0][0] = v[0];
        w[1][0] = v[1];

        Iterator<Vertex> it = v[0].neighboursIterator();
        int k = 1;
        for (int i = 0; i < 3; i++) {
            Vertex u = it.next();
            if (u != v[1]) {
                w[0][k] = u;
                k++;
            }
        }
        it = v[1].neighboursIterator();
        k = 1;
        for (int i = 0; i < 3; i++) {
            Vertex u = it.next();
            if (u != v[0]) {
                w[1][k] = u;
                k++;
            }
        }
        return w;
    }

    /**
     * Does a subtree (defined by v and its neighbour fromV) contain e? Used by
     * sub-tree prune and re-graft.
     */
    protected boolean subtreeContainsEdge(Vertex v, Vertex fromV, Edge e) {
        Vertex w;
        Edge f;
        Iterator<Vertex> it = v.neighboursIterator();
        for (int i = 0; i < v.degree(); i++) {
            w = it.next();
            if (!(w.equals(fromV))) {
                try {
                    f = v.getEdge(w);
                    if (f.equals(e)) {
                        return true;
                    }
                } catch (AlgorithmException anErr) {
                    javax.swing.JOptionPane.showMessageDialog(null, "Serious problem propagating down a subtree: " + anErr.getMessage(), "Algorithm Error", javax.swing.JOptionPane.ERROR_MESSAGE);
                    System.exit(1);
                }
                if (subtreeContainsEdge(w, v, e)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Delete a subtree. The subtree is specified by an edge and a vertex in
     * that edge: the subtree hanging from the vertex is deleted. The edge is
     * assigned length stumpLength, and the vertex passed in is assigned a new
     * taxon name. If stumpName is null then the vertex and edge are also
     * deleted.
     */
    public void deleteSubTree(Vertex rootVertex, Edge theEdge, String stumpName, double stumpLength) throws AlgorithmError {
        if (!rootVertex.getEdges().contains(theEdge)) {
            throw new AlgorithmError("Request to delete a subtree using bad info.");
        }

        Vertex[] w = theEdge.getVertices();
        if (w[0] == rootVertex) {
            recursiveDelete(rootVertex, w[1]);
        } else if (w[1] == rootVertex) {
            recursiveDelete(rootVertex, w[0]);
        } else {
            throw new AlgorithmError("Request to delete a subtree using bad info.");
        }

        if (stumpName != null) {
            rootVertex.label = stumpName;
            theEdge.setLength(stumpLength);
            taxa.add(stumpName);
        } else {
            cut(theEdge);
            vertices.remove(rootVertex);
        }
        generateSplits();
    }

    private void recursiveDelete(Vertex v, Vertex fromV) {
        if (v.degree() == 1) {
            taxa.remove(v.label);
        }
        try {
            Iterator<Vertex> it = v.getNeighbours().iterator();
            while (it.hasNext()) {
                Vertex w = it.next();
                if (w != fromV) {
                    recursiveDelete(w, v);
                    cut(v.getEdge(w));
                    vertices.remove(w);
                }
            }
        } catch (AlgorithmException anErr) {
            System.out.println("Error recursively deleting a subtree. This should not be possible!\n" + anErr.getMessage());
        }
    }

    /* -------------------------------------------------------------------------- */

    /* Additional non-core methods */

    /* Get an ordered list of vertices.
     Util method that should produce the same order for two trees with the same topology. */
    public ArrayList<Vertex> getOrderedListOfVertices() {

        class CompareVerticesByLabel implements Comparator {

            public int compare(Object v1, Object v2) {
                return ((Vertex) v1).label.compareTo(((Vertex) v2).label);
            }
        }

        ArrayList<Vertex> theList = new ArrayList();

        // Get leaves and sort
        Iterator<Vertex> it = getVertexIterator();
        for (int i = 0; i < numVertices(); i++) {
            Vertex v = it.next();
            if (v.degree() == 1) {
                theList.add(v);
            }
        }
        Collections.sort(theList, new CompareVerticesByLabel());

        // Add vertices in order
        Vertex toAdd;
        do {
            // Loop thru pairs in currentVertices
            Vertex v1 = null, v2 = null;
            toAdd = null;

            pairSearch:
            for (int i = 0; i < theList.size(); i++) {
                v1 = theList.get(i);
                for (int j = (i + 1); j < theList.size(); j++) {
                    v2 = theList.get(j);

                    // Are vertices v1 and v2 joined by a mutual neighbour?
                    Iterator<Vertex> itV = v1.neighboursIterator();
                    for (int k = 0; k < v1.degree(); k++) {
                        Vertex neighbour = itV.next();
                        if ((neighbour.isConnectedTo(v2)) && (!theList.contains(neighbour))) {
                            // Mutual neighbour!
                            toAdd = neighbour;
                            break pairSearch;
                        }
                    }

                }
            }

            if (toAdd != null) {
                theList.add(toAdd);
                //System.out.println("Index "+theList.indexOf(toAdd)+ " formed from "+v1.label+" index "+theList.indexOf(v1)+" to "+v2.label+" index "+theList.indexOf(v2));
            }
        } while (toAdd != null);
        // Sanity check:
        if (theList.size() != numVertices()) {
            System.out.println("Algorithm Error: Failed sanity check in Tree.getOrderedListOfVertices.");
        }

        return theList;
    }

    /**
     * A class for ordering vertices. Used by toOrderedString etc.
     */
    protected class CompareVerticesByOrdering implements Comparator {

        private ArrayList<Vertex> orderedList;

        public CompareVerticesByOrdering() {
            orderedList = getOrderedListOfVertices();
        }

        public int compare(Object v1, Object v2) {
            if (!orderedList.contains((Vertex) v1)) {
                System.out.println("AlgorithmError: request to compare two vertices not in tree?");
            }
            if (!orderedList.contains((Vertex) v2)) {
                System.out.println("AlgorithmError: request to compare two vertices not in tree?");
            }
            if (orderedList.indexOf((Vertex) v1) == orderedList.indexOf((Vertex) v2)) {
                return 0;
            } else if (orderedList.indexOf((Vertex) v1) < orderedList.indexOf((Vertex) v2)) {
                return -1;
            } else {
                return 1;
            }
        }
    }

    /**
     * Find a vertex in the middle of the tree. Used for "rooting" for
     * outputting as a string.
     */
    public Vertex pickCentralVertex(boolean topologyFlag) {

        // First seek a real root -- degree 2 vertex
        Iterator<Vertex> it = getVertexIterator();
        for (int i = 0; i < vertices.size(); i++) {
            Vertex v = it.next();
            if (v.degree() == 2) {
                return v;
            }
        }
        // If we're here then there was no deg 2 vertex.

        /* Loop thru internal vertices
         - work out max distance from a leaf
         - choose vertex which minimises this */
        ArrayList<Vertex> intVertices = new ArrayList();
        ArrayList<Vertex> orderedVertices = getOrderedListOfVertices();
        intVertices.addAll(orderedVertices.subList(numTaxa(), vertices.size()));

        Vertex centralVertex = null;
        double minDist = 0.0;
        for (int i = 0; i < intVertices.size(); i++) {
            Vertex v = intVertices.get(i);
            double x = getMaxLeafDistance(v, topologyFlag);
            if ((i == 0) || (x < minDist)) {
                centralVertex = v;
                minDist = x;
            }
        }

        return centralVertex;
    }

    /* -------------------------------------------------------------------------- */
}
