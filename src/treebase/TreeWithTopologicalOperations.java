/*
 * TreeWithTopologicalOperations
   Copyright (C) 2013 Tom M. W. Nye

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

    Contact the author at: <tom.nye@ncl.ac.uk>
 */

package treebase;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * Extends Tree to incorporate topological operations
 * NB: always recalculate splits after any topological operation!
 */

public class TreeWithTopologicalOperations extends Tree implements TopologicalOperations {

    /* Constructors */

    /** No arg constructor for use by extending classes */
    protected TreeWithTopologicalOperations() {
        super();
    }

    /** Constructors based on Newick strings */
    public TreeWithTopologicalOperations(String theString) throws AlgorithmException {
        super(theString);
    }
    public TreeWithTopologicalOperations(java.io.File theFile) throws AlgorithmException {
        super(theFile);
    }

    /** Constructor from a hashmap from splits to lengths */
    public TreeWithTopologicalOperations(HashMap<Split,Double> splitsToLengths) throws AlgorithmError {
        super(splitsToLengths);
    }

    /** Constructor from a Tree  */
    public TreeWithTopologicalOperations(Tree t) {
        super();
        this.buildFromTemplate(t);
    }


    /** Clone  */
    public TreeWithTopologicalOperations clone() {
        return this.clone(null);
    }

    /** Clone a tree but also return the correspondence from old to new vertices.
     Pass in null if this information is not needed, otherwise an empty map */
    public TreeWithTopologicalOperations clone(HashMap<Vertex,Vertex> corresp) {
        TreeWithTopologicalOperations t = new TreeWithTopologicalOperations();
        t.buildFromTemplate(this, corresp);
        return t;
    }


    /** Perform NNI on an edge.
     NB: THIS ASSUMES FULLY RESOLVED TREE!
     Pass in a random bit for choice to get a random choice of the two possibilities.
     If you want a particular choice of NNI then use the method below which takes an
     array of vertices as an argument.
     */
    public Edge performNNI(Edge e, boolean choice, double newLength) throws AlgorithmException {
        Vertex[][] v = getVerticesAdjacentToEdge(e);
        return performNNI(e, v, newLength, choice);
    }

    /** Perform NNI on an edge.
     NB: THIS ASSUMES FULLY RESOLVED TREE!
     Specify the resultant split you'd like.
     */
    public Edge performNNI(Edge e, Split p, double newLength) throws AlgorithmException {
        if (!edges.contains(e)) throw new AlgorithmError("Bad edge sent to request to perform NNI.");
        Vertex[][] v = getVerticesAdjacentToEdge(e);
        // Check how array v matches p
        HashSet<String> a0 = new HashSet<String>();
        getLeaves(v[0][1],v[0][0],a0);
        HashSet<String> b0 = new HashSet<String>();
        getLeaves(v[0][2],v[0][0],b0);
        HashSet<String> a1 = new HashSet<String>();
        getLeaves(v[1][1],v[1][0],a1);
        HashSet<String> b1 = new HashSet<String>();
        getLeaves(v[1][2],v[1][0],b1);

        Split q = makeSplitFromFourSubsets(a0,b0,a1,b1);
        Split r = makeSplitFromFourSubsets(a0,b1,a1,b0);

        if (p.equals(q)) {
            return performNNI(e, v, newLength, true);
        }
        else if(p.equals(r)) {
            return performNNI(e, v, newLength, false);
        }
        else {
            throw new AlgorithmError("Bad split sent to request to perform NNI.");
        }
    }

    /** Perform an NNI on an edge.
     NB: THIS ASSUMES FULLY RESOLVED TREE!
     Pass in the array of vertices obtained from getVerticesAdjacentToEdge.
     This array is modified by the method so it can be passed back to reverseNNI.
     */
    @Override
    public Edge performNNI(Edge e, Vertex[][] v, double newLength, boolean choice) throws AlgorithmException {
        if (!edges.contains(e)) throw new AlgorithmError("Bad edge sent to request to perform NNI.");

        // If the vertex at either end of e is of degree two, throw exception
        if(e.vertex1.degree()==2 || e.vertex2.degree()==2) throw new AlgorithmException("Cannot perform NNI on edge connected to a degree 2 vertex.");
        
        if (splits!=null) {
            splits.remove(e);
        }

        // Cut the edge
        this.cut(e);
        // Add two new vertices
        Vertex w0 = this.addNewVertex("");
        Vertex w1 = this.addNewVertex("");
        Edge newEdge = this.connect(w0, w1, newLength);

        /* Now connect up the vertices */
        Edge f;
        if (choice) {
            f = v[0][0].getEdge(v[0][1]);
            f.changeVertex(v[0][0], w0);
            f = v[1][0].getEdge(v[1][1]);
            f.changeVertex(v[1][0], w0);
            f = v[0][0].getEdge(v[0][2]);
            f.changeVertex(v[0][0], w1);
            f = v[1][0].getEdge(v[1][2]);
            f.changeVertex(v[1][0], w1);
        }
        else {
            f = v[0][0].getEdge(v[0][1]);
            f.changeVertex(v[0][0], w0);
            f = v[1][0].getEdge(v[1][2]);
            f.changeVertex(v[1][0], w0);
            f = v[0][0].getEdge(v[0][2]);
            f.changeVertex(v[0][0], w1);
            f = v[1][0].getEdge(v[1][1]);
            f.changeVertex(v[1][0], w1);
        }
        this.remove(v[0][0]);
        this.remove(v[1][0]);
        v[0][0] = w0;
        v[1][0] = w1;
        if (choice) {
            Graph.Vertex w = v[1][1];
            v[1][1] = v[0][2];
            v[0][2] = w;
        } else {
            Graph.Vertex w = v[1][2];
            v[1][2] = v[0][2];
            v[0][2] = w;
        }
        
        // Update splits
        if (splits!=null) {
            splits.put(newEdge, getSingleSplit(newEdge));
        }

        return newEdge;
    }

    /** Reverse the NNI on an edge.
     NB: THIS ASSUMES FULLY RESOLVED TREE!
     The vertex[][] should be that passed to performNNI (and modified by that method).
     The boolean should be the same as that passed to performNNI in order to reverse the move.
     The edge returned is the new copy of oldEdge.
     */
    @Override
    public Edge reverseNNI(Edge oldEdge, Edge newEdge, Vertex[][] v, boolean choice) throws AlgorithmException {
        if (!edges.contains(newEdge)) throw new AlgorithmError("Bad edge sent to request to reverse NNI.");

        // If the vertex at either end of e is of degree two, throw exception
        if(newEdge.vertex1.degree()==2 || newEdge.vertex2.degree()==2) throw new AlgorithmException("Cannot reverse NNI on edge connected to a degree 2 vertex. This should not be possible.");

        if (splits!=null) {
            splits.remove(newEdge);
        }
 
        // Cut the edge
        this.cut(newEdge);
        // Add two new vertices
        Vertex w0 = this.addNewVertex(""); //when I come to override, addNewVertex can contain the sequences which
        Vertex w1 = this.addNewVertex(""); //I'll be able to get from oldEdge
        Edge oldEdgeCopy = this.connect(w0, w1, oldEdge.getLength()); //when I come to override, connect can contain
                                                   //the subs history which I can get from oldEdge

        // Now connect up the vertices
        Edge f;
        if (choice) {
            f = v[0][0].getEdge(v[0][1]);
            f.changeVertex(v[0][0], w0);
            f = v[1][0].getEdge(v[1][1]);
            f.changeVertex(v[1][0], w0);
            f = v[0][0].getEdge(v[0][2]);
            f.changeVertex(v[0][0], w1);
            f = v[1][0].getEdge(v[1][2]);
            f.changeVertex(v[1][0], w1);
        }
        else {
            f = v[0][0].getEdge(v[0][1]);
            f.changeVertex(v[0][0], w0);
            f = v[1][0].getEdge(v[1][2]);
            f.changeVertex(v[1][0], w0);
            f = v[0][0].getEdge(v[0][2]);
            f.changeVertex(v[0][0], w1);
            f = v[1][0].getEdge(v[1][1]);
            f.changeVertex(v[1][0], w1);
        }
        this.remove(v[0][0]);
        this.remove(v[1][0]);
        
        // Update splits
        if (splits!=null) {
            splits.put(oldEdgeCopy, getSingleSplit(oldEdgeCopy));
        }

        return oldEdgeCopy;
    }

    /** Get splits which are NNI's of a fully resolved edge */
    public Split[] getNNISplits(Edge e) throws AlgorithmException {
        Vertex[][] v = getVerticesAdjacentToEdge(e);
        HashSet<String> a0 = new HashSet<String>();
        getLeaves(v[0][1],v[0][0],a0);
        HashSet<String> b0 = new HashSet<String>();
        getLeaves(v[0][2],v[0][0],b0);
        HashSet<String> a1 = new HashSet<String>();
        getLeaves(v[1][1],v[1][0],a1);
        HashSet<String> b1 = new HashSet<String>();
        getLeaves(v[1][2],v[1][0],b1);

        Split[] p = new Split[2];
        p[0] = makeSplitFromFourSubsets(a0,b0,a1,b1);
        p[1] = makeSplitFromFourSubsets(a0,b1,a1,b0);
        return p;
    }


    /** Carry out an SPR.
     This operates at the level of topology: edge lengths are not handled carefully.
     The graft edge is simply divided in two.
     The method does not report which vertices or edges have been created, so subsequent
     adjustment of edge lengths is hard to achieve.
     Regard this method as a means of generating the tree topology for a given LGT.
     See methods below which deal more carefully with edge lengths and which are
     intended for MCMC schemes for inferring phylogenies.
     */
    public void performTopologicalSPR(Edge pruneEdge, Edge graftEdge) {
        try {
            cut(pruneEdge);

            /* Work out which vertex on the prune edge attaches to the "main tree"
             (ie. the component of the graph containing the graft edge)
             and which vertex is on the pruned sub-tree.
             */

            Vertex[] v = pruneEdge.getVertices();
            Vertex subtreeVertex = null;
            Vertex maintreeVertex = null;
            if (this.subtreeContainsEdge(v[0], v[1], graftEdge)) {
                subtreeVertex = v[1];
                maintreeVertex = v[0];
            }
            else if (this.subtreeContainsEdge(v[1], v[0], graftEdge)) {
                subtreeVertex = v[0];
                maintreeVertex = v[1];
            }
            else {
                javax.swing.JOptionPane.showMessageDialog(null,"Serious problem performing SPR: subtree not found.","Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
                System.exit(1);
            }

            // OK got the vertices

            // Subdivide the graft edge
            Vertex attachVertex = this.divide(graftEdge, "");
            // Connect the new vertex to the pruned tree
            connect(subtreeVertex, attachVertex, pruneEdge.getLength());
            // Absorb degree two vertices if necessary
            if (maintreeVertex.degree()==2) {
                removeDegreeTwoVertex(maintreeVertex);
            }

        }
        catch (AlgorithmException anErr) {
                javax.swing.JOptionPane.showMessageDialog(null,"Unable to perform SPR: "+anErr.getMessage(),"Algorithm Error",javax.swing.JOptionPane.ERROR_MESSAGE);
        }
    }


    /** Carry out a reversible SPR move on a fully resolved tree. You can specify new lengths
        for the four new edges created during the move. Note that regrafting the subtree on
        edge E creates two new edges on either side of the attached subtree. The edge which defines
        the same split as the original edge, E, is given length lengthSplitEdge1. The other is given
        the length lengthSplitEdge2. NOTE: slight inconsistency in the case that the prune and graft
        edges are adjacent. The four new branch lengths are not well defined: lengthMergedEdge and
        lengthSplitEdge2 refer to the same edge. The former length is used.
        The method return an array of edges of dimension 4 x 2. Components 0-3 in the first column contain
        edges from the original tree: the two edges merged into one where the subtree is plucked
        (entry 0 is the edge defining the same split as the merged edge in the new tree, entry 1 is the
        other); the edge which is split in to two (due to regrafting) and the edge which is pruned.
        Components 0-3 in the second column contain edges from the new tree: the new edges on either side
        of the newly attached subtree (in the order above), and the graft edge and the prune edges for the
        reverse move.*/
    @Override
    public Edge[][] performSPR(Edge pruneEdge, Edge graftEdge, double lengthMergedEdge, double lengthSplitEdge1, double lengthSplitEdge2, double lengthPrunedEdge) throws AlgorithmError {
        //performSPR(pruneEdge, graftEdge, something+something, graftEdge.getLength()/2.0, graftEdge.getLength()/2.0, pruneEdge.getLength());
        Edge[][] edgeArray = new Edge[4][2];
        // Record the edges which are plucked off and split in to two
        edgeArray[2][0] = graftEdge;
        edgeArray[3][0] = pruneEdge;
        // Are the edges adjacent?
        boolean areAdjacent = pruneEdge.getNeighbouringEdges().contains(graftEdge);
        try {
            // Work out which vertex on the prune edge attaches to the "main tree"
            // (ie. the component of the graph containing the graft edge)
            // and which vertex is on the pruned sub-tree.
            Vertex[] v = pruneEdge.getVertices();
            Vertex subtreeVertex = null;
            Vertex maintreeVertex = null;
            if (this.subtreeContainsEdge(v[0], v[1], graftEdge)) {
                subtreeVertex = v[1];
                maintreeVertex = v[0];
            }
            else if (this.subtreeContainsEdge(v[1], v[0], graftEdge)) {
                subtreeVertex = v[0];
                maintreeVertex = v[1];
            }
            else {
                throw new AlgorithmError("Serious problem performing SPR: subtree not found.");
            }

            Iterator<Edge> itE = maintreeVertex.getEdges().iterator();
            Edge eM1 = itE.next();
            Edge eM3 = itE.next();
            if(eM1==pruneEdge) eM1 = itE.next();
            if(eM3==pruneEdge) eM3 = itE.next();
            // Find splits on the edges which will be merged and the split defined by the grafting egde
            Split originalSplitM1 = this.getSingleSplit(eM1);
            Split originalSplitM3 = this.getSingleSplit(eM3);
            Split originalSplitS = this.getSingleSplit(graftEdge);

            // Pluck off the subtree
            cut(pruneEdge);

            // OK got the vertices

            // Subdivide the graft edge
            Vertex attachVertex = this.divide(graftEdge, "");
            itE = attachVertex.getEdges().iterator();
            Edge eS1 = itE.next();
            Edge eS3 = itE.next();

            // Connect the new vertex to the pruned tree
            connect(subtreeVertex, attachVertex, lengthPrunedEdge);

            // Record the prune edge for the reverse move:
            edgeArray[3][1] = attachVertex.getEdge(subtreeVertex);

            // Absorb degree two vertex
            Iterator<Vertex> it = maintreeVertex.neighboursIterator();
            Vertex v1 = it.next();
            Vertex v3 = it.next();
            // Remove the degree two vertex and recale the resulting merged edge
            removeDegreeTwoVertex(maintreeVertex);
            v1.getEdge(v3).setLength(lengthMergedEdge);
            // Record the graft edge for the reverse move
            edgeArray[2][1] = v1.getEdge(v3);

            // Now that the new tree has been constructed, find the split on new
            // merged edge and the splits on edges eS1 and eS3
            Split newSplitS1, newSplitS3;
            if(areAdjacent) {
                if(edges.contains(eS1)) newSplitS1 = this.getSingleSplit(eS1);
                else newSplitS1 = this.getSingleSplit(v1.getEdge(v3));
                if(edges.contains(eS3)) newSplitS3 = this.getSingleSplit(eS3);
                else newSplitS3 = this.getSingleSplit(v1.getEdge(v3));
            } else {
                newSplitS1 = this.getSingleSplit(eS1);
                newSplitS3 = this.getSingleSplit(eS3);
            }
            if(originalSplitS.equals(newSplitS1)) {
                // If true, the edge eS1 defines the same split as the original
                // edge prior to being split in to two. Otherwise, this should
                // be true of the edge eS3
                eS1.setLength(lengthSplitEdge1);
                eS3.setLength(lengthSplitEdge2);
                edgeArray[0][1] = eS1;
                edgeArray[1][1] = eS3;
            } else if(originalSplitS.equals(newSplitS3)) {
                eS1.setLength(lengthSplitEdge2);
                eS3.setLength(lengthSplitEdge1); //lengthSplitEdge1 applies to the edge defining the same split
                edgeArray[0][1] = eS3;
                edgeArray[1][1] = eS1;
            } else throw new AlgorithmError("Problem performing SPR: matching splits not found. This should not be possible.");

            Split newSplitM = this.getSingleSplit(v1.getEdge(v3));
            if(originalSplitM1.equals(newSplitM)) {
                // If true, the merged edge on the new tree defines the same
                // split as the original edge eM1. Otherwise,  this is true of
                // the original edge eM3.
                edgeArray[0][0] = eM1; //edgeArray[0] applies to the edge defining the same split
                edgeArray[1][0] = eM3;
            } else if(originalSplitM3.equals(newSplitM)) {
                edgeArray[0][0] = eM3;
                edgeArray[1][0] = eM1;
            } else throw new AlgorithmException("Problem performing SPR: matching splits not found. This should not be possible.");


        }
        catch (AlgorithmException anErr) {
                throw new AlgorithmError("Algorithm Error: Unable to perform SPR: "+anErr.getMessage());
        }

        return edgeArray;
    }

    /* Reverses the SPR move which returns an array of edges. Pass this array into
       reverseSPR. The edges returned by this method are the four edges created during
       the move.*/
    @Override
    public Edge[] reverseSPR(Edge[][] edgeArray) throws AlgorithmError {
        Edge mergedEdge1 = edgeArray[0][0];
        Edge mergedEdge2 = edgeArray[1][0];
        Edge splitEdge = edgeArray[2][0];
        Edge prunedEdge = edgeArray[3][0];
        Edge reverseGraftEdge = edgeArray[2][1];
        Edge reversePruneEdge = edgeArray[3][1];

        Edge[] originalEdges = new Edge[4];

        try {
            /* Work out which vertex on the prune edge attaches to the "main tree"
             (ie. the component of the graph containing the graft edge)
             and which vertex is on the pruned sub-tree.
             */
            Vertex[] v = reversePruneEdge.getVertices();
            Vertex subtreeVertex = null;
            Vertex maintreeVertex = null;
            if (this.subtreeContainsEdge(v[0], v[1], reverseGraftEdge)) {
                subtreeVertex = v[1];
                maintreeVertex = v[0];
            }
            else if (this.subtreeContainsEdge(v[1], v[0], reverseGraftEdge)) {
                subtreeVertex = v[0];
                maintreeVertex = v[1];
            }
            else {
                throw new AlgorithmError("Serious problem reversing SPR: subtree not found.");
            }

            // Find split defined by the grafting egde
            Split originalSplitS = this.getSingleSplit(reverseGraftEdge);

            // Pluck off the subtree
            cut(reversePruneEdge);

            // OK got the vertices

            // Subdivide the graft edge
            Vertex attachVertex = this.divide(reverseGraftEdge, "");
            Iterator<Edge> itE = attachVertex.getEdges().iterator();
            Edge eS1 = itE.next();
            Edge eS3 = itE.next();

            // Connect the new vertex to the pruned tree
            originalEdges[3] = connect(subtreeVertex, attachVertex, prunedEdge.getLength());

            // Now the pruned subtree is attached, find which of the two edges
            // eS1 and eS3 defines the same split as the original edge
            Split newSplitS1 = this.getSingleSplit(eS1);
            Split newSplitS3 = this.getSingleSplit(eS3);
            if(originalSplitS.equals(newSplitS1)) {
                // If true, the edge eS1 defines the same split as the original
                // edge. Otherwise, this should true of the edge eS3
                eS1.setLength(mergedEdge1.getLength());
                eS3.setLength(mergedEdge2.getLength());
                originalEdges[0] = eS1;
                originalEdges[1] = eS3;
            } else if(originalSplitS.equals(newSplitS3)) {
                eS1.setLength(mergedEdge2.getLength());
                eS3.setLength(mergedEdge1.getLength()); //mergedEdge1 applies to the edge defining the same split
                originalEdges[0] = eS3;
                originalEdges[1] = eS1;
            } else throw new AlgorithmError("Problem performing SPR: matching splits not found. This should not be possible.");

            // Absorb degree two vertices if necessary
            Iterator<Vertex> it = maintreeVertex.neighboursIterator();
            Vertex v1 = it.next();
            Vertex v3 = it.next();
            // Remove the degree two vertex and recale the resulting merged edge
            removeDegreeTwoVertex(maintreeVertex);
            originalEdges[2] = v1.getEdge(v3);
            originalEdges[2].setLength(splitEdge.getLength());

        }
        catch (AlgorithmException anErr) {
                throw new AlgorithmError("Algorithm Error: Unable to reverse SPR: "+anErr.getMessage());
        }

        return originalEdges;
    }

    @Override
    public Edge[][] performRootMove(Edge newRootEdge, double lengthMergedEdge, double lengthSplitEdge1, double lengthSplitEdge2) throws AlgorithmError {
        throw new AlgorithmError("Cannot perform root moves on unrooted trees.");
    }

    @Override
    public Edge[] reverseRootMove(Edge[][] edgeArray) throws AlgorithmError {
        throw new AlgorithmError("Cannot perform root moves on unrooted trees.");
    }

    /** Wrapper for split constructor */
    public static Split makeSplitFromFourSubsets(HashSet<String> a1, HashSet<String> a2, HashSet<String> b1, HashSet<String> b2) throws AlgorithmError {
        HashSet<String> setA = new HashSet<String>();
        HashSet<String>setB = new HashSet<String>();
        setA.addAll(a1);
        setA.addAll(a2);
        setB.addAll(b1);
        setB.addAll(b2);
        return new Split(setA,setB);
    }

}
