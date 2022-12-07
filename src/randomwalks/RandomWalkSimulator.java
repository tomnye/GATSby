/*
 * RandomWalkSimulator.java

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

package randomwalks;

/**
 * Set of static methods for performing simulations of random walks
      NB 1: tree must be fully resolved!
      NB 2: much faster to diffuse a true Tree object rather than TreeAsSplits
 */

import simulation.NormalDistribution;
import simulation.GammaDistribution;
import cern.jet.random.tdouble.DoubleUniform;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import treebase.AlgorithmError;
import treebase.AlgorithmException;
import treebase.Graph;
import treebase.Split;
import treebase.Tree;
import treebase.TreeAsSplits;
import treebase.TreeWithTopologicalOperations;
import treedatasets.OperationsOnTreeAsSplits;
//import cern.jet.random.tdouble.engine.DoubleMersenneTwister;

public class RandomWalkSimulator {

    /* Methods based on maps of splits -------------------------------------- */

    /** Diffuse a tree (represented as TreeAsSplits). Internal edges only.
      Do numSteps replacements FOR EACH EDGE.
      NB 1: tree must be fully resolved!
      NB 2: much faster to diffuse a true Tree object rather than TreeAsSplits
     */
    public static void randomWalkTree(TreeAsSplits theTree, double var, int numSteps, boolean randomSweep) {
        NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(var/numSteps));
        DoubleUniform unif = new DoubleUniform(simulation.Random.getEngine());
        RandomWalkSimulator.randomWalkTree(theTree, numSteps, randomSweep, norm, unif);
    }
    public static void randomWalkTree(TreeAsSplits theTree, int numSteps, boolean randomSweep, NormalDistribution norm, DoubleUniform unif) {

        /** If you use this method make sure the normal distribution has variance = total var / numSteps */

        /* First get the non-trivial splits in a fixed order */
        ArrayList<Split> splits = new ArrayList();
        splits.addAll(theTree.getNonTrivialSplits());
        Collections.sort(splits);
        if (splits.size()<theTree.getNumTaxa()-3) {
            System.out.println("Warning: request made to random walk an unresolved tree. Trees must be randomly resolved first.");
        }

        /* Main loop */
        Split q;
        double x, l;
        int ind, topInd;
        for (int i=0; i<numSteps; i++) {
            for (int j=0; j<splits.size(); j++) {

                if (randomSweep)
                    ind = unif.nextIntFromTo(0, splits.size()-1);
                else
                    ind = j;

                x = norm.sample();
                Split p = splits.get(ind);
                l = x+theTree.getSplitLength(p);

                try {
                    if (l>0.0) {
                        theTree.setSplitLength(p, l);
                    }
                    else {
                        topInd = unif.nextIntFromTo(1, 3);
                        if (topInd==1) {
                            theTree.setSplitLength(p, -l);
                        }
                        else {
                            Split[] nniSplit = OperationsOnTreeAsSplits.getNNISplits(theTree, p);
                            q = nniSplit[topInd-2]; // topInd-2 = 0 or 1
                            theTree.remove(p);
                            theTree.add(q, -l);
                            splits.set(ind, q);
                        }
                    }
                }
                catch (AlgorithmException err) {
                    System.out.println("Error setting edge length in a diffusion.");
                }
            } // End loop thru' splits
        } // End main loop

    }

    /** Generate one realization from a single tree.
     numSteps steps PER EDGE each with variance var/numSteps */
    public static TreeAsSplits sampleRandomWalk(TreeAsSplits theTree, double var, int numSteps, boolean randomSweep) {
        NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(var/numSteps));
        DoubleUniform unif = new DoubleUniform(simulation.Random.getEngine());
        return RandomWalkSimulator.sampleRandomWalk(theTree, numSteps, randomSweep, norm, unif);
    }
    public static TreeAsSplits sampleRandomWalk(TreeAsSplits theTree, int numSteps, boolean randomSweep, NormalDistribution norm, DoubleUniform unif) {
        /** If you use this method make sure the normal distribution has variance = total var / numSteps */
        TreeAsSplits t = theTree.clone();
        RandomWalkSimulator.randomWalkTree(t, numSteps, randomSweep, norm, unif);
        return t;
    }

    /** Generate many realizations from a single tree */
    public static TreeAsSplits[] sampleRandomWalk(TreeAsSplits theTree, double var, int numSteps, int numSamples, boolean randomSweep) {
        TreeAsSplits[] res = new TreeAsSplits[numSamples];
        NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(var/numSteps));
        DoubleUniform unif = new DoubleUniform(simulation.Random.getEngine());
        for (int i=0; i<numSamples; i++) {
            TreeAsSplits t = theTree.clone();
            RandomWalkSimulator.randomWalkTree(t, numSteps, randomSweep, norm, unif);
            res[i] = t;
        }
        return res;
    }

    /* Methods based on trees ----------------------------------------------- */

    /** Randomly walk a tree. Internal edges only.
      Do numSteps replacements ON EACH EDGE with variance var/numSteps.
      NB: tree must be fully resolved!  */
    public static void randomWalkTree(TreeWithTopologicalOperations theTree, double var, int numSteps, boolean randomSweep) {
        randomWalkTree(theTree, numSteps, randomSweep, new NormalDistribution(0.0, Math.sqrt(var/numSteps)), new DoubleUniform(simulation.Random.getEngine()));
    }
    public static void randomWalkTree(TreeWithTopologicalOperations theTree, int numSteps, boolean randomSweep, NormalDistribution norm, DoubleUniform unif) {

        /* First get the internal edges in a fixed order */
        ArrayList<Graph.Edge> edges = new ArrayList();
        Iterator<Graph.Edge> itE = theTree.getEdgeIterator();
        for (int i=0; i<theTree.numEdges(); i++) {
            Graph.Edge e = itE.next();
            if (!Tree.isEdgeTerminal(e)) edges.add(e);
        }

        if (edges.size()<theTree.numTaxa()-3) {
            System.out.println("Warning: request made to diffuse an unresolved tree. Trees must be randomly resolved first.");
        }

        /* Main loop */

        Graph.Edge f;
        double x, l;
        int ind, topInd;
        for (int i=0; i<numSteps; i++) {
            for (int j=0; j<edges.size(); j++) {

                if (randomSweep)
                    ind = unif.nextIntFromTo(0, edges.size()-1);
                else
                    ind = j;

                x = norm.sample();
                Graph.Edge e = edges.get(ind);
                l = x+e.getLength();

                try {
                    if (l>0.0) {
                        e.setLength(l);
                    }
                    else {
                        
                        topInd = unif.nextIntFromTo(1, 3);
                        if (topInd==1) {
                            e.setLength(-l);
                        }
                        else {
                            boolean choice = (topInd==2); // choice=true if topInd==2, choice=false if topInd==3
                            f = theTree.performNNI(e, choice, -l);
                            edges.set(ind, f); // Replace e in the list of edges with edge f
                        }

                    }
                }
                catch (AlgorithmException err) {
                    System.out.println("Error setting edge length in a diffusion."+err.getMessage());
                }
            } // End loop thru' splits
        } // End main loop

    }

    /** Generate one realization from a single tree */
    public static Tree sampleRandomWalk(Tree theTree, double var, int numSteps, boolean randomSweep) {
        NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(var/numSteps));
        DoubleUniform unif = new DoubleUniform(simulation.Random.getEngine());
        return RandomWalkSimulator.sampleRandomWalk(theTree, numSteps, randomSweep, norm, unif);
    }
    public static Tree sampleRandomWalk(Tree theTree, int numSteps, boolean randomSweep, NormalDistribution norm, DoubleUniform unif) {
        TreeWithTopologicalOperations t = new TreeWithTopologicalOperations(theTree);
        randomWalkTree(t, numSteps, randomSweep, norm, unif);
        return t;
    }
    /** Generate many realizations from a single tree */
    public static Tree[] sampleRandomWalk(Tree theTree, double var, int numSteps, boolean randomSweep, int numSamples) {
        NormalDistribution norm = new NormalDistribution(0.0, Math.sqrt(var/numSteps));
        DoubleUniform unif = new DoubleUniform(simulation.Random.getEngine());
        Tree[] res = new Tree[numSamples];
        for (int i=0; i<numSamples; i++) {
            TreeWithTopologicalOperations t = new TreeWithTopologicalOperations(theTree);
            randomWalkTree(t, numSteps, randomSweep, norm, unif);
            res[i] = t;
        }
        return res;
    }

    /* Parallel method to sample large numbers of trees */
    public static HashSet<Tree> sampleRandomWalkParallel(Tree theTree, double var, int numSteps, boolean randomSweep, int numSamples, int numProcessors)  {
        int chunk = Math.round(((float)numSamples)/((float)numProcessors));
        int start=0, end=chunk, seed;
        HashSet<Tree> res = new HashSet();

        ParallelRWChunk[] processes = new ParallelRWChunk[numProcessors];
        for (int i=0; i<numProcessors; i++) {
            // Create processes
            seed = simulation.Random.getEngine().nextInt();
            if (i==(numProcessors-1)) end=numSamples;
            ParallelRWChunk p = new ParallelRWChunk(theTree, numSamples, var, numSteps, randomSweep, seed);
            p.start();
            processes[i] = p;
            start = end;
            end = start+chunk;
        }
        // Wait
        try {
             for (int i=0; i<numProcessors; i++) {
                processes[i].join();
                res.addAll(processes[i].getTrees());
            }
        }
        catch (InterruptedException ex) {
           System.out.println("Error: something funny happened with a thread interrupt while projecting trees.");
        }
        return res;

    }

    public static class ParallelRWChunk extends Thread {

        private Tree theTree;
        private int numSamples, numSteps;
        private HashSet<Tree> res;
        private double sigma;
        private NormalDistribution norm;
        private DoubleUniform unif;
        private boolean randomSweep;


        public ParallelRWChunk(Tree t, int nSamples, double v, int nSteps, boolean rs, int seed) {
            theTree = t;
            numSamples = nSamples;
            sigma = Math.sqrt(v/nSteps);
            numSteps = nSteps;
            res = new HashSet();
            DoubleMersenneTwister engine = new DoubleMersenneTwister(seed);
            norm = new NormalDistribution(0.0, sigma, engine);
            unif = new DoubleUniform(engine);
            randomSweep = rs;
        }

       public void run() {
            for (int i=0; i<numSamples; i++) {
                TreeWithTopologicalOperations t = new TreeWithTopologicalOperations(theTree);
                randomWalkTree(t, numSteps, randomSweep, norm, unif);
                res.add(t);
            }
       }
       
       public HashSet<Tree> getTrees() {
           return res;
       }
    }


    /** Use a gamma distribution to perturb pedant edges. */
    public static void perturbPendantEdges(TreeAsSplits theTree, double var) {
        GammaDistribution gamma = new GammaDistribution(1.0, 1.0); // Arguments are shape and scale
        ArrayList<Split> splits = new ArrayList();
        splits.addAll(theTree.getSplits());
        splits.removeAll(theTree.getNonTrivialSplits());
        Collections.sort(splits);

        Iterator<Split> it = splits.iterator();
        double mu, x;
        Split p;
        try {
            while (it.hasNext()) {
                p = it.next();
                mu = theTree.getSplitLength(p);
                if (Math.abs(mu)>1.0E-10) {
                    x = gamma.sample(mu*mu/var, var/mu); // mean = shape*scale = mu; var = shape*scale^2 = sigma^2
                    theTree.setSplitLength(p, x);
                }
            }
        }
        catch (AlgorithmError anErr) {
            System.out.println("Error perturbing pendant edge lengths for a diffusion. This should not be possible. ");
        }
    }

}

