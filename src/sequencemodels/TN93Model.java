/**
    TN93Model
    Implementation of the TN93 model

    Copyright (C) 2011  Tom M. W. Nye

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

package sequencemodels;

import cern.colt.matrix.tdouble.DoubleFactory2D;

/**
 * Implementation of the TN93 model
 * See Ziheng Yang, Computational Molecular Evolution, Section 1.2.3
 */

public class TN93Model implements SubstitutionModel, SequenceDistanceCalculator {
    
    /** Version number for serialization - increment when structural changes are made */
    private static final long serialVersionUID = 1L;

    /* Model parameters */
    protected double piA, piG, piC, piT, piY, piR, alpha1, alpha2, beta;
    private Alphabet theAlphabet = DNAAlphabet.getInstance();
    private double[][] Q;

    /** Constructor: fix rate matrix etc.
     p[] must be an array length 4 containing the stationary distribution in order A, G, C, T*/
    public TN93Model(double[] p, double kappa1, double kappa2) {
        // Fix params
        setParameters(p, kappa1, kappa2);
    }
    
    @Override
    public Alphabet getAlphabet() {
        return theAlphabet;
    }
    
    protected void setParameters(double[] p, double kappa1, double kappa2) {
         // Fix params etc
        alpha1=kappa1; alpha2=kappa2; beta=1.0;
        piA=p[0]; piG=p[1]; piC=p[2]; piT=p[3];
        piR=piA+piG; piY=piC+piT;

        // Fix Q
        Q = new double[4][4];
        Q[0][0] =  -(alpha2*piG+beta*piY); //AA
        Q[0][1] = alpha2*piG; //AG
        Q[0][2] = beta*piC; //AC
        Q[0][3] = beta*piT; //AT
        Q[1][0] = alpha2*piA; //GA
        Q[1][1] = -(alpha2*piA+beta*piY); //GG
        Q[1][2] = beta*piC; //GC
        Q[1][3] = beta*piT; //GT
        Q[2][0] = beta*piA; //CA
        Q[2][1] = beta*piG; //CG
        Q[2][2] = -(alpha1*piT+beta*piR); //CC
        Q[2][3] = alpha1*piT; //CT
        Q[3][0] = beta*piA; //TA
        Q[3][1] = beta*piG; //TG
        Q[3][2] = alpha1*piC; //TC
        Q[3][3] = -(alpha1*piC+beta*piR); //TT


        // Normalize to overall rate 1.0
        double fac = StaticLikelihoodCalculator.normalizeQ(Q, getStationaryDistrib());
        beta = beta / fac;
        alpha1 = alpha1 / fac;
        alpha2 = alpha2 / fac;

    }
    
    

    
    /* For likelihood calculations ------------------------------------------ */
    
    @Override
    public double[][] getTransitionMatrix(double t) {

        double e2 = Math.exp(-beta*t);
        double e3 = Math.exp(-(piR*alpha2+piY*beta)*t);
        double e4 = Math.exp(-(piY*alpha1+piR*beta)*t);

        double yr = piY/piR;
        double ry = 1.0/yr;

        double[][] P = new double[4][4];
        P[0][0] = Math.min(Math.max(0.0,piA+piA*yr*e2+piG/piR*e3),1.0); //AA
        P[0][1] = Math.min(Math.max(0.0,piG+piG*yr*e2-piG/piR*e3),1.0); //AG
        P[0][2] = Math.min(Math.max(0.0,piC*(1-e2)),1.0); //AC
        P[0][3] = Math.min(Math.max(0.0,piT*(1-e2)),1.0); //AT
        P[1][0] = Math.min(Math.max(0.0,piA+piA*yr*e2-piA/piR*e3),1.0); //GA
        P[1][1] = Math.min(Math.max(0.0,piG+piG*yr*e2+piA/piR*e3),1.0); //GG
        P[1][2] = Math.min(Math.max(0.0,piC*(1-e2)),1.0); //GC
        P[1][3] = Math.min(Math.max(0.0,piT*(1-e2)),1.0); //GT
        P[2][0] = Math.min(Math.max(0.0,piA*(1-e2)),1.0); //CA
        P[2][1] = Math.min(Math.max(0.0,piG*(1-e2)),1.0); //CG
        P[2][2] = Math.min(Math.max(0.0,piC+piC*ry*e2+piT/piY*e4),1.0); //CC
        P[2][3] = Math.min(Math.max(0.0,piT+piT*ry*e2-piT/piY*e4),1.0); //CT
        P[3][0] = Math.min(Math.max(0.0,piA*(1-e2)),1.0); //TA
        P[3][1] = Math.min(Math.max(0.0,piG*(1-e2)),1.0); //TG
        P[3][2] = Math.min(Math.max(0.0,piC+piC*ry*e2-piC/piY*e4),1.0); //TC
        P[3][3] = Math.min(Math.max(0.0,piT+piT*ry*e2+piC/piY*e4),1.0); //TT

        return P;
    }
    @Override
    public void getTransitionMatrix(double t, double[][] P) {
        double e2 = Math.exp(-beta*t);
        double e3 = Math.exp(-(piR*alpha2+piY*beta)*t);
        double e4 = Math.exp(-(piY*alpha1+piR*beta)*t);

        double yr = piY/piR;
        double ry = 1.0/yr;

        P[0][0] = Math.min(Math.max(0.0,piA+piA*yr*e2+piG/piR*e3),1.0); //AA
        P[0][1] = Math.min(Math.max(0.0,piG+piG*yr*e2-piG/piR*e3),1.0); //AG
        P[0][2] = Math.min(Math.max(0.0,piC*(1-e2)),1.0); //AC
        P[0][3] = Math.min(Math.max(0.0,piT*(1-e2)),1.0); //AT
        P[1][0] = Math.min(Math.max(0.0,piA+piA*yr*e2-piA/piR*e3),1.0); //GA
        P[1][1] = Math.min(Math.max(0.0,piG+piG*yr*e2+piA/piR*e3),1.0); //GG
        P[1][2] = Math.min(Math.max(0.0,piC*(1-e2)),1.0); //GC
        P[1][3] = Math.min(Math.max(0.0,piT*(1-e2)),1.0); //GT
        P[2][0] = Math.min(Math.max(0.0,piA*(1-e2)),1.0); //CA
        P[2][1] = Math.min(Math.max(0.0,piG*(1-e2)),1.0); //CG
        P[2][2] = Math.min(Math.max(0.0,piC+piC*ry*e2+piT/piY*e4),1.0); //CC
        P[2][3] = Math.min(Math.max(0.0,piT+piT*ry*e2-piT/piY*e4),1.0); //CT
        P[3][0] = Math.min(Math.max(0.0,piA*(1-e2)),1.0); //TA
        P[3][1] = Math.min(Math.max(0.0,piG*(1-e2)),1.0); //TG
        P[3][2] = Math.min(Math.max(0.0,piC+piC*ry*e2-piC/piY*e4),1.0); //TC
        P[3][3] = Math.min(Math.max(0.0,piT+piT*ry*e2+piC/piY*e4),1.0); //TT

    }
    @Override
    public double[] getStationaryDistrib() {
        double[] x = new double[]{piA, piG, piC, piT};
        return x;
    }
    @Override
    public double stationaryProb(int i) {
        if(i==0) return piA;
        else if(i==1) return piG;
        else if(i==2) return piC;
        else if(i==3) return piT;
        else throw new java.lang.ArrayIndexOutOfBoundsException("Index for stationary prob of TN93 model must lie between 0 and 3.");
    }
    /*@Override
    public double getQEntry(int i, int j) {
        return Q[i][j];
    }
    @Override
    public double getRhoEntry(int i, int j) throws treebase.AlgorithmException {
        if(i==j) {
            throw new treebase.AlgorithmException("AlgorithmExcepion: diagonal elements in rho matrix not defined.");
        }
        else if((i+j)==1) return alpha2;
        else if((i+j)==5) return alpha1;
        else return beta;
    }*/
    
    /* Utilities ------------------------------------------------------------ */

    public void changeParameters(double[] p, double kappa1, double kappa2) {
        setParameters(p, kappa1, kappa2);
    }

    /* Distance calcs ------------------------------------------------------- */
    
    @Override
    public double calcDistance(int[] a, int[] b, double infiniteDist) throws treebase.AlgorithmError {
        if (a.length!=b.length) throw new treebase.AlgorithmError("Sequence length mismatch in TN93 distance calculation");
        // Count up where a and b match
        double[] res = countTransitionsTransversions(a, b);
        double v = res[0];
        double s1 = res[1];
        double s2 = res[2];

        double x;

        x = 1.0-0.5*piY*s1/(piT*piC)-0.5*v/piY;
        if (x<=0.0) return infiniteDist;
        double a1 = -Math.log(x);

        x = 1.0-0.5*piR*s2/(piA*piG)-0.5*v/piR;
        if (x<=0.0) return infiniteDist;
        double a2 = -Math.log(x);

        x = 1.0-0.5*v/(piY*piR);
        if (x<=0.0) return infiniteDist;
        double bb = -Math.log(x);

        double d = 2*piT*piC/piY*(a1-piR*bb);
        d += 2*piA*piG/piR*(a2-piY*bb);
        d += 2*piY*piR*bb;

        return d;
    }

    protected static double[] countTransitionsTransversions(int[] a, int[] b) {
        double[] res = new double[3]; // This will contain S_1, S_2, V
        double s1 = 0.0, s2 = 0.0, v = 0.0;

        for (int i=0; i<a.length; i++) {

            if (a[i]<0) {
                if (b[i]<0) {
                    v += 0.5;
                    s1 += 0.125;
                    s2 += 0.125;
                }
                if ((b[i]==0)||(b[i]==1)) {
                    v += 0.5;
                    s2 += 0.25;
                }
                if ((b[i]==2)||(b[i]==3)) {
                    v += 0.5;
                    s1 += 0.25;
                }
            }
            else if (b[i]<0) {
                if ((a[i]==0)||(a[i]==1)) {
                    v += 0.5;
                    s2 += 0.25;
                }
                if ((a[i]==2)||(a[i]==3)) {
                    v += 0.5;
                    s1 += 0.25;
                }
            }
            // OK, a and b both +ve
            else if (((a[i]<=1)&&(b[i]>=2))||((b[i]<=1)&&(a[i]>=2))) {
                // Transversion
                v += 1.0;
            }
            else if ((a[i]==2)&&(b[i]==3)||(a[i]==3)&&(b[i]==2)) {
                // CT or TC
                s1 += 1.0;
            }
            else if ((a[i]==0)&&(b[i]==1)||(a[i]==1)&&(b[i]==0)) {
                // Ag or GA
                s2 += 1.0;
            }
        } // End loop thru sites

        res[0] = v/a.length;
        res[1] = s1/a.length;
        res[2] = s2/a.length;

        return res;
    }

    /** Distance calculation from alignment */
    @Override
    public double[][] calcDistanceMatrix(Alignment a, double infiniteDist) {
        double[][] d = new double[a.getNumTaxa()][a.getNumTaxa()];
        double x;

        try {
            for (int i=0; i<a.getNumTaxa(); i++) {
                d[i][i] = 0.0;
                for (int j=0; j<i; j++) {
                    x = calcDistance(a.getRow(i), a.getRow(j),infiniteDist);
                    d[i][j] = x;
                    d[j][i] = x;
                }
            }
        }
        catch (treebase.AlgorithmError e) {
            javax.swing.JOptionPane.showMessageDialog(null,"Error with alignment -- difference numbers of characters in rows. ","Error",javax.swing.JOptionPane.ERROR_MESSAGE);
        }

        return d;
    }

    
    /* Test area -----------------------------------------------------------  */
    
    public String getQStr() {
        return DoubleFactory2D.dense.make(Q).toString();
    }
    


}
