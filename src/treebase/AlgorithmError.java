package treebase;

/*
 * AlgorithmError.java

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

public class AlgorithmError extends AlgorithmException {
    
    /**
     * Creates a new instance of <code>AlgorithmError</code> without detail message.
     */
    public AlgorithmError() {
        super("An unspecified AlgorithmError has been thrown...");
    }
    
    
    /**
     * Constructs an instance of <code>AlgorithmError</code> with the specified detail message.
     * @param msg the detail message.
     */
    public AlgorithmError(String msg) {
        super(msg);
    }
}
