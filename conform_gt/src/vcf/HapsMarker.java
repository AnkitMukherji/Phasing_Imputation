/*
 * Copyright 2013-2016 Brian L. Browning
 *
 * This file is part of the conform-gt program.
 *
 * Licensed under the Apache License, Version 2.0 (the License);
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package vcf;

/**
 * <p>Interface {@code HapsMarker} represents marker alleles for a 
 * list of haplotype pairs.
 * </p>
 * All instances of {@code HapsMarkers} are required to be
 * immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface HapsMarker extends MarkerContainer {

     /**
     * Returns the allele on the specified haplotype.
     * @param haplotype a haplotype index
     * @return the allele on the specified haplotype
     *
     * @throws IndexOutOfBoundsException if
     * {@code haplotype < 0 || haplotype >= this.nHaps()}
     */
    int allele(int haplotype);

    /**
     * Returns the first allele for the specified haplotype pair.
     * @param hapPair a haplotype pair index
     * @return the first allele for the specified haplotype pair
     *
     * @throws IndexOutOfBoundsException if
     * {@code hapPair < 0 || hapPair >= this.nHapPairs()}
     */
    int allele1(int hapPair);

    /**
     * Returns the second allele for the specified haplotype pair.
     * @param hapPair a haplotype pair index
     * @return the second allele for the specified haplotype pair
     *
     * @throws IndexOutOfBoundsException if
     * {@code hapPair < 0 || hapPair >= this.nHapPairs()}
     */
    int allele2(int hapPair);

    /**
     * Returns the marker.
     * @return the marker
     */
    @Override
    Marker marker();

    /**
     * Returns the number of haplotypes.  The returned value is equal to
     * {@code 2*this.nHapPairs()}.
     * @return the number of haplotypes
     */
    int nHaps();

    /**
     * Returns the number of haplotype pairs.  The returned value is
     * equal to {@code this.nHaps()/2}.
     * @return the number of haplotype pairs
     */
    int nHapPairs();
}
