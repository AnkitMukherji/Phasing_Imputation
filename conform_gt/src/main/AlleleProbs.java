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
package main;

import beagleutil.Samples;
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Interface {@code AlleleProbs} represents per-haplotype allele
 * probabilities for a list of samples.
 * </p>
 * <p>All instances of {@code AlleleProbs} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface AlleleProbs {

    /**
     * Returns the probability that the specified marker allele is
     * present on the first haplotype of the specified sample.
     *
     * @param marker a marker index
     * @param sample a sample index
     * @param allele an allele index
     * @return the probability that the specified marker allele is
     * present on the first haplotype of the specified sample
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.marker(marker).nAlleles()}
     */
    float alProb1(int marker, int sample, int allele);

    /**
     * Returns the probability that the specified marker allele is
     * present on the second haplotype of the specified sample.
     *
     * @param marker a marker index
     * @param sample a sample index
     * @param allele an allele index
     * @return the probability that the specified marker allele is
     * present on the second haplotype of the specified sample.
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.marker(marker).nAlleles()}
     */
    float alProb2(int marker, int sample, int allele);

    /**
     * Returns the phased genotype probability, equal to
     * {@code (this.allele1(marker, sample, allele1)
     *        * this.allele2(marker, sample, allele2))}.
     *
     * @param marker a marker index
     * @param sample a sample index
     * @param allele1 allele index of the allele on the first haplotype
     * @param allele2 allele index of the allele on the second haplotype
     * @return the phased genotype probability equal to
     *       {@code (this.allele1(marker, sample, allele1)
     *              * this.allele2(marker, sample, allele2))}
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code allele1 < 0 || allele1 >= this.marker(marker).nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code allele2 < 0 || allele2 >= this.marker(marker).nAlleles()}
     */
    float gtProb(int marker, int sample, int allele1, int allele2);

    /**
     * Returns the marker allele with maximum probability for the
     * first haplotype of the specified sample. If more than one allele
     * has maximum probability, one of the alleles with maximum
     * probability will be returned.
     * @param marker a marker index
     * @param sample a sample index
     * @return the marker allele with maximum probability for the
     * first haplotype of the specified sample
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    int allele1(int marker, int sample);

    /**
     * Returns the marker allele with maximum probability for the
     * second haplotype of the specified sample. If more than one allele
     * has maximum probability, one of the alleles with maximum
     * probability will be returned.
     * @param marker a marker index
     * @param sample a sample index
     * @return the marker allele with maximum probability for the
     * second haplotype of the specified sample
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    int allele2(int marker, int sample);

    /**
     * Returns the specified marker.
     * @param marker a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    Marker marker(int marker);

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    Markers markers();

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    int nMarkers();

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    int nSamples();

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    Samples samples();
}
