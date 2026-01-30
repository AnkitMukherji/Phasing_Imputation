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

import vcf.Markers;
import vcf.Marker;
import beagleutil.Samples;

/**
 * <p>Interface {@code GenotypeValues} represents a value for each
 * possible genotype for each sample at each marker.
 * </p>
 * <p>All instances of {@code GenotypeValues} are required to be thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface GenotypeValues {

    /**
     * Returns the specified genotype value.
     *
     * @param marker a marker index
     * @param sample a sample index
     * @param genotype a genotype index
     * @return the specified genotype value
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code genotype < 0 || genotype >= this.marker(marker).nGenotypes()}
     */
    float value(int marker, int sample, int genotype);

    /**
     * Adds the specified genotype values to the stored genotype values
     * for the specified sample.  This method is equivalent to
     * <pre>
     * for (m=0; m&lt;this.nMarkers(); ++m) {
     *     offset = this.markers().sumGenotypes(m);
     *     for (gt=0; gt&lt;this.marker(m).nGenotypes(); ++gt) {
     *         this.add(marker, sample, gt, values[offset + gt])
     *     }
     * }
     * </pre>
     *
     * @param sample a sample index
     * @param values an array of length {@code this.markers.sumGenotypes()}
     * containing the genotype values to be added.
     *
     * @throws IndexOutOfBoundsException if
     * if {@code sample < 0 || sample >= this.nSamples()}
     * @throws IllegalArgumentException if
     * {@code values.length != this.markers().sumGenotypes()}
     * @throws NullPointerException if {@code values == null}
     */
    void add(int sample, double[] values);

    /**
     * Adds the specified genotype value to the stored genotype value.
     *
     * @param marker a marker index
     * @param sample a sample index
     * @param genotype a genotype index
     * @param value the value to be added
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code genotype < 0 || genotype >= this.marker(marker).nGenotypes()}
     */
    void add(int marker, int sample, int genotype, double value);

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    int nMarkers();

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    Markers markers();

    /**
     * Returns the specified marker.
     * @param marker a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    Marker marker(int marker);

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

    /**
     * Returns a string representation of {@code this}. The exact details
     * of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString();
}
