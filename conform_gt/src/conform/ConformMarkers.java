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
package conform;

import beagleutil.ChromInterval;
import blbutil.Const;
import blbutil.FileIt;
import blbutil.InputIt;
import blbutil.SampleFileIt;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import vcf.IntervalVcfIt;
import vcf.Marker;
import vcf.VcfEmission;
import vcf.VcfIt;

/**
 * Class {@code ConformMarkers} stores the markers found in a VCF file
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ConformMarkers {

    private static final int SIZE = 100_000;

    private final List<Marker> markerList;
    private final Map<Marker, Integer> markerIndexMap;

    public ConformMarkers(File vcfFile, ChromInterval chromInterval) {
        this.markerList = new ArrayList<>(SIZE);
        this.markerIndexMap = new HashMap<>(SIZE);

        try (FileIt<String> it1 = InputIt.fromGzipFile(vcfFile);
                SampleFileIt<VcfEmission> it2 = VcfIt.create(it1,
                        VcfIt.toBitSetGT);
                SampleFileIt<VcfEmission> it3 = new IntervalVcfIt<>(it2,
                        chromInterval)) {
            while (it3.hasNext()) {
                Marker m = it3.next().marker();
                markerList.add(m);
                Integer prevValue = markerIndexMap.put(m, markerList.size());
                if (prevValue != null) {
                    String s = "Duplicate marker [" + vcfFile + "]: " + m;
                    throw new IllegalArgumentException(s);
                }
            }
        }
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return markerList.size();
    }

    /**
     * Returns the specified marker.
     * @param index a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nMarkers()}
     */
    public Marker marker(int index) {
        return markerList.get(index);
    }

    /**
     * Returns the index of the specified marker, or {@code -1} if the
     * marker is not present in the list of markers represented by {@code this}.
     * @param marker a marker
     * @return the index of the specified marker
     * @throws NullPointerException if {@code marker == null}
     */
    public int indexOf(Marker marker) {
        if (marker==null) {
            throw new NullPointerException(Marker.class.toString());
        }
        Integer i = markerIndexMap.get(marker);
        return i==null ? -1 : i;
    }

    /**
     * Returns a map which maps each VCF record POS field to the
     * list of markers having that position.
     * @return a map which maps each VCF record POS field to the
     * list of markers having that position
     */
    public Map<Integer, List<Marker>> posMap() {
        Map<Integer, List<Marker>> posMap = new HashMap<>(markerList.size());
        for (Marker m : markerList) {
            List<Marker> list = posMap.get(m.pos());
            if (list==null) {
                list = new ArrayList<>(1);
                posMap.put(m.pos(), list);
            }
            list.add(m);
        }
        return posMap;
    }

    /**
     * Returns a map which maps each VCF record ID field to its corresponding
     * marker.
     * @return a map which maps each VCF record ID field it its corresponding
     * marker
     */
    public Map<String, Marker> idMap() {
        Map<String, Marker> idMap = new HashMap<>(markerList.size());
        for (Marker m : markerList) {
            for (int j=0, n=m.nIds(); j<n; ++j) {
                Marker prev = idMap.put(m.id(j), m);
                if (prev!=null) {
                    String s = "Non-unique marker ID: " + m.id(j)
                            + Const.nl + prev + Const.nl + m;
                    throw new IllegalArgumentException(s);
                }
            }
        }
        return idMap;
    }
}
