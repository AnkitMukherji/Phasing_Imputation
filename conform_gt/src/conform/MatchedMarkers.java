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

import beagleutil.Phase;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import vcf.BasicMarker;
import vcf.Marker;

/**
 * Class {@code ConformGtMarkers} finds the sublist of markers in
 * a reference VCF file and a target VCF file that have the same identifier
 * or position in both files and that have consistent alleles
 * when allowing for a strand switch in the target file marker. When mapping
 * target markers to reference markers, a match is first attempted using
 * marker identifiers. If no match is found with marker identifiers, a
 * match is then attempted using marker position.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class MatchedMarkers {

    private final List<Marker> refMarkers;
    private final List<Marker> targMarkers;
    private final List<Phase> strand;

    /**
     * Constructs a {@code ConformGtMarkers} instance
     * Private constructor prevents public instantiation
     * @param ref the list of reference markers
     * @param targ the list of target markers
     * @param matchId {@code true} if reference and target markers are to
     * be matched by identifier and {@code false} if reference and target
     * markers are to be matched by position
     * @param log a {@code PrintWriter} to which log message will be written
     *
     * @throws NullPointerException if
     * {@code ref==null || targ==null || log==null}
     */
    public MatchedMarkers(ConformMarkers ref, ConformMarkers targ,
            boolean matchId, PrintWriter log) {
        this.refMarkers = new ArrayList<>(targ.nMarkers());
        this.targMarkers = new ArrayList<>(targ.nMarkers());
        this.strand = new ArrayList<>(targ.nMarkers());
        Matcher matcher = matchId ? idMatch(ref.idMap()) : posMatch(ref.posMap());
        match(ref, targ, matcher, log);
    }

    /**
     * Returns the number of markers that were matched in the reference
     * and target data.
     * @return the number of markers that were matched in the reference
     * and target data
     */
    public int nMarkers() {
        return refMarkers.size();
    }

    /**
     * Returns the specified reference marker in the list of matched markers.
     * @param index index of a marker in the list of matched markers
     * @return the specified reference marker in the list of matched markers
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nMarkers()}
     */
    public Marker refMarker(int index) {
        return refMarkers.get(index);
    }

    /**
     * Returns the specified target marker in the list of matched markers.
     * @param index index of a marker in the list of matched markers
     * @return the specified target marker in the list of matched markers
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nMarkers()}
     */
    public Marker targetMarker(int index) {
        return targMarkers.get(index);
    }

    /**
     * Returns the relationship of the chromosome strands for the reference
     * allele of the specified marker in the list of matched markers.
     * @param index index of a marker in the list of matched markers
     * @return  the relationship of the chromosome strands for the reference
     * allele of the specified marker in the list of matched markers
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nMarkers()}.
     */
    public Phase strand(int index) {
        return strand.get(index);
    }

    private void match(ConformMarkers ref, ConformMarkers targ, Matcher matcher,
            PrintWriter log) {
        int prevRefIndex = -1;
        for (int j=0; j<targ.nMarkers(); ++j) {
            Marker targMarker = targ.marker(j);
            List<Marker> posMatches = matcher.match(targMarker);
            List<Marker> consistent = consistentMarkers(posMatches, targMarker);
            if (consistent.isEmpty()) {
                ConformMain.printUnmatchedMarkerLog(targMarker,
                        "NOT_IN_REFERENCE", log);
            }
            else if (consistent.size()>1) {
                ConformMain.printUnmatchedMarkerLog(targMarker,
                        "MULTIPLE_REF_MATCHES", log);
            }
            else {
                Marker refMarker = consistent.get(0);
                int refIndex = ref.indexOf(refMarker);
                if (refIndex == prevRefIndex) {
                    ConformMain.printUnmatchedMarkerLog(targMarker,
                            "DUPLICATE_MARKER", log);
                }
                else if (refIndex < prevRefIndex) {
                    ConformMain.printUnmatchedMarkerLog(targMarker,
                            "MARKER_OUT_OF_ORDER", log);
                }
                else {
                    refMarkers.add(refMarker);
                    targMarkers.add(targMarker);
                    strand.add(strand(refMarker, targMarker));
                    if (refIndex > prevRefIndex) {
                        prevRefIndex = refIndex;
                    }
                }
            }
        }
    }

    private static List<Marker> consistentMarkers(List<Marker> refMarkers,
            Marker targMarker) {
        List<Marker> consistent = new ArrayList<>(1);
        for (Marker refMarker : refMarkers) {
            if (strand(refMarker, targMarker)!=Phase.INCONSISTENT) {
                consistent.add(refMarker);
            }
        }
        return consistent;
    }

    private static Phase strand(Marker refMarker, Marker targetMarker) {
        Set<String> refAlleles = toAlleleSet(refMarker);
        Set<String> targetAlleles = toAlleleSet(targetMarker);
        Set<String> flippedTargetAlleles =
                toAlleleSet(BasicMarker.flipStrand(targetMarker));
        boolean sameConsistent = refAlleles.containsAll(targetAlleles);
        boolean oppConsistent = refAlleles.containsAll(flippedTargetAlleles);
        if (sameConsistent && oppConsistent) {
           return Phase.UNKNOWN;
        }
        else if (sameConsistent && oppConsistent==false) {
            return Phase.IDENTICAL;
        }
        else if (sameConsistent==false && oppConsistent) {
            return Phase.OPPOSITE;
        }
        else {
            return Phase.INCONSISTENT;
        }
    }

    private static Set<String> toAlleleSet(Marker marker) {
        Set<String> alleles = new HashSet<>(3);
        for (int j=0, n=marker.nAlleles(); j<n; ++j) {
            alleles.add(marker.allele(j));
        }
        return alleles;
    }

    private static interface Matcher {
        public List<Marker> match(Marker marker);
    }

    private Matcher idMatch(final Map<String, Marker> idMap) {
        return new Matcher() {
            @Override
            public List<Marker> match(Marker marker) {
                List<Marker> matches = new ArrayList<>(1);
                for (int j=0; j<marker.nIds(); ++j) {
                    Marker match = idMap.get(marker.id(j));
                    if (match!=null) {
                        matches.add(match);
                    }
                }
                return matches;
            }
        };
    }

    private Matcher posMatch(final Map<Integer, List<Marker>> posMap) {
        return new Matcher() {
            @Override
            public List<Marker> match(Marker marker) {
                List<Marker> matches = posMap.get(marker.pos());
                return matches==null ? Collections.emptyList() : matches;
            }
        };
    }
}

