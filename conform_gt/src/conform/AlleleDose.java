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

import blbutil.Const;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import vcf.BasicMarker;
import vcf.Marker;
import vcf.VcfEmission;

/**
 * Class {@code AlleleDose} represents reference and target data
 * for a marker.  Class {@code AlleleDose} is immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class AlleleDose {

    private final Marker refMarker;
    private final int refAllele;
    private final VcfEmission filtTarg;
    private final VcfEmission unfiltTarg;

    private final int[] alleleMap;        // ve.allele -> marker.allele

    private final int[] dose;
    private final int[] alleleCnts;

    /**
     * Constructs a {@code VcfRecordPair} instance for two VCF records
     * corresponding to the same marker.
     * @param refMarker the reference marker
     * @param refAllele the reference allele
     * @param filtTarg the target genotype data (after sample filtering)
     * @param unfiltTarg the target genotype data (before sample filtering)
     * @param flipTargAlleles {@code true} if the alleles in the targetRecord
     * should be flipped to the alternate strand
     * @throws IllegalArgumentException if
     * {@code filtTarg.marker().equals(unfiltTarg.marker() == false}
     * @throws IllegalArgumentException if the alleles in the target VCF
     * record marker (after flipping alleles if {@code flipAlleles == true})
     * are not a subset of the alleles in the reference marker
     * @throws IndexOutOfBoundsException if
     * {@code refAllele < 0 || refAllele >= refMarker.nAlleles()}
     * @throws NullPointerException if
     * {@code refMarker == null || filtTargarg == null || unfiltTarg == null}
     */
    public AlleleDose(Marker refMarker, int refAllele, VcfEmission filtTarg,
            VcfEmission unfiltTarg, boolean flipTargAlleles) {
        if (filtTarg.marker().equals(unfiltTarg.marker())==false) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (refAllele<0 || refAllele >= refMarker.nAlleles()) {
            throw new IndexOutOfBoundsException(String.valueOf(refAllele));
        }
        Marker targMarker = filtTarg.marker();
        if (flipTargAlleles) {
            targMarker = BasicMarker.flipStrand(targMarker);
        }
        int[] alMap = alleleMap(targMarker, refMarker);
        if (alMap==null) {
            throw new IllegalArgumentException("inconsistent data");
        }
        int targAllele = indexOf(alMap, refAllele);

        this.refMarker = refMarker;
        this.refAllele = refAllele;
        this.filtTarg = filtTarg;
        this.unfiltTarg = unfiltTarg;
        this.alleleMap = alMap;
        this.dose = dose(filtTarg, targAllele);
        this.alleleCnts = mappedAlleleCnts(filtTarg, alMap);
    }

    private static int indexOf(int[] ia, int value) {
        for (int j=0; j<ia.length; ++j) {
            if (ia[j]==value) {
                return j;
            }
        }
        return -1;
    }

    /* returns null if no allele map from domain alleles to range alleles exists */
    private static int[] alleleMap(Marker domain, Marker range) {
        Map<String, Integer> rangeMap = alleleMap(range);
        int[] alleleMap = new int[domain.nAlleles()];
        for (int j=0; j<alleleMap.length; ++j) {
            Integer targIndex = rangeMap.get(domain.allele(j));
            if (targIndex==null) {
                return null;
            }
            else {
                alleleMap[j] = targIndex;
            }
        }
        return alleleMap;
    }

    private static Map<String, Integer> alleleMap(Marker marker) {
        int nAlleles = marker.nAlleles();
        Map<String, Integer> map = new HashMap<>(nAlleles);
        for (int j=0; j<nAlleles; ++j) {
            map.put(marker.allele(j), j);
        }
        return map;
    }

    private static int[] dose(VcfEmission rec, int allele) {
        int[] dose = new int[rec.nSamples()];
        for (int j=0; j<dose.length; ++j) {
            int a1 = rec.allele1(j);
            int a2 = rec.allele2(j);
            if (a1<0 || a2<0) {
                dose[j] = -1;
            }
            else {
                if (a1==allele) {
                    ++dose[j];
                }
                if (a2==allele) {
                    ++dose[j];
                }
            }
        }
        return dose;
    }


    private static int[] mappedAlleleCnts(VcfEmission ve, int[] alMap) {
        int[] cnts = new int[max(alMap) + 1];
        for (int j=0; j<ve.nHaps(); ++j) {
            int a = ve.allele(j);
            if (a >= 0) {
                ++cnts[alMap[a]];
            }
        }
        return cnts;
    }

    private static int max(int[] ia) {
        int max = Integer.MIN_VALUE;
        for (int i : ia) {
            if (i > max) {
                max = i;
            }
        }
        return max;
    }

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    public int nSamples() {
        return filtTarg.nSamples();
    }

    /**
     * Returns the marker.
     * @return the marker
     */
    public Marker marker() {
        return refMarker;
    }

    /**
     * Returns the reference allele which defines the dose.
     * @return the reference allele which defines the dose
     */
    public int refAllele() {
        return refAllele;
    }

    /**
     * Returns the number of non-missing alleles.
     * @return the number of non-missing alleles
     */
    public int nNonmissingAlleles() {
        return Arrays.stream(alleleCnts).sum();
    }

    /**
     * Returns the number of copies of the specified allele.
     * @param allele an allele index
     * @return the number of copies of the specified allele
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.marker().nAlleles()}
     */
    public int count(int allele) {
        return alleleCnts[allele];
    }

    /**
     * Returns the dose of the specified sample
     * @param sample a sample index
     * @return the dose of the specified sample
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || samples >= dose.length}
     */
    public int dose(int sample) {
        return dose[sample];
    }

    /**
     * Returns the correlation of the specified allele doses. Returns
     * 0.0 if the correlation is undefined.
     * @param x allele doses for a marker
     * @param y allele doses for a marker
     * @return the correlation of the specified allele doses
     * @throws IllegalArgumentException if
     * {@code x.nSamples() != y.nSamples()}
     * @throws NullPointerException if {@code x == null || y == null}
     */
    public static double cor(AlleleDose x, AlleleDose y) {
        if (x==null || y==null) {
            return 0.0;
        }
        if (x.dose.length != y.dose.length) {
            String s = "inconsistent number of alleles";
            throw new IllegalArgumentException(s);
        }
        int cnt = 0;
        int sumX = 0;
        int sumY = 0;
        int sumXX = 0;
        int sumXY = 0;
        int sumYY = 0;
        for(int j=0; j<x.dose.length; ++j) {
            if (x.dose[j]>=0 && y.dose[j]>=0) {
                ++cnt;
                int valX = x.dose[j];
                int valY = y.dose[j];
                sumX += valX;
                sumY += valY;
                sumXX += valX*valX;
                sumXY += valX*valY;
                sumYY += valY*valY;
            }
        }
        if (cnt==0) {
            return 0.0f;
        }
        else if ( (cnt*sumXX == (sumX*sumX)) || (cnt*sumYY == (sumY*sumY)) ) {
            return 0.0f;
        }
        else {
            double n = (double) cnt;
            double meanX = sumX/n;
            double meanY = sumY/n;
            double meanXX = sumXX/n;
            double meanXY = sumXY/n;
            double meanYY = sumYY/n;
            double stdX = Math.sqrt(meanXX - meanX*meanX);
            double stdY = Math.sqrt(meanYY - meanY*meanY);
            double covXY = meanXY - meanX*meanY;
            return  covXY / (stdX*stdY) ;
        }
    }

    /**
     * Prints a VCF record with the target data represented by
     * {@code this} using the specified {@code PrintWriter}.
     * @param out a {@code PrintWriter}.
     */
    public void printVCF(PrintWriter out) {
        out.print(refMarker);
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR);      // QUAL
        out.print(Const.tab);
        out.print("PASS");                  // FILTER
        if (refMarker.end() != -1) {
            out.print(Const.tab);
            out.print("END=");              // INFO
            out.print(refMarker.end());
        }
        else {
            out.print(Const.tab);
            out.print(Const.MISSING_DATA_CHAR);  // INFO
        }
        out.print(Const.tab);
        out.print("GT");                    // FORMAT
        for (int j=0; j<unfiltTarg.nSamples(); ++j) {
            int a1 = unfiltTarg.allele1(j);
            int a2 = unfiltTarg.allele2(j);
            out.print(Const.tab);
            out.print( a1==-1 ? Const.MISSING_DATA_CHAR : String.valueOf(alleleMap[a1]));
            out.print(unfiltTarg.isPhased(j) ? Const.phasedSep : Const.unphasedSep);
            out.print( a2==-1 ? Const.MISSING_DATA_CHAR : String.valueOf(alleleMap[a2]));
        }
        out.println();
    }

    /**
     * Returns a string representation of {@code this}.  The exact details
     * of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(100);
        sb.append("[");
        sb.append(Const.nl);
        sb.append("refMarker: ");
        sb.append(refMarker);
        sb.append(Const.nl);
        sb.append("targMarker: ");
        sb.append(filtTarg.marker());
        sb.append(Const.nl);
        sb.append("refAllele: ");
        sb.append(refAllele);
        sb.append(Const.nl);
        sb.append("alleleCnts: ");
        sb.append(Arrays.toString(alleleCnts));
        sb.append(Const.nl);
        sb.append("alleleMap: ");
        sb.append(Arrays.toString(alleleMap));
        sb.append(Const.nl);
        sb.append("]");
        return sb.toString();
    }
}
