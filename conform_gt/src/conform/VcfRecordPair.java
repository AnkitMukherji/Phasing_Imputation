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
import blbutil.Const;
import java.io.PrintWriter;
import vcf.VcfEmission;

/**
 * Class {@code VcfRecordPair} represents reference and target data
 * for a marker.  Class {@code VcfRecordPair} is immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfRecordPair {

    private final VcfEmission ref;
    private final VcfEmission targ;
    private final Phase relStrand;

    private final AlleleDose refAlleleDose;
    private final AlleleDose targAlleleDose;
    private final AlleleDose flippedTargAlleleDose;

    /**
     * Constructs a {@code VcfRecordPair} instance for two VCF records
     * corresponding to the same marker.
     * @param ref VCF record for reference file
     * @param filtTarg the target genotype data (after sample filtering)
     * @param unfiltTarg the target genotype data (before sample filtering)
     * @param relStrand the possible relative strands for reference and target
     * markers.
     * @throws IllegalArgumentException if
     * {@code filtTarg.marker().equals(unfiltTarg.marker() == false}
     * @throws IllegalArgumentException if the alleles in the target VCF
     * record marker (after flipping alleles if {@code flipAlleles == true})
     * are not a subset of the alleles in the reference marker
     * @throws NullPointerException if
     * {@code ref == null || filtTarg == null || unfiltTarg == null}
     */
    public VcfRecordPair(VcfEmission ref, VcfEmission filtTarg,
            VcfEmission unfiltTarg, Phase relStrand) {
        this.ref = ref;
        this.targ = filtTarg;
        this.relStrand = relStrand;
        this.refAlleleDose = new AlleleDose(ref.marker(), 0, ref, ref, false);

        if (relStrand==Phase.UNKNOWN || relStrand==Phase.IDENTICAL) {
            this.targAlleleDose = new AlleleDose(ref.marker(), 0, filtTarg,
                    unfiltTarg, false);
        }
        else {
            this.targAlleleDose = null;
        }
        if (relStrand==Phase.UNKNOWN || relStrand==Phase.OPPOSITE) {
            this.flippedTargAlleleDose = new AlleleDose(ref.marker(), 0,
                    filtTarg, unfiltTarg, true);
        }
        else {
            this.flippedTargAlleleDose = null;
        }
    }

    private static double absZ(AlleleDose cntA, AlleleDose cntB) {
        if (cntA==null || cntB==null) {
            return Double.MAX_VALUE;
        }
        int xCnt = cntA.count(0);
        int yCnt = cntB.count(0);
        int nx = cntA.nNonmissingAlleles();
        int ny = cntB.nNonmissingAlleles();
        if (nx==0 || ny==0) {
            return Double.MAX_VALUE;
        }
        else if ((xCnt + yCnt)==0 || (xCnt+yCnt)==(nx + ny)) {
            return 0.0f;
        }
        else {
            double px = xCnt / (double) nx;
            double py = yCnt / (double) ny;
            double p = (xCnt + yCnt) / (double) (nx + ny);
            double var = ( (1.0/nx) + (1.0/ny) )* p * (1-p);
            return Math.abs(px - py)/Math.sqrt(var);
        }
    }

    /**
     * Return the reference genotypes.
     * @return the reference genotypes
     */
    public VcfEmission ref() {
        return ref;
    }

    /**
     * Return the target genotypes.
     * @return the target genotypes
     */
    public VcfEmission targ() {
        return targ;
    }

    /**
     * Returns the frequency of {@code this.refAllele()} in the
     * reference data.  The returned value is {@code Double.NaN} if
     * there are no non-missing alleles in the reference data.
     * @param allele the allele index
     * @return the frequency of {@code this.refAllele()} in the
     * reference data
     */
    public double refFreq(int allele) {
        return alleleFreq(refAlleleDose, allele);
    }

    /**
     * Returns the frequency of {@code this.refAllele()} in the
     * target data.  The returned value is {@code Double.NaN} if
     * there are no non-missing alleles in the target data.
     * @param allele a reference marker allele index
     * @return the frequency of {@code this.refAllele()} in the
     * target data
     */
    public double targetFreq(int allele) {
        return alleleFreq(targAlleleDose, allele);
    }

    /**
     * Returns the frequency of {@code this.refAllele()} in the
     * flipped target data.  The returned value is {@code Double.NaN} if
     * there are no non-missing alleles in the target data.
     * @param allele a reference marker allele index
     * @return the frequency of {@code this.refAllele()} in the
     * flipped target data
     */
    public double flippedTargetFreq(int allele) {
        if (targAlleleDose==null) {
            return Double.NaN;
        }
        return alleleFreq(flippedTargAlleleDose, allele);
    }

    private static double alleleFreq(AlleleDose alDose, int allele) {
        if (alDose==null) {
            return Double.NaN;
        }
        int den = alDose.nNonmissingAlleles();
        if (den==0) {
            return Double.NaN;
        }
        else {
            return alDose.count(allele) / (double) den;
        }
    }

    /**
     * Returns the absolute z-statistic of a test for equal reference
     * allele frequencies in the reference and target data.
     * The returned value is {@code Double.NaN} if
     * there are no non-missing alleles in the reference data or
     * in the target data.
     * @return the absolute z-statistic of a test for equal reference
     * allele frequencies in the reference and target data.
     */
    public double absZ() {
        return absZ(refAlleleDose, targAlleleDose);
    }

    /**
     * Returns the absolute z-statistic of a test for equal reference
     * allele frequencies in the reference and strand-flipped target data.
     * The returned value is {@code Double.NaN} if
     * there are no non-missing alleles in the reference data or
     * in the strand-flipped target data.
     * @return the absolute z-statistic of a test for equal reference
     * allele frequencies in the reference and strand-flipped target data
     */
    public double flippedAbsZ() {
        return absZ(refAlleleDose, flippedTargAlleleDose);
    }

    /**
     * Prints a VCF record with the target data represented by
     * {@code this} using the specified {@code PrintWriter}.
     * @param out a {@code PrintWriter}
     */
    public void printTarget(PrintWriter out) {
        targAlleleDose.printVCF(out);
    }

    /**
     * Prints a VCF record with the strand-flipped target data represented by
     * {@code this} using the specified {@code PrintWriter}.
     * @param out a {@code PrintWriter}
     */
    public void printFlippedTarget(PrintWriter out) {
        flippedTargAlleleDose.printVCF(out);
    }

    public Phase allelePhase() {
        return relStrand;
    }

    /**
     * Returns a string representation of {@code this}.  The exact details
     * of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(100);
        sb.append(relStrand);
        sb.append(Const.nl);
        sb.append(refAlleleDose.toString());
        sb.append(targAlleleDose.toString());
        sb.append(flippedTargAlleleDose.toString());
        return sb.toString();
    }

    /**
     * Returns the correlation of the reference allele dosage
     * for the specified two markers for the reference samples.
     * Returns 0.0 if the correlation is undefined.
     * @param x reference and target data for a marker.
     * @param y reference and target data for a marker.
     * @return  the correlation of the reference allele dosage
     * for the specified two markers for the reference samples.
     */
    public static double refCor(VcfRecordPair x, VcfRecordPair y) {
        return AlleleDose.cor(x.refAlleleDose, y.refAlleleDose);
    }

    /**
     * Returns the correlation of the reference allele dosage
     * for the specified two markers for the target samples.  Returns
     * 0.0 if the correlation is undefined.
     * @param x reference and target data for a marker
     * @param flipX {@code true} if strand-flipped target data
     * should be used for {@code x}, and {@code false} otherwise
     * @param y reference and target data for a marker.
     * @param flipY {@code true} if strand-flipped target data
     * should be used for {@code y}, and {@code false} otherwise
     * @return  the correlation of the reference allele dosage
     * for the specified two markers for the target samples
     * @throws IllegalArgumentException if
     * {@code x.targ.nSamples() != y.targ.nSamples()}
     * @throws NullPointerException if {@code x == null || y == null}
     */
    public static double targetCor(VcfRecordPair x, boolean flipX,
            VcfRecordPair y, boolean flipY) {
        if (flipX) {
            if (flipY) {
                return AlleleDose.cor(x.flippedTargAlleleDose,
                        y.flippedTargAlleleDose);
            }
            else {
                return AlleleDose.cor(x.flippedTargAlleleDose, y.targAlleleDose);
            }
        }
        else {
            if (flipY) {
                return AlleleDose.cor(x.targAlleleDose, y.flippedTargAlleleDose);
            }
            else {
                return AlleleDose.cor(x.targAlleleDose, y.targAlleleDose);
            }
        }
    }
}
