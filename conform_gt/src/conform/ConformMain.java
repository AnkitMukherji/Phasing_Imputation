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

import blbutil.SampleFileIt;
import beagleutil.ChromInterval;
import beagleutil.Phase;
import static beagleutil.Phase.IDENTICAL;
import static beagleutil.Phase.INCONSISTENT;
import static beagleutil.Phase.OPPOSITE;
import static beagleutil.Phase.UNKNOWN;
import blbutil.Const;
import blbutil.FileIt;
import blbutil.FileUtil;
import blbutil.Filter;
import blbutil.InputIt;
import blbutil.Utilities;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import vcf.FilterUtil;
import vcf.IntervalVcfIt;
import vcf.Marker;
import vcf.VcfIt;
import vcf.VcfEmission;
import vcf.VcfWriter;

/**
 * <p>Class {@code ConformGtMain} conforms genotype data in a VCF file
 * to match the genotype data in a reference file.
 * </p>
 * <p>
 * Reference: Fisher, R.A. Statistical Methods, Experimental Design, and
 * Scientific Inference.  Oxford: Oxford University Press, 1990.
 * Section 33 (p. 194).
 * </p>
 * <p>
 * Reference:
 *  http://demonstrations.wolfram.com/NullDistributionOfTheCorrelationCoefficient/
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ConformMain {

    private static final String version = "conform-gt.24May16.cee";
    public static final int MIN_NSAMPLES = 20;

    private final int maxInconsistentStrand = 1;
    private final int minStrandDiff = 2;

    private final int windowOverlap = 100;

    private final ConformPar par;
    private final File logFile;
    private final File conformedFile;
    private final Filter<String> sampleFilter;
    private final ChromInterval chromInterval;
    private final MatchedMarkers matchedMarkers;

    private int globalIndex = -1;
    private VcfRecordPair[] window = new VcfRecordPair[0];
    private Phase[] allelePhase = new Phase[0];
    private Phase[] freqPhase = new Phase[0];
    private Phase[] corPhase = new Phase[0];

    /**
     * Main method for {@code ConformGtMain} class.
     * See {@code ConformGtPar.usage()} method for command line
     * syntax.
     * @param args command line arguments
     */
    public static void main(String[] args) {
        if (args.length==0) {
            System.out.println(ConformPar.usage());
            System.exit(0);
        }
        Locale.setDefault(Locale.US);
        ConformPar par = parameters(args);
	ConformMain main = new ConformMain(par);
        main.run();
    }

    public ConformMain(ConformPar par) {
        this.par = par;
        this.logFile = new File(par.out() + ".log");
        this.conformedFile = new File(par.out() + ".vcf.gz");
        this.sampleFilter = FilterUtil.sampleFilter(par.excludesamples());
        this.chromInterval = ChromInterval.parse(par.chrom());

        ConformMarkers refMarkers = new ConformMarkers(par.ref(), chromInterval);
        ConformMarkers targMarkers = new ConformMarkers(par.gt(), chromInterval);

        try (PrintWriter log=FileUtil.printWriter(logFile)) {
            printLogHeader(log);
            this.matchedMarkers = new MatchedMarkers(refMarkers, targMarkers,
                    par.matchID(), log);
        }
    }

    private void run() {
        Filter<String> acceptAll = Filter.acceptAllFilter();
        boolean appendLog = true;
        try (SampleFileIt<VcfEmission> refIt = vcfIt(par.ref(), sampleFilter,
                        chromInterval);
                SampleFileIt<VcfEmission> targIt1 = vcfIt(par.gt(), sampleFilter,
                        chromInterval);
                SampleFileIt<VcfEmission> targIt2 = vcfIt(par.gt(), acceptAll,
                        chromInterval);
                PrintWriter log = FileUtil.printWriter(logFile, appendLog);
                PrintWriter out = FileUtil.bgzipPrintWriter(conformedFile)) {

            VcfWriter.writeMetaLinesGT(targIt2.samples().ids(), version, out);
            while (finished()==false) {
                advanceWindow(refIt, targIt1, targIt2, log, out);
            }
            printMarkers(log, out, 0, window.length);
        }
    }

    private static SampleFileIt<VcfEmission> vcfIt(File file,
            Filter<String> sampleFilter, ChromInterval chromInterval) {
        FileIt<String> it1 = InputIt.fromGzipFile(file);
        SampleFileIt<VcfEmission> it2 = VcfIt.create(it1, sampleFilter,
                        VcfIt.toBitSetGT);
        return new IntervalVcfIt<>(it2, chromInterval);
    }

    private void printMarkers(PrintWriter log, PrintWriter out, int start,
            int end) {
        for (int j=start; j<end; ++j) {
            Phase mergedPhase = allelePhase[j];
            if (par.strict()==true || allelePhase[j]==UNKNOWN) {
                mergedPhase = mergePhase(allelePhase[j], freqPhase[j],
                        corPhase[j]);
            }
            printLog(window[j].targ().marker(), allelePhase[j], freqPhase[j],
                    corPhase[j], mergedPhase, log);
            if (mergedPhase == IDENTICAL) {
                window[j].printTarget(out);
            }
            else if (mergedPhase == OPPOSITE) {
                window[j].printFlippedTarget(out);
            }
        }
    }

    private void advanceWindow(SampleFileIt<VcfEmission> refIt,
            SampleFileIt<VcfEmission> filtTargIt,
            SampleFileIt<VcfEmission> unfiltTargIt,
            PrintWriter log, PrintWriter out) {
        int newWindowSize = 2*windowOverlap;
        List<VcfRecordPair> newWindow = new ArrayList<>(newWindowSize);
        List<Phase> newAllelePhase = new ArrayList<>(newWindowSize);
        List<Phase> newFreqPhase = new ArrayList<>(newWindowSize);
        List<Phase> newCorPhase = new ArrayList<>(newWindowSize);

        int overlap = Math.min(windowOverlap, window.length);
        int overlapStart = window.length - overlap;
        printMarkers(log, out, 0, overlapStart);
        for (int j=overlapStart; j<window.length; ++j) {
            newWindow.add(window[j]);
            newAllelePhase.add(allelePhase[j]);
            newFreqPhase.add(freqPhase[j]);
            newCorPhase.add(corPhase[j]);
        }
        for (int j=overlap; j<newWindowSize && finished()==false; ++j) {
            VcfRecordPair pair = nextConformPair(refIt, filtTargIt, unfiltTargIt);
            newWindow.add(pair);
            newAllelePhase.add(pair.allelePhase());
            newFreqPhase.add(freqPhase(pair));
            newCorPhase.add(UNKNOWN);
        }
        this.window = newWindow.toArray(new VcfRecordPair[0]);
        this.allelePhase = newAllelePhase.toArray(new Phase[0]);
        this.freqPhase = newFreqPhase.toArray(new Phase[0]);
        this.corPhase = newCorPhase.toArray(new Phase[0]);
        updateStrandUsingCorrelation();
    }

    private boolean finished() {
        return (globalIndex+1) >= matchedMarkers.nMarkers();
    }

    private VcfRecordPair nextConformPair(SampleFileIt<VcfEmission> refIt,
            SampleFileIt<VcfEmission> filtTargIt,
            SampleFileIt<VcfEmission> unfiltTargIt) {
        ++globalIndex;
        Marker refMarker = matchedMarkers.refMarker(globalIndex);
        Marker targetMarker = matchedMarkers.targetMarker(globalIndex);
        Phase phase = matchedMarkers.strand(globalIndex);
        VcfEmission refRec = readRec(refIt, refMarker);
        VcfEmission filtTargRec = readRec(filtTargIt, targetMarker);
        VcfEmission unfiltTargRec = readRec(unfiltTargIt, targetMarker);
        if (refRec==null || filtTargRec==null || unfiltTargRec==null) {
            String s = "ERROR: Modification detected to input VCF files"
                    + Const.nl + refMarker
                    + Const.nl + targetMarker;
            throw new IllegalArgumentException(s);
        }
        else {
            return new VcfRecordPair(refRec, filtTargRec, unfiltTargRec, phase);
        }
    }

    private static VcfEmission readRec(SampleFileIt<VcfEmission> it,
            Marker marker) {
        VcfEmission rec = null;
        while (it.hasNext() && rec==null) {
            VcfEmission candidate = it.next();
            if (marker.equals(candidate.marker())) {
                rec = candidate;
            }
        }
        return rec;
    }

    private Phase freqPhase(VcfRecordPair pair) {
        double minZDiff = 4.0;
        double absZ = pair.absZ();
        double flippedAbsZ = pair.flippedAbsZ();
        if (flippedAbsZ >= (absZ + minZDiff)) {
            return IDENTICAL;
        }
        else if (absZ >= (flippedAbsZ + minZDiff)) {
            return OPPOSITE;
        }
        else {
            return UNKNOWN;
         }
    }

    private CorCounts[] getCorCounts() {
        CorCounts[] corCounts = new CorCounts[window.length];
        for (int j=0; j<window.length; ++j) {
            corCounts[j] = getCorCounts(j);
        }
        return corCounts;
    }

    private CorCounts getCorCounts(int index) {
        VcfRecordPair focus = window[index];
        double minAbsRefCor = minAbsCor(focus.refFreq(0), focus.ref().nSamples());
        double minAbsTargetCor = minAbsCor(focus.targetFreq(0), focus.targ().nSamples());
        int sameCnt = 0;
        int oppCnt = 0;
        int informativeCnt = 0;
        for (int j=0; j<window.length; ++j) {
            if (j!=index
                    && (freqPhase[j]==IDENTICAL || freqPhase[j]==OPPOSITE)
                    && (allelePhase[j]==freqPhase[j] || allelePhase[j]==UNKNOWN)) {
                VcfRecordPair anchor = window[j];
                boolean flipAnchor = (freqPhase[j]==OPPOSITE);
                double refCor = VcfRecordPair.refCor(focus, anchor);
                if (Math.abs(refCor) > minAbsRefCor) {
                    ++informativeCnt;
                    double cor = VcfRecordPair.targetCor(focus, false, anchor,
                            flipAnchor);
                    double fCor = VcfRecordPair.targetCor(focus, true, anchor,
                            flipAnchor);
                    if (refCor < -minAbsRefCor) {
                        if (cor < -minAbsTargetCor) {
                            ++sameCnt;
                        }
                        if (fCor < -minAbsTargetCor) {
                            ++oppCnt;
                        }
                    }
                    else if (refCor > minAbsRefCor) {
                        if (cor > minAbsTargetCor) {
                            ++sameCnt;
                        }
                        if (fCor > minAbsTargetCor) {
                            ++oppCnt;
                        }
                    }
                }
            }
        }
        return new CorCounts(sameCnt, oppCnt, informativeCnt);
    }

    /*
     * This method uses the following approximation:
     * for a large sample of N pairs of values,
     * for small or moderate correlations, the sample correlation is
     * distributed normally about the true correlation rho with
     * variance Math.pow((1 - rho*rho), 2)/(N-1).
     *
     * See Section 33 (p. 194) of
     * R.A. Statistical Methods, Experimental Design, and
     * Scientific Inference.  Oxford: Oxford University Press, 1990.
     */
    private static double minAbsCor(double freq, int nSamples) {
        double minHighFreqStdDev = 5.0;
        double minLowFreqStdDev = 7.0;
        double stdDev = 1.0/Math.sqrt(nSamples - 1);
        if (freq > 0.3 && freq < 0.7) {
            return minHighFreqStdDev * stdDev;
        }
        else {
            return minLowFreqStdDev * stdDev;
        }
    }

    private void updateStrandUsingCorrelation() {
        CorCounts[] corCounts = getCorCounts();
        assert corPhase.length == corCounts.length;
        for (int j=0; j<corPhase.length; ++j) {
            CorCounts cc = corCounts[j];
            Phase newCorPhase = strandFromCorCounts(cc);
            if (newCorPhase==INCONSISTENT) {
                corPhase[j] = INCONSISTENT;
            }
            else {
                switch (corPhase[j]) {
                    case UNKNOWN:
                        corPhase[j] = newCorPhase;
                        break;
                    case IDENTICAL:
                        if (newCorPhase==OPPOSITE) {
                            corPhase[j] = INCONSISTENT;
                        }
                        break;
                    case OPPOSITE:
                        if (newCorPhase==IDENTICAL) {
                            corPhase[j] = INCONSISTENT;
                        }
                        break;
                    case INCONSISTENT:
                        break;
                }
            }
        }
    }

    private Phase strandFromCorCounts(CorCounts cc) {
        if (cc.oppCnt()<=maxInconsistentStrand
                && (cc.sameCnt()-cc.oppCnt()) >= minStrandDiff ) {
            return IDENTICAL;
        }
        else if (cc.sameCnt()<=maxInconsistentStrand
                && (cc.oppCnt()-cc.sameCnt()) >= minStrandDiff ) {
            return OPPOSITE;
        }
        else if (cc.sameCnt()>maxInconsistentStrand
                && cc.oppCnt()>maxInconsistentStrand) {
            return INCONSISTENT;
        }
        else {
            return UNKNOWN;
        }
    }

    private static Phase mergePhase(Phase allelePhase, Phase freqPhase,
            Phase corPhase) {
        Phase freqCorPhase = mergePhase(freqPhase, corPhase);
        switch (freqCorPhase) {
            case IDENTICAL:
                if (allelePhase==IDENTICAL || allelePhase==UNKNOWN) {
                    return IDENTICAL;
                }
                else {
                    return INCONSISTENT;
                }
            case OPPOSITE:
                if (allelePhase==OPPOSITE || allelePhase==UNKNOWN) {
                    return OPPOSITE;
                }
                else {
                    return INCONSISTENT;
                }
            case UNKNOWN:
                return UNKNOWN;
            case INCONSISTENT:
                return INCONSISTENT;
            default:
                return INCONSISTENT;
        }
    }

    private static Phase mergePhase(Phase phaseA, Phase phaseB) {
        switch (phaseA) {
            case IDENTICAL:
                if (phaseB==IDENTICAL || phaseB==UNKNOWN) {
                    return IDENTICAL;
                }
                else {
                    return INCONSISTENT;
                }
            case OPPOSITE:
                if (phaseB==OPPOSITE || phaseB==UNKNOWN) {
                    return OPPOSITE;
                }
                else {
                    return INCONSISTENT;
                }
            case UNKNOWN:
                return phaseB;
            case INCONSISTENT:
                return INCONSISTENT;
            default:
                return INCONSISTENT;
        }
    }

    private static void printLogHeader(PrintWriter log) {
        log.print("CHROM");
        log.print(Const.tab);
        log.print("POS");
        log.print(Const.tab);
        log.print("ID");
        log.print(Const.tab);
        log.print("REF");
        log.print(Const.tab);
        log.print("ALT");
        log.print(Const.tab);
        log.print("ALLELE");
        log.print(Const.tab);
        log.print("FREQ");
        log.print(Const.tab);
        log.print("R2");
        log.print(Const.tab);
        log.print("SUMMARY");
        log.print(Const.tab);
        log.println("INFO");
    }

    static void printUnmatchedMarkerLog(Marker marker, String summary,
            PrintWriter log) {
        log.print(marker);
        log.print(Const.tab);
        log.print("NOT_PERFORMED");
        log.print(Const.tab);
        log.print("NOT_PERFORMED");
        log.print(Const.tab);
        log.print("NOT_PERFORMED");
        log.print(Const.tab);
        log.print("REMOVED");
        log.print(Const.tab);
        log.println(summary);
    }

    private static void printLog(Marker marker, Phase allele, Phase freq,
            Phase cor, Phase combined, PrintWriter log) {
        log.print(marker);
        log.print(Const.tab);
        log.print(summary(allele));
        log.print(Const.tab);
        log.print(summary(freq));
        log.print(Const.tab);
        log.print(summary(cor));
        log.print(Const.tab);
        log.print(disposition(combined));
        log.print(Const.tab);
        log.println(summary(combined));
    }

    private static String summary(Phase phase) {
        switch(phase) {
            case IDENTICAL:
                return "SAME_STRAND";
            case OPPOSITE:
                return "OPPOSITE_STRAND";
            case UNKNOWN:
                return "UNKNOWN_STRAND";
            case INCONSISTENT:
                return "INCONSISTENT_STRAND";
            default:
                return "INCONSISTENT_STRAND";
        }
    }

    private static String disposition(Phase phase) {
        switch(phase) {
            case IDENTICAL:
                return "PASS";
            case OPPOSITE:
                return "PASS";
            case UNKNOWN:
                return "FAIL";
            case INCONSISTENT:
                return "FAIL";
            default:
                return "FAIL";
        }
    }

    private static final class CorCounts {
        private int sameCnt;
        private int oppCnt;
        private int informativeCnt;

        public CorCounts(int sameCnt, int oppCnt, int informativeCnt) {
            if (sameCnt < 0) {
                throw new IllegalArgumentException("sameCnt: " + sameCnt);
            }
            if (oppCnt < 0) {
                throw new IllegalArgumentException("oppCnt: " + oppCnt);
            }
            if (informativeCnt < 0) {
                throw new IllegalArgumentException("informativeCnt: "
                        + informativeCnt);
            }
            this.sameCnt = sameCnt;
            this.oppCnt = oppCnt;
            this.informativeCnt = informativeCnt;
        }

        public int sameCnt() {
            return sameCnt;
        }

        public int oppCnt() {
            return oppCnt;
        }

        public int informativeCnt() {
            return informativeCnt;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder(30);
            sb.append("same=");
            sb.append(sameCnt);
            sb.append(" opp=");
            sb.append(oppCnt);
            sb.append(" informative=");
            sb.append(informativeCnt);
            return sb.toString();
        }
    }

    /*
     * Checks that specified parameters are consistent, and prints error
     * message and exits if parameters are inconsistent.
     *
     * @param args the command line arguments.
     */
    private static ConformPar parameters(String[] args) {
        ConformPar par = new ConformPar(args);
        File outPrefix = new File(par.out());
        if (outPrefix.isDirectory()) {
            String s = "ERROR: \"out\" parameter cannot be a directory: \""
                    + par.out() + "\"";
            Utilities.exit(ConformPar.usage() + s);
        }
        File vcfOut = new File(par.out() + ".vcf.gz");
        File logOut = new File(par.out() + ".log");
        if (vcfOut.exists()) {
            String s = "ERROR: VCF output file already exists: " + vcfOut;
            Utilities.exit(ConformPar.usage() + s);
        }
        if (par.ref().equals(vcfOut)) {
            String s = "ERROR: output VCF file equals reference file: " + par.ref();
            Utilities.exit(ConformPar.usage() + s);
        }
        if (par.gt().equals(vcfOut)) {
            String s = "ERROR: output VCF file equals \"gt\" file: " + par.gt();
            Utilities.exit(ConformPar.usage() + s);
        }
        if (par.ref().equals(logOut)) {
            String s = "ERROR: output log file equals reference file: " + par.ref();
            Utilities.exit(ConformPar.usage() + s);
        }
        if (par.gt().equals(logOut)) {
            String s = "ERROR: output log file equals \"gt\" file: " + par.gt();
            Utilities.exit(ConformPar.usage() + s);
        }
        if (ChromInterval.parse(par.chrom())==null) {
            String s = "ERROR: invalid \"chrom\" parameter: \""
                    + par.chrom() + "\"";
            Utilities.exit(ConformPar.usage() + s);
        }
        return par;
    }
}
