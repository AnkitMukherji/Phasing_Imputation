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

import blbutil.Const;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.MathContext;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import main.AlleleProbs;
import main.ConstrainedAlleleProbs;
import main.GenotypeValues;

/**
 * <p>Class {@code VcfWriter} contains static methods for writing data in
 * VCF 4.2 format.
 * </p>
 * <p>Instances of class {@code VcfWriter} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfWriter {

    private static final String PASS = "PASS";
    private static final DecimalFormat df2 = new DecimalFormat("#.##");
    private static final DecimalFormat df3 = new DecimalFormat("#.###");
    private static final DecimalFormat df2_fixed = new DecimalFormat("0.00");
    private static final MathContext mc2 = new MathContext(2);
    private static final BigDecimal ONE = new BigDecimal(1.0);

    private static final String fileformat = "##fileformat=VCFv4.2";

    private static final String afInfo = "##INFO=<ID=AF,Number=A,Type=Float,"
            + "Description=\"Estimated ALT Allele Frequencies\">";
    private static final String ar2Info = "##INFO=<ID=AR2,Number=1,Type=Float,"
            + "Description=\"Allelic R-Squared: estimated squared correlation between "
            + "most probable REF dose and true REF dose\">";
    private static final String dr2Info = "##INFO=<ID=DR2,Number=1,Type=Float,"
            + "Description=\"Dosage R-Squared: estimated squared correlation between "
            + "estimated REF dose [P(RA) + 2*P(RR)] and true REF dose\">";
    private static final String impInfo = "##INFO=<ID=IMP,Number=1,Type=Flag,"
            + "Description=\"Imputed marker\">";

    private static final String gtFormat = "##FORMAT=<ID=GT,Number=1,Type=String,"
            + "Description=\"Genotype\">";
    private static final String dsFormat = "##FORMAT=<ID=DS,Number=1,Type=Float,"
            +"Description=\"estimated ALT dose [P(RA) + P(AA)]\">";
    private static final String glFormat = "##FORMAT=<ID=GL,Number=G,Type=Float,"
            + "Description=\"Log10-scaled Genotype Likelihood\">";
    private static final String gpFormat = "##FORMAT=<ID=GP,Number=G,Type=Float,"
            + "Description=\"Estimated Genotype Probability\">";

    private static final String shortChromPrefix= "#CHROM" + Const.tab + "POS"
            + Const.tab + "ID" + Const.tab + "REF" + Const.tab + "ALT"
            + Const.tab + "QUAL" + Const.tab + "FILTER" + Const.tab + "INFO";

    private static final String longChromPrefix =
            shortChromPrefix + Const.tab + "FORMAT";


    private VcfWriter() {
        // private constructor prevents instantiation
    }

    /**
     * Writes VCF meta-information lines and header line to the specified
     * {@code PrintWriter}. Only one FORMAT subfield, the GT subfield,
     * is described in the meta-information lines.
     * @param sampleIds the sample identifiers
     * @param source a description of the data source, or {@code null} if
     * no description is to be printed
     * @param out the {@code PrintWriter} to which VCF meta-information
     * lines will be written
     * @throws NullPointerException if {@code out == null}
     * @throws NullPointerException if
     * {@code sampleIds == null}, or if {@code sampleIds[j] == null} for any
     * {@code j} satisfying {@code (0 <= j && j < <sampleIds.length)}
     */
    public static void writeMetaLinesGT(String[] sampleIds, String source,
            PrintWriter out) {
        boolean printGT = true;
        boolean printGP = false;
        boolean printGL = false;
        writeMetaLines(sampleIds, source, printGT, printGP, printGL, out);
    }

    /**
     * Writes VCF meta-information lines and header line to the specified
     * {@code PrintWriter}.
     * @param sampleIds the sample identifiers
     * @param source a description of the data source, or {@code null} if
     * no description is to be printed
     * @param printGT {@code true} if the meta-information lines
     * will describe the GT FORMAT subfield and {@code false} otherwise
     * @param printGP {@code true} if the meta-information lines
     * will describe the GP FORMAT subfield and {@code false} otherwise
     * @param printGL {@code true} if the meta-information lines
     * will describe the GL FORMAT subfield and {@code false} otherwise
     * @param out the {@code PrintWriter} to which VCF meta-information lines
     * will be written.
     * @throws NullPointerException if {@code out == null}
     * @throws NullPointerException if
     * {@code sampleIds == null}, or if {@code sampleIds[j] == null} for any
     * {@code j} satisfying {@code (0 <= j && j < sampleIds.length)}
     */
    public static void writeMetaLines(String[] sampleIds, String source,
            boolean printGT, boolean printGP, boolean printGL, PrintWriter out) {
        out.print(fileformat);
        out.print(Const.nl);
        out.print("##filedate=");
        out.print(now());
        out.print(Const.nl);
        if (source != null) {
            out.print("##source=\"");
            out.print(source);
            out.println("\"");
        }
        if (printGP) {
            out.println(afInfo);
            out.println(ar2Info);
            out.println(dr2Info);
            out.println(impInfo);
        }
        if (printGT) {
            out.println(gtFormat);
        }
        if (printGL) {
            out.println(glFormat);
        }
        if (printGP) {
            out.println(dsFormat);
            out.println(gpFormat);
        }
        out.print(longChromPrefix);
        for (String id : sampleIds) {
            if (id==null) {
                throw new NullPointerException("id==null");
            }
            out.print(Const.tab);
            out.print(id);
        }
        out.println();
    }

    private static String now() {
        String dateFormat = "yyyyMMdd";
        Calendar cal = Calendar.getInstance();
        SimpleDateFormat sdf = new SimpleDateFormat(dateFormat);
        return sdf.format(cal.getTime());
    }

    /**
     * Writes the specified genotype data  as VCF records to the specified
     * {@code PrintWriter}.
     * @param gv the scaled sample posterior genotype probabilities
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @param out the {@code PrintWriter} to which VCF records will
     * be written.
     *
     * @throws IllegalArgumentException if
     * {@code haps.markers().equals(gv.markers()) == false}
     * @throws IndexOutOfBoundsException if
     * {@code (start < 0 || start > end || end > haps.nMarkers())}
     * @throws NullPointerException if
     * {@code (gv == null || out == null)}
     */
    public static void appendRecords(GenotypeValues gv, int start, int end,
            PrintWriter out) {
        if (start > end) {
            throw new IllegalArgumentException("start=" + start + " end=" + end);
        }
        for (int marker=start; marker<end; ++marker) {
            printFixedFields(gv, marker, out);
            for (int sample=0, n=gv.nSamples(); sample<n; ++sample) {
                print_GT_DS_GP(gv, marker, sample, out);
            }
            out.println();
        }
    }


    private static void print_GT_DS_GP(GenotypeValues gv, int marker, int sample,
            PrintWriter out) {
        int nAlleles = gv.marker(marker).nAlleles();
        int nGenotypes = gv.marker(marker).nGenotypes();
        float[] dose = new float[nAlleles];
        int bestA1 = -1;
        int bestA2 = -1;
        int gt = 0;
        float sum = 0f;
        float maxGP = 0f;
        for (int a2=0; a2<nAlleles; ++a2) {
            for (int a1=0; a1<=a2; ++a1) {
                float value = gv.value(marker, sample, gt++);
                if (value > maxGP) {
                    bestA1 = a1;
                    bestA2 = a2;
                    maxGP = value;
                }
                dose[a1] += value;
                dose[a2] += value;
                sum += value;
            }
        }
        out.print(Const.tab);
        out.print(bestA1 == -1 ? Const.MISSING_DATA_STRING : bestA1);
        out.print(Const.unphasedSep);
        out.print(bestA2 == -1 ? Const.MISSING_DATA_STRING : bestA2);
        for (int al = 1; al < dose.length; ++al) {
            out.print( (al==1) ? Const.colon : Const.comma);
            out.print(df2.format(dose[al]/sum));
        }
        for (gt=0; gt<nGenotypes; ++gt) {
            out.print(gt==0 ? Const.colon : Const.comma);
            double v = gv.value(marker, sample, gt)/sum;
            out.print(df2.format(v));
        }
    }

    /**
     * Writes the specified genotype data as VCF records to the specified
     * {@code PrintWriter}.
     * @param alProbs the sample haplotype pairs
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @param gprobs {@code true} if the GP field should be printed, and
     * {@code false} otherwise.
     * @param out the {@code PrintWriter} to which VCF records will
     * be written
     *
     * @throws IndexOutOfBoundsException if
     * {@code (start < 0 || start > end || end > alProbs.nMarkers())}
     * @throws NullPointerException if {@code haps == null || out == null}
     */
    public static void appendRecords(ConstrainedAlleleProbs alProbs,
            int start, int end, boolean gprobs, PrintWriter out) {
        if (start > end) {
            throw new IllegalArgumentException("start=" + start + " end=" + end);
        }
        boolean ds = true;
        String format = format(ds, gprobs);
        for (int m=start; m<end; ++m) {
            Marker marker = alProbs.marker(m);
            boolean isImputed = alProbs.isImputed(m);
            GprobsStatistics gps = new GprobsStatistics(alProbs, m);
            String info = info(gps, isImputed);
            printFixedFields(marker, info, format, out);
            for (int sample=0, n=alProbs.nSamples(); sample<n; ++sample) {
                printGTandDose(alProbs, m, sample, ds, out);
                if (gprobs) {
                     printGP(alProbs, m, sample, out);
                }
            }
            out.println();
        }
    }

    /**
     * Writes the specified genotype data as VCF records to the specified
     * {@code PrintWriter}.
     * @param alProbs the sample haplotype pairs
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @param out the {@code PrintWriter} to which VCF records will
     * be written
     *
     * @throws IndexOutOfBoundsException if
     * {@code (start < 0 || start > end || end > alProbs.nMarkers())}
     * @throws NullPointerException if {@code haps == null || out == null}
     */
    public static void appendRecords(AlleleProbs alProbs, int start, int end,
            PrintWriter out) {
        boolean ds = false;
        boolean gp = false;
        String info = Const.MISSING_DATA_STRING;
        String format = format(ds, gp);
        if (start > end) {
            throw new IllegalArgumentException("start=" + start + " end=" + end);
        }
        for (int m=start; m<end; ++m) {
            Marker marker = alProbs.marker(m);
            printFixedFields(marker, info, format, out);
            for (int sample=0, n=alProbs.nSamples(); sample<n; ++sample) {
                printGTandDose(alProbs, m, sample, ds, out);
            }
            out.println();
        }
    }

    private static void printGTandDose(AlleleProbs alProbs, int marker, int
            sample, boolean ds, PrintWriter out) {
        out.print(Const.tab);
        out.print(alProbs.allele1(marker, sample));
        out.append(Const.phasedSep);
        out.print(alProbs.allele2(marker, sample));
        if (ds) {
            int nAlleles = alProbs.marker(marker).nAlleles();
            for (int j = 1; j < nAlleles; ++j) {
                float p1 = alProbs.alProb1(marker, sample, j);
                float p2 = alProbs.alProb2(marker, sample, j);
                out.print( (j==1) ? Const.colon : Const.comma );
                out.print(df2.format(p1 + p2));
            }
        }
    }

    private static void printGP(AlleleProbs alProbs, int marker, int sample,
            PrintWriter out) {
        int nAlleles = alProbs.marker(marker).nAlleles();
        for (int a2=0; a2<nAlleles; ++a2) {
            for (int a1=0; a1<=a2; ++a1) {
                out.print((a2 == 0 && a1 == 0) ? Const.colon : Const.comma);
                float gtProb = alProbs.gtProb(marker, sample, a1, a2);
                if (a1 != a2) {
                    gtProb += alProbs.gtProb(marker, sample, a2, a1);
                }
                out.print(df2.format(gtProb));
            }
        }
    }

    /**
     * Prints the first 9 VCF record fields for the specified marker to
     * the specified {@code PrintWriter}.  Only one VCF FORMAT subfield,
     * the GT subfield, is printed.
     *
     * @param marker a marker
     * @param out the {@code PrintWriter} to which the first 9 VCF record
     * fields will be written
     *
     * @throws NullPointerException if {@code marker == null || out == null}
     */
    public static void printFixedFieldsGT(Marker marker, PrintWriter out) {
        out.print(marker);
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // QUAL
        out.print(Const.tab);
        out.print(PASS);                    // FILTER
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // INFO
        out.print(Const.tab);
        out.print("GT");
    }

    private static void printFixedFields(GenotypeValues gv, int marker,
            PrintWriter out) {
        GprobsStatistics gpm = new GprobsStatistics(gv, marker);
        double[] alleleFreq = gpm.alleleFreq();
        out.print(gv.marker(marker));
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // QUAL
        out.print(Const.tab);
        out.print(PASS);                    // FILTER
        out.print(Const.tab);
        out.print("AR2=");                  // INFO
        out.print(df2_fixed.format(gpm.allelicR2()));
        out.print(";DR2=");
        out.print(df2_fixed.format(gpm.doseR2()));
        for (int j=1; j<alleleFreq.length; ++j) {
            out.print( (j==1) ? ";AF=" : Const.comma);
            out.print(formatProb(alleleFreq[j]));
        }
        out.print(Const.tab);
        out.print("GT:DS:GP");
    }

    private static String info(GprobsStatistics gps, boolean isImputed) {
        double[] alleleFreq = gps.alleleFreq();
        StringBuilder sb = new StringBuilder(20);
        sb.append("AR2=");                  // INFO
        sb.append(df2_fixed.format(gps.allelicR2()));
        sb.append(";DR2=");
        sb.append(df2_fixed.format(gps.doseR2()));
        for (int j=1; j<alleleFreq.length; ++j) {
            sb.append( (j==1) ? ";AF=" : Const.comma);
            sb.append(formatProb(alleleFreq[j]));
        }
        if (isImputed) {
            sb.append(";IMP");
        }
        return sb.toString();
    }

    private static String format(boolean ds, boolean gp) {
        if (ds) {
            if (gp) {
                return "GT:DS:GP";
            }
            else {
                return "GT:DS";
            }
        }
        else {
            return "GT";
        }
    }

    private static void printFixedFields(Marker marker, String info,
            String format, PrintWriter out) {
        out.print(marker);
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // QUAL
        out.print(Const.tab);
        out.print(PASS);                    // FILTER
        out.print(Const.tab);
        out.print(info);
        out.print(Const.tab);
        out.print(format);
    }

    private static String formatProb(double d) {
        if (d>=0 && d <= 0.5) {
            return new BigDecimal(d).round(mc2).toPlainString();
        }
        else if (d <= 1.0) {
            return new BigDecimal(d-1.0).round(mc2).add(ONE).toString();
        }
        else {
            throw new IllegalArgumentException(String.valueOf(d));
        }
    }
}
