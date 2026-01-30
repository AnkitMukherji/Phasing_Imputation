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
import blbutil.Validate;
import java.io.File;
import java.util.Map;

/**
 * Class {@code ConformGtPar} represents the command line
 * parameters and default parameters used by the {@code ConformGtMain}
 * class.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ConformPar {

    private final File ref;
    private final File gt;
    private final String chrom;
    private final String out;
    private final File excludesamples;
    private final boolean strict;
    private final boolean matchID;

    public ConformPar(String[] args) {
        Map<String, String> argsMap = Validate.argsToMap(args, '=');

        ref = Validate.getFile(Validate.stringArg("ref", argsMap, true, null, null));
        gt= Validate.getFile(Validate.stringArg("gt", argsMap, true, null, null));
        chrom = Validate.stringArg("chrom", argsMap, true, null, null);
        out = Validate.stringArg("out", argsMap, true, null, null);

        strict = Validate.booleanArg("strict", argsMap, false, false);
        String match = Validate.stringArg("match", argsMap, false, "ID",
                new String[] {"ID", "POS", "id", "pos"});
        matchID = (match.equalsIgnoreCase("ID"));
        excludesamples = Validate.getFile(Validate.stringArg("excludesamples",
                argsMap, false, null, null));
        Validate.confirmEmptyMap(argsMap);
    }

    public static String usage() {
        String nl = Const.nl;
        return  nl + "usage: java -jar conform-gt.24May16.cee.jar [arguments]" + nl + nl
                + "where [arguments] have the format" + nl
                + "  ref=<reference VCF file with GT data>                         (required)" + nl
                + "  gt=<target VCF file with GT data>                             (required)" + nl
                + "  chrom=<[chrom] or [chrom]:[start]-[end]>                      (required)" + nl
                + "  out=<output file prefix>                                      (required)" + nl
                + "  match=<ID or POS (field for matching VCF records)>            (default: ID)" + nl
                + "  strict=<true if strand alignment requires MAF or R2 evidence> (default: false)" + nl
                + "  excludesamples=<file with 1 sample ID per line>               (optional)" + nl + nl

                + "Two output files are created:" + nl
                + "  <out prefix>.vcf.gz - reference-matched target data." + nl
                + "  <out prefix>.log    - summary of result for each target marker." + nl;
    }

    /**
     * Returns the ref parameter.
     * @return the ref parameter
     */
    public File ref() {
        return ref;
    }

    /**
     * Returns the gt parameter.
     * @return the gt parameter
     */
    public File gt() {
        return gt;
    }

    /**
     * Returns the chrom parameter.
     *
     * @return the chrom parameter
     */
    public String chrom() {
        return chrom;
    }

    /**
     * Returns the out parameter.
     * @return the out parameter
     */
    public String out() {
        return out;
    }

    /**
     * Returns {@code true} if the argument "match=ID" is specified, or
     * if of the match argument is absent, and returns {@code false} if
     * "match=POS" is specified.
     * @return {@code true} if the argument "match=ID" is specified, or
     * if of the match argument is absent, and returns {@code false} if
     * "match=POS" is specified
     */
    public boolean matchID() {
        return matchID;
    }

    /**
     * Returns the strict parameter.
     * @return the strict parameter
     */
    public boolean strict() {
        return strict;
    }

    /**
     * Returns the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified.
     *
     * @return the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified
     */
    public File excludesamples() {
        return excludesamples;
    }
}
