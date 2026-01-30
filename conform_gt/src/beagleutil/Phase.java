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
package beagleutil;

/**
 * <p>Class {@code Phase} represents the equivalence of two phased genotypes
 * for a marker or for a set of markers.  Genotype equivalence is defined
 * in terms of allele equivalence.  Two alleles are equivalent if
 * either allele is missing or if both alleles are non-missing and equal.
 * </p>
 * <p>
 * For the case of a single marker with phased (i.e. ordered) genotypes
 * ({@code a1}, {@code a2})
 * and ({@code b1}, {@code b2}), then
 * <br>
 * 1) the genotypes have IDENTICAL phase if a) alleles {@code a1}
 * and {@code b1} are equivalent, b) alleles {@code a2} and
 * {@code b2} are equivalent, and c) either alleles {@code a1}
 * and {@code b2} are not equivalent or alleles
 * {@code a2} and {@code b1} are not equivalent.
 * <br>
 * 2) the genotypes have OPPOSITE phase if a) alleles {@code a1}
 * and {@code b2} are equivalent, b) alleles {@code a2} and
 * {@code b1} are equivalent, and c) either alleles {@code a1}
 * and {@code b1} are not equivalent or alleles {@code a2} and
 * {@code b2} are not equivalent.
 * <br>
 * 3) the genotypes have UNKOWN phase if a) alleles {@code a1}
 * and {@code b1} are equivalent, b) alleles {@code a2} and
 * {@code b2} are equivalent, c) alleles {@code a1} and
 * {@code b2} are equivalent, and d) alleles {@code a2} and
 * {@code b1} are equivalent.
 * <br>
 * 4) the genotypes have INCONSISTENT phase if a) either alleles
 * {@code a1} and {@code b1} are not equivalent or alleles
 * {@code a2} and {@code b2} are not equivalent, and
 * b) either alleles {@code a1} and {@code b2} are not equivalent
 * or alleles {@code a2} and {@code b1} are not equivalent.
 * </p>
 * For the case of two sets of phased genotypes for the same markers,
 * the two sets have
 * <br>
 * 1) IDENTICAL phase if the phase is IDENTICAL for at least
 * one marker and is either IDENTICAL or UNKNOWN for all markers.
 * <br>
 * 2) OPPOSITE phase if the if the phase is OPPOSITE for at
 * lease one marker and is either OPPOSITE or UNKNOWN for all markers.
 * <br>
 * 3) UNKNOWN phase if the phase is UNKNOWN for all markers.
 * <br>
 * 4) INCONSISTENT phase if a) the phase is INCONSISTENT for at least one
 * marker or if b) the relative phase is IDENTICAL for at least one marker and
 * OPPOSITE for at least one marker.
 *
 * @author Brian L. Browning
 */
public enum Phase {
        IDENTICAL,
        OPPOSITE,
        UNKNOWN,
        INCONSISTENT
}
