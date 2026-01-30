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

import java.util.Arrays;

/**
 * <p>Class {@code ChromIds} is a singleton class that represents a
 * list of chromosome identifiers.
 * </p>
 * The singleton instance of {@code ChromIds} is thread-safe.
 *
 * @author Brian L. Browning
 */
public final class ChromIds {

    private static final ChromIds chromIds = new ChromIds();

    private final ThreadSafeIndexer<String> instance;

    private ChromIds() {
        // private constructor to restrict instantiation.
        int initCapacity = 4;
        this.instance = new ThreadSafeIndexer<>(initCapacity);
    }

    /**
     * Returns the singleton {@code ChromIds} instance.
     * @return the singleton {@code ChromIds} instance
     */
    public static ChromIds instance() {
        return chromIds;
    }

    /**
     * Returns the index of the specified chromosome identifier.  If
     * the chromosome identifiers is not yet indexed, the chromosome identifier
     * will be indexed. Chromosome identifier indices are assigned in
     * consecutive order beginning with 0.
     * @param id a chromosome identifier
     * @return the index of the specified chromosome identifier
     * @throws IllegalArgumentException if {@code id.isEmpty()}
     * @throws NullPointerException if {@code id == null}
     */
    public int getIndex(String id) {
        if (id.isEmpty()) {
            throw new IllegalArgumentException("id.isEmpty()");
        }
        return instance.getIndex(id);
    }

    /**
     * Returns the index of the specified chromosome identifier, or returns
     * {@code -1} if the specified chromosome identifier is not indexed.
     *
     * @param id a chromosome identifier.
     * @return the index of the specified chromosome identifier, or
     * {@code -1} if the specified chromosome identifier is not indexed.
     *
     * @throws IllegalArgumentException if {@code id.isEmpty()}
     * @throws NullPointerException if {@code id == null}
     */
    public int getIndexIfIndexed(String id) {
        if (id.isEmpty()) {
            throw new IllegalArgumentException("id.isEmpty()");
        }
        return instance.getIndexIfIndexed(id);
    }

    /**
     * Returns the number of indexed chromosomes identifiers.
     * @return the number of indexed chromosomes identifiers
     */
    public int size() {
        return instance.size();
    }

    /**
     * Returns the chromosome identifier with the specified index.
     * @param index a chromosome identifier index.
     * @return the specified chromosome identifier.
     * @throws IndexOutOfBoundsException if
     * {@code  index < 0 || index >= this.size()}
     */
    public String id(int index) {
        return instance.item(index);
    }

    /**
     * Returns the list of chromosome identifiers as an array.
     * The returned array will have length {@code this.size()}, and
     * it will satisfy {@code this.ids()[k].equals(this.id(k)) == true}
     * for {@code  0 <= k < this.size()}.
     *
     * @return an array of chromosome identifiers
     */
    public String[] ids() {
        return instance.items().toArray(new String[0]);
    }

    /**
     * Returns  {@code java.util.Arrays.toString(this.ids())}.
     *
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        return Arrays.toString(this.ids());
    }
}
