/*
 * Copyright 2020, Yun (Leo) Zhang <imzhangyun@gmail.com>
 *
 * This file is part of HISAT-3N.
 *
 * HISAT-3N is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT-3N is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT-3N.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef UTILITY_3N_TABLE_H
#define UTILITY_3N_TABLE_H

#include <algorithm>
#include <functional>
#include <ostream>
#include <iostream>
#include <queue>
#include <sstream>
#include <thread>

using namespace std;

const long long int loadingBlockSize = 60000;
const char convertFrom = 'C';
const char convertTo = 'T';
const char convertFromComplement = 'G';
const char convertToComplement = 'A';


/**
 * the simple data structure to bind quality score and position (on reference)
 * together.
 */
class PosQuality {
  public:
    int readPos; // 0-based
    int refPos;  // 0-based
    char qual;
    bool converted;
    bool remove;

    PosQuality(int &inputPos) {
        readPos = inputPos;
        refPos = inputPos;
        remove = true;
    }

    void setQual(char &inputQual, bool inputConverted) {
        qual = inputQual;
        converted = inputConverted;
        remove = false;
    }
};

/**
 * the base class for string we need to search.
 */
class string_search {
  public:
    int start;
    string s;
    int stringLen;

    void initialize() {
        start = 0;
        stringLen = 0;
        s.clear();
    }

    void loadString(string intputString) {
        s = intputString;
        stringLen = s.size();
        start = 0;
    }
};

/**
 * to store CIGAR string and search segments in it.
 */
class CIGAR : public string_search {
  public:
    bool getNextSegment(int &len, char &symbol) {
        if (start == stringLen) {
            return false;
        }
        len = 0;
        int currentIndex = start;
        while (true) {
            if (isalpha(s[currentIndex])) {
                len = stoi(s.substr(start, currentIndex - start));
                symbol = s[currentIndex];
                start = currentIndex + 1;
                return true;
            }
            currentIndex++;
        }
    }
};

/**
 * to store MD tag and search segments in it.
 */
class MD_tag : public string_search {
  public:
    bool getNextSegment(string &seg) {
        if (start >= stringLen) {
            return false;
        }
        seg.clear();
        int currentIndex = start;
        bool deletion = false;

        while (true) {
            if (currentIndex >= stringLen) {
                start = currentIndex + 1;
                return !seg.empty();
            }
            if (seg.empty() && s[currentIndex] == '0') {
                currentIndex++;
                continue;
            }
            if (isalpha(s[currentIndex])) {
                if (seg.empty()) {
                    seg = s[currentIndex];
                    start = currentIndex + 1;
                    return true;
                } else {
                    if (deletion) {
                        seg += s[currentIndex];
                        // currentIndex++;
                    } else {
                        start = currentIndex;
                        return true;
                    }
                }
            } else if (s[currentIndex] == '^') {
                if (seg.empty()) {
                    seg = s[currentIndex];
                    deletion = true;
                } else {
                    start = currentIndex;
                    return true;
                }
            } else { // number
                if (seg.empty()) {
                    seg = s[currentIndex];
                } else {
                    if (deletion || isalpha(seg.back())) {
                        start = currentIndex;
                        return true;
                    } else {
                        seg += s[currentIndex];
                    }
                }
            }
            currentIndex++;
        }
    }
};

/**
 * simple safe queue
 */
template <typename T> class UnsafeQueue {
  private:
    queue<T> queue_;

    string getReadName(string *line) {
        int startPosition = 0;
        int endPosition;

        endPosition = line->find("\t", startPosition);
        string readName =
            line->substr(startPosition, endPosition - startPosition);
        return readName;
    }

  public:
    int size() {

        int s = queue_.size();
        return s;
    }

    /**
     * return true if the queue is not empty and pop front and get value.
     * return false if the queue is empty.
     */
    bool popFront(T &value) {

        bool isEmpty = queue_.empty();
        if (!isEmpty) {
            value = queue_.front();
            queue_.pop();
        }
        return !isEmpty;
    }

    void push(T value) { queue_.push(value); }

    bool empty() {

        bool check = queue_.empty();

        return check;
    }
};

/**
 * store one chromosome and it's stream position
 */
class ChromosomeFilePosition {
  public:
    string chromosome;
    streampos linePos;
    ChromosomeFilePosition(string inputChromosome, streampos inputPos) {
        chromosome = inputChromosome;
        linePos = inputPos;
    }

    bool operator<(const ChromosomeFilePosition &in) const {
        return chromosome < in.chromosome;
    }
};

/**
 * store all chromosome and it's stream position
 */
class ChromosomeFilePositions {
  public:
    vector<ChromosomeFilePosition> pos;

    /**
     * input the chromosome name and it's streamPos, if it is not in pos, add
     * it.
     */
    void append(string &chromosome, streampos &linePos) {
        pos.push_back(ChromosomeFilePosition(chromosome, linePos));
    }

    const string &getChromesomeString(int index) {
        return pos[index].chromosome;
    }
    /**
     * make binary search on pos for target chromosome name
     */
    int findChromosome(string &targetChromosome, int start, int end) {
        if (start <= end) {
            int middle = (start + end) / 2;
            if (pos[middle].chromosome == targetChromosome) {
                return middle;
            }
            if (pos[middle].chromosome > targetChromosome) {
                return findChromosome(targetChromosome, start, middle - 1);
            }
            return findChromosome(targetChromosome, middle + 1, end);
        } else {
            // cannot find the chromosome! throw!
            cerr << "Cannot find the chromosome: " << targetChromosome
                 << " in reference file." << endl;
            throw 1;
        }
    }

    /**
     * given targetChromosome name, return its streampos
     */
    streampos getChromosomePosInRefFile(string &targetChromosome) {
        int index = findChromosome(targetChromosome, 0, pos.size() - 1);
        // assert(pos[index].chromosome == targetChromosome);
        return pos[index].linePos;
    }

    /**
     * sort the pos by chromosome name
     */
    void sort() { std::sort(pos.begin(), pos.end()); }
};

#endif // UTILITY_3N_TABLE_H
