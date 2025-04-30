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

#ifndef POSITION_3N_TABLE_H
#define POSITION_3N_TABLE_H

#include "alignment_3n_table.h"
#include <cassert>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

using namespace std;

/**
 * basic class to store reference position information
 */
class Position {
  public:
    short chromosomeId;
    long long int location; // 1-based position
    char strand;            // +(REF) or -(REF-RC)
    unsigned short convertedCount = 0;
    unsigned short unconvertedCount = 0;
    bool empty = true;

    void initialize() {
        // chromosome.clear();
        chromosomeId = -1;
        location = -1;
        strand = '?';
        convertedCount = 0;
        unconvertedCount = 0;
        empty = true;
    }

    Position() { initialize(); };

    /**
     * return true if there is mapping information in this reference position.
     */
    inline bool isEmpty() { return empty; }

    /**
     * set the chromosome, location (position), and strand information.
     */

    inline void set(int InChromosomeId, long long int inputLoc) {
        chromosomeId = InChromosomeId;
        location = inputLoc + 1;
    }

    inline void set(char inputStrand) { strand = inputStrand; }

    /**
     * append the SAM information into this position.
     */
    void appendBase(PosQuality &input, Alignment &a) {
        empty = false;
        if (input.converted) {
            convertedCount++;
        } else {
            unconvertedCount++;
        }
    }
};

/**
 * store all reference position in this class.
 */
class Positions {
  public:
    Position refPositions[loadingBlockSize*2+123];

    string chromosome; // current reference chromosome name.'
    int curChromosomeId;
    int refPosStartPtr, refPosEndPtr;
    long long int
        location; // current location (position) in reference chromosome.
    long long int refCoveredPosition; // this is the last position in reference
                                      // chromosome we loaded in refPositions.
    ifstream refFile;
    ChromosomeFilePositions
        chromosomePos; // store the chromosome name and it's streamPos. To
                       // quickly find new chromosome in file.

    Alignment tmpAlignment;

    Positions(string inputRefFileName) {
        refFile.open(inputRefFileName, ios_base::in);
        LoadChromosomeNamesPos();
        refPosStartPtr = refPosEndPtr = location = refCoveredPosition = 0;
        chromosome = "";
    }

    ~Positions() {
        refFile.close();
    }

    inline int Mod(int x) { return x >= 2*loadingBlockSize+67 ? x - (2*loadingBlockSize+67) : x; }

    void startOutput(bool final_ = false) {
        int start_id = refPosStartPtr;
        int end_id = final_ ? refPosEndPtr : Mod(refPosStartPtr + loadingBlockSize);
        for (int i = start_id; i != end_id; i = Mod(i+1)) {
            Position &pos = refPositions[i];
            if (!(pos.isEmpty() || pos.strand == '?')) {
                const string &chr =
                    chromosomePos.getChromesomeString(pos.chromosomeId);
                cout << chr << '\t' << pos.location << '\t'
                        << pos.strand << '\t' << pos.convertedCount << '\t'
                        << pos.unconvertedCount << '\n';
            }
        }
    }

    /**
     * given the target Position output the corresponding position index in
     * refPositions.
     */
    int getIndex(long long int &targetPos) {
        int firstPos = refPositions[refPosStartPtr].location;
        int ret = Mod(targetPos - firstPos + refPosStartPtr);
        return ret;
    }

    /**
     * given reference line (start with '>'), extract the chromosome
     * information. this is important when there is space in chromosome name.
     * the SAM information only contain the first word.
     */
    string getChrName(string &inputLine) {
        string name;
        for (int i = 1; i < inputLine.size(); i++) {
            char c = inputLine[i];
            if (isspace(c)) {
                break;
            }
            name += c;
        }

        return name;
    }

    /**
     * Scan the reference file. Record each chromosome and its position in file.
     */
    void LoadChromosomeNamesPos() {
        string line;
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // this line is chromosome name
                chromosome = getChrName(line);
                streampos currentPos = refFile.tellg();
                chromosomePos.append(chromosome, currentPos);
            }
        }
        chromosomePos.sort();
        chromosome.clear();
    }

    /**
     * get a fasta line (not header), append the bases to positions.
     */
    inline void appendRefPosition(string &line, int &cur) {

        // check the base one by one
        int len = line.size();

#pragma unroll(60)
        for (int i = 0; i < len; i++) {
            refPositions[Mod(cur + i)].initialize();
            refPositions[Mod(cur + i)].set(curChromosomeId, location + i);
            char b = line[i];
            if (b == convertFrom) {
                refPositions[Mod(cur + i)].set('+');
            } else if (b == convertFromComplement) {
                refPositions[Mod(cur + i)].set('-');
            }
        }
        location += len;
        cur = Mod(cur + len);
    }

    /**
     * initially load reference sequence for 2 million bp
     */
    void loadNewChromosome(string targetChromosome, int &meetNext) {
        meetNext = 0;
        refFile.clear();
        // find the start position in file based on chromosome name.
        streampos startPos =
            chromosomePos.getChromosomePosInRefFile(targetChromosome);
        chromosome = targetChromosome;
        curChromosomeId = chromosomePos.findChromosome(
            targetChromosome, 0, chromosomePos.pos.size() - 1);
        refFile.seekg(startPos, ios::beg);
        refCoveredPosition = 2 * loadingBlockSize;
        refPosStartPtr = 0;

        string line;
        location = 0;
        refPosEndPtr = 0;
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // this line is chromosome name
                meetNext = 1;
                break;                // meet next chromosome, return it.
            } else {
                if (line.empty()) {
                    continue;
                }
                appendRefPosition(line, refPosEndPtr);
                if (location >= refCoveredPosition) {
                    break;
                }
            }
        }
    }

    bool flag_ = false;
    /**
     * load more Position (loadingBlockSize bp) to positions
     * if we meet next chromosome, return false. Else, return ture.
     */
    void loadMore(int &meetNext) {
        meetNext = 0;
        refCoveredPosition += loadingBlockSize;
        string line;
        refPosStartPtr = Mod(refPosStartPtr + loadingBlockSize);
        while (refFile.good()) {
            getline(refFile, line);
            if (line.front() == '>') { // meet next chromosome, return.
                meetNext = 1;
                break;
            } else {
                if (line.empty()) {
                    continue;
                }
                appendRefPosition(line, refPosEndPtr);
                if (location >= refCoveredPosition) {
                    break;
                }
            }
        }
    }

    /**
     * add position information from Alignment into ref position.
     */
    void appendPositions(Alignment &newAlignment) {
        if (!newAlignment.mapped || newAlignment.bases.empty()) {
            return;
        }
        long long int startPos = newAlignment.location; // 1-based position
        // find the first reference position in pool.
        int index = getIndex(newAlignment.location);

        for (int i = 0; i < newAlignment.sequence.size(); i++) {
            PosQuality *b = &newAlignment.bases[i];
            if (b->remove) {
                continue;
            }

            Position &pos = refPositions[Mod(index + b->refPos)];
            if (pos.location != startPos + b->refPos) {
                cerr << "Error: position mismatch. pos.location is " << pos.location << " which is refPositions+" << Mod(index + b->refPos) << ", but startPos is " << startPos << ", and b->refPos is " << b->refPos << endl;
                cerr << "newAlignment.location = " << newAlignment.location <<  ", index = " << index << ", b->refPos = " << b->refPos << endl;
                cerr << "refPositions[refPosStartPtr].location = " << refPositions[refPosStartPtr].location << ", refPosStartPtr = " << refPosStartPtr << endl;
                cerr << "refPositions[refPosEndPtr-1].location = " << refPositions[refPosEndPtr-1].location << ", refPosEndPtr = " << refPosEndPtr << endl;
                exit(-1);
            }
            // assert(pos.location == startPos + b->refPos);
            // assert(0 <= b->refPos && b->refPos <= loadingBlockSize);

            if (pos.strand == '?') {
                // this is for CG-only mode. read has a 'C' or 'G' but not 'CG'.
                continue;
            }
            pos.appendBase(newAlignment.bases[i], newAlignment);
        }
    }

    void appendSync(string line) {
        tmpAlignment.parse(line);
        appendPositions(tmpAlignment);
    }
};

#endif // POSITION_3N_TABLE_H