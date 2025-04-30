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

#include "position_3n_table.h"
#include <getopt.h>

using namespace std;

const int inf = 1234567890;
string refFileName;
bool uniqueOnly = false;
bool multipleOnly = false;


void printHelp(const char *s) {
    printf("Usage: %s u|m <alignment file>\n", s);
    printf("example: %s u /mnt/ramdisk/rna/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa\n", s);
    exit(-1);
}

bool fileExist(string &filename) {
    ifstream file(filename);
    return file.good();
}

void parseOptions(int argc, const char **argv) {
    // ./hisat-3n-table u /mnt/ramdisk/rna/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
    if (argc != 3) printHelp(argv[0]);
    uniqueOnly = argv[1][0] == 'u';
    multipleOnly = argv[1][0] == 'm';
    if (!uniqueOnly && !multipleOnly) printHelp(argv[0]);
    refFileName = argv[2];
    if (!fileExist(refFileName))
        cerr << "reference (FASTA) file is not exist." << endl, throw(1);
}

Positions *positions;


/**
 * give a SAM line, extract the chromosome and position information.
 * return true if the SAM line is mapped. return false if SAM line is not maped.
 */
bool getSAMChromosomePos(string line, string &chr, long long int &pos) {
    int startPosition = 0;
    int endPosition = 0;
    int count = 0;

    while ((endPosition = line.find("\t", startPosition)) != string::npos) {
        if (count == 2) {
            chr = line.substr(startPosition, endPosition - startPosition);
        } else if (count == 3) {
            pos =
                stoll(line.substr(startPosition, endPosition - startPosition));
            if (chr == "*") {
                return false;
            } else {
                return true;
            }
        }
        startPosition = endPosition + 1;
        count++;
    }
    return false;
}

int hisat_3n_table() {
    Positions positions(refFileName);

    // main function, initially 2 load loadingBlockSize (2,000,000) bp of
    // reference, set reloadPos to 1 loadingBlockSize, then load SAM data. when
    // the samPos larger than the reloadPos load 1 loadingBlockSize bp of
    // reference. when the samChromosome is different to current chromosome,
    // finish all sam position and output all.
    istream *alignmentFile = &cin;

    string *line;         // temporary string to get SAM line.
    string samChromosome; // the chromosome name of current SAM line.
    long long int samPos; // the position of current SAM line.
    long long int
        reloadPos; // the position in reference that we need to reload.
    long long int lastPos = 0; // the position on last SAM line. compare lastPos
                               // with samPos to make sure the SAM is sorted.

    static char buff[1000007];
    while (true) {
        if (fgets(buff, sizeof(buff), stdin) == NULL) break;
        string line(buff);

        if (line.empty() || line.front() == '@') {
            continue;
        }
        // if the SAM line is empty or unmapped, get the next SAM line.
        if (!getSAMChromosomePos(line, samChromosome, samPos)) continue;
        // if the samChromosome is different than current positions' chromosome,
        // finish all SAM line. then load a new reference chromosome.
        if (samChromosome != positions.chromosome) {
            cerr << "chromosome changed from " << positions.chromosome << " to "
                 << samChromosome << endl << flush;
            positions.startOutput(true);

            int meetNext;
            positions.loadNewChromosome(samChromosome, meetNext);
            reloadPos = meetNext ? inf : loadingBlockSize;
            lastPos = 0;
        }
        // if the samPos is larger than reloadPos, load 1 loadingBlockSize bp in
        // from reference.
        while (samPos > reloadPos) {
            positions.startOutput();
            int meetNext;
            positions.loadMore(meetNext);
            reloadPos += meetNext ? inf : loadingBlockSize;
        }
        if (lastPos > samPos) {
            cerr << "The input alignment file is not sorted. Please use sorted "
                    "SAM "
                    "file as alignment file."
                 << endl;
            throw 1;
        }
        positions.appendSync(line);
        lastPos = samPos;
    }

    // prepare to close everything.

    // move all position to outputPool
    positions.startOutput(true);
    return 0;
}

int main(int argc, const char **argv) {
    ios::sync_with_stdio(false);
    int ret = 0;

    try {
        parseOptions(argc, argv);
        ret = hisat_3n_table();
    } catch (std::exception &e) {
        cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
        cerr << "Command: ";
        for (int i = 0; i < argc; i++)
            cerr << argv[i] << " ";
        cerr << endl;
        return 1;
    } catch (int e) {
        if (e != 0) {
            cerr << "Error: Encountered internal HISAT-3N exception (#" << e
                 << ")" << endl;
            cerr << "Command: ";
            for (int i = 0; i < argc; i++)
                cerr << argv[i] << " ";
            cerr << endl;
        }
        return e;
    }

    return ret;
}
