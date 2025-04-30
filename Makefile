
all: hisat-3n-table

hisat-3n-table:
	g++ -O3 -flto -msse2 -funroll-loops -g3 -std=c++11 -DPOPCNT_CAPABILITY -pthread -o hisat-3n-table hisat_3n_table.cpp

clean:
	rm -f hisat-3n-table
