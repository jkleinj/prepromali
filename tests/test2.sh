#! /bin/sh

valgrind --track-origins=yes --leak-check=full ../src/prepromali --maliIn 'test0.fastaln' --query_i 'F0SZF4_SYNGF/215-324' --firstSeqResidue_i 12 --firstStrResidue_i 14 --maliJoin 'test1.fastaln' --query_j 'O42614_CYPCA/1-208' --firstSeqResidue_j 26 --firstStrResidue_j 28 || exit 1

