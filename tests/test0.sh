#! /bin/sh

valgrind --track-origins=yes --leak-check=full ../src/prepromali --maliIn 'test0.fastaln' --query_i 'E0P250_9FIRM/221-330' || exit 1

