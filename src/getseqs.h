/*===============================================================================
getseqs.h : Read alignment of FASTA sequences
(C) 2004-2013 Jens Kleinjung and John Romein
Read the COPYING file for license information.
================================================================================*/

#ifndef GETSEQS_H
#define GETSEQS_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arg.h"
#include "safe.h"
#include "seq.h"

void read_sequences(FILE *seqfile, Mali *mali, Arg *arg);

#endif
