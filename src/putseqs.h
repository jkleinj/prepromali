/*===============================================================================
putseqs.h : write alignment
(C) 2006-2013 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef PUTSEQS_H
#define PUTSEQS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arg.h"
#include "seq.h"

void write_malij_fasta(FILE *malijOutFile, Mali *malijproc, Arg *arg);
void write_malij_headless(FILE *headlessijOutFile, Mali *malijproc, Arg *arg);
void write_mali_gapless(FILE *gaplessiOutFile, Mali *maligapless, Arg *arg);
void write_mali_join(FILE *maliJoinOutFile, Mali *malijoin, Arg *arg);
void write_mali_join_headless(FILE *headlessJoinOutFile, Mali *maliproc, Arg *arg);
void write_query_annotated(FILE *annotatedOutFile, char *p0, int j, Arg *arg, int type, int cont);

#endif
