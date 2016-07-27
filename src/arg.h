/*==============================================================================
arg.h : parse command line arguments
(C) 2013 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARG_H
#define ARG_H

#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*____________________________________________________________________________*/
/* structures */

/* variables for commmand line arguments */
typedef struct  
{
    FILE *maliInFile;
	char *maliInFileName;
	char *query_i;
	int noquery_i;
	int firstSeqResidue_i;
	int firstStrResidue_i;
	char *query_j;
	int noquery_j;
	int firstSeqResidue_j;
	int firstStrResidue_j;
    FILE *maljFile;
    char *maljFileName;
	int join;
    FILE *maliOutFile;
	char *maliOutFileName;
    FILE *maljOutFile;
	char *maljOutFileName;
    FILE *headlessiOutFile;
	char *headlessiOutFileName;
    FILE *gaplessiOutFile;
	char *gaplessiOutFileName;
    FILE *annotatediOutFile;
	char *annotatediOutFileName;
    FILE *headlessjOutFile;
	char *headlessjOutFileName;
    FILE *annotatedjOutFile;
	char *annotatedjOutFileName;
    FILE *maliJoinOutFile;
	char *maliJoinOutFileName;
	char *headlessJoinOutFileName;
    FILE *annotatedJoinOutFile;
	char *annotatedJoinOutFileName;
    FILE *headlessJoinOutFile;
	int silent;
} Arg;

/*____________________________________________________________________________*/
/* prototypes */
int parse_args(int argc, char **argv, Arg *arg);

#endif
