/*==============================================================================
putseqs.c : Routines for printing sequence alignments
(C) 2004-2013 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include <stdio.h>
#include <ctype.h>

#include "putseqs.h"

/*___________________________________________________________________________*/
/* print processed FASTA malij */
void write_malij_fasta(FILE *malijOutFile, Mali *malijproc, Arg *arg)
{
	unsigned int i;

	if (! arg->silent) fprintf(stdout, "\t\t%d sequences\n\t\t%d length\n",
						malijproc->nSeq, malijproc->length);

	for (i = 0; i < malijproc->nSeq; ++ i) {
		fprintf(malijOutFile, "%s%s\n",
			malijproc->sequence[i].name, malijproc->sequence[i].residue);
	}
}

/*___________________________________________________________________________*/
/* print headless malij */
void write_malij_headless(FILE *headlessijOutFile, Mali *malijproc, Arg *arg)
{
	unsigned int i;

	if (! arg->silent) fprintf(stdout, "\t\t%d sequences\n\t\t%d length\n",
						malijproc->nSeq, malijproc->length);

	for (i = 0; i < malijproc->nSeq; ++ i) {
		fprintf(headlessijOutFile, "%s\n",
			malijproc->sequence[i].residue);
	}
}

/*___________________________________________________________________________*/
/* print gapless mali */
void write_mali_gapless(FILE *gaplessiOutFile, Mali *maligapless, Arg *arg)
{
	unsigned int i;

	if (! arg->silent) fprintf(stdout, "\t\t%d sequences\n\t\t%d length\n",
						maligapless->nSeq, maligapless->length);

	for (i = 0; i < maligapless->nSeq; ++ i) {
		fprintf(gaplessiOutFile, "%s%s\n",
			 maligapless->sequence[i].name, maligapless->sequence[i].residue);
	}
}

/*___________________________________________________________________________*/
/* print joined FASTA mali */
void write_mali_join(FILE *maliJoinOutFile, Mali *malijoin, Arg *arg)
{
	unsigned int i;

	if (! arg->silent) fprintf(stdout, "\t\t%d sequences\n\t\t%d length\n",
						malijoin->nSeq, malijoin->length);

	for (i = 0; i < malijoin->nSeq; ++ i) {
		fprintf(maliJoinOutFile, "%s%s\n",
			malijoin->sequence[i].name, malijoin->sequence[i].residue);
	}
}

/*___________________________________________________________________________*/
/* print joined headless mali */
void write_mali_join_headless(FILE *headlessJoinOutFile, Mali *malijoin, Arg *arg)
{
	unsigned int i;

	if (! arg->silent) fprintf(stdout, "\t\t%d sequences\n\t\t%d length\n",
						malijoin->nSeq, malijoin->length);

	for (i = 0; i < malijoin->nSeq; ++ i) {
		fprintf(headlessJoinOutFile, "%s\n",
			malijoin->sequence[i].residue);
	}
}

/*___________________________________________________________________________*/
/* print annotated query sequence */
/* query types: 0 = mali, 1 = malj */
void write_query_annotated(FILE *annotatedOutFile, char *p0, int j, Arg *arg, int type, int cont)
{
    unsigned int i; /* local index */
	int c = 100; /* length of printed segments */
	unsigned int spaces = 1; /* print offset spaces until first annotation */
	int firstResidue = 0;

	/* query sequence */
    for (i = 0; i < c && p0[i] != '\0'; ++ i) {
		fputc(*(p0 + i), annotatedOutFile);
	}
	fputc('\n', annotatedOutFile);

	/* internal numbering, incremented by 'cont' if continued (joined sequence) */
	firstResidue = 0;
    for (i = 0, spaces = 1; i < c && p0[i] != '\0'; ++ i) {
		if ((j + i + cont) % 10 == 0) {
			fprintf(annotatedOutFile, "%-10d", (j + i + cont));
			spaces = 0;
		} else {
			if (spaces) {
				fputc(' ', annotatedOutFile);
			}
		}
	}
	fputc('\n', annotatedOutFile);

	/* UniProt numbering */
	if (type == 0) {
		firstResidue = arg->firstSeqResidue_i;
	} else if (type == 1) {
		firstResidue = arg->firstSeqResidue_j;
	} else {
		assert(1 == 2 && "This should not happen");
	}

    for (i = 0, spaces = 1; i < c && p0[i] != '\0'; ++ i) {
		if ((j + i + firstResidue) % 10 == 0) {
			fprintf(annotatedOutFile, "%-10d", (j + i + firstResidue));
			spaces = 0;
		} else {
			if (spaces) {
				fputc(' ', annotatedOutFile);
			}
		}
	}

	fputc('\n', annotatedOutFile);

	/* PDB numbering */
	if (type == 0) {
		firstResidue = arg->firstStrResidue_i;
	} else if (type == 1) {
		firstResidue = arg->firstStrResidue_j;
	} else {
		assert(1 == 2 && "This should not happen");
	}

    for (i = 0, spaces = 1; i < c && p0[i] != '\0'; ++ i) {
		if ((j + i + firstResidue) % 10 == 0) {
			fprintf(annotatedOutFile, "%-10d", (j + i + firstResidue));
			spaces = 0;
		} else {
			if (spaces) {
				fputc(' ', annotatedOutFile);
			}
		}
	}

	fputc('\n', annotatedOutFile);
	fputc('\n', annotatedOutFile);

	/* update 'global' index j */
	j += i;

    if ((i == c) && (p0[c] != '\0'))
        write_query_annotated(annotatedOutFile, p0 + c, j, arg, type, cont);
}

