/*===============================================================================
seq.h : sequence data structure
(C) 2008-2013 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
================================================================================*/

#ifndef SEQ_H
#define SEQ_H

/*___________________________________________________________________________*/
/** data structures */

/* sequence */
typedef struct  
{
    char *name; /* sequence name */
    char *residue; /* array of residues = sequence */
    int length; /* length of sequence */
} Seq;

/* multiple sequence alignment */
typedef struct  
{
	Seq *sequence; /* mali constituting sequences */
	int nSeq; /* number of sequences */
	int maxnSeq; /* largest number of sequnces */
    int length; /* length of sequence */
	int maxNameLength; /* longest sequence name */
	char *query; /* query (header) of this mali */
	int noquery;
} Mali;

#endif

