/*==============================================================================
getseqs.c : Read FASTA multiple sequence alignment
(C) 2013 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "getseqs.h"

/*___________________________________________________________________________*/
/* string to upper case */
char *strupr(char *str)
{
	char *p = str;

	while ((*p = toupper(*p)))
		p ++;

	return str;
}

/*___________________________________________________________________________*/
void read_sequences(FILE *file, Mali *mali, Arg *arg)
{
	unsigned int sequenceAllocated = 64;
	unsigned int residueAllocated = 1;
	unsigned int lLine = 1024;
	char line[1024];
	int nameLength = 0;

    mali->nSeq = -1;
    mali->length = 0;
    mali->maxNameLength = 0;
    mali->sequence = safe_malloc(sequenceAllocated * sizeof(Seq));

	while(fgets(line, 1024, file) != 0) { 
		assert((strlen(line) < lLine) && "input line must be not longer than 1024 characters");
		/* sequence header */
		if (line[0] == '>') {
			++ mali->nSeq;
			/* allocate space for more sequences if needed */
			if (mali->nSeq >= sequenceAllocated) {
				sequenceAllocated += 64;
				mali->sequence = safe_realloc(mali->sequence, sequenceAllocated * sizeof(Seq));
			}

			/* store sequence name */
			nameLength = strlen(line) + 1;
			mali->sequence[mali->nSeq].name = safe_malloc(nameLength * sizeof(char));
			strcpy(mali->sequence[mali->nSeq].name, line);
			/* longest sequence name in alignment */
			mali->maxNameLength = nameLength > mali->maxNameLength ? nameLength : mali->maxNameLength;
			/* initialise aa sequence (to be continued in 'else' section) */
			residueAllocated = 1;
			mali->sequence[mali->nSeq].residue = safe_malloc(residueAllocated * sizeof(char));
			mali->sequence[mali->nSeq].residue[0] = '\0';
		} else {
			assert(mali->nSeq > -1 && "input sequence format must be FASTA");
			/* concatenate (without NEWLINE) and store aa sequence */
			residueAllocated += strlen(line);
			mali->sequence[mali->nSeq].residue = safe_realloc(mali->sequence[mali->nSeq].residue, residueAllocated);
			strncat(mali->sequence[mali->nSeq].residue, strupr(line), (strlen(line) - 1));
		}
    }

	++ mali->nSeq;
	mali->length = strlen(mali->sequence[0].residue);

    if (mali->nSeq < 2) {
		printf("Exiting: Need at least 2 sequences for multiple alignment\n");
		exit(1);
	}

	if (! arg->silent) fprintf(stdout, "\t%d sequences\n\t%d length\n",
						mali->nSeq, mali->length);
}

