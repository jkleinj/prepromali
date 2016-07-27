/*=============================================================================
prepromali: pre-process multiple alignment 
(C) 2013 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "prepromali.h"

/*____________________________________________________________________________*/
/** create gapless mali (stacked sequences) */
void gapless_mali(Mali *mali, Mali *maligapless)
{
	unsigned int i, j, k;

	maligapless->nSeq = mali->nSeq;
	maligapless->sequence = safe_malloc(maligapless->nSeq * sizeof(Seq));

	/*____________________________________________________________________________*/
	/* assign gapless sequence */
	for (i = 0; i < maligapless->nSeq; ++ i) {
		maligapless->sequence[i].name = safe_malloc((mali->maxNameLength + 1) * sizeof(char));
		maligapless->sequence[i].residue = safe_malloc((mali->length + 1) * sizeof(char));
		
		strcpy(maligapless->sequence[i].name, mali->sequence[i].name);

		for (j = 0, k = 0; j < mali->length; ++ j) {
			if (isalpha(mali->sequence[i].residue[j])) {
				maligapless->sequence[i].residue[k] = mali->sequence[i].residue[j];
				++ k;
			}
		}
		maligapless->sequence[i].residue[k] = '\0';
	}

	maligapless->length = strlen(maligapless->sequence[0].residue);
}

/*____________________________________________________________________________*/
/** process mali or malj, depending on the passed pointer */
void process_mali(Arg *arg, Mali *malij, Mali *malijproc)
{
	unsigned int i, j, k, n;
	int query = -1;
	int found = 0;

	/* determine query position in mali */
	for (i = 0; i < malij->nSeq; ++ i) {
		if (arg->join) {
			/* perform comparison with abbreviated query name */
			n = strcspn(malij->query, "/");
			if (strncmp(&(malij->sequence[i].name[1]), malij->query, n) == 0) {
				query = i;
				++ found;
				break;
			}
		} else {
			/* perform comparison with full query name */
			if (strncmp(&(malij->sequence[i].name[1]), malij->query, strlen(malij->query)) == 0) {
				query = i;
				++ found;
				break;
			}
		}
	}

	/* if not found, set first sequence as default query sequence */
	assert(found && "Query sequence not found!");

	if (! arg->silent) fprintf(stdout, "\tquery is input mali sequence %d\n", query);

	/*____________________________________________________________________________*/
	/* allocate memory to maliproc */
	malijproc->maxNameLength = malij->maxNameLength;
	malijproc->sequence = safe_malloc(malij->nSeq * sizeof(Seq));
	for (i = 0; i < malij->nSeq; ++ i) {
		malijproc->sequence[i].name = safe_malloc((malijproc->maxNameLength + 1) * sizeof(char));
		malijproc->sequence[i].residue = safe_malloc((malij->length + 1) * sizeof(char));
	}

	/*____________________________________________________________________________*/
	/* process names : place query in first position */
	/* see comments below in 'process sequences' */
	strcpy(malijproc->sequence[0].name, malij->sequence[query].name);
	strcpy(malijproc->sequence[query].name, malij->sequence[0].name);
	for (i = 1; i < malij->nSeq; ++ i) {
		if (i != query) {
			strcpy(malijproc->sequence[i].name, malij->sequence[i].name);
		}
	}

	/*____________________________________________________________________________*/
	/* process sequences : create a gapless maliproc with respect to the
		query sequence and place query seqence at the top */
	for (j = 0, k = 0; j < malij->length; ++ j) {
		if (isalpha(malij->sequence[query].residue[j])) {
			/* copy malij[query] sequence to malijproc[0] sequence */
			malijproc->sequence[0].residue[k] = malij->sequence[query].residue[j];
			for (i = 1; i < malij->nSeq; ++ i) {
				/* copy mali[0] sequence to maliproc[query] sequence to complete swap */
				if (i == query) {

					malijproc->sequence[i].residue[k] = malij->sequence[0].residue[j];
				/* otherwise copy malij[i] sequence to malijproc[i] sequence */
				} else {
					malijproc->sequence[i].residue[k] = malij->sequence[i].residue[j];
				}
			}
			++ k;
		}
	}

	/* terminate sequence strings */
	for (i = 0; i < malij->nSeq; ++ i) {
		malijproc->sequence[i].residue[k] = '\0';
	}

	malijproc->nSeq = malij->nSeq;
	malijproc->length = strlen(malijproc->sequence[0].residue);
	malijproc->query = malij->query;
	malijproc->noquery = malij->noquery;
}

/*____________________________________________________________________________*/
/** join malis */
/* note that malijoin.nSeq may be smaller than maliproc.nSeq,
    depending on how many headers can be matched */
void join_mali(Arg *arg, Mali *maliproc, Mali *maljproc, Mali *malijoin)
{
	unsigned int i, j, n;
	char strbuf[128] = "";

	malijoin->nSeq = 0;
	malijoin->length = maliproc->length + maljproc->length;
	malijoin->maxNameLength = maliproc->maxNameLength > maljproc->maxNameLength ? maliproc->maxNameLength : maljproc->maxNameLength;

	malijoin->maxnSeq = maliproc->nSeq > maljproc->nSeq ? maliproc->nSeq : maljproc->nSeq;
	fprintf(stderr, "==> %d (%d %d)\n", malijoin->maxnSeq, maliproc->nSeq, maljproc->nSeq);
	malijoin->sequence = safe_malloc(malijoin->maxnSeq * sizeof(Seq));
	for (i = 0; i < malijoin->maxnSeq; ++ i) {
		/* first sequence contains concatenated headers plus a few characters */
		malijoin->sequence[i].name = safe_malloc(2 * (malijoin->maxNameLength + 4) * sizeof(char));
		malijoin->sequence[i].residue = safe_malloc((malijoin->length + 1) * sizeof(char));
	}

	/* for first sequences */
	/* joined sequence header */
	sprintf(malijoin->sequence[0].name, ">%s + %s\n", maliproc->query, maljproc->query);
	/* joined sequence residues */
	strcpy(malijoin->sequence[0].residue, maliproc->sequence[0].residue);
	strcat(malijoin->sequence[0].residue, maljproc->sequence[0].residue); 
	++ malijoin->nSeq;

	/* for all (except first) sequences of N-terminal alignment */
	for (i = 1; i < maliproc->nSeq; ++ i) {
		fprintf(stderr, "===>%d %d\n", i, malijoin->maxnSeq);
		n = strcspn(maliproc->sequence[i].name, "/");
		/* for all (except first) sequences of C-terminal alignment */
		for (j = 1; j < maljproc->nSeq; ++ j) {
			if (strncmp(maliproc->sequence[i].name, maljproc->sequence[j].name, n) == 0) {
				strncpy(&(strbuf[0]), maliproc->sequence[i].name, n);
				if (! arg->silent)
					fprintf(stdout, "\tjoining sequences %d:%d with common header '%s'\n", i, j, &(strbuf[1]));
				/* joined sequence header */
				/*sprintf(malijoin->sequence[malijoin->nSeq].name, "%s\n", &(strbuf[0]));*/
				/* joined sequence residues */
				/*
				strcpy(malijoin->sequence[malijoin->nSeq].residue, maliproc->sequence[i].residue);
				strcat(malijoin->sequence[malijoin->nSeq].residue, maljproc->sequence[j].residue); 
				*/
				++  malijoin->nSeq;
				break;
			}
		}
	}
}

/*___________________________________________________________________________*/
int main(int argc, char *argv[])
{
	unsigned int i;

	Arg arg; /* command line arguments */
    Mali mali; /* multiple alignment */
    Mali maliproc; /* processed multiple alignment */
    Mali maligapless; /* gapless stacked sequences */
    Mali malj; /* multiple alignment to join */
    Mali maljproc; /* processed multiple alignment to join */
    Mali malijoin; /* joined multiple alignment */

    /*____________________________________________________________________________*/
    /** parse command line arguments */
	parse_args(argc, &(argv[0]), &arg);
	mali.query = arg.query_i;
	mali.noquery = arg.noquery_i;
	malj.query = arg.query_j;
	malj.noquery = arg.noquery_j;

    /*____________________________________________________________________________*/
    /** read input mali */
	if (! arg.silent) fprintf(stdout, "\nInput multiple sequence alignment\n");
	arg.maliInFile = safe_open(arg.maliInFileName, "r");
	read_sequences(arg.maliInFile, &mali, &arg);
	fclose(arg.maliInFile);

    /*____________________________________________________________________________*/
    /** process mali : gapless alignment */
	if (! arg.silent) fprintf(stdout, "\nGapless multiple sequence alignment\n");
	gapless_mali(&mali, &maligapless);

    /*____________________________________________________________________________*/
    /** process mali : gapless alignment with respect to query and
		query sequence in top position */
	if (! arg.silent) fprintf(stdout, "\nProcessing multiple sequence alignment\n");
	process_mali(&arg, &mali, &maliproc);

    /*____________________________________________________________________________*/
    /** join mali : join second alignment to first if requested */
	if (arg.join) {
		if (! arg.silent) fprintf(stdout, "\nInput multiple sequence alignment to join\n");
		arg.maljFile = safe_open(arg.maljFileName, "r");
		read_sequences(arg.maljFile, &malj, &arg);
		fclose(arg.maljFile);
		if (! arg.silent) fprintf(stdout, "\nProcessing multiple sequence alignment to join\n");
		process_mali(&arg, &malj, &maljproc);

		join_mali(&arg, &maliproc, &maljproc, &malijoin);
	}

    /*____________________________________________________________________________*/
    /** write processed mali */
	if (! arg.silent) fprintf(stdout, "\nOutput multiple sequence alignments\n");

	/* UNsorted mali in gapless format */
	if (! arg.silent) fprintf(stdout, "\tprinting unsorted gapless alignment i to '%s'\n", arg.gaplessiOutFileName);
	arg.gaplessiOutFile = safe_open(arg.gaplessiOutFileName, "w");
	write_mali_gapless(arg.gaplessiOutFile, &maligapless, &arg);
	fclose(arg.gaplessiOutFile);

	/* sorted/processed mali in FASTA format */
	if (! arg.silent) fprintf(stdout, "\tprinting sorted/processed FASTA alignment i to '%s'\n", arg.maliOutFileName);
	arg.maliOutFile = safe_open(arg.maliOutFileName, "w");
	write_malij_fasta(arg.maliOutFile, &maliproc, &arg);
	fclose(arg.maliOutFile);

	/* sorted/processed mali in headless format */
	if (! arg.silent) fprintf(stdout, "\tprinting sorted/processed headless alignment i to '%s'\n", arg.headlessiOutFileName);
	arg.headlessiOutFile = safe_open(arg.headlessiOutFileName, "w");
	write_malij_headless(arg.headlessiOutFile, &maliproc, &arg);
	fclose(arg.headlessiOutFile);


	if (arg.join) {
		/* sorted/processed malj in FASTA format */
		if (! arg.silent) fprintf(stdout, "\tprinting sorted/processed FASTA alignment j to '%s'\n", arg.maljOutFileName);
		arg.maljOutFile = safe_open(arg.maljOutFileName, "w");
		write_malij_fasta(arg.maljOutFile, &maljproc, &arg);
		fclose(arg.maljOutFile);

		/* processed malj in headless format */
		if (! arg.silent) fprintf(stdout, "\tprinting sorted/processed headless alignment j to '%s'\n", arg.headlessjOutFileName);
		arg.headlessjOutFile = safe_open(arg.headlessjOutFileName, "w");
		write_malij_headless(arg.headlessjOutFile, &maljproc, &arg);
		fclose(arg.headlessjOutFile);

		/* joined malij in FASTA format */
		if (! arg.silent) fprintf(stdout, "\tprinting joined FASTA alignment to '%s'\n", arg.maliJoinOutFileName);
		arg.maliJoinOutFile = safe_open(arg.maliJoinOutFileName, "w");
		write_mali_join(arg.maliJoinOutFile, &malijoin, &arg);
		fclose(arg.maliJoinOutFile);

		/* joined malij in headless format */
		if (! arg.silent) fprintf(stdout, "\tprinting joined headless alignment to '%s'\n", arg.headlessJoinOutFileName);
		arg.headlessJoinOutFile = safe_open(arg.headlessJoinOutFileName, "w");
		write_mali_join_headless(arg.headlessJoinOutFile, &malijoin, &arg);
		fclose(arg.headlessJoinOutFile);
	}

    /*____________________________________________________________________________*/
    /** write residue number-annotated query sequence */
	if (! arg.silent) fprintf(stdout, "\nOutput residue number-annotated query sequence(s)\n");

	/* query of mali in annotated format */
	if (! arg.silent) fprintf(stdout, "\tprinting annotated query sequence of mali to '%s'\n", arg.annotatediOutFileName);
	arg.annotatediOutFile = safe_open(arg.annotatediOutFileName, "w");
	fprintf(arg.annotatediOutFile, "1. Sequence\n2. Internal numbering\n3. UniProt numbering\n4. PDB numbering\n\n");
	write_query_annotated(arg.annotatediOutFile, &(maliproc.sequence[0].residue[0]), 0, &arg, 0, 0);
	fclose(arg.annotatediOutFile);

	if (arg.join) {
		/* query of malj in annotated format */
		if (! arg.silent) fprintf(stdout, "\tprinting annotated query sequence of malj to '%s'\n", arg.annotatedjOutFileName);
		arg.annotatedjOutFile = safe_open(arg.annotatedjOutFileName, "w");
		fprintf( arg.annotatedjOutFile, "1. Sequence\n2. Internal numbering\n3. UniProt numbering\n4. PDB numbering\n\n");
		write_query_annotated(arg.annotatedjOutFile, &(maljproc.sequence[0].residue[0]), 0, &arg, 1, 0);
		fclose(arg.annotatedjOutFile);

		/* query of joined malij in annotated format */
		if (! arg.silent) fprintf(stdout, "\tprinting annotated query sequence of malij to '%s'\n", arg.annotatedJoinOutFileName);
		arg.annotatedJoinOutFile = safe_open(arg.annotatedJoinOutFileName, "w");
		fprintf( arg.annotatedJoinOutFile, "1. Sequence\n2. Internal numbering\n3. UniProt numbering\n4. PDB numbering\n\n");
		write_query_annotated(arg.annotatedJoinOutFile, &(maliproc.sequence[0].residue[0]), 0, &arg, 0, 0);
		write_query_annotated(arg.annotatedJoinOutFile, &(maljproc.sequence[0].residue[0]), 0, &arg, 1, maliproc.length);
		fclose(arg.annotatedJoinOutFile);
	}

    /*____________________________________________________________________________*/
	/* free memory */
	for (i = 0; i < mali.nSeq; ++ i) {
		free(mali.sequence[i].name);
		free(maligapless.sequence[i].name);
		free(maliproc.sequence[i].name);
		free(mali.sequence[i].residue);
		free(maligapless.sequence[i].residue);
		free(maliproc.sequence[i].residue);
	}
	free(mali.sequence);
	free(maligapless.sequence);
	free(maliproc.sequence);

	if (arg.join) {
		for (i = 0; i < malj.nSeq; ++ i) {
			free(malj.sequence[i].name);
			free(maljproc.sequence[i].name);
			free(malj.sequence[i].residue);
			free(maljproc.sequence[i].residue);
		}
		free(malj.sequence);
		free(maljproc.sequence);

		for (i = 0; i < maliproc.nSeq; ++ i) {
			free(malijoin.sequence[i].name);
			free(malijoin.sequence[i].residue);
		}
		free(malijoin.sequence);
		
	}

    /*____________________________________________________________________________*/
    fprintf(stdout, "\nClean termination\n\n");
	return 0;
}

