/*==============================================================================
arg.c : parse command line arguments
(C) 2013 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "config.h"
#include "arg.h"

/*____________________________________________________________________________*/
/** print version */
static void print_version()
{
	fprintf(stdout, "\n%s %s\n", PACKAGE, VERSION);
}

/*____________________________________________________________________________*/
/** print header */
static void print_header()
{
    fprintf(stdout, "\nprepromali : pre-process multiple alignment (for PSICOV)\n");
}

/*____________________________________________________________________________*/
/** print license */
static void print_license()
{
    fprintf(stdout, "\nCopyright (C) 2013 Jens Kleinjung\n"
			"prepromali is free software and comes with ABSOLUTELY NO WARRANTY.\n"
			"You are welcome to redistribute it under certain conditions.\n"
			"Read the COPYING file for distribution details.\n\n");
}

/*____________________________________________________________________________*/
/** print citation */
static void print_citation()
{
    fprintf(stdout, "\nNone.\n\n");
}

/*____________________________________________________________________________*/
/** set defaults */
static void set_defaults(Arg *arg)
{
	arg->maliInFileName = ""; /* multiple alignment input file */ 
    arg->query_i = ""; /* query i header */
	arg->noquery_i = 0; /* no query i file */ 
    arg->firstSeqResidue_i = 0; /* UniProt sequence of query_i starts at this residue number */
    arg->firstStrResidue_i = 0; /* PDB structure of query_i  starts at this residue number */
    arg->query_j = ""; /* query j header */
	arg->noquery_j = 0; /* no query i file */ 
    arg->firstSeqResidue_j = 0; /* UniProt sequence of query_j starts at this residue number */
    arg->firstStrResidue_j = 0; /* PDB structure of query_j  starts at this residue number */
	arg->maljFileName = ""; /* multiple alignment input file to join */ 
	arg->join = 0;
    arg->maliOutFileName = "prepromali.i.fasta"; /* mali output file */
    arg->maljOutFileName = "prepromali.j.fasta"; /* mali output file */
	arg->headlessiOutFileName = "prepromali.i.headless"; /* headless mali output file */ 
	arg->annotatediOutFileName = "prepromali.i.annotated"; /*  annotated queryi output file */ 
	arg->gaplessiOutFileName = "prepromali.i.gapless"; /* gapless N-terminal mali output file */ 
	arg->headlessjOutFileName = "prepromali.j.headless"; /* headless malj output file */ 
	arg->annotatedjOutFileName = "prepromali.j.annotated"; /*  annotated terminal queryj output file */ 
	arg->maliJoinOutFileName = "prepromali.join.fasta"; /* joined mali output file */ 
	arg->headlessJoinOutFileName = "prepromali.join.headless"; /* joined headless mali output file */ 
	arg->annotatedJoinOutFileName = "prepromali.join.annotated"; /* annotated joined malij output file */ 
	arg->silent = 0;
}

/*____________________________________________________________________________*/
/** check input */
static void check_input(Arg *arg)
{
	assert((strlen(arg->maliInFileName) > 0) && "Use '--maliIn <alignment input>' and consult README");
	if (strcmp(arg->query_i, "") == 0)
		 arg->noquery_i = 1;
	if (strcmp(arg->query_j, "") == 0)
		 arg->noquery_j = 1;
	assert(arg->noquery_i == 0 || arg->noquery_i == 1);
	assert(arg->noquery_j == 0 || arg->noquery_j == 1);
	assert(arg->join == 0 || arg->join == 1);
	assert((strlen(arg->maliOutFileName) > 0) && "Use '--maliOut <FASTA output>' and consult README");
	assert((strlen(arg->headlessiOutFileName) > 0) && "Use '--headlessOut <headless output>' and consult README");
	assert((strlen(arg->gaplessiOutFileName) > 0) && "Use '--gaplessOut <gapless output>' and consult README");
	assert((strlen(arg->headlessjOutFileName) > 0) && "Use '--headlessOut <headless output>' and consult README");
	assert((strlen(arg->maliJoinOutFileName) > 0) && "Use '--maliJoinOut <mali output>' and consult README");
	assert((strlen(arg->headlessJoinOutFileName) > 0) && "Use '--headlessJoinOut <mali output>' and consult README");
}

/*____________________________________________________________________________*/
/** print command line arguments */
static void print_args(Arg *arg)
{
    time_t now;
    time(&now);

    fprintf(stdout, "\tdate: %s", ctime(&now));
    fprintf(stdout, \
                    "\tinput FASTA alignment: %s\n"
                    "\tquery_i: %s\n"
                    "\tfirstSeqResidue_i: %d\n"
                    "\tfirstStrResidue_i: %d\n"
                    "\tinput FASTA alignment to join: %s\n"
                    "\tquery_j: %s\n"
                    "\tfirstSeqResidue_j: %d\n"
                    "\tfirstStrResidue_j: %d\n"
                    "\toutput FASTA alignment: %s\n"
                    "\toutput headless alignment i: %s\n"
                    "\toutput gapless alignment i: %s\n"
                    "\toutput headless alignment j: %s\n",
        arg->maliInFileName, arg->query_i,
		arg->firstSeqResidue_i, arg->firstStrResidue_i,
		arg->maljFileName, arg->query_j,
		arg->firstSeqResidue_j, arg->firstStrResidue_j,
		arg->maliOutFileName, 
		arg->headlessiOutFileName, arg->gaplessiOutFileName,
		arg->headlessjOutFileName);

	if (arg->join)
		fprintf(stdout, "\toutput joined FASTA alignment: %s\n"
						"\toutput joined headless alignment: %s\n",
			arg->maliJoinOutFileName, arg->headlessJoinOutFileName);

    fflush(stdout);
}

/*____________________________________________________________________________*/
/** parse command line long_options */
int parse_args(int argc, char **argv, Arg *arg)
{
	int c;
	const char usage[] = "\tprepromali [--maliIn ...] [--query_i] [OPTIONS ...]\n\
	prepromali [--maliIn ...] [--query_i] [--maliJoin] [--query_j] [OPTIONS ...]\n\
	   --maliIn <alignment input>\t\t\t(mode: mandatory, type: char  , default: void)\n\
	   --query_i <query i name input>\t\t(mode: optional , type: char  , default: void)\n\
	   --noquery_i\t\t\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --firstSeqResidue_i <seq residue number>\t(mode: optional , type: int   , default: 0)\n\
	   --firstStrResidue_i <str residue number>\t(mode: optional , type: int   , default: 0)\n\
	   --maliJoin <alignment input>\t\t\t(mode: optional , type: char  , default: void)\n\
	   --query_j <query j name input>\t\t(mode: optional , type: char  , default: void)\n\
	   --noquery_j\t\t\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --firstSeqResidue_j <seq residue number>\t(mode: optional , type: int   , default: 0)\n\
	   --firstStrResidue_j <str residue number>\t(mode: optional , type: int   , default: 0)\n\
	   --maliOut <processed FASTA output>\t\t(mode: optional , type: char  , default: void)\n\
	   --headlessiOut <processed headless i output>\t(mode: optional , type: no_arg, default: off)\n\
	   --gaplessiOut <processed gapless i output>\t(mode: optional , type: no_arg, default: off)\n\
	   --maljOut <joined alignment output>\t(mode: optional , type: char  , default: void)\n\
	   --headlessjOut <processed headless j output>\t(mode: optional , type: no_arg, default: off)\n\
	   --maliJoinOut <joined alignment output>\t(mode: optional , type: char  , default: void)\n\
	   --headlessJoinOut <joined headless output>\t(mode: optional , type: char  , default: void)\n\
	   --silent\t\t\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --cite\t\t\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --version\t\t\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --help\n";

    set_defaults(arg);

    if (argc < 3) {
		print_header();
        fprintf(stderr, "%s", usage);
		print_license();
        exit(1);
    }

    /** long option definition */
    static struct option long_options[] = {
        {"maliIn", required_argument, 0, 1},
        {"query_i", required_argument, 0, 2},
        {"noquery_i", no_argument, 0, 3},
        {"firstSeqResidue_i", required_argument, 0, 4},
        {"firstStrResidue_i", required_argument, 0, 5},
        {"maliJoin", required_argument, 0, 6},
        {"query_j", required_argument, 0, 7},
        {"noquery_j", no_argument, 0, 8},
        {"firstSeqResidue_j", required_argument, 0, 9},
        {"firstStrResidue_j", required_argument, 0, 10},
        {"maliOut", required_argument, 0, 11},
        {"headlessiOut", required_argument, 0, 12},
        {"gaplessiOut", required_argument, 0, 13},
        {"maljOut", required_argument, 0, 14},
        {"headlessjOut", required_argument, 0, 15},
        {"maliJoinOut", required_argument, 0, 16},
        {"headlessJoinOut", required_argument, 0, 17},
        {"silent", no_argument, 0, 21},
        {"cite", no_argument, 0, 22},
        {"version", no_argument, 0, 23},
        {"help", no_argument, 0, 24},
        {0, 0, 0, 0}
    };

    /** assign parameters to long options */
    while ((c = getopt_long(argc, argv, "1:2:3 4:5:6:7:8:9:10:11:12:13:14:15:16:17:21 22 23 24", long_options, NULL)) != -1) {
        switch(c) {
            case 1:
                arg->maliInFileName = optarg;
                break;
            case 2:
                arg->query_i = optarg;
                break;
            case 3:
                arg->noquery_i = 1;
                break;
            case 4:
                arg->firstSeqResidue_i = atoi(optarg);
                break;
            case 5:
                arg->firstStrResidue_i = atoi(optarg);
                break;
            case 6:
                arg->maljFileName = optarg;
				++ arg->join;
                break;
            case 7:
                arg->query_j = optarg;
                break;
            case 8:
                arg->noquery_j = 1;
                break;
            case 9:
                arg->firstSeqResidue_j = atoi(optarg);
                break;
            case 10:
                arg->firstStrResidue_j = atoi(optarg);
                break;
            case 11:
                arg->maliOutFileName = optarg;
                break;
            case 12:
                arg->headlessiOutFileName = optarg;
                break;
            case 13:
                arg->gaplessiOutFileName = optarg;
                break;
            case 14:
                arg->maljOutFileName = optarg;
                break;
            case 15:
                arg->headlessjOutFileName = optarg;
                break;
            case 16:
                arg->maliJoinOutFileName = optarg;
                break;
            case 17:
                arg->headlessJoinOutFileName = optarg;
                break;
            case 21:
                arg->silent = 1;
                break;
            case 22:
                print_citation();
                exit(0);
            case 23:
				print_version();
				print_license();
                exit(0);
            case 24:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(0);
            default:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(1);
        }
    }

	check_input(arg);
    print_header();
    print_args(arg);

    return 0;
}

