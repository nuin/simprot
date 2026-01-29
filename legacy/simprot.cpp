/*	Simprot  (c) Copyright 2005-2012 by the University Health Network written by Elisabeth Tillier.
	Permission is granted to copy and use this program provided no fee is charged for it
	and provided that this copyright notice is not removed.	    */

#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include <stdio.h>
#include <string.h>
#include <popt.h>
#include <math.h>
#include <float.h>
#include <sys/types.h>
//#include <windows.h> //Windows only
#include <unistd.h> //Linux only
#include <time.h>
#include <stdlib.h> //needed on Linux and not Windows
#include <errno.h>
#include <stdint.h>


/** THIS HEADER CONTAINS THE EIGEN-DECOMPOSED MATRICES **/
#include "eigen.h"
#include "random.h"

using namespace std;


// Constants.
#define INSERTION 1
#define DELETION  0

#define ALPHABET_SIZE 27

#define MAX_TREE_SIZE                  65536
#define BUFFER_SIZE                    65536
#define MAX_TOKEN_SIZE                 65536
#define MAX_SEQUENCE_NAME_LENGTH         128
#define MAX_PATH_NAME                   1024
#define MAX_SEQUENCE_LENGTH            5000 
#define INIT_MAX_GAMMA_RATE_ITERATIONS 10000
#define MAX_GAMMA_THRESHOLD            0.9999999999
/*
#define CR 13            /* Decimal code of Carriage Return char
#define LF 10            /* Decimal code of Line Feed char
#define EOF_MARKER 26    /* Decimal code of DOS end-of-file marker
#define MAX_REC_LEN 1024 /* Maximum size of input buffer
*/

// This constant represents the spaces making up gaps
#define SPACE '-'

// These constants are used to remember where we came
// from when recursively traversing the tree.
#define LEFT  'L'
#define RIGHT 'R'

/**************************** FILE NAMES **************************************/
// Input files                              File contains:
char *TreeFileName ;                        // evolutionary tree
char *RootSequenceFileName ;                 // user root sequence 
char *GapDistFileName ;                 // user root sequence 
char *CorrelFileName ;

// Output files                             File will contain:
/*@null@*/ 
static char *TrueAlignmentFileName = NULL;  // The name of the output file
static char *IndelFileName=NULL; //to outout indel sizes as created;
/*@null@*/ 
static char *FastaFileName = NULL;          // generated sequences in FASTA format
static char *PhylipFileName = NULL;         // (PN Mar05), output in Phylip format
/******************************************************************************/

/**************************** FILES *******************************************/
static FILE *TrueAlignmentFile;
static FILE *FastaFile;
static FILE *PhylipFile;
static FILE *pFile;
/******************************************************************************/

/**************************** LOOKUP TABLES ***********************************/
static char itor[21];    // Convert an {i}nteger to a  {r}esidue
static int  rtoi[26];    // Convert a  {r}esidue to an {i}nteger

static double SymbolCumulativeDensity[21];
static double matrix[20][20][4];
static double row_totals[20][4];

typedef struct correl {
    int pos;
    double value;
};

static correl Correlation[MAX_SEQUENCE_LENGTH];

static double *IndelCumulativeDensity;
static double *IndelLengthCumulativeDensity;
/******************************************************************************/

// Some values
static int max_indel_length     = 2048;
static int root_sequence_length = 50;
static double IndelLengthFreq[2048];

static int debug_mode     = 0;
static int internal_count = 0;

// Parameters for the Gamma distribution
static double alpha = 1.0;
static double beta  = 1.0;

// Gap initiation and extension penatlties
static double GapExtensionPenalty = 1;
static double GapOpenProbability;

double INDEL_FREQ   = 0.03;
int subModel        = 2;
int interleave      = 0;    // (PN Aug05)
int benner          = 0;    // (PN Aug05)
int variablegamma   = 0;    // (PN Aug05)
int bennerk         = -2;   // (PN Sep05)
double indelWeight  = 0;    // (JA Aug06)
double extra_terminal_indels = 0; // (JA Sep06)

double *eigmat;                /* eig matrix variable                  */
double **probmat;              /* prob matrix variable                 */
double freqaa[20];             /* amino acid frequencies               */

double TreeBranchScale = 1.0;  /* Scaling the branch length of the tree file */
double EvolScale       = 3.0;  /* scale factor of distance when determining the length of Indel */
double previousGamma;          /* (PN Sep05) debug test */

double VariableBranch = 0;       /* flag that sets the branch length scale multiplier to random or not */
double BranchExtinction = 0;  /*PN Dec 07, value of the probability of branch being "extinct", or have a negative
                              branch length*/
double InDelRatio = 0.5;
int IndelDistLength= 200;
/****************************************************************************
 *   				STRUCT DEFINITIONS  
 ****************************************************************************/

/**** THE STRUCT REPRESENTING THE TREE NODES ****/
typedef struct TreeNode_ {
    // The {left,right}_profile_self contain exactly the same residue list as
    // "sequence". the {left,right}_profile_self contain exactly the same
    // residue list as "sequence" in the left or right child. However:
    // left_profile_{child,self} are aligned such that the two sequences have
    // their relative gaps marked with SPACE characters.
    // right_profile_{child,self} are aligned in the same way.

    char   *left_profile_child, *left_profile_self;
    char   *right_profile_child, *right_profile_self;

    /* The name, sequence and rate vector of the node */
    char   *name;       // as defined in the tree input file
    char   *sequence;   // the sequence at this node, with no gaps
    double *rate;
    char   *extinct;
    /* The left and right children of the node */
    struct TreeNode_ *left, *right;

    /* The distance on the edge above the node */
    double distance;    // branch length to the parent node
} TreeNode;


/*************** THE STRUCTS REPRESENTING ALIGNMENTS ***************
* These were originally structs of pointers, but no provision was
* made for dynamic memory management. As a consequence, memory
* usage was out of control for large trees, as the maximum possible
* memory allocation was used to guard against buffer overflow.
* The recursive nature of the algorithm meant exponential growth
* of memory usage. String's memory management obviates this problem.
*
* Further optimizations in speed and memory usage might be made
* elsewhere in the program as well, but this is not urgent.
*******************************************************************/
struct Row {
	std::string name;
	std::string sequence;
	Row() { }
	Row(const char* theName, const char* theSequence) : name(theName), sequence(theSequence) { }
};

struct Alignment {
	std::vector<Row> r;
	explicit Alignment(TreeNode* currentNode); // was GetAlignment()
	Alignment(const Alignment& child_aln, TreeNode* current_node, char child_type_flag); // was AlignProfile()
	Alignment(const Alignment& a, const Alignment& b); // was AlignAlignments()
};

/*********************************************************************
 * Modified by RLC from:
 * static Alignment *GetAlignment(TreeNode *current_node);
 *
 * RECURSIVE CONSTRUCTOR TO OBTAIN A TRUE ALIGNMENT OF SEQUENCES AT THE
 * LEAVES OF THE TREE.
 *
 * ARGUMENTS:
 * current_node  (TreeNode *) The root of the (sub) tree for which an
 *                            alignment of the sequences is to be constructed
**********************************************************************/
Alignment::Alignment(TreeNode* current_node)
{
	// If the present node is a leaf, then we just make an alignment
	// of a single sequence (i.e. no indels).
	// Remember: if there is no left child, there will be no right
	// child, and the node is a leaf.
	if (!current_node->left && current_node->sequence != NULL) { // if the node is a leaf
//		printf("%s\n%s\n", current_node->name, current_node->sequence);
		r.resize(1, Row(current_node->name, current_node->sequence)); // This obviates the need for MakeAlignment()
	} else {
		// Here we are getting the alignment for a non-leaf node.
		// Recurse left, then align the present sequence to the alignment from the left subtree:
		const Alignment left_profile(Alignment(current_node->left), current_node, LEFT);
		// Recurse right, then align the present sequence to the alignment from the right subtree:
		const Alignment right_profile(Alignment(current_node->right), current_node, RIGHT);
		// Align the left+present and right+present alignments, and assign it to this object:
		Alignment lr(left_profile, right_profile);
		r.swap(lr.r);
	}
}

/***********************************************************************
 * Modified by RLC from:
 * static Alignment *AlignProfile(Alignment *child_aln, TreeNode *current_node, char child_type_flag);
 *
 * FUNCTION TO ALIGN A SINGLE SEQUENCE TO AN EXISTING ALIGNMENT, USING
 * A PROFILE THAT CONTAINS THE POSITIONS WHERE THE SINGLE ALIGNS WITH
 * THE FIRST SEQUENCE OF THE EXISTING ALIGNMENT.
 *
 * ARGUMENTS:
 * child_aln   : Alignment of sequences in the left or right subtree
 * current_node: (TreeNode *) Node whose sequence is to be added to 
 *		 the 'child_aln'
 * child_type_flag: (char *) Flag indicating whether it is the left or right 
 * 		    child alignment with which the current_node's sequence 
 *		    is to be aligned
 *
 * RESULTS:
 * Constructs a 'true' alignment of all sequences in the
 * (left/right) subtree and the sequence at the 'current_node'
**********************************************************************/
Alignment::Alignment(const Alignment& child_aln, TreeNode* current_node, char child_type_flag) :
	r(child_aln.r.size() + 1) // Make a new alignment
{
	// Copy the name of the sequence to be added:
	r[0].name = current_node->name;

	// Copy the names of the sequences already there:
	const std::size_t m = child_aln.r.size();
	std::size_t z;
	for (z = 0; z < m; ++z) {
		r[z + 1].name = child_aln.r[z].name;
	}

	// Point the profile pointers at the appropriate child sequences:
	char* profile[2];
	if (child_type_flag == LEFT) {
		profile[0] = current_node->left_profile_child;
		profile[1] = current_node->left_profile_self;
	} else {
		profile[0] = current_node->right_profile_child;
		profile[1] = current_node->right_profile_self;
	}

	// Iterate over positions of the new alignment, filling in each column:
	const std::size_t n = std::strlen(profile[0]);
	const std::size_t child_aln_n = child_aln.r[0].sequence.length();
	std::size_t j = 0, k = 0;
	while (j < child_aln_n && k < n) {
		// If the j-th symbol in the first row in the child alignment is equal to 
		// the k-th symbol in the present node's profile 
		if (child_aln.r[0].sequence[j] == profile[0][k]) {
			r[0].sequence += profile[1][k];
			for (z = 0; z < m; ++z) {
				r[z + 1].sequence += child_aln.r[z].sequence[j];
			}
			++j;
			++k;
		} else if (child_aln.r[0].sequence[j] == SPACE) {
			r[0].sequence += SPACE;
			for (z = 0; z < m; ++z) {
				r[z + 1].sequence += child_aln.r[z].sequence[j];
			}
			++j;
		} else {
			r[0].sequence += profile[1][k];
			for (z = 0; z < m; ++z) {
				r[z + 1].sequence += SPACE;
			}
			++k;
		}
	}
}

/********************************************************************
 * Modified by RLC from:
 * static Alignment *AlignAlignments(Alignment *a, Alignment *b);
 *
 * CONSTRUCTOR BUILDS AN ALIGNMENT FROM TWO ALIGNMENTS.
 * THE FIRST ROW OF EACH ALIGNMENT
 * IS CONSIDERED TO BE A PROFILE, AND MUST CONTAIN IDENTICAL
 * CHARACTERS, WITH SPACES AT POSITIONS THAT ARE NOT NECESSARILY
 * IDENTICAL.
 *
 * ARGUMENTS:
 * a     One of the alignments to align
 * b     The other alignment to align
**********************************************************************/
Alignment::Alignment(const Alignment& a, const Alignment& b) :
	r(a.r.size() + b.r.size() - 1) // Make a new alignment
{
	// Copy 'a' into the new alignment:
	const std::size_t a_m = a.r.size();
	std::size_t z;
	for (z = 0; z < a_m; ++z) {
		r[z].name = a.r[z].name;
	}
	
	// Copy 'b' into the new alignment:
	const std::size_t aln_m = r.size();
	for (; z < aln_m; ++z) {
		r[z].name = b.r[z - a_m + 1].name;
	}
	
	// Iterate over the positions of the new alignment, filling in each column:
	const std::size_t a_n = a.r[0].sequence.length();
	const std::size_t b_n = b.r[0].sequence.length();
	std::size_t j = 0, k = 0;
	while (j < a_n && k < b_n) {
		if ((a.r[0].sequence[j] != SPACE) && (a.r[0].sequence[j] == b.r[0].sequence[k])) {
			// present position in 'a' and 'b' is same and not a blank
			for (z = 0; z < a_m; ++z) {
				r[z].sequence += a.r[z].sequence[j];
			}
			for (; z < aln_m; ++z) {
				r[z].sequence += b.r[z - a_m + 1].sequence[k];
			}
			++j;
			++k;
		} else if (a.r[0].sequence[j] == SPACE) {
			// 'a' has a space in the profile row, so we must insert spaces
			// for the rows of 'b'.
			for (z = 0; z < a_m; ++z) {
				r[z].sequence += a.r[z].sequence[j];
			}
			for (; z < aln_m; ++z) {
				r[z].sequence += SPACE;
			}
			++j;
		} else {
			// 'a' and 'b' differ in the profile row,
			// we must select one to make the present
			// column in the new alignment, and we select 'b'.
			for (z = 0; z < a_m; ++z) {
				r[z].sequence += SPACE;
			}
			for (; z < aln_m; ++z) {
				r[z].sequence += b.r[z - a_m + 1].sequence[k];
			}
			++k;
		}
	}
}

/**** THE STRUCTS REPRESENTING SEQUENCES IN LIST FORM. USED BECAUSE IT
      IS EASIER TO DO INDELS IN A LIST THAN AN ARRAY ****/
typedef struct SequenceNode_ {

    char   c;         // Residue at this site
    double rate;      // Evolutionary rate at this site
    int    mark;      // Marker indicating the site beside an indel

    struct SequenceNode_ *prev; // Backward link
    struct SequenceNode_ *next; // Forward link

} SequenceNode;

typedef struct SequenceList_ {
    SequenceNode *head;
    SequenceNode *tail;
    int size;
} SequenceList;

void MemoryRequestFailure(char *function_name) {
  fprintf(stderr, "ERROR %s: memory request failure\n", function_name);
  exit(EXIT_FAILURE);
}


//////////////////////////////////////////////////////////////////// 
///   FUNCTIONS FOR READING AND CONSTRUCTING THE TREE   

/******************************************************************
 * FUNCTION TO GET THE NEXT TOKEN FROM THE STRING REPRESENTING THE
 * TREE.
 * 
 * ARGUMENTS:
 * tree_string: (char **) the string (being parsed) representing the tree 
 * delimiters : (char *) the string representing the set of delimiters
 * 
 * RETURN VALUE: * (char *) a character array containing the token, 
 *		 which may be a delimiter or a delimited string (e.g. a name)
 ******************************************************************/
static char *GetNextToken(char **tree_string, const char *delimiters) {
  int i, j, delimiters_length;
  char buffer[MAX_TOKEN_SIZE];
  char *token;

  /* Get the size of the set of delimiters */
  delimiters_length = (int)strlen(delimiters);

  /* Check to see if the token is a delimiter (a single character) */
  for (j = 0; j < delimiters_length; j++)
    if (**tree_string == delimiters[j]) {
      token = (char *)malloc(sizeof(char));
      if (token != NULL) { // This is for "splint"
	*token = **tree_string;
	(*tree_string)++; // Eat the input
	return token;
      }
      else MemoryRequestFailure("GetNextToken()");
    }

  /* Since the token is not a delimiter, fill the buffer with all characters
     up to, but not including, the next delimiter. */
  memset(buffer, 0, MAX_TOKEN_SIZE * sizeof(char));
  for (i = 0; **tree_string != '\0'; (*tree_string)++) {
    for (j = 0; j < delimiters_length && (**tree_string != delimiters[j]); j++);
    if (j != delimiters_length) break;
    else buffer[i++] = **tree_string;
  }

  /* Allocte some memory to return, and copy the contents of 
     the buffer into it */
  token = (char *)malloc((strlen(buffer) + 1) * sizeof(char));
  if (token != NULL) {
    strcpy(token, buffer);
  }
  else MemoryRequestFailure("GetNextToken()");

  // Pass back the token
  return token;
}

/****************************************************************
 * FUNCTION TO ALLOCATE AND INITIALIZE TREE NODES.
 *
 * RETURN VALUE:
 *(TreeNode *) this function returns a pointer to a newly allocated
 * and Null-initialized tree node.
 ****************************************************************/
static TreeNode *MakeTreeNode() {
  TreeNode *temp;

  // Allocate the space
  temp = (TreeNode *)malloc(sizeof(TreeNode));
  if (temp != NULL) {
    // Zero or NULL everything
    temp->distance = 0.0;
    temp->extinct = "N";
    temp->name = NULL;
    temp->sequence = NULL;
    temp->rate = NULL;
    temp->left = NULL;
    temp->right = NULL;
    temp->left_profile_child = NULL;
    temp->left_profile_self = NULL;
    temp->right_profile_child = NULL;
    temp->right_profile_self = NULL;
  }
  else MemoryRequestFailure("MakeTreeNode()");
  return temp;
}

/****************************************************************
 * RECURSIVE FUNCTION TO FREE A TREE.
 *
 * ARGUMENTS:
 * current_node   (TreeNode *) this is a pointer to the root of the
 *                             (sub)tree being freed
 ****************************************************************/
static void FreeTree(TreeNode *current_node) {
  // Can't free it if it ain't been got...
  if (current_node != NULL) {
    // Recursively free the left subtree
    if (current_node->left) {
      FreeTree(current_node->left);
    }
    // Recursively free the right subtree
    if (current_node->right) {
      FreeTree(current_node->right);
    }
    // Free all the allocated fields
    if (current_node->name != NULL) 
      free(current_node->name);
    if (current_node->sequence != NULL) 
      free(current_node->sequence);
    if (current_node->left_profile_child != NULL) 
      free(current_node->left_profile_child);
    if (current_node->left_profile_self != NULL) 
      free(current_node->left_profile_self);
    if (current_node->right_profile_child != NULL) 
      free(current_node->right_profile_child);
    if (current_node->right_profile_self != NULL) 
      free(current_node->right_profile_self);
    free(current_node);
  }
}

/***********************************************************************
 * RECURSIVE function that flags the nodes that are extinct in the case of
 * internal node extinction
 *
 * ARGUMENTS
 * tree_node
PN Dec 2007
 ***********************************************************************/
void FlagTree(TreeNode *current_node){

  if(current_node->left) {
      FlagTree(current_node->left);
      FlagTree(current_node->right);
  }

  if (current_node->left && strstr(current_node->name, "Neg") != NULL){
      sprintf(current_node->left->name, "%s %s", current_node->left->name, "Neg");
      sprintf(current_node->left->name, "%s %s", current_node->left->name, "Neg");
  }
}


/***********************************************************************
 * RECURSIVE PARSING FUNCTION TO BUILD THE EVOLUTIONARY TREE FROM ITS
 * STRING REPRESENTATION.
 *
 * ARGUMENTS
 * tree_string       (char **) a pointer to the character array holding
 *                             the string representation of the tree
 ***********************************************************************/
static TreeNode *RecursiveParse(char **tree_string, char *ext) {
  int label_flag = 1, distance_flag = 0;
  double randomGamma, x; //PN Dec 2007
  TreeNode * temp;
  char * buffer;

  // This will be returned as the root
  temp = MakeTreeNode();

  while ( ( buffer=GetNextToken(tree_string, "(),:;") ) ) {

    if (*buffer == '(') {
      // Parsing a descendant_list, push onto the stack
      temp->left = RecursiveParse(tree_string, temp->extinct);
      temp->right = RecursiveParse(tree_string, temp->extinct);
//      label_flag = 1;
    } else if ( (*buffer == ',') || (*buffer == ';') || (*buffer == ')') ) {
      // If it was a ',' then we just finished a subtree, if it was a
      // ';' then we just finished a tree. In both cases, we must return.
      // Popping off the stack, finished parsing a descendant_list, subtree, 
      // or tree

      // ZHUOZHI
      {
	// Special case, not covered in ":", the right node of the root,
	// the right node of the right node of the root......

	if ( temp->name == NULL ) {
	  temp->name = (char *)malloc(MAX_SEQUENCE_NAME_LENGTH * sizeof(char));
          if (temp->name != NULL) {
	          sprintf(temp->name, "%s_%2.2d", "Internal", internal_count);
                  internal_count++;
	  } else {
	    MemoryRequestFailure ( "RecursiveParse()" );
	  }
	}
      }

      // !!!!!!!!!!!!!!!!!!!!! WHAT DOES THIS ACTUALLY DO

      free(buffer);
      return temp;
    } else if (*buffer == ':') {
      // This token indicates that the next token will be a branch length.
      // So we must have already processed a leaf label, a subtree, or an
      // internal node label.
      if (label_flag != 0) {
	// We haven't just processed a label, so we must be at an internal
	// node, and there must not be a label. So we make one up.
	temp->name = (char *)malloc(MAX_SEQUENCE_NAME_LENGTH * sizeof(char));
	if (temp->name != NULL) {
	  sprintf(temp->name, "%s_%2.2d", "Internal", internal_count);
	  // Increment the count of internal nodes
	  internal_count++;
	  label_flag = 0;
	}
	else MemoryRequestFailure("RecursiveParse()");
      } else {
	;
     }
      // Set the flag to indicate that we expect a number as the next token
      distance_flag = 1;
    }
    else {
      // If we are here, the token should be a label or distance
      if (label_flag == 1) {
	// This should be a label
	temp->name = (char *)malloc(MAX_SEQUENCE_NAME_LENGTH * sizeof(char));
	if (temp->name != NULL) {
	  strcpy(temp->name,buffer);
	  label_flag = 0;
	} else MemoryRequestFailure("RecursiveParse()");
      } else {
	//PN Dec 2007
        x = rndu();
	if(VariableBranch == 0){
	    // This should be a distance
            if(x <= BranchExtinction){
                  /* read in the branch length multiplied by a scale factor (program parameter)*/
                  temp->distance = 0;
                  temp->extinct = "Y";
                  sprintf(temp->name, "%s %s", temp->name, "Neg");
            } else {
                  temp->distance = strtod(buffer, NULL) * TreeBranchScale;
            }
        } else {
               randomGamma = rndgamma(VariableBranch);
         //        printf("%f\n",randomGamma);

               if(x <= BranchExtinction){
                   temp->distance =  0;
                   temp->extinct = "Y"; //not fully implemented
                   sprintf(temp->name, "%s %s", temp->name, "Neg");
               } else {
                 temp->distance = strtod(buffer, NULL) * TreeBranchScale * randomGamma;
               }
       }
       distance_flag = 0;
      }
    }
  }

  // If we get here, something's wrong
  fprintf(stderr, "ERROR RecursiveParse(): tree format problem\n");
  exit(EXIT_FAILURE);
}

////////////////////////////////////////////////////////////////////////
///   INITIALIZATION ROUTINES
/************************************************************
 * FUNCTION TO INITIALIZE THE PROTEIN MATRICES FROM PHYLIP
 ************************************************************/
void InitProtMats() {
  int i, j, temp;

  eigmat = (double *)malloc(20 * sizeof(double));
  for (i = 0; i < 20; i++)
    if (subModel == 0) eigmat[i] = pameigmat[i];
    else if (subModel == 1) eigmat[i] = jtteigmat[i];
    else eigmat[i] = pmbeigmat[i];

  probmat = (double **)malloc(20 * sizeof(double *));
  for (i = 0; i < 20; i++)
    probmat[i] = (double *)malloc(20 * sizeof(double));

  for (i = 0; i < 20; i++)
    for (j= 0; j < 20; j++)
      if (subModel == 0) probmat[i][j] = pamprobmat[i][j];
      else if (subModel == 1) probmat[i][j] = jttprobmat[i][j];
      else probmat[i][j] = pmbprobmat[i][j]; 

  for(i = 0; i < 20; i++)	{
    temp = (int)(aaSequence[i]-'A');
    rtoi[temp] = i;
    itor[i] = aaSequence[i];
  }
}




/*********************************************************************
 * FUNCTION TO CALCULATE AMINO ACID FREQUENCIES BASED ON EIGMAT
 *********************************************************************/
void MakeProtFreqs() {
  int i, maxeig = 0;

  for (i = 0; i < 20; i++)
    if (eigmat[i] > eigmat[maxeig]) 
      maxeig = i;

  memcpy(freqaa, probmat[maxeig], 20 * sizeof(double));

  for (i = 0; i < 20; i++) 
    freqaa[i] = fabs(freqaa[i]);
}


/*********************************************************************
 * FUNCTION TO INITIALIZE THE TREE. OPENS THE TREE FILE, READS IN THE
 * STRING REPRESENTATION, PARSES IT AND CONSTRUCTS THE TREE.
 *
 * RETURN VALUE:
 * (TreeNode *) a pointer to the root of the initialized tree
 *********************************************************************/
static TreeNode *InitTree() {
  char *TreeString, *buffer, *temp;
  int i, offset = 0, buffer_len;
  FILE *TreeFile;
  TreeNode *t;
  
  if (!(TreeFile = fopen(TreeFileName, "r"))) {
    fprintf ( stderr,
		"ERROR InitTree(): cannot open tree file \"%s\"\n", 
		TreeFileName);
    exit(EXIT_FAILURE);
  }

  // Allocate the buffer, which is temporary storage for individual lines
  // in the tree file.
  buffer = (char *)malloc(BUFFER_SIZE * sizeof(char));

  // Allocate the tree string, which will contain the contents of the tree file
  // without newlines or the EOF
  TreeString = (char *)malloc(MAX_TREE_SIZE * sizeof(char));

  if (buffer != NULL && TreeString != NULL) {

    /* Grab input until the end of the file */
    while ( fscanf(TreeFile, "%s\n", buffer ) != EOF) {
      buffer_len = (int)strlen(buffer);
      for (i = 0; i < buffer_len; i++)
         TreeString[offset + i] = buffer[i];
      offset += i;
    }

    // Free the buffer
    free(buffer);
    // Parse the tree
    temp = TreeString; /* Need to do this because the RecursiveParse() 
                       modifies the location where "tree" points */
    // This is the recursive descent parser
    t = RecursiveParse(&TreeString, NULL);

    // Free the tree string before returning
    free(temp);
    {
	    // Set name for the root of the Tree, by ZHUOZHI
	    char * pe;
	    pe = (char *) malloc ( sizeof(char) * MAX_SEQUENCE_NAME_LENGTH );
	    if ( pe == NULL ) {
	        MemoryRequestFailure ( "InitTree()" );
	    } else {
	        strcpy ( pe, "InternalRoot" );
	        t->name = pe;
	    }
    }
    return t;
  }
  else MemoryRequestFailure("InitTree()");
}

/******************************************
 *
 *
 ******************************************/
static void InitTrueAlignmentFile() {
  if (TrueAlignmentFileName) TrueAlignmentFile = fopen(TrueAlignmentFileName, "w");
  else TrueAlignmentFile = stdout;
}

/**********************************************************************
 * THIS FUNCTION INITIALIZES THE FASTA OUTPUT FILE. THE INITIALIZATION
 * IS DONE HERE, AND NOT IN THE ONE FUNCTION THAT USES THAT FILE,
 * BECAUSE THAT ONE FUNCTION IS RECURSIVE
 *
 * ARGUMENTS:
 * temp_fasta_file_name  (char *) The name of the file for output of the
 *                                sequences in fasta format
 **********************************************************************/
static void InitFastaFile(char *temp_fasta_file_name) {
  if (FastaFileName != NULL) {
    if ((FastaFile = fopen(temp_fasta_file_name, "w")) == NULL) {
      fprintf(stderr, "ERROR InitFastaFile(): cannot open file \"%s\"\n", temp_fasta_file_name);
      exit(EXIT_FAILURE);
    }
  }
  else FastaFile = stdout;
}

//PN March 2005
/**********************************************************************
 * THIS FUNCTION INITIALIZES THE Phylip OUTPUT FILE. THE INITIALIZATION
 * IS DONE HERE, AND NOT IN THE ONE FUNCTION THAT USES THAT FILE,
 * BECAUSE THAT ONE FUNCTION IS RECURSIVE
 *
 * ARGUMENTS:
 * temp_phylip_file_name  (char *) The name of the file for output of the
 *                                sequences in fasta format
 **********************************************************************/
static void InitPhylipFile(char *temp_phylip_file_name) {
  if (PhylipFileName != NULL) {
    if ((PhylipFile = fopen(temp_phylip_file_name, "w")) == NULL) {
      fprintf(stderr, "ERROR InitFastaFile(): cannot open file \"%s\"\n", temp_phylip_file_name);
      exit(EXIT_FAILURE);
    }
  }
  //else PhylipFile = stdout;
}


/*****************************************************************
 * INITIALIZES A VECTOR CONTAINING THE CUMULATIVE DENSITY FOR THE
 * ALPHABET. THIS VECTOR IS USED TO CONVERT A UNIFORMLY RANDOM 
 * VALUE IN [0,1] INTO A SYMBOL THAT IS RANDOM ACCORDING TO THE 
 * BACKGROUND FREQUENCY
 ***************************************************************/
static void InitSymbolCumulativeDensity() {
  int i;
  SymbolCumulativeDensity[0] = 0.0;

  for (i = 1; i <= 19; i++)	{
    SymbolCumulativeDensity[i] = 
      SymbolCumulativeDensity[i-1] + freqaa[i-1];
  }
  SymbolCumulativeDensity[20] = 1.0;
}

/////////////////////////////
/////////////////////////////
///                       ///
///   MUTATION ROUTINES   ///
///                       ///
/////////////////////////////
/////////////////////////////


/********************************************************************
 * WRAPPER FUNCTION TO GET A RANDOM GAMMA RATE.
 *
 * ARGUMENTS:
 * a (double) The alpha parameter to the Gamma distribution 
 * b (double) The beta parameter to the Gamma distribution 
 * a controls more
 *
 * RETURN VALUE:
 * (double) A value chosen at random from the Gamma distribution with
 * parameters 'a' and 'b'.
 ********************************************************************/
static double RandomGammaRate(double a, double b) {
    if (a<0) {
        return 1.0;
    }
  return rndgamma(a);
}

/*********************************************************************
 * FUNCTION TO OBTAIN A RANDOM SYMBOL ACCORDING TO THE BACKGROUND
 * DISTRIBUTION. THE FUNCTION DOES A BINARY SEARCH ON THE VECTOR OF
 * CUMULATIVE DENSITIES FOR EACH SYMBOL
 *
 * RETURN VALUE:
 * (char) A symbol chosen at random from the background distribution
 *********************************************************************/
static char RandomSymbol() {
  int mid, low = 0, high = 20;
  double x;
  int y;

  // Get a uniformly random number in [0,1]
  x = rndu();

  // Do a binary search in the cumulative density vector
  // to find the symbol corresponding to the random number
  while (1) {
    mid = (low+high)/2;
    if (SymbolCumulativeDensity[mid] <= x)
      if (SymbolCumulativeDensity[mid+1] > x) 
         return itor[mid];
      else low = mid;
    else high = mid;
  }

  return '\0'; /* This will also never happen */
}
/// ET May 2007
 static char *ReadRootSequence() {

   char *RootString, *buffer, *temp;
  int i, offset = 0, buffer_len;
  FILE *RootSequenceFile;

  int size = 0;  
  if (!(RootSequenceFile = fopen(RootSequenceFileName, "r"))) {
    fprintf ( stderr,
		"ERROR ReadRoot(): cannot open root sequence file \"%s\"\n", 
		RootSequenceFileName);
    exit(EXIT_FAILURE);
  }
  // Allocate the buffer, which is temporary storage for individual lines
  buffer = (char *)malloc(BUFFER_SIZE * sizeof(char));
  // Allocate the root string, which will contain the contents of the tree file
  // without newlines or the EOF
  RootString = (char *)malloc(MAX_TREE_SIZE * sizeof(char));
  if (buffer != NULL && RootString != NULL) {

    /* Grab input until the end of the file */
    while ( fscanf(RootSequenceFile, "%s\n", buffer ) != EOF) {
      buffer_len = (int)strlen(buffer);
      for (i = 0; i < buffer_len; i++)
	RootString[offset + i] = buffer[i];
      offset += i;
    }

    size=(int)strlen(buffer);
    RootString[size] = '\0';
    return RootString;
  }
  else MemoryRequestFailure("RootSequence()");
}



// RLC's implementation of ReadGapDist(const char* distfile):
void ReadGapDist(const char* distfile)
{
	for (int i = 0; i < max_indel_length; ++i) {
		IndelLengthFreq[i] = 0.0;
	}
	std::ifstream inputFile(distfile);
	if (inputFile.is_open()) {
		double val;
		int count = 0;
		inputFile >> val;
		while (inputFile && count < max_indel_length) {
			IndelLengthFreq[count++] = val;
			inputFile >> val;
		}
		inputFile.close();
	} else {
		printf("Error opening %s\n", distfile);
		std::exit(1);
	}
}

void ReadCorrelatedPairs(const char* CorrelFile)
{
	std::ifstream inputFile(CorrelFile);
    string tempLine;

	if (inputFile.is_open()) {
       while (getline(inputFile,tempLine,'\n')) {
            
            stringstream ss(tempLine);
            string temp;
            int l=0;
           int i=0;
           int j=0;
           double c=0;
            while (getline(ss, temp, '\t')){
                stringstream tmpstring(temp);
                if (l==0) {  tmpstring >> i ;}
                else if (l==1){ tmpstring >> j ;}
                else if (l==2) { tmpstring >> c ; }
                l++;
            }   
           Correlation[i-1].pos=j-1;
           Correlation[j-1].pos=i-1;
           Correlation[i-1].value=c;
           Correlation[j-1].value=c;
           
		}
		inputFile.close();
	} else {
		printf("Error opening %s\n", CorrelFileName);
		std::exit(1);
	}
}




/****************************************************************
 * GET A RANDOM SEQUENCE ACCORDING TO THE BACKGROUND DISTRIBUTION.
 * THIS FUNCTION ALLOCATES MEMORY, SO REMEMBER TO FREE IT.
 *
 * ARGUMENTS:
 * size    (int) The size of the desired sequence
 * 
 * RETURN VALUE:
 * (char *) A pointer to the newly allocated character array
 * containing the sequence.
 ****************************************************************/
static char *RandomSequence(int size) {
  char *s;
  int i;

  // Allocate space for the new sequence
  s = (char *)malloc (size * sizeof(char) + 1);

  if (s != NULL) {
    // Fill the space with random symbols
    for (i = 0; i < size; i++)	{
      s[i] = RandomSymbol();
    }

    s[size] = '\0';
    return s;
  }
  else MemoryRequestFailure("RandomSequence()");
}

/****************************************************************
 * (i) original character state 
 * (k) lambda counter
 ****************************************************************/
static char GetSubstitution(char c, double t) {
 double sum, p, x, exp_eigmat[20];
 int i, j, k;

 i = rtoi[c-'A'];
 x = rndu();

  for (k = 0; k < 20; k++)
    exp_eigmat[k] = exp(t * eigmat[k]);

  sum = 0.0;
  for (j = 0; j < 20 && sum < x; j++) {
    for (p = 0.0, k = 0; k < 20; k++)
      p += probmat[k][j] * probmat[k][i] * exp_eigmat[k];
    sum += p/freqaa[j];
  }
  return itor[j-1];
}
/****************************************************************
 * (i) original character state 
 * (k) lambda counter
 ****************************************************************/
static char GetContact(char c) {
    double sum, p, x, exp_eigmat[20];
    int i, j, k;
    
    i = rtoi[c-'A'];
    x = rndu();
        
    sum = 0.0;
    for (j = 0; j < 20 && sum < x; j++) {
            sum += contactmat[i][j];
    }
    return itor[j-1];
}

/*******************************************************************
 * FUNCTION TO TURN A SEQUENCE, GIVEN IN THE FORM OF A CHARACTER
 * ARRAY, INTO A DOUBLY LINKED LIST REPRESENTATION
 *
 * ARGUMENTS:
 * sequence     (char *) The sequence to be made into a list
 * rate       (double *) The vector of Gamma rates for the sites in the sequence
 *
 * RETURN VALUE:
 * (SequenceList *) A pointer to the newly constructed sequence list
 ********************************************************************/
static SequenceList *MakeSequenceList(char *sequence, double *rate) {
    int i, list_length;
    SequenceList *l;
    SequenceNode *temp1, *temp2 = NULL;

    // Allocate space for the list
    l = (SequenceList *)malloc(sizeof(SequenceList));
    if (l != NULL) {

        // Get the length of the sequence
        list_length = (int)strlen(sequence);

        // Iterate over all positions in the sequence
        l->size = list_length;
        for (i = 0; i < list_length; i++) {

            // Allocate a node for this position, and set the node's variables
            temp1 = (SequenceNode *)malloc(sizeof(SequenceNode));
            if (temp1 != NULL) {

	            // Get a Gamma rate for this site, or use an old one if possible
	            temp1->rate = (rate) ? rate[i] : RandomGammaRate(alpha, beta);
	            temp1->next = NULL;
	            temp1->prev = NULL;
	            temp1->c = sequence[i];
                temp1->mark = 0;

	            if (temp2) {
	                temp1->prev = temp2;
	                temp2->next = temp1;
	            } else l->head = temp1;

	            temp2 = temp1; 
            }
            else MemoryRequestFailure("MakeSequenceList()");
        }
        l->tail = temp1;
        return l;
    }
    else MemoryRequestFailure("MakeSequenceList()");
}


/*****************************************************************
 * FUNCTION TO FREE SPACE USED BY A SEQUENCE LIST.
 *
 * ARGUMENTS:
 * l (SequenceList *) The list whose memory we want to free
 ***************************************************************/
static void FreeSequenceList(SequenceList *l) {
  SequenceNode *temp;

  // Iterate over all nodes, free each
  while (l->head->next) {
    temp = l->head->next;
    l->head->next = l->head->next->next;
    free(temp);
  }

  // Free the head node
  free(l->head);

  // Free the list
  free(l);
}

/***********************************************************************
 * THIS FUNCTION TAKES A SEQUENCE IN THE FORM OF A DOUBLY LINKED LIST,
 * A NODE REPRESENTING A SITE IN THE SEQUENCE, AND AN INTEGER
 * REPRESENTING THE SIZE OF THE INSERTION. THE FUNCTION INSERTS A
 * RANDOM SEQUENCE STARTING IMMEDIATELY BEFORE THE NODE PASSED AS
 * PARAMETER. THIS IS SO THAT AN INSERTION CAN HAPPEN AS A PREFIX OF
 * THE SEQUENCE. THE CASE OF AN INSERTION HAPPENING AS A SUFFIX IS
 * HANDLED BY PASSING A NULL NODE.
 *
 * ARGUMENTS:
 * l    (SequenceList *)  The list into which a segment will be inserted
 * a    (SequenceNode *)  The node prior to which the segment will be inserted
 * size            (int)  The size of the insertion
 *
 * RETURN VALUE:
 * (SequenceNode *) The first node of the inserted segment
 ************************************************************************/
static SequenceNode *Insertion(SequenceList *l, SequenceNode *a, int size) {
  char *insertion_sequence;
  SequenceList *temp_list = NULL;
  SequenceNode *temp_node = NULL;

  // Creates a string representing a sequence, with random amino acids chosen
  // based on the background distribution.
  insertion_sequence = RandomSequence(size);

  // Transform this string (char*) into a doubly-linked list. Pass NULL for the
  // rates, which will result in randomly chosen rates according to the gamma
  // distribution.
  temp_list = MakeSequenceList(insertion_sequence, NULL);

  // We no longer need the char* version of the insertion sequence, as the
  // doubly-linked list version will be used instead.
  free(insertion_sequence);

  // Insert the list representation of the insertion into the list
  // representation of the original sequence, just prior to the node
  temp_node = temp_list->head;
  if (a)
    if (a->prev) {
      a->prev->next = temp_list->head;
      a->prev->next->prev = a->prev;   // used to be a->prev->next->prev = a->prev->next;
      a->prev = temp_list->tail;
      a->prev->next = a;
    }
    else {
      l->head = temp_list->head;
      a->prev = temp_list->tail;
      a->prev->next = a;
    }
  else {
    l->tail->next = temp_list->head;
    l->tail->next->prev = l->tail;   // used to be l->tail->next->prev =  l->tail->next;
    l->tail = temp_list->tail;
  }

  // Free the list (not its nodes) representing the insertion
  free(temp_list);

  // Set the new size of the sequence
  l->size += size;

  // Return the first node of the insertion
  return temp_node;
}

/***********************************************************************
 * THIS FUNCTION TAKES A SEQUENCE IN THE FORM OF A DOUBLY LINKED LIST,
 * A NODE REPRESENTING A SITE IN THE SEQUENCE, AND AN INTEGER
 * REPRESENTING THE SIZE OF THE DELETION.  THE FUNCTION DELETES 'size'
 * NODES, STARTING WITH 'a' AND CONTINUING TOWARD THE TAIL OF THE
 * LIST.
 *
 * ARGUMENTS:
 * 
 * l     (SequenceList *) The list representation of the sequence 
 *                        from which part will be deleted
 * a     (SequenceNode *) The first node to be deleted 
 * size             (int) The size of the deletion, number of nodes to delete
 *
 * RETURN VALUE:
 * (SequenceNode *) The node representing the first site following the deletion
 ******************************************************************************/
static SequenceNode *Deletion(SequenceList *l, SequenceNode *a, int size) {
  SequenceNode *b, *temp;
  int i, flag = 0, extinct =1;

  // Set a flag if the first 
  if (l->head == a) flag = 1;
  else b = a->prev;

  // Iterate over the sequence, deleting nodes.  Stop when enough
  // sites have been deleted, or there are no more to delete
  for (i = 0; i < size && a; i++) {
    temp = a;
    a = a->next;
    if (a) a->prev = temp->prev;
    free(temp);
  }

  // Reset the size of the list
  l->size -= i;
  if (l->size <= 0) {
//    fprintf(stderr, "ERROR Deletion(): Extinction!\n");
//    exit(EXIT_FAILURE);
  extinct=1;
  return NULL;
  }

  // If the deleted segment was a suffix of the list, set the 'tail' pointer
  if (a && !a->next) l->tail = a;

  // If the deleted segment was a pre re are no more nodes in the list
  if (!a) {
    l->tail = b;
  }

  // If the flag is set, then a prefix was deleted, and the list head
  // should be set.
  if (flag) l->head = a;

  // Otherwise, connect the list prior to the deletion with the list
  // following the deletion
  else b->next = a;

  // Return the node just after the deletion
  return a;
}


/***********************************************************************
 * TAKES A TREE NODE AND A SEQUENCE IN THE FORM OF A LIST, AND WRITES
 * THE SEQUENCE REPRESENTED BY THE LIST INTO THE ARRAY THAT HOLDS THE
 * NODE'S SEQUENCE. THIS FUNCTION DOES NOT FREE THE LIST, YOU MUST DO
 * THAT WHEREVER THE FUNCTION IS CALLED.
 *
 * ARGUMENTS:
 * n       (TreeNode *) The node to which a sequence is to be assigned
 * l   (SequenceList *) The sequence to assign to the node
 ***********************************************************************/
static inline void AssignSequence(TreeNode *n, SequenceList *l) {
  SequenceNode *temp; // An iterator for traversing l
  int i;              // A  counter  for traversing l

  // If there's already a sequence there, free it
  if (n->sequence) free(n->sequence);

  // Allocate some space for the sequence and rates
  n->sequence = (char *)malloc((l->size + 1) * sizeof(char));
  n->rate = (double *)malloc(l->size * sizeof(double));

  // Iterate over the list, assigning the characters
  // at the nodes to successive positions in the sequence
  for (temp = l->head, i = 0; i < l->size; i++) {
    n->sequence[i] = temp->c;    // c is the residue at that node
    n->rate[i]     = temp->rate;
    temp           = temp->next;
  }
  n->sequence[i] = '\0';
}


/**********************************************************************
 * FUNCTION TO RETURN THE MUTATION TYPE, WHICH IS A VALUE OF 0 OR 1 TO
 * INDICATE 'deletion' or 'insertion', RESPECTIVELY.
 *
 * RETURN VALUE:
 * (int) A value of 0 or 1 to indicate the type of indel
 
 **********************************************************************/
static inline int GetIndelType() { //ET Jan 2009 allows for different insertion and deletion probs
  double x;
  x=rndu();
  if (x<InDelRatio){
  return INSERTION;
  }
  else {
  return DELETION;
  }
  
//  (int)(rndu() * 1000) % 2; 
  
}

/**********************************************************************
 * FUNCTION TO RETURN THE NUMBER OF INSERTIONS OR DELETIONS GIVEN THE
 * TREE BRANCH LENGTH AND THE LENGTH OF THE SEQUENCE
 *
 * RETURN VALUE:
 * (int) A value of the number of indels
 *
 * distance: distance from the parent, ie branch length
 * length:   length of the parent sequence, in amino acids
 *
 * this is explained in the Simprot paper, in the 
 * "Number of Indels" section
 **********************************************************************/
static int GetNumIndels(double distance, int length) {
    double INDEL_RATE, f;
    int numIndel = 0, i = 0;
    //INDEL_RATE = -1 / EvolScale * log(1 - INDEL_FREQ);
    
    if (benner==3) {
        return 0;
    }
    
    
	INDEL_RATE = -log(1 - INDEL_FREQ); //PN March 2005
    //f = (1 - exp(-INDEL_RATE * distance));
	if (benner > 0)
	{
		f = 1 - exp(-INDEL_RATE * distance);
	}	
	else
	{
		f = 1 - exp((-INDEL_RATE * distance)/EvolScale);
	}

    for (i = 0; i < length; i++)	{
        if(rndu() < f)	{
            numIndel++;
        }
    }

    return numIndel;
}
 
/***********************************************************************
 * FUNCTION THAT CALCULATES GAP PENALTY FOR A GIVEN SET OF PARAMETERS
 * BASED ON THE PREVIOUSLY CALCULATED GAP EXTENSION
 *
 *
 * It uses the cumulative indel density to generate the gap extension
 * and uses the gap extension value to assign a gap penalty
 * PN March 2005
 ***********************************************************************/
double GapPenalty(double gapext)
{
	double gpenal, indel_rate;

	//indel_rate = -log(1 - INDEL_FREQ);
	gpenal = log(INDEL_FREQ/(1-exp(gapext)))+2*gapext;
	GapExtensionPenalty = gpenal;
	//fprintf(stdout, "Gap Extension Penalty: %lf - Gap Insertion Probability: %lf\n\n", gapext, gpenal);
	return gpenal;
}

 /***********************************************************************
 * FUNCTION THAT CALCULATES GAP EXTENSION FOR A GIVEN SET OF PARAMETERS
 *
 *
 * It uses the cumulative indel density to generate the gap extension
 * and uses the gap extension value to assign a gap penalty
 * PN March 2005
 ***********************************************************************/
double GapExtension()
{
	double gext;
	int i;
	double term1, term2, term3, term4, sum, * gqg;

	gqg = (double *) malloc (root_sequence_length*10 * sizeof(double));

	for (i = 1; i < root_sequence_length*10; i++) {
             term1 = 1.027 * 0.01 * exp( (i) / (-0.96));
             term2 = 3.031 * 0.001 * exp( (i) / (-3.13));
             term3 = 6.141 * 0.0001 * exp( (i) / (-14.3));
             term4 = 2.090 * 0.00001 * exp( (i) / (-81.7));
             gqg[i] = term1 + term2 + term3 + term4;
             sum += gqg[i];
        }
	gext = (gqg[3]-gqg[1])/(2*gqg[2]);
	
	GapOpenProbability = gext;
	GapPenalty(gext);
	return gext;
}

 /***********************************************************************
 * FUNCTION THAT CALCULATES GAP EXTENSION FOR A GIVEN SET OF PARAMETERS
 * in the Benner/Zipfian distribution
 *
 * It uses the cumulative indel density to generate the gap extension
 * and uses the gap extension value to assign a gap penalty
 * PN September 2005
 ***********************************************************************/
double GapExtensionBenner()
{
	double gext;
	int i;
	double term1, sum, * gqg;

	gqg = (double *) malloc (root_sequence_length*10 * sizeof(double));

	//temp_max_indel_length = (int)(.05 * root_sequence_length);
	for (i = 1; i < root_sequence_length*10; i++) {
	    term1 = pow(i, bennerk);
            gqg[i] = term1;
            sum += gqg[i];
         }

	gext = (gqg[3]-gqg[1])/(2*gqg[2]);
	GapOpenProbability = gext;
	GapPenalty(gext);
	return gext;
}

double GapExtensionOther()
{

}


 /***********************************************************************
 * FUNCTION TO INITIALIZE THE CUMULATIVE DISTRIBUTION VECTOR FOR
 * INDELS.  THE VECTOR OF THE CUMULATIVE DENSITY FOR INDELS IS USED TO
 * SELECT A RANDOM INDEL POSITION FROM THE DISTRIBUTION BY OBTAINING A
 * UNIFORMLY RANDOM NUMBER IN [0,1], AND CHECKING WHERE THAT FALLS IN
 * THE CUMULATIVE DISTRIBUTION FOR INDELS.
 *
 * Arguments:
 * head (SequenceNode *): The beginning of the sequence
 * length (int): The length of the sequence
 * type (int): The type of indel: insertion or deletion
 *
 * RETURN VALUE:
 * (int) The size of the array
 ***********************************************************************/
static int InitCumulativeIndels(SequenceNode *head, int length, int type)	{
    double sum = 0.0;
    int i, arraySize;
    SequenceNode *current;
    current = head;
    IndelCumulativeDensity[0] = 0.0;
    current = head;

    if (type == INSERTION) {
        // Insertion: position i in the vector IndelCumulativeDensity stores the
        // cumulative density that an indel begins at position i. 
        // It is possible that the insertion sequence is the prefix, and such probability
        // equals to that of inserting at position 1 of the sequence.
        sum += current->rate;
        IndelCumulativeDensity[1] = current->rate;
        arraySize = length + 2;
        i = 2;
    } else {
        // Deletion: position i in IndelCumulativeDensity stores the cumulative 
        // probability that the amino acid at spot i will be deleted.
        // Since there are only 'length' number of amino acids to choose from,
        // the arraySize is 'length + 1'
        // ? why length + 1 rather than just length ?
        arraySize = length + 1;
        i = 1;
    }

    // Each node in the sequence has a rate associated with it. Copy these rates
    // into the IndelCumulativeDensity vector, normalize by the total sum, then
    // add into each cell the sum of all the cells before it: this converts the
    // rates to a cumulative distribution.

    // Copy the rates.
    for (current = head; current && i < arraySize; current = current->next, i++)	{
        IndelCumulativeDensity[i] = current->rate;
        sum += current->rate;
    }  
    if (i != arraySize) {
        fprintf(stderr, "ERROR: Improper initialization of IndelCumulativeDensity with type = %i\n", type);
        exit(EXIT_FAILURE);
    }

    // Normalize by the sum.
    for(i = 0; i < arraySize; i++)	{
        IndelCumulativeDensity[i] /= sum;
    }

    // Add in the sum of all cells before the current one.
    for(i = 1; i < arraySize; i++)	{
        IndelCumulativeDensity[i] += IndelCumulativeDensity[i-1];
    }

    return arraySize;
}

/// Start adding a new function for indel weighting

 /***********************************************************************
 * TO CHANGE:
 * This is the function that will be modified to add indel weighting.
 * 
 * FUNCTION TO INITIALIZE THE CUMULATIVE DISTRIBUTION VECTOR FOR
 * INDELS.  THE VECTOR OF THE CUMULATIVE DENSITY FOR INDELS IS USED TO
 * SELECT A RANDOM INDEL POSITION FROM THE DISTRIBUTION BY OBTAINING A
 * UNIFORMLY RANDOM NUMBER IN [0,1], AND CHECKING WHERE THAT FALLS IN
 * THE CUMULATIVE DISTRIBUTION FOR INDELS.
 *
 * Arguments:
 * head (SequenceNode *): The beginning of the sequence
 * length (int): The length of the sequence
 * type (int): The type of indel: insertion or deletion
 *
 * RETURN VALUE:
 * (int) The size of the array
 ***********************************************************************/
static int InitCumulativeIndelsWeighted(SequenceNode *head, int length, int type) {
    double sum = 0.0;
    int i, arraySize;
    SequenceNode *current;
    current = head;
    IndelCumulativeDensity[0] = 0.0;
    current = head;

    if (type == INSERTION) {
        // Insertion: position i in the vector IndelCumulativeDensity stores the
        // cumulative density that an indel begins at position i. 
        // It is possible that the insertion sequence is the prefix, and such probability
        // equals to that of inserting at position 1 of the sequence.
        sum += current->rate;
        IndelCumulativeDensity[1] = current->rate;
        arraySize = length + 2;
        i = 2;
    } else {
        // Deletion: position i in IndelCumulativeDensity stores the cumulative 
        // probability that the amino acid at spot i will be deleted.
        // Since there are only 'length' number of amino acids to choose from,
        // the arraySize is 'length + 1'
        // ? why length + 1 rather than just length ?
        arraySize = length + 1;
        i = 1;
    }

    // Each node in the sequence has a rate associated with it. Copy these rates
    // into the IndelCumulativeDensity vector, normalize by the total sum, then
    // add into each cell the sum of all the cells before it: this converts the
    // rates to a cumulative distribution.
    // Copy the rates.
    for (current = head; current && i < arraySize; current = current->next, i++)	{
        IndelCumulativeDensity[i] = current->rate;
        sum += current->rate;
    }
    if (i != arraySize) {
        fprintf(stderr, "ERROR: Improper initialization of IndelCumulativeDensity with type = %i\n", type);
        exit(EXIT_FAILURE);
    }

    // Normalize by the sum.
    for(i = 0; i < arraySize; i++)	{
        IndelCumulativeDensity[i] /= sum;
    }

    // Add in the sum of all cells before the current one.
    for(i = 1; i < arraySize; i++)	{
        IndelCumulativeDensity[i] += IndelCumulativeDensity[i-1];
    }

    return arraySize;
}

/// End adding new function for indelWeighting


/********************************************************************
 * FUNCTION TO OBTAIN A RANDOM POSITION FOR AN INDEL. THE FUNCTION
 * SEARCHES A VECTOR CONTAINING THE CUMULATIVE DENSITY OF 
 * AVERAGE GAMMA RATES.
 *
 * ARGUMENTS:
 * length: (int) The length of the array
 *
 * RETURN VALUE:
 * (int) an integer representing the position of the amino acid
 ********************************************************************/
static int GetIndelPosition(int length) {
    int low = 0, mid, high;
    double x;

    // Get a uniformly random number in [0,1]
    x = rndu();
    high = length;

    // Do a binary search in the cumulative density vector
    // to find the indel length corresponding to the random number
    while (1) {
        mid = (low+high)/2;
        if (mid == 0) return 0;
        if (mid == length) return (length - 1);
        if (mid == low) return mid;
        if (IndelCumulativeDensity[mid] <= x) {
        if (IndelCumulativeDensity[mid+1] > x) return mid;
            low = mid;
        }
        else high = mid;
    }
    return -1; /* This will never happen... gcc likes it though... */
}


/********************************************************************
 * JA Sep 06
 * FUNCTION TO OBTAIN A RANDOM POSITION FOR AN INDEL, but only for
 * extra terminal indels. This is based on GetIndelPosition, but
 * instead of searching the cumulative distribution, it just randomly
 * chooses the first or last position.
 *
 * ARGUMENTS:
 * length: (int) The length of the array
 *
 * RETURN VALUE:
 * (int) an integer representing the position of the amino acid
 ********************************************************************/
static int GetIndelPosition_for_extra_terminal_indels(int length) {
    int low = 0, mid, high;
    double x;

    // Get a uniformly random number in [0,1]
    x = rndu();

    if (x > 0.5) { return 0; }
    else         { return length - 1; }
}

/*********************************************************************
 * FUNCTION TO obtain the gamma rate at the starting position of indel
 * ARGUMENTS:
 * head: (SequenceNode *) The beginning of the sequence
 * indelPos: (int)The starting position of indel.
 * type: (int) The type of indel: Insertion or Deletion
 *********************************************************************/
static double GetIndelGammaRate(SequenceNode *head, int indelPos, int type)	{
    int i;

    //PN from Elisabeth
    //SequenceNode *current;
    SequenceNode *current = head;
    // Remember for insertion, the actual index is 1 smaller (in order to accomodate
    // inserting at the beginning of a sequence.

    if (type == INSERTION && indelPos == 0)	{
        return current->rate;
    }
    if (type == INSERTION)	{
        indelPos--;
    }

    for(i = 0, current = head; i < indelPos && current; i++, current = current->next);
    if (!current)	{
        fprintf(stderr, "ERROR: Cannot locate the start position gamma rate for indelPos = %i\n", indelPos);
        exit(EXIT_FAILURE);    
    }
    return current->rate;  
}

/***********************************************************************
 * FUNCTION TO INITIALIZE THE CUMULATIVE DISTRIBUTION VECTOR FOR
 * INDELS LENGTH.  THE VECTOR OF THE CUMULATIVE DENSITY FOR INDELS LENGTH
 * IS USED TO SELECT A RANDOM INDEL SIZE FROM THE DISTRIBUTION BY OBTAINING A
 * UNIFORMLY RANDOM NUMBER IN [0,1], AND CHECKING WHERE THAT FALLS IN
 * THE CUMULATIVE DISTRIBUTION FOR INDELS LENGTH.  The distance factor is
 * 3.  This is a modified implementation of the Qian and Goldstein formula.
 * Argument:
 * length: (int) The length of the sequence
 * distance: (double) The length of the tree branch
 * rate; (double) The gamma rate at the beginning position of the indel
 * RETURN VALUE;
 * int : The size of the array
 ***********************************************************************/
static int InitCumulativeIndelLength(int length, double distance, double rate) {
    double term1, term2, term3, term4, normalDistance, sum = 0.0;
    int i, n, temp_max_indel_length;

    // distance = 1.0;   // depends if we want to take the branch length into consideration
    if (EvolScale == 0)	EvolScale = 3;
    normalDistance = rate * distance / EvolScale;  /* EvolScale is a program parameter */

    IndelLengthCumulativeDensity[0] = 1.0;

    // Usually a gap will have length of at most 5% of the length of sequence,
    // so set this as the maximum gap length.
    temp_max_indel_length = (int)(.05 * length);
    if (temp_max_indel_length > max_indel_length) {
        temp_max_indel_length = max_indel_length;
    }

    for (i = 1; i < temp_max_indel_length && IndelLengthCumulativeDensity[i-1] > DBL_EPSILON; i++) {
        term1 = 1.027 * pow(10, -2) * exp( (i) / (-0.96 * normalDistance) );
        term2 = 3.031 * pow(10, -3) * exp( (i) / (-3.13 * normalDistance) );
        term3 = 6.141 * pow(10, -4) * exp( (i) / (-14.3 * normalDistance) );
        term4 = 2.090 * pow(10, -5) * exp( (i) / (-81.7 * normalDistance) );
        IndelLengthCumulativeDensity[i] = term1 + term2 + term3 + term4;
        sum += IndelLengthCumulativeDensity[i];
    }
    IndelLengthCumulativeDensity[0] = 0.0;

    if (i > max_indel_length) {
        fprintf(stderr, "ERROR: cumulative distribution of indel lengths does not converge\n");
        exit(EXIT_FAILURE);
    }

    n = i;
    for (i = 1; i < n; i++)	{
        IndelLengthCumulativeDensity[i] /= sum;
    }
    for (i = 1; i < n && IndelLengthCumulativeDensity[i-1] < 1.0 - DBL_EPSILON; i++) 	{
        IndelLengthCumulativeDensity[i] += IndelLengthCumulativeDensity[i-1];
    }

    // Returns the size of the array
    return i;
}

//PN August 2005
/***********************************************************************
 * FUNCTION TO INITIALIZE THE CUMULATIVE DISTRIBUTION VECTOR FOR
 * INDELS LENGTH.  THE VECTOR OF THE CUMULATIVE DENSITY FOR INDELS LENGTH
 * IS USED TO SELECT A RANDOM INDEL SIZE FROM THE DISTRIBUTION BY OBTAINING A
 * UNIFORMLY RANDOM NUMBER IN [0,1], AND CHECKING WHERE THAT FALLS IN
 * THE CUMULATIVE DISTRIBUTION FOR INDELS LENGTH. Benner
 *
  RETURN VALUE;
 * int : The size of the array
 ***********************************************************************/
static int InitCumulativeIndelLengthBenner(int length, double distance) {
  double term1, normalDistance, sum = 0.0;
  int i, n, temp_max_indel_length;

  IndelLengthCumulativeDensity[0] = 1.0;

  // Usually a gap will have length of at most 5% of the length of sequence
  temp_max_indel_length = (int)(.05 * length);
  if (temp_max_indel_length > max_indel_length)	{
    temp_max_indel_length = max_indel_length;
  }

  for (i = 1; i < temp_max_indel_length && IndelLengthCumulativeDensity[i-1] > DBL_EPSILON; i++) {
    term1 = pow(i, bennerk); //PN Sept 2005
    IndelLengthCumulativeDensity[i] = term1;
    sum += IndelLengthCumulativeDensity[i];
  }
  IndelLengthCumulativeDensity[0] = 0.0;

  if (i > max_indel_length) {
    fprintf(stderr, "ERROR: cumulative distribution of indel lengths does not converge\n");
    exit(EXIT_FAILURE);
  }

  n = i;
  for (i = 1; i < n; i++)	{
    IndelLengthCumulativeDensity[i] /= sum;
  }
  for (i = 1; i < n && IndelLengthCumulativeDensity[i-1] < 1.0 - DBL_EPSILON; i++) {
    IndelLengthCumulativeDensity[i] += IndelLengthCumulativeDensity[i-1];
  }

  // Returns the size of the array
  return i; 
}

static int InitCumulativeIndelLengthOther(int length, double distance) {
  double term1, normalDistance, sum = 0.0;
  int i, n, temp_max_indel_length;

  IndelLengthCumulativeDensity[0] = 1.0;

  // can define up to max_indel_length
  temp_max_indel_length = IndelDistLength;

  for (i = 1; i < temp_max_indel_length; ++i) {//&& IndelLengthCumulativeDensity[i-1] > DBL_EPSILON; i++) {
    term1 = IndelLengthFreq[i-1]; //ET Jan 2009
    IndelLengthCumulativeDensity[i] = term1;
    sum += IndelLengthCumulativeDensity[i];
  }
  IndelLengthCumulativeDensity[0] = 0.0;

  if (i > max_indel_length) {
    fprintf(stderr, "ERROR: cumulative distribution of indel lengths does not converge\n");
    exit(EXIT_FAILURE);
  }

  n = i;
  for (i = 1; i < n; i++)	{
    IndelLengthCumulativeDensity[i] /= sum;
  }
  for (i = 1; i < n && IndelLengthCumulativeDensity[i-1] < 1.0 - DBL_EPSILON; i++) {
    IndelLengthCumulativeDensity[i] += IndelLengthCumulativeDensity[i-1];
  }

  // Returns the size of the array
  return i; 
}


/********************************************************************
 * FUNCTION TO OBTAIN A RANDOM LENGTH FOR AN INDEL. THE FUNCTION
 * SEARCHES A VECTOR CONTAINING THE CUMULATIVE DENSITY OF 
 * AVERAGE GAP LENGTH PROBABILITY.
 *
 * ARGUMENTS:
 * tempMaxIndellength: (int) The maximum indel length for the sequence
 *
 * RETURN VALUE:
 * (int) an integer representing the length of the indel
 ********************************************************************/
static int GetIndelLength(int tempMaxIndelLength) {
  int low = 0, mid, high;
  double x;

  // Get a uniformly random number in [0,1]
  x = rndu();
  high = tempMaxIndelLength + 1;

  // Do a binary search in the cumulative density vector
  // to find the indel length corresponding to the random number
  while (1) {
    mid = (low+high)/2;
    if (mid == low) {return (mid + 1);}
    if (IndelLengthCumulativeDensity[mid] <= x) {
      if (IndelLengthCumulativeDensity[mid+1] > x) {return (mid + 1);}
      low = mid;
    }
    else high = mid;
  }
  return -1; /* This will never happen... gcc likes it though... */
}


/*********************************************************************
 * FUNCTION TO MARK THE AMINO ACID(S) AFFECTED BY INDELS
 * ARGUMENTS:
 * head: (SequenceNode *) The beginning amino acid of the sequence
 * length: (int) The length of the sequence
 * indelPos : (int) The beginning of the indel, indexed from 0
 * indelLength: (int) The length of the indel
 * type: (int) The type of indel: Insertion or Deletion
 *********************************************************************/
static void MarkPositions(SequenceNode *head, int length, int indelPos, int indelLength, int type)	{
    int i;
    int count = 0;
    int newIndelPos;
    SequenceNode * current;

    // If the indelLength is zero, do nothing
    if (indelLength != 0)	{
        // Make sure that the position of insertion/deletion is valid
        if (indelPos < 0 || indelPos > length)	{
            fprintf(stderr, "ERROR: invalid indel position of %i\n", indelPos);
            exit(EXIT_FAILURE);
        }

        if (type == INSERTION)	{
            // Remember that the indelPos indicates the position of which the inserting seqeunce
            // is to be inserted into (ie. mark the amino acid before indelPos; if indelPos is 0,
            // mark only the first amino acid, if indelPos is length, then mark only the last amino acid)
            for (i = 0, current = head; i < indelPos - 1 && current; i++, current = current->next, count++);
            if (!current)	{
                fprintf(stderr, "ERROR: Cannot locate the start position for insertion, i = %i, indelPos = %i\n", i, indelPos);  
                exit(EXIT_FAILURE);
            }
            current->mark = 1;

            if (indelPos != 0 && current->next)	{
                current->next->mark = 1;
            }
        } else {
            // Remember that the indelPos indicates the first position of the
            // to-be-deleted amino acid.
            // Mark every amino acid starting at indelPos to indelLength number of residue or to the end of the sequence.
            // If the deletion is to reach beyond the end of the sequence, then we choose to delete the amino acids
            // before the indelPos.
            if (indelPos + indelLength > length) {
                indelPos = -indelLength + length;
            }
            if (indelPos < 0 || indelPos > length)	{
                fprintf(stderr, "ERROR: invalid indel position in backward deletion with indelPos = %i\n", indelPos);
	            exit(EXIT_FAILURE);
            }
            for (i = 0, current = head; i < indelPos && current; i++, current = current->next, count++);
            if (!current) {
                fprintf(stderr, "ERROR: Cannot locate the start position for deletion, i = %i, indelPos = %i\n", i, indelPos);
                exit(EXIT_FAILURE);
            }

            // Mark all amino acids up to indelLength or the end of the sequence
            for (i = 0; i < indelLength && current; i++, current = current->next) {
                current->mark = 1;
            }
        }
    }
}


/*********************************************************************
 * FUNCTION TO DETERMINE THE INDEL POSITION.  IT FIRST DETERMINES 
 * THE POSITION OF THE INDEL AND THEN DETERMINES ITS LENGTH.
 * IT FINALLY MARKS THE POSITION OF THE INDEL
 *
 * ARGUMENTS:
 * head: (SequenceNode *) The beginning of the sequence
 * length : (int) The length of the sequence
 * distance : (double) The branch length of the tree
 * type : (int) The type of indel: insertion or deletion
 * RETURN VALUE:
 * int : The length of the indel sequence
 *********************************************************************/
static int InitIndel(SequenceNode *head, int length, double distance, 
                                                    int type, int extra)	{
    int indelPos, indelLength, IndelLengthCumulativeDensitySize, IndelCumulativeDensitySize;
    int i = 0;
    double indelGammaRate;

    // Allocate the cumulative density vectors
    IndelCumulativeDensity = (double *)malloc((MAX_SEQUENCE_LENGTH + 1) * sizeof(double));
    IndelLengthCumulativeDensity = (double *)malloc((max_indel_length + 1) * sizeof(double));

    if (IndelCumulativeDensity != NULL && IndelLengthCumulativeDensity != NULL)	{
        IndelCumulativeDensitySize = InitCumulativeIndels(head, length, type);

        // JA Sep 06: Originally was only the indelPos = GetIndelPosition(...), added
        // the "if" and the option for extra terminal indels.
        if (extra == 0) {
            indelPos = GetIndelPosition(IndelCumulativeDensitySize - 1);
        } else {
            indelPos =
                GetIndelPosition_for_extra_terminal_indels(IndelCumulativeDensitySize - 1);
        }

        indelGammaRate = GetIndelGammaRate(head, indelPos, type);
        previousGamma = indelGammaRate; //PN September 2005

        //**********************************************************      
        //PN August 2005
        if(benner == 1) {
          IndelLengthCumulativeDensitySize = InitCumulativeIndelLengthBenner(length, distance);
        } else if(benner == 0) {
           IndelLengthCumulativeDensitySize = InitCumulativeIndelLength(length, distance, indelGammaRate);
        } else if(benner == 2) {
         IndelLengthCumulativeDensitySize = InitCumulativeIndelLengthOther(length, distance);
		}
        indelLength = GetIndelLength(IndelLengthCumulativeDensitySize - 1);

        //PN August 2005
        //IndelLengthCumulativeDensitySize = InitCumulativeIndelLength(length, distance, indelGammaRate);
        //indelLength = GetIndelLength(IndelLengthCumulativeDensitySize - 1);

        //indelLength = GetIndelLength(IndelLengthCumulativeDensitySize - 1);	
        //IndelLengthCumulativeDensitySize = InitCumulativeIndelLength(length, distance, indelGammaRate);
        //indelLength = GetIndelLength(IndelLengthCumulativeDensitySize - 1);
        //********************************************************

        MarkPositions(head, length, indelPos, indelLength, type);  

        free(IndelCumulativeDensity);
        free(IndelLengthCumulativeDensity);
        return indelLength;

    } else {
        MemoryRequestFailure("InitIndel");
  }
}

/*********************************************************************
 * FUNCTION TO CLEAR THE MARKS OF THE SEQUENCE.
 * ARGUMENTS:
 * head: (SequenceNode *) The beginning of the sequence
 * type: (int) The type of indel: Insertion or Deletion
 *********************************************************************/
static void ClearMarks(SequenceNode *head, int type)	{
  int valid = 1, count = 0;
  SequenceNode *current;
  for (current = head; current; current = current->next)	{
    if (current->mark)	{
      current->mark = 0;
      count++;
    }
  }
  if ( (type == INSERTION && count > 2 )|| (type == DELETION && count != 0) )	{
    fprintf(stderr, "ERROR: error in INSERTION == %i with count == %i", type, count);  
    exit(EXIT_FAILURE);
  }
}

/***********************************************************************
 * THIS FUNCTION RESET THE TEMP_SEQUENCE FOR THE MUTATE FUNCTION
 * SO THAT INDELS ARE CREATED BASED ON A CURRENT SEQUENCE
 *
 * Arguments:
 * l (SequenceList *) The most up-to-date SequenceList after each indel
 * temp_sequence (char *) The sequence to be updated.
 ***********************************************************************/
static char * ResetSequence(SequenceList *l, char *temp_sequence)	{
  SequenceNode *current;
  int i = 0;

  // If there's already a sequence there, free it
  if (temp_sequence) free(temp_sequence);

  // Allocate some space for the temp_sequence
  temp_sequence = (char *)malloc((l->size + 1) * sizeof(char));

  for (current = l->head, i = 0; current; current= current->next, i++)	{
    temp_sequence[i] = current->c;
  }

  temp_sequence[i] = '\0';
  return temp_sequence;
}

/***********************************************************************
 * THIS FUNCTION PERFORMS INSERTIONS AND DELETIONS. IT KEEPS TRACK OF 
 * THE NEWLY GENERATED SEQUENCES AS WELL AS THE ALIGNMENT PROFILES
 *
 * Arguments:
 * type (int) The type of indel: insertion or deletion
 *            Note that  INSERTION == 1, DELETION == 0 are #defined up top
 * l (SequenceList *) The most up-to-date SequenceList after each indel
 * temp_profile_child (char *) The previous alignment profile, which will be updated
 * temp_profile_self (char *) The previous alignment profile, which will be updated
 * temp_sequence (char *) The sequence to be updated.
 * indelSize (int) The size of the indel
 ***********************************************************************/
static void PerformIndel(int type, 
                         SequenceList *l, 
                         char *temp_profile_child, char *temp_profile_self, 
                         char *temp_sequence, 
                         int temp_length, int indelSize) {

    int i = 0, j = 0, k = 0, deleted = 0;
    int prev_length, marked = 0, theCount = 0;
    int delDummy = 0;

    //int dummy = 0; // Commented out instances of "dummy" since it isn't used
    SequenceNode *current;
    SequenceNode *temp_current;

    char *prev_profile_self;
    char *prev_profile_child;

    // Allocate space for the temporary sequences and profiles, 
    // assuming temp_profile_self and temp_profile_child have the same length
    prev_length = (int)strlen(temp_profile_self);
    prev_profile_self = (char *)malloc ((MAX_SEQUENCE_LENGTH + 1)* sizeof(char));
    prev_profile_child = (char *)malloc ((MAX_SEQUENCE_LENGTH + 1)* sizeof(char));
    if (prev_profile_child != NULL && prev_profile_self != NULL) {
        // copy parent->sequence
        strcpy(prev_profile_child, temp_profile_child);
        strcpy(prev_profile_self, temp_profile_self);
    } else MemoryRequestFailure("Indel(), previous sequences");

    if (type == INSERTION)	{

        // Set SequenceNode *current to point to the node BEFORE which the
        // insertion should be performed.
        for (   current = l->head,       i = 0,  j = 0; 
                current && !current->mark &&     j < prev_length;
                current = current->next, i++,    j++               ) {

            // Skip over residues that have been deleted previously
            while (j < prev_length && prev_profile_child[j] == SPACE)	{
                temp_profile_child[i] = prev_profile_child[j];
                temp_profile_self[i]  = prev_profile_self[j];
                i++; j++;
                //dummy = 1; // Never used.
            }

            if (prev_profile_child[j] != current->c)	{
                fprintf(stderr, "ERROR: In PerformIndel insertion, before indel site, profile aa=%c not equal to sequence aa=%c\n",
                    prev_profile_child[j], current->c);
                    exit(EXIT_FAILURE);
            }
            temp_profile_child[i] = prev_profile_child[j];
            temp_profile_self[i]  = prev_profile_self[j];
        }

        // Reached the insertion site (ie, current correctly points to the
        // SequenceNode representing the residue following the new insertion).
        // Now: Update the sequence list and the profiles.
        if (current) { 

            // The insertion should be added as a prefix to the sequence
            if (current == l->head && current->next && !current->next->mark) {
                current = Insertion(l, current, indelSize);

            // The insertion follows at least one other residue. 
            } else {

                // Skip over the residues that have been deleted previously
                // Recall that SPACE is #defined to be "-"
                while (j < prev_length && prev_profile_child[j] == SPACE) {
                    temp_profile_child[i] = prev_profile_child[j];
                    temp_profile_self[i]  = prev_profile_self[j];
                    i++; j++;
                }

                // Ensure we're at the same place in both the profile and the
                // sequence
                if (prev_profile_child[j] != current->c) {
	                fprintf(stderr, 
                            "ERROR: In PerformIndel, the aa before insertion: profile aa = %c not equal to sequence aa = %c\n",
	                        prev_profile_child[j],current->c);
                }

                // Add the residue prior to insertion
                temp_profile_child[i] = prev_profile_child[j];
                temp_profile_self[i]  = prev_profile_self[j];
                i++; j++;

                // Insert in the middle or suffix
                current = Insertion(l, current->next, indelSize);
            }

            // Update the profiles in the middle of the insertion
            // So, the "child" profile reflects the actual insertion: it
            // lists the residues included in the indel itself. The "self", on
            // the other hand, has gaps included where there is an insertion,
            // because it now has a gap relative to its child. For insertions,
            // this is the only place where a SPACE is actually added.
            for (k = 0;
                    (k < indelSize) && (current);
                    current = current->next, k++) {
                temp_profile_child[i+k] = current->c;
                temp_profile_self[i+k]  = SPACE;
            }

            i += indelSize; // Jump over the indel

            // copies over the rest of the profiles
            for (; j < prev_length; i++, j++) {
                temp_profile_child[i] = prev_profile_child[j];
                temp_profile_self[i] = prev_profile_self[j];
            }
            temp_profile_child[i] = '\0';
            temp_profile_self[i] = '\0';

            ClearMarks(l->head, type);
        }
    }// End type == INSERTION

    // Begin type == DELETION
    else {

        // Traverses to the deletion sites and updates the profiles
        for(    i = 0, j = 0, current = l->head; 
                current && !current->mark && j < prev_length; 
                current = current->next, i++, j++   ) {

            // Skip over deleted residues
            while (j < prev_length && prev_profile_child[j] == SPACE)	{
                temp_profile_child[i] = prev_profile_child[j];
                temp_profile_self[i] = prev_profile_self[j];
                i++; j++;
                //dummy = 1; // Never used.
            }
            if (prev_profile_child[j] != current->c)	{
                fprintf(stderr, "ERROR: In PerformIndel Deletion, before indel site, profile aa=%c not equal to sequence aa=%c\n",
                            prev_profile_child[j], current->c);
                exit(EXIT_FAILURE);
            }
        }

        // Reached the deletion site, now iterates through the residues and updates the profiles
        if (current) {
            temp_current = current;
            while (temp_current && temp_current->mark && j < prev_length && deleted < indelSize) {

                // Delete a residue that existed in original parent, insert space to temp_profile_child
                if (temp_current->c == prev_profile_self[j] && temp_current->c == prev_profile_child[j]) {
                    temp_profile_child[i] = SPACE;
                    temp_profile_self[i] = temp_current->c;
                    temp_current = temp_current->next;
                    i++; j++; deleted++;
                    continue;
                }

                // Delete a residue that was added in previous iterations of indel, advance in temp_profiles and temp_current
                if (prev_profile_self[j] == SPACE && temp_current->c == prev_profile_child[j])	{
	                temp_current = temp_current->next;
                    j++; deleted++;
	                continue;
	            }
// Delete residues that was deleted in previous iterations of indel, do not move temp_current and skip over in prev_profiles
                while (j < prev_length && temp_profile_child[j] == SPACE)	{
                    temp_profile_child[i] = prev_profile_child[j];
                    temp_profile_self[i] = prev_profile_self[j];
                    i++;
                    j++;
                    delDummy = 1;
                }

                if (delDummy) {
                    delDummy = 0;
                    continue;
                }

                // Should not reach here
                fprintf(stderr, "\ndie!!!!!!!!!!\n");
            }

            // Performs deletion on the sequence
            current = Deletion(l, current, indelSize);
			
            // Copies the remaining profiles
            for (; j < prev_length; i++, j++)	{
                temp_profile_child[i] = prev_profile_child[j];
                temp_profile_self[i] = prev_profile_self[j];
            }
            temp_profile_child[i] = '\0';
            temp_profile_self[i] = '\0';

            ClearMarks(l->head, type);
        }
    }// end deletion

    // Free the allocated space for previous profiles
    free(prev_profile_self);
    free(prev_profile_child);
}

/******************************************************************************
 * FUNCTION TO MUTATE A SEQUENCE BY AN AMOUNT SPECIFIED BY THE
 * DISTANCE TO THE CHILD. THE MUTATIONS ARE INSERTIONS, DELETIONS OR
 * SUBSITUTIONS. WHEN THE MUTATION IS FINISHED, THE FUNCTION CALLS
 * ITSELF RECURSIVELY ON THE NODE CORRESPONDING TO THE SEQUENCE THAT
 * RESULTED FROM THE MUTATIONS.
 *
 * Note that the comment above is note quite accurate: Mutate is not a recursive
 * function. Rather, the Evolve function is a (recursive) preorder tree
 * traversal which calls Mutate once on each node. (JA Aug06)
 *
 * ARGUMENTS:
 * parent: (TreeNode *) The tree node whose sequence is to be mutated
 * child : (TreeNode *) The tree node that will be assigned the mutated sequence
 * child_type_flag: (int) The flag indicating whether the child is a left
 *                        or right child of the parent.
 *
 * ie. _distance from the parent_ is used to _mutate the child's sequence_
 *****************************************************************************/
static void Mutate (TreeNode *parent, TreeNode *child, int child_type_flag) {

    int type, indelSize, temp_length, i = 0, j = 0, k = 0, m = 0; // (PN Sep05)
    int numIndels, valid, count;
    double sumrate = 0; // (PN Aug05)
    double r;
    
    char temp_c;
    char * temp_sequence;
    char * temp_profile_child;
    char * temp_profile_self;

    SequenceList *l;
    SequenceNode *current;

    // Allocate space for the temporary sequences and profiles
    temp_length        = (int)strlen(parent->sequence);
    temp_sequence      = (char *)malloc ((temp_length + 1)* sizeof(char));
    temp_profile_child = (char *)malloc( (int)((MAX_SEQUENCE_LENGTH + 1) * sizeof(char)) );
    temp_profile_self  = (char *)malloc( (int)((MAX_SEQUENCE_LENGTH + 1) * sizeof(char)) );
    if (temp_sequence != NULL && temp_profile_child != NULL && temp_profile_self != NULL) {
        strcpy(temp_sequence,      parent->sequence);
        strcpy(temp_profile_child, parent->sequence);
        strcpy(temp_profile_self,  parent->sequence);
    } else MemoryRequestFailure("Mutate(), temporary sequences");

    // Convert the non-gapped parent sequence (*char) 
    // to a doubly-linked list (SequenceList)
    l = MakeSequenceList(temp_sequence, parent->rate);

    // (PN Aug05) Normalize the rates at each amino acids in the sequence.
    for (current = l->head, m = 0; current && i < temp_length; current = current->next, m++) {
          sumrate += current->rate;
    }
    for (current = l->head, m = 0; current && i < temp_length; current = current->next, m++) {
         current->rate /= sumrate;
         current->rate *= temp_length;
    }
    // (PN Aug05)
    // If the branch length ("distance") is zero, then the child is the same as
    // the parent and there is nothing to do. Otherwise, continue mutating.
    // This section is the main work of the mutate function.
    if (child->distance > 0.0) {

        // Choose the number of indels to be inserted into the sequence. This is
        // based on: branch length, parent sequence length, Benner vs QG, and
        // some global/user-defined constants.
        numIndels = GetNumIndels(child->distance, temp_length);
        // Perform each indel: 
        // (1) decide on an insertion or deletion
        // (2) run InitIndel: decide on the size and position of the indel, and
        // mark the appropriate positions in the sequence.
        // (3) run PerformIndel to actually do the indel
        //     - choosing where to put the indel is done inside PerformIndel

       if (IndelFileName !=NULL) {
               fprintf (pFile, ">distance %f\n",child->distance);
       }
 

        for (count = 0; count < numIndels; count++)	{
            type = GetIndelType(); // choose "insertion" vs "deletion" randomly
				if (IndelFileName !=NULL) {
					if (type==INSERTION) {
						fprintf (pFile, "Ins %d\n",indelSize);
					}
					else {
						fprintf (pFile, "Del %d\n",indelSize);
					}
				}

            // indelSize is actually the length of the indel
            // 0 indicates that this is a normal indel, not an "extra terminal
            // indel", ie that the CDF vector should be used rather than
            // restricting the insertion point to the first or last positions.
            indelSize = InitIndel(l->head, temp_length, child->distance, type, 0);
            if (indelSize != 0)	{
                PerformIndel(type, l, temp_profile_child, temp_profile_self, temp_sequence, temp_length, indelSize); 
                // Reset temp_sequence and temp_length for the next iteration
                temp_sequence = ResetSequence(l, temp_sequence);
                temp_length = (int)strlen(temp_sequence);
            }
        }

        ////////////////////////////////////////////////////////////////////////////
        // JA Sep 06: If we want extra indels, rerun the code from the loop
        // above extra times, this time restricting the choice of indel
        // position. Note that this code is copied exactly (thought it was
        // easier to read this way than if the loop above was changed.)
        if (extra_terminal_indels > 0) {
            numIndels = extra_terminal_indels * GetNumIndels(child->distance, temp_length);
            for (count = 0; count < numIndels; count++) {
                type = GetIndelType();
                indelSize = InitIndel(l->head, temp_length, child->distance, type, 1);
                if (indelSize != 0)	{
                    PerformIndel(type, l, temp_profile_child, temp_profile_self, temp_sequence, temp_length, indelSize);

                    // Reset temp_sequence and temp_length for the next iteration
                    temp_sequence = ResetSequence(l, temp_sequence);
                    temp_length = (int)strlen(temp_sequence);
                }
            }
        }
        ////////////////////////////////////////////////////////////////////////////

        // Substitution is done on a per-residue basis, excpt if a correlation
        // corelated mutations?
               
        current = l->head;
        i = 0;
        temp_length = (int)strlen(temp_profile_child);
               
               
                      
        for (current = l->head, i = 0; current && i < temp_length; current = current->next, i++)	{
            // Skip those deleted residues
            while (temp_profile_child[i] == SPACE)	{
                i++;
            }
            if (subModel == 2)	{
            // PMB matrix
                temp_c = GetSubstitution(current->c, child->distance * current->rate);
            } else {
            // PAM and JTT matrix
                temp_c = GetSubstitution (current->c, 100.0 * child->distance * current->rate);
            }
            
            // corelated mutations?
            if (Correlation[i].pos <i) {
                r= rndu() ;
                if (r<Correlation[i].value) {
                    temp_c=temp_profile_child[Correlation[i].pos];
                    //temp_c=GetContact(temp_profile_child[Correlation[i].pos]);
                    // good idea but the contact matrix is highly variable, I don't think correlations will be maintained,
                    
                }
            
            }
            
            
            if (temp_c != current->c) {
                current->c = temp_c;
            }
            temp_profile_child[i] = current->c;
        }

        // Convert the linked-list l to a char* representing the sequence, 
        // and store it as "child->sequence" (this has no gaps)
        AssignSequence(child, l);
        FreeSequenceList(l);       //l is no longer needed; free the memory

    } else { 
        // child->distance = 0.0, so the child is the same as the parent.
        // Thus, no actual mutation is necessary: 
        //      just copy the parent node to the child node exactly.
        child->sequence = (char *)malloc(MAX_SEQUENCE_LENGTH);
        strcpy(child->sequence, parent->sequence);
        strcpy(temp_profile_child, child->sequence);
        strcpy(temp_profile_self, child->sequence);
    }

    // Copy temp_profile_self and temp_profile_child to permanent places in the
    // TreeNodes child and self.
    if (child_type_flag) { // This node is a left child of its parent.
        parent->left_profile_child = (char *)malloc(strlen(temp_profile_child)+1);
        strcpy (parent->left_profile_child, temp_profile_child);
        parent->left_profile_self = (char *)malloc(strlen(temp_profile_self)+1);
        strcpy (parent->left_profile_self, temp_profile_self);
    }
    else { // This node is a right child of its parent.
        parent->right_profile_child = (char *)malloc(strlen(temp_profile_child)+1);
        strcpy(parent->right_profile_child, temp_profile_child);
        parent->right_profile_self = (char *)malloc(strlen(temp_profile_self)+1);
        strcpy(parent->right_profile_self, temp_profile_self);
    }
    free(temp_profile_child);
    free(temp_profile_self);

} // end Mutate function

/***********************************************************************
 * THIS RECURSIVE FUNCTION EVOLVES THE ROOT SEQUENCE BY TRAVERSING THE
 * TREE. THE 'Mutate()' METHOD IS CALLED AT EACH NODE IN THE TREE TO
 * OBTAIN THE NEW SEQUENCES.
 *
 * ALIGNMENT:
 * current_node (TreeNode *) The root of the (sub)tree for which
 *                           sequences are to be evolved
 ***********************************************************************/
static void Evolve(TreeNode *current_node) {
    char * name; //PN Dec 2007

    // Make sure the current node is not NULL
    name = current_node->name;
    if (current_node->sequence) {
        // Evolve the LEFT subtree
        if (current_node->left ) {
            Mutate (current_node, current_node->left, 1);
            Evolve (current_node->left);
        }
        // Evolve the RIGHT subtree
        if (current_node->right) {
            Mutate (current_node, current_node->right, 0);
            Evolve (current_node->right);
        }
    } else {
        fprintf(stderr, "ERROR Evolve(): trying to evolve a node with no sequence.\n");
    }
 }


/***********************************************************************
 * USED TO REMOVE ALL INFORMATION STORED AT THE NODES OF THE TREE IN
 * CASE MULTIPLE SIMULATIONS OF THE SAME TREE ARE DESIRED. THIS IS AN
 * ALTERNATIVE TO REREADING THE TREE FROM THE FILE.
 *
 * ARGUMENTS:
 * current_node (TreeNode *) The root of the (sub)tree to be cleared.
 ***********************************************************************/
static void ClearTree(TreeNode *current_node) {

  if (current_node->left) {

    // Recurse left
    ClearTree(current_node->left);

    // Free and null the left sequence data
    free(current_node->left_profile_child);
    current_node->left_profile_child = NULL;
    free(current_node->left_profile_self);
    current_node->left_profile_self = NULL;
  }

  if (current_node->right) {
    // Recurse right
    ClearTree(current_node->right);
    // Free and null the right sequence data
    free(current_node->right_profile_child);
    current_node->right_profile_child = NULL;
    free(current_node->right_profile_self);
    current_node->left_profile_child = NULL;
  }

  // Free and null the sequence
  free(current_node->sequence);
  current_node->sequence = NULL;
}


/*   OUTPUT ROUTINES   */
/**************************************************************
 * THIS ONE'S JUST FOR DEBUGGING. IT PRINTS THE SEQUENCE LISTS.
 *
 * ARGUMENTS:
 * the_list    (SequenceList *) Pointer to the list to be printed
 ****************************************************************/
static void PrintSequenceList(SequenceList *the_list) {
  SequenceNode *temp;
  for (temp = the_list->head; temp; temp = temp->next)
    fprintf(stderr, "%c", temp->c);
  fprintf(stderr, "\n");
}

/**********************************************************************
 * RECURSIVE FUNCTION TO PRINT THE SEQUENCES AT EACH LEAF OF THE TREE.
 * THE SEQUENCES ARE PRINTED IN FASTA FILE, WITH THE NAMES COMING FROM
 * THE NAMES OF THE NODE AT EACH LEAF.
 *
 * ARGUMENTS:
 * current_node  (TreeNode *) The root of the (sub) tree for which the
 *                            sequences are to be printed (in .faa format)
 **********************************************************************/
static void PrintFastaFormat(TreeNode *current_node) {
  int i,count,seqLength,flag = 0;

  char * name;
  name = current_node->name;

  if (current_node->left) {
//     printf("%s\n", current_node->name);
    PrintFastaFormat(current_node->left);
    flag = 1;
  }
  if (current_node->right ) {
//     printf("%s\n", current_node->name);
    PrintFastaFormat(current_node->right);
    flag = 1;
  }

  if (!flag) {
    if(current_node->sequence != NULL && strstr(name, "Neg") == NULL){
	fprintf(FastaFile, ">%s\n", current_node->name);
	seqLength = (int)strlen(current_node->sequence);
	count = 0;
	for (i = 0; i < seqLength; i++)	{
	fprintf(FastaFile, "%c", current_node->sequence[i]);
	count++;
	    if (count == 79) {
		 fprintf(FastaFile, "\n");
		 count = 0;
	    }
        }
        fprintf(FastaFile, "\n");
     }
  }
}

/**********************************************************************
 * FUNCTION TO PRINT THE SEQUENCES AT EACH LEAF OF THE TREE.
 * THE SEQUENCES ARE PRINTED IN FASTA FILE, WITH THE NAMES COMING FROM
 * THE NAMES OF THE NODE AT EACH LEAF.
 *
 * ARGUMENTS:
 * 
 **********************************************************************/
static void PrintPhylipFormat(const Alignment& aln, char *temp_phylip_file_name, char *GapColumn) {
  int i, j, count, nodes = 0, seq_size = 0;

  if (PhylipFileName == NULL) PhylipFile = stdout;
  else if ((PhylipFile = fopen(temp_phylip_file_name, "w"))  == NULL) {
    fprintf(stderr, "ERROR PrintPhylip(): cannot open file \"%s\"\n", temp_phylip_file_name);
    exit(EXIT_FAILURE);
  }

//checking the number of terminal nodes in the tree. The aln variable has all nodes including internal
  const int aln_m = aln.r.size();
  for (i = 0; i < aln_m; i++) {
    if (strncmp(aln.r[i].name.c_str(), "Internal", 8) != 0 && strstr(aln.r[i].name.c_str(), "Neg") == NULL) {
      nodes++;
    }
   }

//checking sequence size. The variable that holds the sequence has only the aas not the indels
  const int aln_n = aln.r[0].sequence.length();
  for (j = 0; j < aln_n; j++) {
    if (j < aln_n && aln.r[0].sequence[j] != '\0' && GapColumn[j] == '0') {
      seq_size++;
    }
  }
  fprintf(PhylipFile, "%d", nodes);
  fprintf(PhylipFile, "  %d\n", seq_size);

  for (i = 0; i < aln_m; i++) {
    if (strncmp(aln.r[i].name.c_str(), "Internal", 8) != 0 && strstr(aln.r[i].name.c_str(), "Neg") == NULL) {
      count = 0;
      char temp_label[11];
      strncpy (temp_label, aln.r[i].name.c_str(), 10);
      //fprintf(PhylipFile, "%-10s", aln.r[i].name.c_str());
      fprintf(PhylipFile, "%-10s", temp_label);
      fprintf(PhylipFile, "\t");
      for (j = 0; j < aln_n; j++) {
        // Skip over the gap column
        while (j < aln_n && aln.r[i].sequence[j] != '\0' && GapColumn[j] == '1') {
          j++;
        }
        if (j < aln_n && aln.r[i].sequence[j] != '\0' && GapColumn[j] == '0') {
          fprintf(PhylipFile, "%c", aln.r[i].sequence[j]);
          count++;
        }
      }
      fprintf(PhylipFile, "\n");
    }
  }
  if (PhylipFileName) fclose(PhylipFile);
}

/***********************************************************************
 * FUNCTION TO PRINT THE TREE, INCLUDING ALL INFORMATION AVAILIBLE FOR
 * EACH NODE. THIS FUNCTION IS RECURSIVE, AND ONLY PRINTS THOSE SEQUENCES
 * WITH NAMES GIVEN IN THE TREE FILE.
 *
 * ARGUMENTS:
 * current_node  (TreeNode *) The root of the (sub) tree to be printed.
 ***********************************************************************/
static void PrintTree(TreeNode *current_node) {
  if (current_node->left) PrintTree(current_node->left);
  fprintf(stdout, "[%lld]\t%s\t%s\t%lf\t<%lld><%lld>\n", 
	 ((uint64_t)current_node>>12)%1024, current_node->name, 
	  current_node->sequence, current_node->distance, 
	  ((uint64_t)current_node->left>>12)%1024, ((uint64_t)current_node->right>>12)%1024);
  printf("\nSelf  %s   ", current_node->left_profile_self);
  printf("       %s\n", current_node->right_profile_self);
  printf("Left  %s   ", current_node->left_profile_child);
  printf("Right  %s\n\n", current_node->right_profile_child);
  if (current_node->right) PrintTree(current_node->right);
}

/***********************************************************************
 * FUNCTION TO DETECT GAP COLUMNS.  IT KEEPS A RECORD IN A BINARY STRING
 * WHICH WILL BE USED LATER ON IN PRINT ALIGNMENT
 *
 * ARGUMENTS:
 * aln  (TreeNode *) The alignment
 *
 * RETURN VALUE:
 * (char *) The binary string storing the existence of gap column
 ***********************************************************************/
static char *FindGapColumn(const Alignment& aln) {
  int x = 0, y = 0, rowI = 0, columnI = 0;
  int allGap, isLeaf;
  char *GapColumn;
  const int aln_n = aln.r[0].sequence.length();
  GapColumn = (char *)malloc((int)( (aln_n) * sizeof(char) ));


  if (GapColumn != NULL) {
    // Find the first non-internal node sequence
	const int aln_m = aln.r.size();
    for (rowI = 0; rowI < aln_m && strncmp(aln.r[rowI].name.c_str(), "Internal", 8) == 0; rowI++)	{
    }
    x = 0;
    while (rowI < aln_m && x < aln_n && aln.r[rowI].sequence[x] != '\0')	{
      // Skip over non-gap columns
      for (allGap = 0; x < aln_n && aln.r[rowI].sequence[x] != '\0' &&
           aln.r[rowI].sequence[x] != SPACE; x++) {
        GapColumn[x] = (char)(allGap + 48);
      }

      // a Gap is found at pos x, traverse down the column to check whether the entire column has gaps
      if (x < aln_n && aln.r[rowI].sequence[x] != '\0' && aln.r[rowI].sequence[x] == SPACE)	{
        for (y = rowI + 1, allGap = 1; allGap && y < aln_m; y++)	{
          // Skip over internal node sequences
          while ( y < aln_m && strncmp(aln.r[y].name.c_str(), "Internal", 8) == 0 )	{
            y++;
          }
          // Determine if the character at column x of row y is a gap
          if (y < aln_m)	{
            allGap &= (aln.r[y].sequence[x] == SPACE);
          }
        }
        // Update the GapColumn string
        GapColumn[x] = (char)(allGap + 48);
        x++;
      }
    }
    GapColumn[x] = '\0';
  } else {
    MemoryRequestFailure("FindGapColumn()");
  }

  return GapColumn;
}


/***********************************************************************
 * FUNCTION TO PRINT AN ALIGNMENT. USED TO PRINT THE 'TRUE' ALIGNMENT
 * OF THE GENERATED SEQUENCES.
 *
 * ARGUMENTS:
 * aln                            The Alignment to be printed
 * temp_true_alignment_file_name       (char *) The name of the file in which to print the alignment
 * GapColumn			       (char *) The string indicating whether the entire column is gap
 ****************************************************************************************************/
static void PrintAlignment(const Alignment& aln, char *temp_true_alignment_file_name, char *GapColumn) {
  int i, j, count;

  if (TrueAlignmentFileName == NULL) TrueAlignmentFile = stdout;
  else if ((TrueAlignmentFile = fopen(temp_true_alignment_file_name, "w")) == NULL) {
    fprintf(stderr, "ERROR PrintAlignment(): cannot open file \"%s\"\n", temp_true_alignment_file_name);
    exit(EXIT_FAILURE);
  }

  const int aln_m = aln.r.size();
  const int aln_n = aln.r[0].sequence.length();
  for (i = 0; i < aln_m; i++) {
    if (strncmp(aln.r[i].name.c_str(), "Internal", 8) != 0 && strstr(aln.r[i].name.c_str(), "Neg") == NULL) {
      count = 0;
      fprintf(TrueAlignmentFile, ">");
      fprintf(TrueAlignmentFile, "%-20s", aln.r[i].name.c_str());
      fprintf(TrueAlignmentFile, "\n");
      for (j = 0; j < aln_n; j++)	{
        // Skip over the gap column
        while (j < aln_n && aln.r[i].sequence[j] != '\0' && GapColumn[j] == '1') {
          j++;
        }
        if (j < aln_n && aln.r[i].sequence[j] != '\0' && GapColumn[j] == '0') {
          fprintf(TrueAlignmentFile, "%c", aln.r[i].sequence[j]);
          count++;
          if (count == 79) {
            fprintf(TrueAlignmentFile, "\n");
            count = 0;
          }
        }
      }
      fprintf(TrueAlignmentFile, "\n");
    }
  }
  if (TrueAlignmentFileName) fclose(TrueAlignmentFile);
}


/////////////////////////////////
/////////////////////////////////
///                           ///
///   COMMAND LINE OPTIONS    ///
///                           ///
/////////////////////////////////
/////////////////////////////////

// Currently Taken:
// a b c d e f g i j k l m p r s t v w x y z u q o h
// Currently Available:
// n  

static struct poptOption optionsTable[] = {
  {
    "alignment",
    'a',
    POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
    &TrueAlignmentFileName,
    0,
    "name of output alignment file, in fasta format",
    "<string>"
  },
  {
    "branch",
    'b',
    POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
    &TreeBranchScale,
    0,
    "branch length scale multiplier",
    "<double>",
  },
  {
    "eFactor",
    'c',
    POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
    &EvolScale,
    0,
    "Evolutionary Scale factor for the distribution of indel lengths",
    "<double>",
  },
  {
    "debug",
    'd',
    POPT_ARG_NONE,
    &debug_mode,
    0,
    "debug mode",
    "<debug mode>",
  },
  {
    "tree",
    'f',
    POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
    &TreeFileName,
    0,
    "name of tree file, only bifurcations are allowed",
    "<string>"
  },
  {
    "indelFrequency",
    'g',
    POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
    &INDEL_FREQ,
    0,
    "the indel frequency for evolutionary time c",
    "<double>",
  },
  {
    "maxIndel",
    'l',
    POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
    &max_indel_length,
    0,
    "the maximum insertion/deletion length",
    "<int>"
  },
  {
    "subModel",
    'p',
    POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
    &subModel,
    0,
    "substitution model: 0 for PAM, 1 for JTT, 2 for PMB",
    "<int>"
  },
  {
    "rootLength",
    'r',
    POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
    &root_sequence_length,
    0,
    "root sequence length",
    "<int>"
  },
  {
    "sequence",
    's',
    POPT_ARG_STRING | POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
    &FastaFileName,
    0,
    "name of output sequence file",
    "<string>"
  },
  {
     "phylip",
     'j',
     POPT_ARG_STRING | POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
     &PhylipFileName,
     0,
     "name of output phylip file",
     "<string>"
  },
  {
    "alpha",
    'x',
    POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
    &alpha,
    0,
    "gamma alpha. Set to -1 for equal evolutionary rates",
    "<double>",
  },
  {//PN August 2005
    "interleaved",
    'i',
    POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
    &interleave,
    0,
    "interleaved output",
    "<int>",
  },
    {//PN August 2005
    "benner",
    'y',
    POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
    &benner,
    0,
    "benner",
    "<int>",
  },
  {//PN September 2005
    "variableGamma",
    'v',
    POPT_ARG_INT | POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
    &variablegamma,
    0,
    "variable gamma",
    "<int>"
  },
  {//PN Sept 2005
    "bennerk",
    'k',
    POPT_ARG_INT | POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
    &bennerk,
    0,
    "benner k factor",
    "<int>"
  },
  {//JA Aug 2006
    "indelWeight",
    'w',
    POPT_ARG_DOUBLE | POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
    &indelWeight,
    0,
    "indel weighting parameter",
    "<double>"
  },
  {//JA Sep 2006
    "extraTerminalIndels",
     't',
     POPT_ARG_DOUBLE | POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
     &extra_terminal_indels,
     0,
     "extra terminal indels parameter",
     "<double>"
  },
  {// ET May 2007
    "rootSequence",
    'z',
     POPT_ARG_STRING | POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
     &RootSequenceFileName,
     0,
     "name of root sequence file",
     "<string>"
  },
   {// ET Jan 2009
    "indel distribution",
    'u',
     POPT_ARG_STRING | POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
     &GapDistFileName,
     0,
     "name of indel distribution file",
     "<string>"
  },

  {// ET Jan 2009
    "indel output",
    'o',
     POPT_ARG_STRING | POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
     &IndelFileName,
     0, 
     "file to output length of indels as created",
     "<string>"
  },

  {// ET Jan 2012 //brings back 2003 memories though ...
        "correlation file",
        'h',
        POPT_ARG_STRING | POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
        &CorrelFileName,
        0,
        "Correlated site pairs: i    j   correlation, only one site per. Indels are disabled for this option.",
        "<string>"
    },


  {// PN Dec 2007
    "variableBranch",
    'e',
    POPT_ARG_DOUBLE| POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
    &VariableBranch,
    0,
    "variable branch length scale multiplier",
    "<double>"
  },
  {// PN Dec 2007
    "branchExtinction",
    'm',
    POPT_ARG_DOUBLE | POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
    &BranchExtinction,
    0,
    "branch extinction probability",
    "<double>"
  },

  {// PN Dec 2007
    "InsDelRatio",
    'q',
    POPT_ARG_DOUBLE | POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT,
    &InDelRatio,
    0,
    "prob of insertion (1 - prob of deletion) (would be 0.5 for equal frequencies)",
    "<double>"
  },

  POPT_AUTOHELP // No comma
  POPT_TABLEEND
};

// ZHUOZHI
void printTree ( TreeNode * t, char * p )
{
    char * q;
    if ( t != NULL ) {
       q = (char*) malloc ( sizeof(char) * (strlen(p)+3) );
       strcpy (q, p);
       strcat (q, "  ");
       printTree ( t->left,  q );
       printTree ( t->right, q );
       printf ( "%s%s\n", p, t->name );
    }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main (int argc, char **argv) {
    int i;
    char c, *temp_fasta_file_name = NULL, *temp_true_alignment_file_name = NULL, *temp_phylip_file_name;
    poptContext optCon;
    Alignment *aln;     // A "true" alignment
    TreeNode *t;        // The tree.
    char *GapColumn;    // A string storing 0 or 1: 
                        //    1 means a column of gaps is found in the alignment;
                        //    0 otherwise

    //FILE * seedfile;
    //int fogo = GetCurrentProcessId();
    //seedfile = fopen("seed.txt", "at");	
    //fprintf(seedfile, "test %d\n", fogo);
    //fclose (seedfile);

    printf ( "WARNING: The maximum length of the sequence is %d\n", 
                                                MAX_SEQUENCE_LENGTH );
    ///PN
    // Initialize the random number generators
    //SetSeed(time(NULL));
    SetSeed(getpid()); //Linux uses getpid() to set the seed value - changed due to implementation of multiple runs
    //Windows uses GetCurrentProcessId()
    //SetSeed(GetCurrentProcessId());

    ///PN
    /*********** GET COMMAND LINE ARGUMENTS ***********************************/
    optCon=poptGetContext("simprot", argc, (const char **)argv, optionsTable, 0);
    if ((int)(c = poptGetNextOpt(optCon)) >= 0) {
        fprintf(stderr, "ERROR bad option ?");
        return EXIT_FAILURE;
    }

    if ((int)(c = poptGetNextOpt(optCon)) < -1) {
        fprintf(stderr, "Simprot: bad argument %s: %s\n", 
        poptBadOption(optCon, POPT_BADOPTION_NOALIAS), poptStrerror(c));
        poptPrintUsage(optCon, stdout, 0);
        return EXIT_FAILURE;
    }
    /**************************************************************************/

    /*******  INITIALIZATION STUFF  *******************************************/
    t = InitTree();
FlagTree(t);
    /********* MATRIX INITIALIZATION STUFF ************************************/
    InitProtMats();
    MakeProtFreqs();
    InitSymbolCumulativeDensity();

    /**************** DEBUG MODE **********************************************/  
    if (debug_mode == 1) {
        printf("Options:\n\n");
        printf("Inputs:\n%s\n\n", TreeFileName);
        printf("Outputs:\n%s\n%s\n\n", TrueAlignmentFileName, FastaFileName);
        printf("Others:\nroot sequence length:\t%d\nmax indel length:\t%d\n\n",
	                                        root_sequence_length, max_indel_length);
    }

    /* Prepare the output file name buffers */
    if (FastaFileName != NULL)
        temp_fasta_file_name = (char *)malloc((strlen(FastaFileName) + 
						10 /* MAGIC */) * sizeof(char));
    if (TrueAlignmentFileName != NULL)
        temp_true_alignment_file_name = (char *)malloc((
	    strlen(TrueAlignmentFileName) + 10 /* MORE MAGIC */) * sizeof(char));
    if (PhylipFileName != NULL) //PN
	    temp_phylip_file_name = (char *)malloc((strlen(PhylipFileName) + 10) /* no magic here, just plain code */ * sizeof(char)); //PN
   
   if (GapDistFileName != NULL) {
    // read in gap distribution into IndelLengthFreq, and clear the rest of the array:
		ReadGapDist(GapDistFileName);
		benner=2;
	} 
    
    if (CorrelFileName != NULL) {
        // read in coorrelated pairs and value:
		ReadCorrelatedPairs(CorrelFileName);
		benner=3;
	}     
    
   if (IndelFileName != NULL) {
    //print in gap distribution
      pFile = fopen (IndelFileName,"w");

	}


    /***** EVOLVE *************************************************************/
    if (RootSequenceFileName != NULL){
        // read in root sequence
        t->sequence = ReadRootSequence();
        root_sequence_length=(int)strlen(t->sequence);
    }
    else {
        // Generate the root sequence
        t->sequence = RandomSequence(root_sequence_length);
    }
    // Evolve the root sequence
    Evolve(t);

    // JA Aug 06 - print the tree for testing
//      debug_mode = 1;

    // If we are debugging, print the tree
    if (debug_mode == 1) PrintTree(t);

    //PN September 2005
    if (benner == 1) {
        GapExtensionBenner();
    }
    else if (benner == 0) {
		GapExtension(); //PN March 2005
    }
	else if (benner == 2) { // other distribution from file
		GapExtensionOther(); //ET Jan 2009
	}
	
    //FILE * gapsfile;
    //gapsfile = fopen("gaps.txt", "w");
    fprintf(stdout, "Gap Extension Penalty: %lf \nGap Insertion Penalty: %lf\n\n", GapOpenProbability,GapExtensionPenalty); //PN March 2005
    //fclose (gapsfile);
    /**** Print the sequences in Fasta format *****/
    if (FastaFileName)
        sprintf(temp_fasta_file_name, "%s", FastaFileName);
    if (PhylipFileName)
        sprintf(temp_phylip_file_name, "%s", PhylipFileName);
    // Because PrintFastaFormat() is recursive
    InitFastaFile(temp_fasta_file_name);
    PrintFastaFormat(t);
    // Because PrintFastaFormat() is recursive
    if (FastaFileName) fclose(FastaFile);
   //if (PhylipFileName) fclose(PhylipFile);

    /***** Print the true alignment *****/
	aln = NULL;
	GapColumn = NULL;
	if (TrueAlignmentFileName != NULL || PhylipFileName != NULL) {
		aln = new Alignment(t);
		GapColumn = FindGapColumn(*aln);
	} // GetAlignment(t) and FindGapColumn were originally replicated in both sections below.
    if (TrueAlignmentFileName != NULL) {
        sprintf(temp_true_alignment_file_name, "%s", TrueAlignmentFileName);
        PrintAlignment(*aln, temp_true_alignment_file_name, GapColumn);
    }
    // (PN Mar05)
    if (PhylipFileName != NULL) {
         sprintf(temp_phylip_file_name, "%s", PhylipFileName);
        //InitPhylipFile(temp_phylip_file_name); //PN March 2005  
        PrintPhylipFormat(*aln, temp_phylip_file_name, GapColumn);	//PN March 200
    }
	delete GapColumn;
	delete aln;

    // Clear the tree to go again
    ClearTree(t);
	
	   if (IndelFileName != NULL) {
    //print in gap distribution
     fclose (pFile);

	}

	
	
    return EXIT_SUCCESS;
}

