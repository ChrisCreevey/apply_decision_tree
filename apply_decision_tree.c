/*
 *  Apply_Decision_tree.c
 *  
 *
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <time.h>
#include <math.h>
#include <signal.h>
#include <unistd.h> 

/****** Define  ********/

#define DESCRIPTION_LENGTH 100000
#define MAXNUMPHENOTYPES 100

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TEMP
#define TEMP 2
#endif



/******** Structure definitions *************/

struct node {
	
	int name;					/* This holds the node number */
	char *label;				/* This holds the node label (genefamily name or result) */

	float left_value; 				/* this holds the cutoff value to decide true/false */
	float right_value; 				/* this holds the cutoff value to decide true/false */

	char left_condition[DESCRIPTION_LENGTH]; 		/* this holds the text version of the condition ">","<", ">=", etc... */
	char right_condition[DESCRIPTION_LENGTH]; 		/* this holds the text version of the condition ">","<", ">=", etc... */

	char phenotype[100];

	int traversal_count; /* counter to record how many times we visited each noe fo the tree during assessment of genomes */
	float phenotype_count; /* Records the model infomration for the leaves of the tree on how many genomes ended up at this point of the tree */

	int arff_line_num;	/* holds the line number that the information */
	float *freqs;		/* hold the frequencies of this gene familily across all the genomes in the arff file */

	struct node *left;   	/* This points to a daughter node if the condition is true */
	struct node *right;   	/* This points to a daughter node if the condition is false */

	struct node *parent;		/* This points to the parent node, however its only set on the first sibling at each level */
	} node_type;

struct node **node_array = NULL;
int numfields=0, nodecount=0, edgecount=0, num_genomes=0, useGenomeNames = TRUE, numphenotypes=0, numgenes=0; 
float **gene_pos_contrib_to_pheno = NULL, **gene_neg_contrib_to_pheno = NULL;
char **genome_names=NULL, **genomes=NULL, **phenotypes, **genenames = NULL, file_stem[1000], outfile_name[10000];
FILE *outfile = NULL, *path_outfile = NULL, *pos_contrib_file = NULL, *neg_contrib_file = NULL;




/**** function definitions ***/

void print_tree_details(struct node *position);
void read_tree (FILE * treefile);
void read_arff (FILE * arff_file); 
void read_genome_names (FILE * genome_names_file);
int is_genefam_in_tree(struct node *position, char *famname);
void asses_genomes(struct node *position, int genome);
void print_tree_traversal_counts(struct node *position);
void print_newick_treefile(struct node *position, FILE *newickfile);
void calculate_gene_pheno_contribtion(int gene_num);
void count_all_phenos_pos(struct node *position, int gene_num);
void count_all_phenos_neg(struct node *position, int gene_num);
void print_tree_paths(struct node *position, char *path);
void print_gene_contib_phenos(void);




	
 int main(int argc, char const *argv[]){

	 FILE *treefile = NULL, *arff_file = NULL, *newickfile = NULL, *genomenamefile = NULL;
	 int i, j, found=FALSE;
	 char command[10000] = "grep \"^N[0-9]\" ", path[1000000], pathoutfile_name[10000], node_phenotype[100];

	 path[0] = '\0'; outfile_name[0] = '\0'; pathoutfile_name[0] = '\0'; file_stem[0] ='\0';

	 if(argc < 3)
	 	{	
	 	printf("apply_decision_tree\n\n\tUsage:\n\tapply_decision_tree <DOT-formmated tree file> <ARFF-formmated data file> <OPTIONAL genome-names-file>\n\n");
	 	printf("Output predictions for all genomes will be written to a file called <genome-names-file>.predictions.txt\n" );
	 	printf("All paths though the decision tree will be written to a file called <genome-names-file>.paths.txt\n\n" );
	 	exit(0);
	 	}
	 if(argc == 3)
	 	useGenomeNames = FALSE;

	 strcpy(file_stem, argv[1]);

	 strcat(command, argv[1]);
	 strcat(command, " | sed 's/\\[label=\"//g' | sed 's/\" *\\]//g'| sed 's/\\\".*\\]$//g' | sed \"s/\\(N[0-9]*\\) /\\1	/g\" | sed '/->/s/^/EDGE	/g' | sed '/^N[0-9]/s/^/NODE	/g' | awk '{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"#\"}' > formatted_tree_file.txt");
 	/* code */
	system(command); /* reformat the tree file using a shell script */

	/* STEP 1: read the treefile and build the decisiontree in memory */
	
	treefile = fopen("formatted_tree_file.txt", "r");
	printf("Starting reading decision tree file\n");
	read_tree (treefile); 
	fclose(treefile);

 	/* STEP 2: read in the genefamily profile infromation */
 	arff_file = fopen(argv[2], "r");
 	printf("Starting reading AFF file\n");
	read_arff (arff_file); /* read the arff_file and build the decisiontree in memory */
 	fclose(arff_file);

    /* STEP 2.5 Read in the genome_names file */
 	if(useGenomeNames){
	 	genomenamefile = fopen(argv[3], "r");

	 	printf("Starting to read genome name file\n");
	 	read_genome_names (genomenamefile);
	 	fclose(genomenamefile);
	 }

 	sprintf(outfile_name, "%s.predictions.txt", file_stem);
 	outfile = fopen(outfile_name, "w");

 	/* STEP 3: assess all the genome information based on the decision tree */
 	printf("Printing Genome predictions to %s\n", outfile_name);
 	fprintf(outfile, "Genome_name\tGenome_number\tNode_number\tNode_label\n");
 	for(i=0; i<num_genomes; i++)
 		{
 		asses_genomes(node_array[0], i);
 		}
 	fclose(outfile);
 	/* STEP 4: print out node traversal counts across the tree for from the genome assessments */
 	printf("\nDecision Tree traversal counts:\n");
 	print_tree_traversal_counts(node_array[0]);

 	/* STEP 5: print Newick formatted output of decision tree with traversal stats included */
 	printf("\nCreating Newick formatted version of decision tree\n");
 	sprintf(outfile_name, "%s.newick_output.txt", file_stem);
 	newickfile = fopen(outfile_name, "w");
 	fprintf(newickfile, "(");
 	print_newick_treefile(node_array[0], newickfile);
 	fprintf(newickfile, ");\n");
 	fclose(newickfile);
 	printf("\tNewick tree written to file newick_output.txt\n");

 	/* STEP 5.5 calculate contribution of each gene to each phenotype on the tree */

 	genenames=malloc(nodecount*sizeof(char *)); /* Allocate memory to hold the names of all the genes */
 	for(i=0; i<nodecount; i++){
 		genenames[i] = malloc(1000*sizeof(char));
 		genenames[i][0] = '\0';
 		}


 	phenotypes=malloc(MAXNUMPHENOTYPES*sizeof(char *));
 	for(i=0; i<=MAXNUMPHENOTYPES; i++)
 		{
 		phenotypes[i]=malloc(100*sizeof(char));
 		phenotypes[i][0] = '\0';
 		}
 		/* Count how many phenotypes are in the tree */
 	for(i=0; i<nodecount; i++){
 		if(node_array[i]->left == NULL && node_array[i]->right == NULL) /* this is a terminal node (leaf) */
 			{
 			strcpy(node_array[i]->phenotype, strtok(node_array[i]->label, " ")); /* Get the phenotype */
 			node_array[i]->phenotype_count = atof(strtok(NULL, "( /)")); /* capture the phenotype count on this leaf */
 			found=FALSE;
 			for(j=0; j<numphenotypes; j++) /* Look to see if this is the first time we've seen this phenotype */
 				{
 				if(strcmp(node_array[i]->phenotype, phenotypes[j])== 0)
 						found=TRUE;
 				}
 			if(!found){ /* If this is the first time we've seen this phenotype */
 				strcpy(phenotypes[numphenotypes], node_array[i]->phenotype);
 				numphenotypes++;
 				}
			}
		else  /* This is an internal node, so capture the gene name */
			{
			found = FALSE;
			for(j=0; j<numgenes; j++) /* Look to see if this is the first time we've seen this gene */
 				{
 				if(strcmp(node_array[i]->label, genenames[j])== 0)
 						found=TRUE;
 				}
 			if(!found){ /* If this is the first time we've seen this gene */
 				strcpy(genenames[numgenes], node_array[i]->label);
 				numgenes++;
 				}
			}
 		}

 	printf("Number of Phenotypes = %d\n", numphenotypes);
 	printf("phenos found:\t");
 	for(j=0; j<numphenotypes; j++) printf("%s\t", phenotypes[j]);
 	printf("\n");
 	printf("Genes found:\t");
 	for(j=0; j<numgenes; j++) printf("%s\t", genenames[j]);
 	printf("\n");
 
 	/* Assign the global arrays that record the values contribution of each gene to each phenotype */
 	gene_pos_contrib_to_pheno = malloc(numgenes*sizeof(float*));
 	gene_neg_contrib_to_pheno = malloc(numgenes*sizeof(float*));

 	for(i=0; i<numgenes; i++){
 		gene_pos_contrib_to_pheno[i] = malloc(numphenotypes*sizeof(float));
 		gene_neg_contrib_to_pheno[i] = malloc(numphenotypes*sizeof(float));
 		for(j=0; j<numphenotypes; j++){
 			gene_pos_contrib_to_pheno[i][j]=0;
 			gene_neg_contrib_to_pheno[i][j]=0; /* assign all values to zero to start */
 			}
 		}

 	for(j=0; j<numgenes; j++){  /* Now for each unique gene name, find all instnaces of it in the tree and calculate the toal contribution to each phenotype */
	 	calculate_gene_pheno_contribtion(j);
	 }

	 print_gene_contib_phenos();

 	/* STEP 6 print all paths through trees (iE the paths to resistance) */

 	sprintf(pathoutfile_name, "%s.paths.txt", file_stem);
 	printf("Printing paths through Decision trees to %s\n", pathoutfile_name);
 	path_outfile = fopen(pathoutfile_name, "w");

 	fprintf(path_outfile, "Path_terminal_node\tLabel\tTime_traversed\tPath\n");
 	print_tree_paths(node_array[0], path);

 		/* test that the tree and data was read in correctly */
/* 	printf("Starting Tree traversal\n"); */
 /*	print_tree_details(node_array[0]); */ /* node 0 is always the top if the tree */
 /*	printf("Done tree traversal\n");  */

 	/* Clean up allocated memory before quitting */
 	for(i=0; i<nodecount; i++)
 		{
 		free(node_array[i]->freqs);
 		free(node_array[i]->label);
 		free(node_array[i]);
 		}
 	free(node_array);

 	return 0;
	}


void print_gene_contib_phenos(void){

	int i=0, j=0;
	 /* Print out the resultsing counts for contribution of each gene to each phenotype */
	 sprintf(outfile_name, "%s.Gene_presense_contrib_to_phenotype.txt", file_stem);
	 pos_contrib_file = fopen(outfile_name, "w");
	 sprintf(outfile_name, "%s.Gene_absense_contrib_to_phenotype.txt", file_stem);
	 neg_contrib_file = fopen(outfile_name, "w");

	 fprintf(pos_contrib_file, "Gene_name");  /* Print the header */
	 fprintf(neg_contrib_file, "Gene_name");  /* Print the header */
	 for(i=0; i<numphenotypes; i++) {
	 	fprintf(pos_contrib_file, "\t%s", phenotypes[i]);
	 	fprintf(neg_contrib_file, "\t%s", phenotypes[i]);
	 	}
	 fprintf(pos_contrib_file, "\n");
	 fprintf(neg_contrib_file, "\n");

	 /* For each gene print the total contributions to each phenotype */
	 for(i=0; i<numgenes; i++){
	 	fprintf(pos_contrib_file, "%s", genenames[i]);
	 	fprintf(neg_contrib_file, "%s", genenames[i]);
	 	for(j=0; j<numphenotypes; j++){
	 		fprintf(pos_contrib_file, "\t%f", gene_pos_contrib_to_pheno[i][j]);
	 		fprintf(neg_contrib_file, "\t%f", gene_neg_contrib_to_pheno[i][j]);
	 		}
		fprintf(pos_contrib_file, "\n");
		fprintf(neg_contrib_file, "\n");
	 	}
	 fclose(pos_contrib_file);
	 fclose(neg_contrib_file);

	}

void print_tree_details(struct node *position)
	{
	struct node *start = position;
	int i;
	
	while(position != NULL)
		{
		printf("Node number %d\n", position->name);
		if(position->left != NULL){
			printf("If %s is %s %f then goto node %d\n", position->label, position->left_condition, position->left_value, position->left->name);
			printf("If %s is %s %f then goto node %d\n", position->label, position->right_condition, position->right_value, position->right->name);
			printf("ARFFLINENUM= %d\n", position->arff_line_num);
			printf("Data from ARFF file for node %s (numfields=%d):\n", position->label, numfields);
			if(position->freqs!=NULL){
				for(i=0; i<num_genomes; i++){
					printf("%d;", (int)position->freqs[i]);
					}
				printf("\n"); 
				printf("FOUND IN ARFF FILE!!\n");
			}else{
				printf("NOT FOUND IN ARFF FILE\n");
			}
		}else{
			printf("RESULT: %s\n", position->label);
			}

		if(position->left != NULL) print_tree_details(position->left);
		if(position->right != NULL) print_tree_details(position->right);
		position=NULL;
		}
	}


void read_tree (FILE * treefile)
	{
	char nodename[DESCRIPTION_LENGTH], type[DESCRIPTION_LENGTH], nodelabel[DESCRIPTION_LENGTH], comparison[DESCRIPTION_LENGTH];
	int i, j, node_num, from_node, to_node, edge_num, l;
	float value;

	/* read the decision tree data to capture the information on the nodes */
	/* count the numbver of nodes */
	while(!feof(treefile)){	
 		fscanf(treefile, "%s\t%s\t%s\t%s\n", type, nodename, comparison, nodelabel);
 		printf("Type = %s, Nodename = %s, Nodelabel = %s" "%s\n", type, nodename, comparison, nodelabel); 
 		if(strcmp(type, "NODE")==0) nodecount++;
 		if(strcmp(type, "EDGE")==0) edgecount++;
 		}
 	printf("\t%d nodes and %d edges found in decision tree\n", nodecount, edgecount);

 	/* assign array of structures for storing the decision tree information */
 	node_array=malloc(nodecount*sizeof(node_type));
 	for(i=0; i<nodecount; i++){ /* define the nodes of the tree */
 		node_array[i]=malloc(sizeof(node_type));
 		node_array[i]->name=i;
 		node_array[i]->label=malloc(DESCRIPTION_LENGTH*sizeof(char));
 		node_array[i]->label[0]='\0';
 		node_array[i]->phenotype[0]='\0';
 		node_array[i]->left_value=0;
  		node_array[i]->right_value=0;
 		node_array[i]->left_condition[0]='\0';
 		node_array[i]->right_condition[0]='\0';
 		node_array[i]->traversal_count = 0;
 
 		node_array[i]->arff_line_num=-1;
 		node_array[i]->freqs=malloc(DESCRIPTION_LENGTH*sizeof(float));

 		node_array[i]->left=NULL;
 		node_array[i]->right=NULL;
 		node_array[i]->parent=NULL;
 		}


 	/*re-read in the treefile capturing the information in the tree and building it in memory */

	rewind(treefile); /* go back to the start of the file */

	while(!feof(treefile)){	
 		fscanf(treefile, "%s\t%s\t%s\t%s#\n", type, nodename, comparison, nodelabel); 
 		if(strcmp(type, "NODE")==0) {
			/* if this represents NODE information */
	 		
	 		sscanf(nodename, "N%d", &node_num); /* extract the number from the node ID */
	 		/*printf("NODE: %d, %s->%s\n", node_num, comparison, nodelabel); */
   			
			node_array[node_num]->name = node_num;
			strcpy(node_array[node_num]->label, comparison); /* copy the node name */
			if(strcmp(nodelabel, "#")!= 0){	
				strcat(node_array[node_num]->label, " ");
				strcat(node_array[node_num]->label, nodelabel);
				}

	 		}else{
 			if(strcmp(type, "EDGE")==0){

				/* if this represents EDGE information */
	 		
	 			sscanf(nodename, "N%d->N%d", &from_node, &to_node); /* extract the number from the node ID */
	 			sscanf(nodelabel, "%f#", &value);
	 			/*strcpy(node_array[node_num], nodelabel); */
	 			/*printf("EDGE: from node %d, to node %d, if value is %s %f\n", from_node, to_node, comparison, value);*/

	 			if(node_array[from_node]->left == NULL){ /* if we haven;t assigned the left edge yet */
	 				strcpy(node_array[from_node]->left_condition, comparison);	
	 				node_array[from_node]->left_value = value;
	 				node_array[from_node]->left = node_array[to_node]; /*assign the left pointer to the to_node */
	 				node_array[to_node]->parent=node_array[from_node];

	 			}else{ /* assign the right edge */
	 				if(node_array[from_node]->right == NULL){
		 				strcpy(node_array[from_node]->right_condition, comparison);	
		 				node_array[from_node]->right_value = value;
		 				node_array[from_node]->right = node_array[to_node]; /*assign the right pointer to the to_node */
		 				node_array[to_node]->parent=node_array[from_node];

		 			}else{
		 				printf("ERROR: more than two child edges assigned to Node %d in input file\n", from_node);
		 				exit(1);
		 				}

	 				}

	 		}else{
	 				printf("Error type %s not recognised, please re-format the input tree file\n", type);
	 				exit(1);
 				}	
 			}
 		}
	}



void read_arff (FILE * arff_file)
	{
	char tmptext[DESCRIPTION_LENGTH], tmptext2[DESCRIPTION_LENGTH], *fam_name=NULL, c;
	int num=0, l, i;
	int foundnode=-1, linenum=0;
	const char delims[3] = " \n";
	const char delims2[3] = ",\n";
    char *token;
    float *genome_freq = NULL;


    genome_names=malloc(DESCRIPTION_LENGTH*sizeof(char*));
    for(i=0; i<DESCRIPTION_LENGTH; i++)
    	{
    	genome_names[i]=malloc(1000*sizeof(char));
    	genome_names[i][0] = '\0';
    	}

	fam_name=malloc(DESCRIPTION_LENGTH*sizeof(char));
	fam_name[0]='\0';
	tmptext[0]='\0';
	tmptext2[0]='\0';

	fgets(tmptext2, DESCRIPTION_LENGTH, arff_file);
	while(tmptext2[0] == '%' && !feof(arff_file)) fgets(tmptext2, DESCRIPTION_LENGTH, arff_file); /* Skip comment lines */
	if(!feof(arff_file)){

		token=strtok(tmptext2, delims); /* get first token (everything up to the first space) */
		/*printf("token=%s\n", token);*/
	 	if(strcmp(token, "@RELATION")!=0){
	 		printf("\tERROR: File not in ARFF format\n");	
		 
		}else{
	 		/* Read in the attribute section */
	 		fgets(tmptext2, DESCRIPTION_LENGTH, arff_file);
	 		while(tmptext2[0] == '%' && !feof(arff_file)) fgets(tmptext2, DESCRIPTION_LENGTH, arff_file); /* Skip comment lines */
	 		if(!feof(arff_file)){
		 	/*	printf("line=%s\n", tmptext2); */
		 		token=strtok(tmptext2, delims); /* get first token (everything up to the first space) */
		 		strcpy(tmptext, token);
		 		token= strtok(NULL, delims); /* Get next token (gene family name) */
		 		if(token!=NULL) strcpy(fam_name, token);
		 	/*	printf("tmptext=%s fam_name=%s\n",tmptext, fam_name); */

		 		while(strcmp(tmptext, "@ATTRIBUTE") == 0 && !feof(arff_file)){
		 			foundnode=is_genefam_in_tree(node_array[0], fam_name);  /* for each genefamily in the arff file, search the decision tree to see if it is there */
		 			linenum++;
		 			if(foundnode != -1){
		 				node_array[foundnode]->arff_line_num=linenum-1; /* If the gene family is in the decision tree, then foundnode will be equal to the node number that contains the genefamily */
		 			/*	printf("found attribute %s at node %d in decision tree! - linenum=%d\n", fam_name, foundnode, linenum-1); */
		 				}
					if(!feof(arff_file))fgets(tmptext2, DESCRIPTION_LENGTH, arff_file);
					while(tmptext2[0] == '%' && !feof(arff_file)) fgets(tmptext2, DESCRIPTION_LENGTH, arff_file); /* Skip comment lines */
					if(!feof(arff_file)){
				 		token=strtok(tmptext2, delims); /* get first token (everything up to the first space) */
				 		if(token!=NULL) strcpy(tmptext, token);

				 		token= strtok(NULL, delims); /* Get next token (gene family name) */
				 		if(token!=NULL) strcpy(fam_name, token);
				 		}
		 			}

		 		printf("\tNum attributes found = %d\n", linenum-1);
		 		numfields=linenum-1; /* number of attributes (fields) in arff file */
		 		if(strcmp(tmptext2, "@DATA")==0){ /* start reading the data */
		 			linenum=0;
		 			
		 			genome_freq = malloc((numfields+1)*sizeof(float)); /* this will hold the tokenised text for each data line */
		 			while(!feof(arff_file)){

		 				/* for each line of data, read in the numbers, tokenise them and then pass the array to add the data to the tree */
						fgets(tmptext, DESCRIPTION_LENGTH, arff_file); /* read the line of data*/
						while(tmptext[0] == '%' && !feof(arff_file)) fgets(tmptext, DESCRIPTION_LENGTH, arff_file); /* Skip comment lines */
		 				if(!feof(arff_file)){
			 				i=0;
			 				token=strtok(tmptext, delims2); /* get first token (everything up to the first comma) from the data line */
							while(token != NULL ){
								genome_freq[i]=atof(token); /* this assigns the values from the data lines to the array (converted to double/float) */
								strcpy(genome_names[linenum], token);
								i++; 
								token = strtok(NULL, delims2); /* Get next token (value) */

								}	

							/* this line represents one genome, the nodes on the decision tree represent one gene family in the dataset */
							/* For each node we need to add the information of its frequency from this genome */
			 				for(i=0; i<nodecount; i++) /* for each node of the tree */
			 					{
			 					if(node_array[i]->arff_line_num >= 0){
			 						node_array[i]->freqs[linenum] = genome_freq[node_array[i]->arff_line_num]; 	
			 						}else{
			 						node_array[i]->freqs[linenum] = 0;  /* if the gene family is not arff file, then set the frequency to zero */
			 						}				
			 					}

			 				linenum++;
			 				}
		 				} /* end while */
		 			num_genomes=linenum;
		 			printf("\tnumber of genomes = %d\n", num_genomes);
		 			}	
			 	}
			}	
		}
	 free(fam_name);
	}

/* Search tree to see if a specified gene family is contained in it */
int is_genefam_in_tree(struct node *position, char *famname)
	{
	int foundnode = -1;
	if(strcmp(position->label, famname)==0){
		foundnode=position->name;
	}else{
	if(position->left != NULL) foundnode=is_genefam_in_tree(position->left, famname);
	if(foundnode == -1){
		if(position->right != NULL) foundnode=is_genefam_in_tree(position->right, famname);
		}

	}
	return(foundnode);
}

void asses_genomes(struct node *position, int genome)
	{

	position->traversal_count++; /* increment this counter in this node to reflect taht we visited here in the assessment */

	/* Check Left codition */
	if(position->left != NULL){
		if(strcmp(position->left_condition, "<")==0)
			{	
			if(position->freqs[genome] < position->left_value)
				asses_genomes(position->left, genome);
			}else{
			if(strcmp(position->left_condition, "<=")==0)
				{
				if(position->freqs[genome] <= position->left_value)
					asses_genomes(position->left, genome);
				}else{
					if(strcmp(position->left_condition, ">")==0)
						{	
						if(position->freqs[genome] > position->left_value)
							asses_genomes(position->left, genome);
						}else{
							if(strcmp(position->left_condition, ">=")==0)
								{	
								if(position->freqs[genome] >= position->left_value)
									asses_genomes(position->left, genome);
								}else{
									if(strcmp(position->left_condition, "=")==0)
										{	
										if(position->freqs[genome] == position->left_value)
											asses_genomes(position->left, genome);
										}
								}
						}
				}	

			}
		}

	/* Check right codition */
	if(position->right != NULL){
		if(strcmp(position->right_condition, "<")==0)
			{	
			if(position->freqs[genome] < position->right_value)
				asses_genomes(position->right, genome);
			}else{
			if(strcmp(position->right_condition, "<=")==0)
				{
				if(position->freqs[genome] <= position->right_value)
					asses_genomes(position->right, genome);
				}else{
					if(strcmp(position->right_condition, ">")==0)
						{	
						if(position->freqs[genome] > position->right_value)
							asses_genomes(position->right, genome);
						}else{
							if(strcmp(position->right_condition, ">=")==0)
								{	
								if(position->freqs[genome] >= position->right_value)
									asses_genomes(position->right, genome);
								}else{
									if(strcmp(position->right_condition, "=")==0)
										{	
										if(position->freqs[genome] == position->right_value)
											asses_genomes(position->right, genome);
										}
									else{
										printf("ERROR:Unable to process genome %d through decision tree\n", genome);
										}
								}
						}
				}	

			}
		}

	if(position->left == NULL && position->right == NULL)
		{
		if(useGenomeNames)	
			fprintf(outfile, "%s\t%d\t%d\t%s\n", genomes[genome], genome, position->name, position->label);
		else
			fprintf(outfile, "%d\t%d\t%s\n", genome, position->name, position->label);
		}

	}

void print_tree_traversal_counts(struct node *position)
	{	
		printf("Node number %d\tLabel:%s\tVisit Count:%d\n", position->name, position->label, position->traversal_count);
		if(position->left != NULL) print_tree_traversal_counts(position->left);
		if(position->right != NULL) print_tree_traversal_counts(position->right);
	}


/* This traverses the tree to identify all terminal paths - IE the paths to resistance */
void print_tree_paths(struct node *position, char *path)
	{
		char left_path[100000], right_path[100000];	

		left_path[0] = '\0'; right_path[0]= '\0';

		if(position->left != NULL) 
			{
				sprintf(left_path, "%s; %s %s %f", path, position->label, position->left_condition, position->left_value);
				print_tree_paths(position->left, left_path);

			}
		if(position->right != NULL)
			{
				sprintf(right_path, "%s; %s %s %f", path, position->label, position->right_condition, position->right_value );
				print_tree_paths(position->right, right_path);
			}
		if(position->right == NULL && position->left == NULL) {
			/*printf("Node number %d\tLabel:%s\tVisit Count:%d\n", position->name, position->label, position->traversal_count);*/
			fprintf(path_outfile,"%d\t%s\t%d\t%s\n", position->name, position->label, position->traversal_count, path);
		}
	}



void print_newick_treefile(struct node *position, FILE *newickfile)
	{
	char *newlabel, tmpstring[DESCRIPTION_LENGTH];
	int i;
	
	newlabel=malloc(DESCRIPTION_LENGTH*sizeof(char));
	newlabel[0]='\0';
	tmpstring[0] = '\0';

	/* get rid of parentheses in labels as they will conflict with the newick format */
	for(i=0; i<DESCRIPTION_LENGTH; i++){
		switch(position->label[i]){
			case '(':
				newlabel[i]='{';
				break;
			case ')':
				newlabel[i]='}';
				break;
			case '\0':
				newlabel[i]='\0';
				i=DESCRIPTION_LENGTH;
				break;
			default:
				newlabel[i]=position->label[i];
				break;
			}

		}
	if(position->left != NULL)
		{
		fprintf(newickfile, "(");
	 	print_newick_treefile(position->left, newickfile);
	 	fprintf(newickfile, ",");
		}
	if(position->parent != NULL)  /* Go up one level in the tree to get the condition which results in a genome coming to this node */
		{	
		if(position->parent->right == position){
				sprintf(tmpstring, "%s%f", position->parent->right_condition, position->parent->right_value);
				
			}
		else{
				sprintf(tmpstring, "%s%f", position->parent->left_condition, position->parent->left_value);
				
			}
		}
	if(position->right != NULL) 
		{
		print_newick_treefile(position->right, newickfile);
		fprintf(newickfile, ")Prev%s-N%d-%s-%d",tmpstring, position->name, newlabel, position->traversal_count);
		}
	if(position->left == NULL && position->right == NULL)
		{	
		fprintf(newickfile, "Prev%s-N%d-%s-%d",tmpstring, position->name, newlabel, position->traversal_count);
		}

	free(newlabel);

	}


/* This fuinction is provided a gene name, then it does two steps. */
/* 1) find all nodes (one by one) that have this gene name in them */
/* 2) for each of the nodes calculates the Contributions to each phenotype and adde to the global arrays */


void calculate_gene_pheno_contribtion(int gene_num)  
	{
	int i;

	for(i=0; i<nodecount; i++){
		if(strcmp(node_array[i]->label, genenames[gene_num]) == 0 ){
			if(strchr(node_array[i]->left_condition, '<') != NULL){
				count_all_phenos_neg(node_array[i]->left, gene_num);
			}else{
				if(strchr(node_array[i]->left_condition, '>') != NULL){
					count_all_phenos_pos(node_array[i]->left, gene_num);
					}
				else{
					printf("warning: INTERNAL NODE %s DOESN'T HAVE AN OBVIOUS '>' OR '<' LEFT CONDITION, this will not be counted\n", genenames[gene_num]);
					}
				}
			if(strchr(node_array[i]->right_condition, '<') != NULL){
				count_all_phenos_neg(node_array[i]->right, gene_num);
			}else{
				if(strchr(node_array[i]->right_condition, '>') != NULL){
					count_all_phenos_pos(node_array[i]->right, gene_num);
					}
				else{
					printf("warning: INTERNAL NODE %s DOESN'T HAVE AN OBVIOUS '>' OR '<' RIGHT CONDITION, this will not be counted\n", genenames[gene_num]);
					}
				}
			
			}

		}

	
	}

/* This is called by calculate_gene_pheno_contribtion and sums all phenotype contributions */	
void count_all_phenos_pos(struct node *position, int gene_num)
	{
	int i=0;
	if(position->right != NULL) 
			count_all_phenos_pos(position->right, gene_num);
		
	if(position->left != NULL)
		count_all_phenos_pos(position->left, gene_num);

	if(position->left == NULL && position->right == NULL){
		i=0;
		while(strcmp(position->phenotype, phenotypes[i]) != 0)i++;
		gene_pos_contrib_to_pheno[gene_num][i]+=position->phenotype_count;
		}
	}


/* This is called by calculate_gene_pheno_contribtion and sums all phenotype contributions */	
void count_all_phenos_neg(struct node *position, int gene_num)
	{
	int i=0;
	if(position->right != NULL) 
		count_all_phenos_neg(position->right, gene_num);
		
	if(position->left != NULL)
		count_all_phenos_neg(position->left, gene_num);

	if(position->left == NULL && position->right == NULL){
		i=0;
		while(strcmp(position->phenotype, phenotypes[i]) != 0)i++;
		gene_neg_contrib_to_pheno[gene_num][i]+=position->phenotype_count;
		}
	}


void read_genome_names (FILE * genome_names_file)
	{
	int i;
	i=0;

	genomes=malloc(num_genomes*sizeof(char*));
	for(i=0; i<num_genomes; i++)
		{
		genomes[i]=malloc(1000*sizeof(char));
		genomes[i][0] = '\0';
		}
	i=0;
	while(!feof(genome_names_file)){
		if(i >= num_genomes) {
			printf("WARNING! The number of names in the genome names file doesn;t match the number of genomes lines in the arff file!\n");
			}
		else{
			fscanf(genome_names_file, "%s\n", genomes[i]);
			}
		i++;
		}
	}
