#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <boost/regex.hpp>
#include "util.h"
#include "proc_definition.h"
#include "pdb.h"
#include "dcd.h"

using namespace boost;

struct atom_def
{
	const char *atName;
	const char *resName;
	const char *segName;
	const char *resNum;		
};

struct hbonder
{
	const char *segName1;
	const char *segName2;

	const char *resName1;
	const char *resName2;

	const char *res_off1;
	const char *res_off2;
	
	const char *at1;
	const char *at2;

	int unit_at1; // 0 (any), 1, 2
	int unit_at2; 
	struct hbonder *next;
};

struct ion_def
{
	const char *seg_name;
	const char *res_off;	
	const char *res_name; 
	double max_dist;
};

struct pair_definition
{
	int seg[2];
	int sym;	
	int info;
	atom_def *lipid1;
	atom_def *lipid2;
	int nat1;	
	int nat2;	
	hbonder *bonders;
	int v2;

	ion_def *ions;
	int nions;
	int force_hbonds;
	int allow_no_ions;

	double cutoff;
	double binary_benefit;
	double binary_penalty;
	double hbond0_penalty;
	char *psfFile;

};


static struct pair_definition pdef;

// these functionalities are only available if you specify a PSF:

struct pdef_site_record
{
	int index;
	char *atName;
	char *resName;
	char *segName;
	int res;
};

struct pdef_site_record *allSites = NULL;
int sitesN = 0;
int sitesNSpace = 0;

int allowNoIons( void )
{
	return pdef.allow_no_ions;
}

double getCutoff( void )
{
	return pdef.cutoff;
}

int getPDefIndex( const char *atName, const char *resName, const char *segName, int res )
{
	for( int t = 0; t < sitesN; t++ )
	{
		if( res == allSites[t].res && 
			!strcasecmp( atName, allSites[t].atName) &&
			!strcasecmp( resName, allSites[t].resName) &&
			!strcasecmp( segName, allSites[t].segName) )
		{
			return allSites[t].index;
		}  
	}	
	
	return -1;
}

void clearPDefPSF( void )
{
	for( int i = 0; i < sitesN; i++ )
	{
		free(allSites[sitesN].atName );
		free(allSites[sitesN].resName );
		free(allSites[sitesN].segName );
	}
	free(allSites);
	allSites = NULL;
}

int loadPDefPSF( void )
{
	if( allSites )
	{
		printf("Error. Parameter definition PSF already loaded.\n");
		exit(1);
	}
	
	sitesNSpace = 10;
	allSites = (pdef_site_record *)malloc( sizeof(pdef_site_record) * sitesNSpace );

	if( pdef.psfFile )
	{
		FILE *theFile = fopen(pdef.psfFile,"r");

		if( !theFile )	
		{
			printf("Couldn't open pdef/PSF file '%s'.\n", pdef.psfFile );
			exit(1);
		}

		struct atom_rec *at = NULL;
	
		if( strlen( pdef.psfFile ) > 3 && !strncasecmp( pdef.psfFile + strlen(pdef.psfFile)-3, "pdb", 3 ) )
		{
			loadPSFfromPDB( theFile );
			rewind(theFile);
			at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );
			loadPDB( theFile, at ) ;
		}
		else
		{
			loadPSF( theFile );
			at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );
			loadDummy( at );
		}

		int pres = -1;	
		int cur_hex = 0;

		char psegid[256];

		int *mapped = (int *)malloc( sizeof(int) * curNAtoms() );
		int *linked = (int *)malloc( sizeof(int) * curNAtoms() );
		
		for( int a = 0; a < curNAtoms(); a++ )
		{
			mapped[a] = -1;
			linked[a] = -1;
		}

		int cur_unit = -1, cur_map = -1;
		for( int pass = 0; pass < 2; pass++ )
		{
			pres = -1;
			psegid[0] = '\0';

			int op = 0;
		
			if( is_sym() && pass == 1 ) break;

			for( int a = 0; a < curNAtoms(); a++ )
			{
				int is_new = 0;
				if( is_seg(pass) )
				{
					if( strcasecmp( at[a].segid, psegid) )
					{
						op = 0;
						is_new = 1;
					}
				}	
				else
				{
					if( at[a].res != pres )
					{
						op = 0;
						is_new = 1;
					}
				}
		
				if( use_atom( USE_ANYWHERE, pass, at[a].res, at[a].atname, at[a].resname, at[a].segid ) )
				{
					int used_unit = cur_unit;
					if( is_new )
					{
						if( pass == 1 && mapped[a] >= 0 )
						{
							cur_map = mapped[a];
					//		linked[cur_unit] = mapped[a];
						}
						else
						{
							cur_unit++;
							cur_map=-1;
							used_unit = cur_unit;
						}
						if( pass == 0 )
							mapped[a] = cur_unit;
						op=1;
					}
					
					if( cur_map >= 0 )
						used_unit = cur_map;
			
					if( sitesNSpace == sitesN )
					{
						sitesNSpace *= 2;
						allSites = (pdef_site_record *)realloc( allSites, sizeof(struct pdef_site_record) * sitesNSpace ); 
					}
	
					allSites[sitesN].atName = (char *)malloc( sizeof(char) * (1+strlen(at[a].atname) ) );
					allSites[sitesN].resName = (char *)malloc( sizeof(char) * (1+strlen(at[a].resname) ) );
					allSites[sitesN].segName = (char *)malloc( sizeof(char) * (1+strlen(at[a].segid) ) );
		
					strcpy( allSites[sitesN].atName, at[a].atname );
					strcpy( allSites[sitesN].resName, at[a].resname );
					strcpy( allSites[sitesN].segName, at[a].segid );
					allSites[sitesN].res = at[a].res;
	
//					printf("Setting index to %d\n", cur_unit );
					allSites[sitesN].index = used_unit; 
					sitesN++;
				}	
			
				pres = at[a].res;	
				strcpy( psegid, at[a].segid );
			}		
		}

		for( int a = 0; a < curNAtoms(); a++ )
			at[a].zap();
				
		
	}
	else
		return -1;		

	return 1;
}

int countTrackedSites( void )
{
	if( !pdef.psfFile )
		return -1;

}

int info_mode( void )
{
	return pdef.info;
}

int is_seg( int pair )
{
	if( is_sym ) return pdef.seg[0];

	return pdef.seg[pair];
}

int is_sym( void)
{
	return pdef.sym;
}

void load_pair_definition( const char *pairDefinitionFileName )
{
	pdef.cutoff=14;
	pdef.info=0;
	pdef.sym = 0;
	pdef.seg[0] = 0;
	pdef.seg[1] = 0;
	pdef.bonders = NULL;
	pdef.force_hbonds = 0;
	pdef.allow_no_ions = 0;
	pdef.v2 = 0;

	pdef.binary_benefit = 0;
	pdef.binary_penalty = 0;
	pdef.hbond0_penalty = 0;
	pdef.psfFile = NULL;

	FILE *theFile = fopen(pairDefinitionFileName, "r");

	if( !theFile )
	{
		printf("Couldn't open pair definition file '%s'.\n", pairDefinitionFileName );
		exit(1);
	}

	char *buffer = (char *)malloc( sizeof(char) * 1000000 );

	for( int pass = 0; pass < 2; pass++ )
	{
		pdef.nat1=0;
		pdef.nat2=0;
		pdef.nions=0;
		rewind(theFile);
		while( !feof(theFile) )
		{
			getLine( theFile, buffer );
			if( feof(theFile) ) break;
			
			char *tuse = buffer;
			
			while( *tuse && isspace(*tuse) ) tuse+=1; 
	
			if( *tuse == '#' ) continue;
			if( !*tuse ) continue;
	
			if( !strncasecmp( tuse, "INFO", 4 ) )
				pdef.info = 1;		
			else if( !strncasecmp( tuse, "SEG1", 4 ) )
				pdef.seg[0] = 1;		
			else if( !strncasecmp( tuse, "SEG2", 4 ) )
				pdef.seg[1] = 1;		
			else if( !strcasecmp( tuse, "V2" ) )
				pdef.v2 = 1;		
			else if( !strncasecmp( tuse, "SYM", 3 ) )
				pdef.sym = 1;	
			else if( !strncasecmp( tuse, "ALLOW_NO_IONS", 13 ) )
				pdef.allow_no_ions = 1;	
			else if( !strncasecmp( tuse, "FORCE_HBONDS", 11 ) )
				pdef.force_hbonds = 1;	
			else if( !strncasecmp( tuse, "CUTOFF", 6 ) )
			{
				double cutoff;
				int nr = sscanf( tuse, "CUTOFF %lf", &cutoff );
				if( nr == 1 )
					pdef.cutoff = cutoff;
				else
				{
					printf("Failed to read cutoff from line '%s'.\n", buffer );
					exit(1);
				}
			}
			else if( !strncasecmp( tuse, "BENEFIT", 7 ) )
			{
				double benefit;
				int nr = sscanf( tuse, "BENEFIT %lf", &benefit );
				if( nr == 1 )
					pdef.binary_benefit = benefit;
				else
				{
					printf("Failed to read benefit from line '%s'.\n", buffer );
					exit(1);
				}
			}
			else if( !strncasecmp( tuse, "PENALTY", 7 ) )
			{
				double penalty;
				int nr = sscanf( tuse, "PENALTY %lf", &penalty );
				if( nr == 1 )
					pdef.binary_penalty = penalty;
				else
				{
					printf("Failed to read penalty from line '%s'.\n", buffer );
					exit(1);
				}
			}
/*			// not using this,, much too confusing
 * else if( !strncasecmp( tuse, "HBOND0", 7 ) )
			{
				double HBOND0;
				int nr = sscanf( tuse, "HBOND0 %lf", &HBOND0 );
				if( nr == 1 )
					pdef.hbond0_penalty = HBOND0;
				else
				{
					printf("Failed to read HBOND0 from line '%s'.\n", buffer );
					exit(1);
				}
			}
*/
			else if( !strncasecmp( tuse, "ION", 3 ) )
			{
				if( pass == 1 )
				{
					char *res_read = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					char *seg_read = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					char *res_num = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					
					double max_dist;	
					int nr = sscanf(tuse, "ION %s %s %s %lf",
						res_num, res_read, seg_read, &max_dist ); 
		
					if( nr != 4 )
					{
						printf("Failed to read ion record from line '%s'. Was expecting \"ION pair# res#-regex res-regex seg-regex max_distance\"\n", buffer );
						exit(1); 
					}
					
					res_read = (char *)realloc( res_read, sizeof(char)*(1+strlen(res_read)) );
					seg_read = (char *)realloc( seg_read, sizeof(char)*(1+strlen(seg_read)) );
					res_num = (char *)realloc( res_num, sizeof(char)*(1+strlen(res_num)) );
				
					pdef.ions[pdef.nions].res_name = res_read;
					pdef.ions[pdef.nions].seg_name = seg_read;
					pdef.ions[pdef.nions].res_off = res_num;
					pdef.ions[pdef.nions].max_dist = max_dist;

					pdef.nions++;

				}
				else if( pass == 0 )
					pdef.nions++;
			}
			else if( !strncasecmp( tuse, "HBOND", 5) )
			{
				if( pass == 1 )
				{
					char *res_read1 = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					char *at_read1 = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					char *seg_read1 = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					char *res_num1 = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					char *res_read2 = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					char *at_read2 = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					char *seg_read2 = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					char *res_num2 = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
		
					char unit_code1 = '*', unit_code2 = '*';
					int pair_number = 0;
					if( pdef.v2 )
					{
						int nr = sscanf(tuse, "HBOND %s %s %s %c %s %s %s %s %c %s", 
							res_num1, at_read1, res_read1, &unit_code1, seg_read1,
							res_num2, at_read2, res_read2, &unit_code2, seg_read2 );
		
						if( nr != 10 )
						{
							printf("Failed to read hbonding record from line '%s'. Was expecting \"HBOND res#-regex at-regex res-regex {*|1|2} seg-regex ...\"\n", buffer );
							exit(1); 
						}
					}
					else
					{
						int nr = sscanf(tuse, "HBOND %s %s %s %s %s %s %s %s", 
							res_num1, at_read1, res_read1, seg_read1,
							res_num2, at_read2, res_read2, seg_read2 );
		
						if( nr != 8 )
						{
							printf("Failed to read hbonding record from line '%s'. Was expecting \"HBOND res#-regex at-regex res-regex seg-regex\"\n", buffer );
							exit(1); 
						}
					}
					at_read1  = (char *)realloc( at_read1, sizeof(char)*(1+strlen(at_read1)) );
					res_read1 = (char *)realloc( res_read1, sizeof(char)*(1+strlen(res_read1)) );
					seg_read1 = (char *)realloc( seg_read1, sizeof(char)*(1+strlen(seg_read1)) );
					res_num1 = (char *)realloc( res_num1, sizeof(char)*(1+strlen(res_num1)) );
					
					at_read2  = (char *)realloc( at_read2, sizeof(char)*(1+strlen(at_read2)) );
					res_read2 = (char *)realloc( res_read2, sizeof(char)*(1+strlen(res_read2)) );
					seg_read2 = (char *)realloc( seg_read2, sizeof(char)*(1+strlen(seg_read2)) );
					res_num2 = (char *)realloc( res_num2, sizeof(char)*(1+strlen(res_num2)) );

					
					hbonder *hbond = (struct hbonder *)malloc( sizeof(struct hbonder) );
					hbond->at1 = at_read1;
					hbond->resName1 = res_read1;
					hbond->segName1 = seg_read1;
					hbond->res_off1 = res_num1;
					hbond->at2 = at_read2;
					hbond->resName2 = res_read2;
					hbond->segName2 = seg_read2;
					hbond->res_off2 = res_num2;
					switch( unit_code1 )
					{
						case '*':
							hbond->unit_at1 = 0; // could be any
							break;
						case '1':
							hbond->unit_at1 = 1; 
							break;
						case '2':
							hbond->unit_at1 = 2; 
							break;
						default:
							printf("Can't interpret hydrogen bond code '%c'.\n", unit_code1 );
							exit(1);
							break;
					}
					switch( unit_code2 )
					{
						case '*':
							hbond->unit_at2 = 0; // could be any
							break;
						case '1':
							hbond->unit_at2 = 1; 
							break;
						case '2':
							hbond->unit_at2 = 2; 
							break;
						default:
							printf("Can't interpret hydrogen bond code '%c'.\n", unit_code2 );
							exit(1);
							break;
					}

					hbond->next = pdef.bonders;
					pdef.bonders = hbond;
				}

			}	
			else if( !strncasecmp( tuse, "ATOM", 4) )
			{
				if( pass == 1 )
				{
					char *res_read = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					char *at_read = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					char *seg_read = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
					char *res_num = (char*)malloc( sizeof(char) * (1+strlen(tuse) ) );
		
					int pair_number = 0;
					int nr = sscanf(tuse, "ATOM %d %s %s %s %s", &pair_number, res_num, at_read, res_read, seg_read );
		
					if( nr != 5 )
					{
						printf("Failed to read atom record from line '%s'. Was expecting \"ATOM sub_pair_# res# at-regex res-regex seg-regex\"\n", buffer );
						exit(1); 
					}
	
					at_read  = (char *)realloc( at_read, sizeof(char)*(1+strlen(at_read)) );
					res_read = (char *)realloc( res_read, sizeof(char)*(1+strlen(res_read)) );
					seg_read = (char *)realloc( seg_read, sizeof(char)*(1+strlen(seg_read)) );
					res_num = (char *)realloc( res_num, sizeof(char)*(1+strlen(res_num)) );

					if( pair_number == 1 )
					{
						pdef.lipid1[pdef.nat1].atName = at_read;
						pdef.lipid1[pdef.nat1].resName = res_read;
						pdef.lipid1[pdef.nat1].segName = seg_read;
						pdef.lipid1[pdef.nat1].resNum = res_num;

						pdef.nat1++;
					}
					else if( pair_number == 2 )
					{
						if( pdef.sym )
						{
							printf("ERROR. Doing a symmetric lipid interaction but an atom in sub-pair 2 was used. Only do sub-pair 1.\n");
							exit(1);
						}
						pdef.lipid2[pdef.nat2].atName = at_read;
						pdef.lipid2[pdef.nat2].resName = res_read;
						pdef.lipid2[pdef.nat2].segName = seg_read;
						pdef.lipid2[pdef.nat2].resNum = res_num;

						pdef.nat2++;
					}
					else
					{
						if( pdef.sym == 0 )
							printf("Pair# should be either 1 or 2.\n");
						else
							printf("Pair# should be 1.\n");
						exit(1);
					}
				}	
				else
				{				
					int res_off = 0;
					int pair_number = 0;
					int nr = sscanf(tuse, "ATOM %d %d", &pair_number, &res_off );

					if( pair_number == 1 )
						pdef.nat1++;
					else if( pair_number == 2 )
						pdef.nat2++;
				}
			}
			else if( !strncasecmp( tuse, "PSF", 3) )
			{
				if( pass == 1 )
				{
					pdef.psfFile = (char *)malloc( sizeof(char) * ( 1 + strlen(tuse) ) );
					int nr = sscanf( tuse, "PSF %s", pdef.psfFile );	

					if( nr != 1 )
					{
						printf("Failed to read the filename from PSF directive '%s'.\n", tuse );
						exit(1);
					}
				}
			}	
		}
		
		if( pass == 0 )
		{
			pdef.lipid1 = (atom_def *)malloc( sizeof(atom_def) * pdef.nat1 );
			pdef.lipid2 = (atom_def *)malloc( sizeof(atom_def) * pdef.nat2 );
			pdef.ions = (ion_def *)malloc( sizeof(ion_def) * pdef.nions );
		}
	}	
	
	// resort the list in reverse so that they show up in the order of the definition file.
	hbonder *new_list = NULL;
	hbonder *next = NULL;

	for( hbonder *hbond = pdef.bonders; hbond; hbond = next )
	{
		next = hbond->next;

		hbond->next = new_list;
		new_list = hbond;
	}
	
	pdef.bonders = new_list;

	free(buffer);
	fclose(theFile);
}

int force_hbonds( void )
{
	return pdef.force_hbonds;
}

int nextra( void )
{
	int nex = pdef.nions;

	return nex;
} 

int used_as_ion( int res, const char *resname, const char *segname, double *check_r )
{
	int matches = -1;
	char resNumBuffer[1024];
	sprintf(resNumBuffer, "%d", res );

	double check_rval = 0;

	for( int i = 0; i < pdef.nions; i++ )
	{
		regex wild_seg(pdef.ions[i].seg_name);
		regex wild_res(pdef.ions[i].res_name);
		regex wild_off(pdef.ions[i].res_off);
		
		int ok = 1;

		if( !regex_match(segname, wild_seg ) )
			ok = 0;
		if( !regex_match(resname, wild_res ) )
			ok = 0;
		if( !regex_match(resNumBuffer, wild_off ) )
			ok = 0;

		if( ok )
		{
			matches = i;
			if( pdef.ions[i].max_dist > check_rval )
			 	check_rval = pdef.ions[i].max_dist;
		}
	}

	if( check_r && matches >= 0 )
		*check_r = check_rval;

	return matches;
}

int use_atom( int use_type, int seg, int res, const char *atname, const char *resname, const char *segname )
{
	if( pdef.sym ) seg = 0;


	int use_it = 0;

	char resNumBuffer[1024];

	sprintf(resNumBuffer, "%d", res );

	if( use_type == USE_ANYWHERE )
	{
		if( pdef.seg[seg]  )
		{
			// just has to match the segment

			if( seg == 0  && pdef.nat1 > 0 )
			{
				regex wild_seg(pdef.lipid1[0].segName);
				
				if( regex_match(segname, wild_seg ) )
					use_it = 1;
			}
			else if( seg == 1  && pdef.nat2 > 0 )
			{
				regex wild_seg(pdef.lipid2[0].segName);
				
				if( regex_match(segname, wild_seg ) )
					use_it = 1;
			}
		}
		else
		{
			// has to match segment and residue name.
			if( seg == 0  && pdef.nat1 > 0 )
			{
				regex wild_res(pdef.lipid1[0].resName);
				regex wild_seg(pdef.lipid1[0].segName);
				
				if( regex_match( resname, wild_res) && regex_match(segname, wild_seg ) )
				{
					/*printf("Setting %s based on matches: %s and %s.\n", resname, 
								(regex_match( resname, wild_res) ? "true" : "false"),
								(regex_match( segname, wild_seg) ? "true" : "false") );
				*/	use_it = 1;
				}
			}
			else if( seg == 1  && pdef.nat2 > 0 )
			{
				regex wild_res(pdef.lipid2[0].resName);
				regex wild_seg(pdef.lipid2[0].segName);
				
				if( regex_match( resname, wild_res) && regex_match(segname, wild_seg ) )
					use_it = 1;
			}
		}
	}
	else
	{
		if( seg == 0 )
		{
			for( int i = 0; i < pdef.nat1; i++ )
			{ 
				regex wild_at(pdef.lipid1[i].atName);
				regex wild_res(pdef.lipid1[i].resName);
				regex wild_seg(pdef.lipid1[i].segName);
				regex wild_resn(pdef.lipid1[i].resNum);
		
				if( regex_match( atname, wild_at) && regex_match( resname, wild_res) && regex_match(segname, wild_seg ) && regex_match( resNumBuffer, wild_resn) )
					use_it = 1;
			}
		}
		else
		{
			for( int i = 0; i < pdef.nat2; i++ )
			{ 

				regex wild_at(pdef.lipid2[i].atName);
				regex wild_res(pdef.lipid2[i].resName);
				regex wild_seg(pdef.lipid2[i].segName);
				regex wild_resn(pdef.lipid2[i].resNum);
		
				if( regex_match( atname, wild_at) && regex_match( resname, wild_res) && regex_match(segname, wild_seg ) && regex_match( resNumBuffer, wild_resn))
					use_it = 1;
			}
		}
	}

	return use_it;
}

int binaryHBonder( struct atom_rec *at1, struct atom_rec *at2) // is this a binary hbond motif?
{
	char resNumBuffer1[256];
	char resNumBuffer2[256];

	sprintf(resNumBuffer1, "%d", at1->res );
	sprintf(resNumBuffer2, "%d", at2->res );

	int index = 0;
	for( hbonder *hbond = pdef.bonders; hbond; hbond = hbond->next )
	{
		regex wild_at1(hbond->at1); 
		regex wild_res1(hbond->resName1);
		regex wild_seg1(hbond->segName1);
		regex wild_resn1(hbond->res_off1);
			
		if( regex_match( at1->atname, wild_at1) && regex_match( at1->resname, wild_res1) && regex_match(at1->segid, wild_seg1 ) && regex_match(resNumBuffer1, wild_resn1 ) )
		{
			regex wild_at2(hbond->at2);
			regex wild_res2(hbond->resName2);
			regex wild_seg2(hbond->segName2);
			regex wild_resn2(hbond->res_off2);
			
			if( regex_match( at2->atname, wild_at2) && regex_match( at2->resname, wild_res2) && regex_match(at2->segid, wild_seg2 ) && regex_match(resNumBuffer2, wild_resn2 ) )
			{
				return index;

			}
		}

		if( pdef.sym )
		{
				// swap, look at reverse.
				
			regex wild_at1(hbond->at2); 
			regex wild_res1(hbond->resName2);
			regex wild_seg1(hbond->segName2);
			regex wild_resn1(hbond->res_off2);
			if( regex_match( at1->atname, wild_at1) && regex_match( at1->resname, wild_res1) && regex_match(at1->segid, wild_seg1 ) && regex_match(resNumBuffer1, wild_resn1 ) )
			{
				regex wild_at2(hbond->at1);
				regex wild_res2(hbond->resName1);
				regex wild_seg2(hbond->segName1);
				regex wild_resn2(hbond->res_off1);
			
				if( regex_match( at2->atname, wild_at2) && regex_match( at2->resname, wild_res2) && regex_match(at2->segid, wild_seg2 ) && regex_match(resNumBuffer2, wild_resn2 ) )
				{
					return index;
	
				}
			}
		}

		index++;
	}

	return -1;	
}

double pair_binary_penalty(void)
{
	return pdef.binary_penalty; 
}
double pair_binary_benefit(void)
{
	return pdef.binary_benefit; 
}
double pair_hbond0_penalty(void)
{
	return pdef.binary_penalty; 
}
