#define __spec_comparisonc__
#include "pdb.h"
#include "comparison.h"
// comparison routines for GM1.

int useAtomInChi2( struct atom_rec *at )
{
	if( !strncasecmp( at->resname, "CER180",4 ) )
	{
		if( !strcasecmp( at->atname, "HNF") ) 
			return 1;
		if( !strcasecmp( at->atname, "OF") ) 
			return 1;
		if( !strcasecmp( at->atname, "HO3") ) 
			return 1;
		if( !strcasecmp( at->atname, "O3") ) 
			return 1;
	}
			
	if( !strcasecmp( at->atname, "C1") && !strcasecmp( at->resname, "BGLC") && at->res == 2 ) 
		return 1;

	return 0;
}


// 64 most frequent hydrogen bonders.
bond_data allBonders[] = 
{
{"OF", "HNF", 1, 1}, // ca. 5.5k out of 38k structures @ 10% GM1
{"O12", "HO7", 6, 6},
{"O11", "HO7", 6, 6},
{"O12", "HN", 6, 6},
{"O11", "HN", 6, 6},
{"O2", "HO6", 3, 3},
{"OF", "H4S", 1, 1},
{"OF", "HO6", 1, 2},
{"OF", "H2G", 1, 1},
{"O1", "H1S", 2, 1}, //10
{"O2", "HN", 2, 4},
{"O3", "HNF", 1, 1},
{"O4", "HO6", 2, 3},
{"O6", "HO4", 4, 6},
{"O6", "HO6", 2, 2},
{"O11", "HO8", 6, 6},
{"OF", "H5S", 1, 1},
{"O12", "HO8", 6, 6},
{"O8", "HO9", 6, 6},
{"O5", "H1S", 2, 1}, //20
{"O6", "HO2", 3, 3},
{"O11", "HO2", 6, 5},
{"O6", "HO3", 2, 2},
{"O6", "H1", 2, 3},
{"O12", "HO2", 6, 5},
{"O9", "HO9", 6, 6},
{"O8", "HO8", 6, 6},
{"O4", "HO6", 6, 4},
{"O4", "HO2", 6, 5},
{"O4", "HN", 6, 6}, //30
{"O11", "HO9", 6, 6},
//{"O12", "HO9", 6, 6}, //32
//{"O3", "HO2", 2, 2},
//{"O6", "HO2", 4, 5},
//{"O6", "HN", 4, 6},
//{"O3", "HO6", 2, 2},
//{"O3", "HN", 5, 6},
//{"O11", "HO3", 6, 5},
//{"O2", "HO2", 2, 2},
//{"O6", "HO9", 2, 6},
//{"O2", "HN", 5, 6},
//{"O6", "HO7", 3, 6},
//{"O12", "HO3", 6, 5},
//{"O3", "HO2", 1, 2},
//{"O6", "HO4", 3, 6},
//{"O2", "HO4", 5, 6},
//{"O6", "HN", 4, 4},
//{"O4", "HO4", 6, 4},
//{"O4", "HO6", 4, 4},
//{"O2", "HO6", 5, 4},
//{"O6", "HN", 3, 4},
//{"O9", "HN", 6, 4},
//{"OF", "HO2", 1, 2},
//{"O3", "HO6", 1, 2},
//{"OF", "H3G", 1, 1},
//{"O3", "HO6", 2, 3},
//{"O8", "HO7", 6, 6},
//{"O4", "H61", 2, 3},
//{"O6", "HO4", 3, 4},
//{"OF", "H3S", 1, 1},
//{"O7", "HN", 6, 4},
//{"O12", "H7", 6, 6},
//{"O6", "H3", 2, 2},
//{"O4", "HO4", 4, 4}
};

int binaryHBonder( struct atom_rec *at1, struct atom_rec *at2) // is this a binary hbond motif?
{
	int nbonders = sizeof(allBonders) / sizeof(bond_data);

	int bx_use = -1;

	for( int b = 0; b < nbonders; b++ )
	{
		if(!( !strcasecmp( at1->atname, allBonders[b].atName1) && !strcasecmp( at2->atname, allBonders[b].atName2)) && 
		!( !strcasecmp( at1->atname, allBonders[b].atName2) && !strcasecmp( at2->atname, allBonders[b].atName1))) continue;

		if( allBonders[b].res1 != -1 )
		{
			if( at1->res != allBonders[b].res1 )
				continue;
		}
		
		if( allBonders[b].res2 != -1 )
		{
			if( at1->res != allBonders[b].res2 )
				continue;
		}

		bx_use = b;
	}	

	return bx_use;
}
















