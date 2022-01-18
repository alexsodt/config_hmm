//#include <OpenCL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "util.h"
//#include "proc_definition.h"

int matches_atomselect( const char *select, const char *target );

int main( int argc, char **argv )
{
	char *buffer = (char *)malloc( sizeof(char) * 1000000 );
	char *buffer2 = (char *)malloc( sizeof(char) * 1000000 );
		
		const char *command_string = 
		
		"Commands: list_modes list_lipids av_curvature mode_time_series\n"
		"\tlist_modes [q_cutoff=any]\n"
		"\tlist_lipids [lipidName]\n"
		"\tav_gauss [lipid_index=any] [lipid_name=any] [seg_name=any] [q_cutoff=any]\n"
		"\tav_curvature [lipid_index=any] [lipid_name=any] [seg_name=any] [q_cutoff=any] [start=t1] [stop=t2]\n"
		"\tmode_curvature [lipid_name=any] [atom_name=any] [definition.inp state_assignment] [n_select=any]\n"
		"\tmode_thick [lipid_name=any] [atom_name=any] [definition.inp state_assignment] [n_select=any]\n"
		"\tmode_curvature_err [lipid_name=any] [atom_name=any] [start=t1] [stop=t2] [definition.inp state_assignment] [n_select=any]\n"
		"\tmode_thick_err [lipid_name=any] [atom_name=any] [start=t1] [stop=t2] [definition.inp state_assignment] [n_select=any]\n"
		"\ttrack_mode mode_index [time_between_frames=unset]\n"
		"\tkc\n"
		"\ttrack_rho mode_index res_name [time_between_frames=unset]\n";

	if( argc < 4 )
	{
		printf("Syntax: fourierExtract.exe hq.txt ldata.txt COMMAND ...\n");


		printf("%s", command_string );
		return -1;
	}

	FILE *hqFile = fopen(argv[1],"r");
	if( !hqFile )
	{
		printf("Couldn't open file '%s'.\n", argv[1] );
		return -1;
	}

	FILE *ldataFile = fopen(argv[2],"r");
	if( !ldataFile )
	{
		printf("Couldn't open file '%s'.\n", argv[2] );
		return -1;
	}

	getLine( hqFile, buffer );

	int bilayer_mode = 0;
	int line_skip = 0;
	if( !strncasecmp(buffer, "mode=bilayer",strlen("mode=bilayer")))
		bilayer_mode = 1;
	else
		rewind(hqFile);

	const char *cmd = argv[3];

	struct lipid
	{
		char *resName;
		char *atomName;
		char *segid;
		int leaflet;
		int res;
		int pdef_index; // if a definition file is included, this links into the state sequence.
		double avc;
		double navc;
	};

	int nlipids = 0;
	int nlipidsSpace = 10;	

	lipid *lipidInfo = (lipid *)malloc( sizeof(lipid) * nlipidsSpace );

	getLine( ldataFile, buffer );

	char *parse = buffer + 1; // skip the comment line.
		
	while( *parse == ' ' || *parse == '\t' ) parse += 1;
	while( *parse != ' ' && *parse != '\t' ) parse += 1;

	int nsphingo = 0;

	while( *parse )
	{	
		while( *parse == ' ' || *parse == '\t' ) parse += 1;
		
		char resName[256];
		char atomName[256];
		int leaflet;

		int tx = 0;
		while( *parse && *parse != ',' ) { resName[tx] = *parse; resName[tx+1] = '\0'; parse += 1; tx++; }
		if( !*parse )
		{
			printf("LDATA read error.\n");
			exit(1);
		}
		parse += 1;
		tx = 0;
		while( *parse && *parse != ',' ) { atomName[tx] = *parse; atomName[tx+1] = '\0'; parse += 1; tx++; }
		if( !*parse )
		{
			printf("LDATA read error.\n");
			exit(1);
		}
		parse += 1;
		leaflet = atoi(parse);
		while( *parse && *parse != ',' ) parse+=1;
		parse += 1;
		int res = atoi(parse);
		while( *parse && *parse != ',' ) parse+=1;
		parse += 1;

		char segid[256];
		tx = 0;
		while( *parse && *parse != ' ' && *parse != '\t' ) { segid[tx] = *parse; segid[tx+1] = '\0'; parse += 1; tx++; }
		while( *parse && *parse != ' ' && *parse != '\t' ) parse+=1;
		while( *parse && *parse == ' ' || *parse == '\t' ) parse+=1;

//		printf("READ %s %s %d %d %s\n", resName, atomName, leaflet, res, segid );

		if( nlipids == nlipidsSpace )
		{
			nlipidsSpace *= 2;
			lipidInfo = (lipid *)realloc( lipidInfo, sizeof(lipid) * nlipidsSpace );
		}

		lipidInfo[nlipids].resName = (char*)malloc( sizeof(char)*(1+strlen(resName)) );
		strcpy( lipidInfo[nlipids].resName, resName );
		lipidInfo[nlipids].atomName = (char*)malloc( sizeof(char)*(1+strlen(atomName)) );
		strcpy( lipidInfo[nlipids].atomName, atomName );
		lipidInfo[nlipids].segid = (char*)malloc( sizeof(char)*(1+strlen(segid)) );
		strcpy( lipidInfo[nlipids].segid, segid );
		lipidInfo[nlipids].leaflet = leaflet;
		lipidInfo[nlipids].res = res;
		lipidInfo[nlipids].avc = 0;
		lipidInfo[nlipids].navc = 0;

		lipidInfo[nlipids].pdef_index = -1;

		nlipids++;
	}		
	
	int nmodes = 0;
	int nmodesSpace = 10;
	double *modes = (double *)malloc( sizeof(double) * 2 * nmodesSpace );
	getLine( hqFile, buffer );


	parse = buffer + 1; // skip the comment line.
		
	int N_QX, N_QY;
	while( *parse == ' ' || *parse == '\t' ) parse += 1;
	N_QX = atoi(parse);
	while( *parse != ' ' && *parse != '\t' ) parse += 1;
	while( *parse == ' ' || *parse == '\t' ) parse += 1;
	N_QY = atoi(parse);
	while( *parse != ' ' && *parse != '\t' ) parse += 1;
	while( *parse == ' ' || *parse == '\t' ) parse += 1;
	// DCD
	while( *parse != ' ' && *parse != '\t' ) parse += 1;
	while( *parse == ' ' || *parse == '\t' ) parse += 1;

	while( *parse )
	{	
		while( *parse == ' ' || *parse == '\t' ) parse += 1;
	
		double qx,qy;

		qx = atof(parse);
		while( *parse && *parse != ',' ) parse += 1;

		if( !*parse ) break;

		parse += 1; // comma
		qy = atof(parse);
		while( *parse && *parse != ' ' && *parse != '\t' ) parse += 1;
		while( *parse == ' ' || *parse == '\t' ) parse += 1;

		if( nmodes == nmodesSpace )
		{
			nmodesSpace *= 2;
			modes = (double *)realloc( modes, sizeof(double) * 2 * nmodesSpace );
		}

		modes[2*nmodes+0] = qx;
		modes[2*nmodes+1] = qy;

		nmodes++;
	}		
	

	double *mode_sq = (double *)malloc( sizeof(double) *  nmodes );			
	memset( mode_sq, 0, sizeof(double) *  nmodes );
	double *n_mode_sq = (double *)malloc( sizeof(double) * nmodes );
	memset( n_mode_sq, 0, sizeof(double) * nmodes );
	double *hq = (double *)malloc( sizeof(double) * 2 * N_QX * N_QY *2 );
	double *lxy = (double *)malloc( sizeof(double) * 3 * nlipids );

	int do_gauss = 0;
	int do_curv = 0;
	int do_thick = 0;
	int do_config = 0;
	int do_multistate = 0;
	int mode_curv = 0;
	int mode_thick = 0;
	int n_select = -1;
	int for_kc = 0;
	int track_rho = 0;
	int track_mode = 0;
	int lipid_select = -1;
	int mode_select = -1;
	char res_select[256];	
	res_select[0] = '\0';
	char atomname_select[256];	
	atomname_select[0] = '\0';
	char seg_select[256];	
	seg_select[0] = '\0';
	double q_cutoff = -1;
	double time_between_frames = 0;
	int start_frame=-1;
	int stop_frame=-1;

	FILE *pair_assignment = NULL;
	const char *pdef_fileName = NULL;

	if( !strcasecmp( cmd, "list_modes") ) 
	{

		if( argc > 4 && strcasecmp( argv[4], "any") )
			q_cutoff = atof(argv[4]);
		printf("MODES\n");
		printf("------------\n");
		for( int m = 0; m < nmodes; m++ )
		{
			double q = sqrt(modes[2*m+0]*modes[2*m+0] + modes[2*m+1]*modes[2*m+1]);

			if( q_cutoff < 0 || q < q_cutoff )
				printf("index %d qx %4.5lf qy %4.5lf qtot %4.5lf\n", m, modes[2*m+0], modes[2*m+1], q );
		}
		
		return -1;
	}
	else if( !strcasecmp( cmd, "list_lipids") ) 
	{
		printf("LIPIDS\n");
		printf("------------\n");
		for( int l = 0; l < nlipids; l++ )
			printf("index %d segid %s resname %s resid %d atom %s leaflet: %s\n", l, lipidInfo[l].segid, lipidInfo[l].resName, lipidInfo[l].res, lipidInfo[l].atomName, (lipidInfo[l].leaflet < 0 ? "lower" : "upper" ) );	

		return -1;
	}
	else if( !strcasecmp( cmd, "mode_curvature") || !strcasecmp( cmd, "config_curvature") || !strcasecmp( cmd, "multistate_curvature") || !strcasecmp( cmd, "mode_curvature_err")  
		|| !strcasecmp( cmd, "mode_thick") || !strcasecmp( cmd, "config_thick") || !strcasecmp( cmd, "multistate_thick") || !strcasecmp( cmd, "mode_thick_err")  ) 

	{
		if( !strcasecmp( cmd, "config_curvature") ) 
			do_config = 1;	
		if( !strcasecmp( cmd, "multistate_curvature") ) 
			do_multistate = 1;	

		
		if( !strcasecmp( cmd, "mode_thick") || !strcasecmp( cmd, "config_thick") || !strcasecmp( cmd, "multistate_thick") || !strcasecmp( cmd, "mode_thick_err")  )
			do_thick = 1; 

		mode_curv = 1;
		do_curv = 1;	
		//"\tav_curvature [lipid_index=any] [lipid_name=any] [seg_name=any] [q_cutoff=any]\n"

		// lipid_index
		if( argc > 4 && strcasecmp( argv[4], "any") )
		{
			if( strlen(argv[4]) < 200 )
				strcpy( res_select, argv[4] );
		}

		if( argc > 5 && strcasecmp( argv[5], "any") )
		{
			if( strlen(argv[5]) < 200 )
				strcpy( atomname_select, argv[5] );
		}

		if( !strcasecmp( cmd, "mode_curvature_err" ) || !strcasecmp( cmd, "mode_thick_err" ))
		{
			if( argc > 6 )
				start_frame = atoi(argv[6]);
			if( argc > 7 )
				stop_frame = atoi(argv[7]);	
			
			if( argc > 9 )
			{
				pdef_fileName = argv[8];

				pair_assignment = fopen(argv[9], "r");
	
				if( !pair_assignment ) 
				{
					printf("Couldn't open sphingolipid assignment file '%s'.\n", argv[9]);
					exit(1);
				}
			}
			if( argc > 10 && strcasecmp( argv[10], "any") )
			{
				n_select = atoi(argv[10]);		
			}
		} 
		else if( argc > 6 )
		{
			if( argc < 8 )
			{
				printf("State-assignment breakdown of curvature requires a definition.inp and state-assignment file.\n");
				exit(1);
			}

			pdef_fileName = argv[6];


			pair_assignment = fopen(argv[7], "r");

			if( !pair_assignment ) 
			{
				printf("Couldn't open state assignment file '%s'.\n", argv[7]);
				exit(1);
			}
			if( argc > 8 && strcasecmp( argv[8], "any") )
			{
				n_select = atoi(argv[8]);		
			}
			if( argc > 9 )
				start_frame = atoi(argv[9]); 
			if( argc > 10 )
				stop_frame = atoi(argv[10]); 
		}	

	}
	else if( !strcasecmp( cmd, "av_curvature") || !strcasecmp( cmd, "av_gauss") ) 
	{	
		if( !strcasecmp( cmd, "av_gauss" ) )
			do_gauss = 1;
		
		if( !strcasecmp( cmd, "av_thick" ) )
			do_thick = 1;

		do_curv = 1;	
		//"\tav_curvature [lipid_index=any] [lipid_name=any] [seg_name=any] [q_cutoff=any]\n"

		// lipid_index
		if( argc > 4 && strcasecmp( argv[4], "any") )
		{
			if( strlen(argv[4]) < 200 )
				strcpy( atomname_select, argv[4] );
		}		

//		if( argc > 4 && strcasecmp( argv[4], "any") )
//			lipid_select = atoi(argv[4]);	
		if( argc > 5 && strcasecmp( argv[5], "any") )
		{
			if( strlen(argv[5]) < 200 )
				strcpy( res_select, argv[5] );
		}		
		if( argc > 6 && strcasecmp( argv[6], "any") )
		{
			if( strlen(argv[6]) < 200 )
				strcpy( seg_select, argv[6] );
		}
		if( argc > 7 && strcasecmp( argv[7], "any" ) )
		{
			q_cutoff = atof(argv[7]);	

			if( q_cutoff >= 1 || q_cutoff <= 1e-10 )
			{
				mode_select = atoi(argv[7]);
				q_cutoff = -1;
			}
		}
		
		if( argc > 8 )
			start_frame = atoi(argv[8]); 
		if( argc > 9 )
			stop_frame = atoi(argv[9]); 
				
	}
	else if( !strcasecmp( cmd, "track_mode") ) 
	{
		track_mode = 1;

		if( argc < 5 )
		{
			printf("Track mode requires an index to track.\n");
			printf("%s", command_string );
			return -1;
		}
		
		mode_select = atoi(argv[4]);
	
		if( mode_select < 0 || mode_select >= nmodes )
		{
			printf("Invalid mode index.\n");
			return -1;
		}

		if( argc >= 6 )
		{
			time_between_frames = atof(argv[5]);
		}
	}
	else if( !strcasecmp( cmd, "track_rho" ) )
	{
		track_rho = 1;
		
		if( argc < 6 )
		{
			printf("Track rho requires an mode index and residue to track.\n");
			printf("%s", command_string );
			return -1;
		}
		
		mode_select = atoi(argv[4]);
		
		if( strlen(argv[5]) < 200 )
			strcpy( res_select, argv[5] );
		else
		{
			printf("Input too long.\n");
			printf("%s", command_string );
			return -1;
		}
		
		if( argc >= 7 )
		{
			time_between_frames = atof(argv[6]);
		}
	}
	else if( !strcasecmp( cmd, "kc" ) )
	{
		for_kc = 1;
	}
	else
	{
		printf("Unknown command '%s'.\n", cmd );
		printf("%s", command_string );
		return -1;
	}
	
#if 0		
	if( pdef_fileName )
	{
		load_pair_definition( pdef_fileName );

		int has_psf = loadPDefPSF() == 1 ;

		if( has_psf )
		{
			for( int l = 0; l < nlipids; l++ )
			{
				lipidInfo[l].pdef_index = getPDefIndex( lipidInfo[l].atomName, lipidInfo[l].resName, lipidInfo[l].segid, lipidInfo[l].res );
	
			}	
		}
		else
		{ 
			int index_offset = 0;
			for( int pass = 0; pass < 2; pass++ )
			{
				if( is_sym() && pass > 0 ) break;
		
				for( int l = 0; l < nlipids; l++ )
				{
					if( !(atomname_select[0]) || matches_atomselect( atomname_select, lipidInfo[l].atomName) ) // || !strcasecmp( lipidInfo[l].atomName, atomname_select) )
					{
						if( use_atom( USE_ANYWHERE, pass, lipidInfo[l].res, lipidInfo[l].atomName, lipidInfo[l].resName, lipidInfo[l].segid ) )
						{
							lipidInfo[l].pdef_index = index_offset;
							index_offset++;
						}
					}
				}
			}
		}
	}
#endif

	int npasses = 1;

//	if( pair_assignment )
//		npasses = 2;

	int nstateSpace = 1;
	int nstates = 1;
	int max_time = 20000;
	int cur_time = 0;

	// plan to to A-Z,a-z+0	
	nstateSpace = 53;
	nstates = 53;
	int max_nstates = 1; // the states we actually get, could be less than 53 (which is the space for configs, A-Z,a-z

	if( do_config) max_nstates = 53;
	if( do_multistate ) max_nstates = 11; // monomer, 0-9

	double *tot_av_c_m = (double *)malloc( sizeof(double) * nstateSpace * nmodes  );
	double *tot_av_c_t_m = (double *)malloc( sizeof(double) * nstateSpace * max_time * nmodes  );
	double *ntot_t_m = (double *)malloc( sizeof(double) * nstateSpace * max_time * nmodes  );
	double *ntot_m = (double *)malloc( sizeof(double) * nstateSpace * nmodes  );
	memset( tot_av_c_m, 0, sizeof(double) * nstateSpace * nmodes );
	memset( ntot_m, 0, sizeof(double) * nstateSpace * nmodes );
	memset( tot_av_c_t_m, 0, sizeof(double) * nstateSpace * max_time * nmodes );
	memset( ntot_t_m, 0, sizeof(double) * nstateSpace * max_time * nmodes );

	double *tot_av_c_t = (double *)malloc( sizeof(double) * nstateSpace * max_time  );
	double *ntot_t = (double *)malloc( sizeof(double) * nstateSpace * max_time  );
	double *tot_av_c = (double *)malloc( sizeof(double) * nstateSpace  );
	double *ntot = (double *)malloc( sizeof(double) * nstateSpace  );
	memset( tot_av_c, 0, sizeof(double) * nstateSpace );
	memset( ntot, 0, sizeof(double) * nstateSpace );
	memset( tot_av_c_t, 0, sizeof(double) * nstateSpace * max_time );
	memset( ntot_t, 0, sizeof(double) * nstateSpace * max_time );

	int nmodes_used = 0;

	char *used_string = (char *)malloc( sizeof(char) * 1024 );
	int atoms_use[30];
	int atoms_use_len[30];
	int natom_use=0;
	int use_init = 0;

	double avK=0,navK=0;
	double avJ=0,navJ=0;
	int g_nmodes_used = 0;
	sprintf( used_string, "#Used: ");

	for( int pass = 0; pass < npasses; pass++ )
	{
		rewind(hqFile);

		int abort = 0;

		// re-read initial line for possible second pass.
		getLine(hqFile, buffer );
		
		if( !strncasecmp(buffer, "mode", 4 ) )
			getLine(hqFile, buffer );
			 
		
		rewind(ldataFile);
		// re-read initial line for possible second pass.
		getLine(ldataFile, buffer );
			
		if( pair_assignment )
			rewind(pair_assignment);

		int frame = 0;
		while( !feof(hqFile) && !abort )
		{
			getLine( hqFile, buffer );
			if( feof(hqFile) ) break;
			parse = buffer + 1;
			int mode = 0;
			while( *parse == ' ' || *parse == '\t' ) parse += 1;
			double La = atof(parse);
			while( *parse && *parse != ' ' && *parse != '\t' ) parse += 1;
			while( *parse == ' ' || *parse == '\t' ) parse += 1;
			double Lb = atof(parse);
			while( *parse && *parse != ' ' && *parse != '\t' ) parse += 1;
			while( *parse == ' ' || *parse == '\t' ) parse += 1;
			double LC = atof(parse);
			while( *parse && *parse != ' ' && *parse != '\t' ) parse += 1;
			while( *parse == ' ' || *parse == '\t' ) parse += 1;
		
			while( *parse && mode < nmodes )
			{
				double hqx,hqy,l_hqx,l_hqy;
	
				hqx = atof(parse);
				while( *parse && *parse != ',' ) parse += 1;
				parse += 1;
				hqy = atof(parse);		
	
				if( !bilayer_mode )
				{
					while( *parse && *parse != ',' ) parse += 1;
					parse += 1;
					l_hqx = atof(parse);		
					while( *parse && *parse != ',' ) parse += 1;
					parse += 1;
					l_hqy = atof(parse);		
				}
				while( *parse && *parse != ' ' && *parse != '\t' ) parse += 1;
				while( *parse == ' ' || *parse == '\t' ) parse += 1;
					
				mode_sq[mode] += hqx*hqx+hqy*hqy; 
				n_mode_sq[mode] += 2;

				if( bilayer_mode )
				{
					hq[2*mode+0] = hqx;
					hq[2*mode+1] = hqy;

				}
				else
				{
					hq[4*mode+0] = hqx;
					hq[4*mode+1] = hqy;
					hq[4*mode+2] = l_hqx;
					hq[4*mode+3] = l_hqy;
				}
				mode++;	
			} 
			
			getLine( ldataFile, buffer );
			if( feof(ldataFile) ) break;
			parse = buffer + 1;

			int lip = 0;
	
			while( *parse && lip < nlipids )
			{
				double x,y,z;
	
				x = atof(parse);
				while( *parse && *parse != ',' ) parse += 1;
				parse += 1;
				y = atof(parse);		
				while( *parse && *parse != ',' ) parse += 1;
				parse += 1;
				z = atof(parse);
				while( *parse && *parse != ' ' && *parse != '\t' && *parse != ',' ) parse += 1;
				if( *parse == ',' )
				{
					parse += 1;
					int new_leaf = 1;
					if( *parse == '-' )
						new_leaf = -1;
					parse += 1;
					while( *parse && *parse != ' ' && *parse != '\t' ) parse += 1;

					// V2 format. re-assign leaflet.
			
					if( new_leaf != lipidInfo[lip].leaflet  )
					{
//						printf("%s_%d flipped at frame=%d.\n", lipidInfo[lip].resName, lipidInfo[lip].res, frame );
					}
					lipidInfo[lip].leaflet = new_leaf;
		
				}
				while( *parse == ' ' || *parse == '\t' ) parse += 1;
				
				lxy[3*lip+0] = x;
				lxy[3*lip+1] = y;
				lxy[3*lip+2] = z;
	
				

				lip++;	
			} 
				
			double rho_q_i_0 = 0;	
			double rho_q_i_1 = 0;	
			double rho_q_o_0 = 0;	
			double rho_q_o_1 = 0;	
	

			if( do_curv || track_rho )
			{
				if( cur_time >= max_time )
				{
					printf("ERROR: increase max_time (%d).\n", max_time );
					exit(1);
				}
				if( pair_assignment )
				{
					getLine( pair_assignment, buffer2 );
					if( feof(pair_assignment ) ) 
					{
						abort = 1;
						break;
					}
				}
			}
			
			if( start_frame >= 0 && frame < start_frame ) 	
			{
				frame++;
				continue;
			}
			if( stop_frame >= 0 && frame >= stop_frame ) 	
			{
				frame++;
				continue;
			}
			if( do_curv || track_rho )
			{
				double avc = 0;
				double navc = 0;
	
				double leaflet_n[2] = {0,0};
	
				for( int lip = 0; lip < nlipids; lip++ )
				{
					if( lipidInfo[lip].leaflet > 0 )
						leaflet_n[0] += 1;
					else
						leaflet_n[1] += 1;
				}
	
				double av_dz = 0;
				double av_dz_c = 0;
				double nav_dz = 0;
				
				int max_simul=10;
				int nhmm[nlipids];
				int cur_obs[nlipids*max_simul];
				memset( nhmm, 0, sizeof(int) * nlipids );
				memset( cur_obs, 0, sizeof(int) * nlipids * max_simul );
				
				if( pair_assignment )
				{
					if( do_config || do_multistate  )
					{
						char *tread = buffer2;	
	
						int  curp = 0;
						while( *tread )
						{
							int nconf = *tread - '0';
	
							tread += 1;
	
							if( curp >= nlipids )
							{
								printf("There are too many lipids listed in the state file. Expected no more than %d.\n",
										nlipids );
								exit(1);
							}

							if( nconf == 0 )
							{ // if there are no similarity centers listed, it's a monomer; assign to state 0.
								cur_obs[curp*max_simul+nhmm[curp]] = 0;
								nhmm[curp] = 1;
							}
							else
							{
								for( int x = 0; x < nconf; x++ )
								{	
									char obs = *tread;
									int iobs = 0;

									if( do_config )
									{
										if( obs >= 'A' && obs <= 'Z' )
											iobs = 1 + obs - 'A'; // state 0 is the monomer, everything else are dimer configs.
										else if( obs >= 'a' && obs <= 'z' )
											iobs = 1 + 26 + obs - 'a';
										else
										{
											printf("Similarity center read error.\n");
											exit(1);
										}
									}
									else if( do_multistate )
									{
										if( obs >= '6' )
										{
											printf("Check this.\n");
										}

										if( obs >= '0' && obs <= '9' )
											iobs = 1 + obs - '0'; // state 0 is the monomer, everything else are dimer configs.
										else
										{
											printf("Similarity center read error.\n");
											exit(1);
										}
									}
									// we keep a record of each state that the lipid is in.
									// max_simul is the maximum number of simultaneous states a lipid can be involved in.
									cur_obs[curp*max_simul+nhmm[curp]] = iobs;
									// increment the number of states we have stored here.
									nhmm[curp]++;
									tread += 1;
								}
							}
			
							curp++;
						}
					}
					else		
					{	
						char *tread = buffer2;	
	
						int  curp = 0;
						while( *tread )
						{
							char obs = *tread;
							int nconf = 1;
							int iobs = 0;

							switch ( obs )
							{
								case 'X':
									iobs = 0;
									break;
								case '0':
									iobs = 1;
									break;
								case '1':
								case '2':
								case '3':
								case '4':
								case '5':
								case '6':
								case '7':
								case '8':
								case '9':
									iobs = 2 + obs - '1';
									break;
								default:
									printf("Cannot interpret state code '%c'.\n", obs );
									exit(1);
									break;
							}
							
							if( curp >= nlipids )
							{
								printf("There are too many lipids listed in the state file. Expected no more than %d.\n",
										nlipids );
								exit(1);
							}

	
							if( iobs >= max_nstates ) {
								max_nstates = iobs+1;
							}
							cur_obs[curp*max_simul+nhmm[curp]] = iobs;
							nhmm[curp] = 1;

							tread += 1;
							curp++;
						}
					}
				}
					
				double my_use_hq[2] = {0,0};

	

				for( int lip = 0; lip < nlipids; lip++ )
				{
					int cur_state = 0;
					int n_states_to_write = 1;
					int states[1+max_simul];
					states[0] = cur_state;

					if( pair_assignment )
					{
						int indx;

						indx = lipidInfo[lip].pdef_index;
						if( indx < 0 ) continue;

						// put the current state(s) into the state array.
						n_states_to_write = nhmm[indx];
						memcpy( states, cur_obs + indx* max_simul, sizeof(int) * nhmm[indx] ); 
					}
					
					if( seg_select[0] && strcasecmp( lipidInfo[lip].segid, seg_select) ) continue;	
					if( res_select[0] && strcasecmp( lipidInfo[lip].resName, res_select) ) continue;	
//					if( atomname_select[0] && strcasecmp( lipidInfo[lip].atomName, atomname_select) ) continue;	
					if( atomname_select[0] && !matches_atomselect( atomname_select, lipidInfo[lip].atomName) )
						continue;
	
 // || !strcasecmp( lipidInfo[l].atomName, atomname_select) )
					if( lipid_select >= 0 && lip != lipid_select ) continue;
						

					double x = lxy[3*lip+0];
					double y = lxy[3*lip+1];
					double z = lxy[3*lip+2];
	
					double lipid_c = 0;
					double lipid_c_0 = 0;
					double lipid_c_1 = 0;
					nmodes_used = 0;
				
					if( !use_init )
					{
						char used[256];
						sprintf(used, "%s:%s", lipidInfo[lip].resName, lipidInfo[lip].atomName );
						int repeated = 0;

						for( int ax = 0; ax < natom_use; ax++ )
						{
							if( strlen(used) == atoms_use_len[ax] && !strncasecmp( used, used_string + atoms_use[ax], strlen(used) ) )
								repeated = 1; 
						}

						if( natom_use < 29 && !repeated && strlen(used) + strlen(used_string) < 1023 )
						{
							atoms_use[natom_use] = strlen(used_string);
							atoms_use_len[natom_use] = strlen(used);
							strcpy( used_string+strlen(used_string), used );
							strcpy( used_string+strlen(used_string), " " );
							natom_use++;
						}
					}
	
					if( do_gauss )
					{
						int nmodes_used = 0;
						double l_avK = 0;

						double c_mat[4] = {0,0,0,0};

						for( int mode1 = 0; mode1 < nmodes; mode1++ )
						{
							double q1x = modes[2*mode1+0];
							double q1y = modes[2*mode1+1];
							double q1 = sqrt(q1x*q1x+q1y*q1y);
							

							if( fabs(q1) < 1e-10 ) continue;
							if( q_cutoff >= 0 && q1 > q_cutoff ) continue;
		
							double l_ntot = leaflet_n[0] + leaflet_n[1];
							double use_hq1[2] = {0,0};
						
							if( bilayer_mode )
							{
								use_hq1[0] = hq[2*mode1+0];
								use_hq1[1] = hq[2*mode1+1];
							}
							else
							{
								if( lipidInfo[lip].leaflet > 0 )
								{
									use_hq1[0] = hq[4*mode1+0];
									use_hq1[1] = hq[4*mode1+1];
								}
								else
								{
									use_hq1[0] = hq[4*mode1+2];
									use_hq1[1] = hq[4*mode1+3];
								}
							}
		
							double c_x = 0, c_y = 0;


							c_mat[0] += lipidInfo[lip].leaflet * q1x*q1x * cos( q1x * x + q1y * y ) * use_hq1[0]; // xx
							c_mat[1] += lipidInfo[lip].leaflet * q1x*q1y * cos( q1x * x + q1y * y ) * use_hq1[0]; // xy
							c_mat[2] += lipidInfo[lip].leaflet * q1x*q1y * cos( q1x * x + q1y * y ) * use_hq1[0]; // yx	 
							c_mat[3] += lipidInfo[lip].leaflet * q1y*q1y * cos( q1x * x + q1y * y ) * use_hq1[0]; // yy
							
							c_mat[0] += -lipidInfo[lip].leaflet * q1x*q1x * sin( q1x * x + q1y * y ) * use_hq1[1]; // xx
							c_mat[1] += -lipidInfo[lip].leaflet * q1x*q1y * sin( q1x * x + q1y * y ) * use_hq1[1];	
							c_mat[2] += -lipidInfo[lip].leaflet * q1x*q1y * sin( q1x * x + q1y * y ) * use_hq1[1];	
							c_mat[3] += -lipidInfo[lip].leaflet * q1y*q1y * sin( q1x * x + q1y * y ) * use_hq1[1];	
							
							nmodes_used += 1;
						}

						double cxx = c_mat[0];
						double cxy = c_mat[1];
						double cyx = c_mat[2];
						double cyy = c_mat[3];

						double c_1 = 0.5 * (cxx+cyy-sqrt(cxx*cxx+4*cxy*cxy-2*cxx*cyy+cyy*cyy));							
						double c_2 = 0.5 * (cxx+cyy+sqrt(cxx*cxx+4*cxy*cxy-2*cxx*cyy+cyy*cyy));							

						avK += c_1*c_2;	
						avJ += c_1+c_2;	
						if( g_nmodes_used == 0 ) g_nmodes_used = nmodes_used;
						navK += 1;
						navJ += 1;
					}
					else
					{
						for( int mode = 0; mode < nmodes; mode++ )
						{
							double qx = modes[2*mode+0];
							double qy = modes[2*mode+1];
							double q = sqrt(qx*qx+qy*qy);
	
							if( track_rho && mode != mode_select ) continue;
		
							if( fabs(q) < 1e-10 ) continue;
							if( q_cutoff >= 0 && q > q_cutoff) continue;
							if( mode_select != -1 && mode != mode_select ) continue;	
							nmodes_used++;				
		
							double l_ntot = leaflet_n[0] + leaflet_n[1];
							double use_hq[2] = {0,0};
						
							if( bilayer_mode )
							{
								use_hq[0] = hq[2*mode+0];
								use_hq[1] = hq[2*mode+1];
							}
							else
							{
								if( lipidInfo[lip].leaflet > 0 )
								{
									use_hq[0] = hq[4*mode+0];
									use_hq[1] = hq[4*mode+1];
								}
								else
								{
									use_hq[0] = hq[4*mode+2];
									use_hq[1] = hq[4*mode+3];
								}
							}
		
							if( fabs(qx-0.029920) < 1e-5 && fabs(qy) < 1e-8 )
							{
								my_use_hq[0] = use_hq[0];	
								my_use_hq[1] = use_hq[1];	
							}
	
							if( lipidInfo[lip].leaflet < 0 )
							{
								rho_q_i_0 += cos(qx * x + qy * y) / (La*Lb); 
								rho_q_i_1 -= sin(qx * x + qy * y) / (La*Lb); 
							}
							else
							{
								rho_q_o_0 += cos(qx * x + qy * y) / (La*Lb); 
								rho_q_o_1 -= sin(qx * x + qy * y) / (La*Lb); 
							}
							double cr_c=0,cr_s=0;
	
							if( do_thick )
							{
								if( lipidInfo[lip].leaflet > 0 )
								{
									cr_c   =  cos( qx * x + qy * y ) * use_hq[0]; 
									cr_c  +=  cos( qx * x + qy * y ) * use_hq[0]; 
									cr_s   =- sin( qx * x + qy * y ) * use_hq[1];
									cr_s  +=- sin( qx * x + qy * y ) * use_hq[1];
		
								}
								else
								{
									cr_c   = cos( qx * x + qy * y ) * use_hq[0]; 
									cr_c  += cos( qx * x + qy * y ) * use_hq[0]; 
									cr_s   =-sin( qx * x + qy * y ) * use_hq[1];
									cr_s  +=-sin( qx * x + qy * y ) * use_hq[1];
								}
							}
							else
							{	
								if( lipidInfo[lip].leaflet > 0 )
								{
									cr_c   = lipidInfo[lip].leaflet * qx*qx * cos( qx * x + qy * y ) * use_hq[0]; 
									cr_c  += lipidInfo[lip].leaflet * qy*qy * cos( qx * x + qy * y ) * use_hq[0]; 
									cr_s   =-lipidInfo[lip].leaflet * qx*qx * sin( qx * x + qy * y ) * use_hq[1];
									cr_s  +=-lipidInfo[lip].leaflet * qy*qy * sin( qx * x + qy * y ) * use_hq[1];
		
								}
								else
								{
									cr_c   = lipidInfo[lip].leaflet * qx*qx * cos( qx * x + qy * y ) * use_hq[0]; 
									cr_c  += lipidInfo[lip].leaflet * qy*qy * cos( qx * x + qy * y ) * use_hq[0]; 
									cr_s   =-lipidInfo[lip].leaflet * qx*qx * sin( qx * x + qy * y ) * use_hq[1];
									cr_s  +=-lipidInfo[lip].leaflet * qy*qy * sin( qx * x + qy * y ) * use_hq[1];
								}
							}
						
							if( fabs(qy) < 1e-9 ) {
	
							
							if( fabs(qx+2.992000e-02) < 1e-5 ) 
							{
							}
	
							}
	
							
	//						printf("x: %le y: %le cr_c: %le cr_s: %le qx %le qy %le\n", x, y, cr_c, cr_s, qx, qy );
		
							av_dz_c += lipidInfo[lip].leaflet * qx*qx *cos( qx * x + qy * y ) * use_hq[0];
							av_dz_c += lipidInfo[lip].leaflet * qx*qy *cos( qx * x + qy * y ) * use_hq[0];
							av_dz_c += -lipidInfo[lip].leaflet * qx*qx *sin( qx * x + qy * y ) * use_hq[1];
							av_dz_c += -lipidInfo[lip].leaflet * qy*qy *sin( qx * x + qy * y ) * use_hq[1];
	//						av_dz_c += lipidInfo[lip].leaflet * sin( qx * x + qy * y ) * use_hq[1];
							av_dz   += z;
							nav_dz  += 1;
	
							lipid_c_0 += cr_c;
							lipid_c_1 += cr_s;
	
	
							lipid_c += cr_c+cr_s;
							lipidInfo[lip].avc += cr_c+cr_s;
				
	
							if( n_select == -1 || n_select == n_states_to_write )
							{
								for( int sx = 0; sx < n_states_to_write; sx++ )
								{	
									int cur_state_write = states[sx];
		
									tot_av_c_t_m[cur_state_write*max_time*nmodes+mode*max_time+cur_time] += cr_c + cr_s;
									tot_av_c_m[cur_state_write*nmodes+mode] += cr_c + cr_s;
		
									ntot_m[cur_state_write*nmodes+mode] += 1;
									ntot_t_m[cur_state_write*nmodes*max_time+mode*max_time+cur_time] += 1;
		
									tot_av_c_t[cur_state_write*max_time+cur_time] += cr_c + cr_s;
									tot_av_c[cur_state_write] += cr_c + cr_s;
								}
							}
	
							avc += cr_c + cr_s;
	//						if( frame == 225 && (!strcasecmp(lipidInfo[lip].segid, "GLPA47") || !strcasecmp(lipidInfo[lip].segid, "GLPA49")) )
	//							printf("cr: %le cr_s: %le\n", cr_c, cr_s );
						}
					}
					
//					printf("%d %s %lf %lf c %lf %lf\n", lipidInfo[lip].leaflet, lipidInfo[lip].resName, x, y, lipid_c_0, lipid_c_1 );
			

//					if( seg_select[0] )
//						printf("lipid %s frame %d c %le\n", lipidInfo[lip].segid, frame, lipid_c); 
	
					lipidInfo[lip].navc += 1;
			
					if( n_select == -1 || n_select == n_states_to_write )
					{
						for( int sx = 0; sx < n_states_to_write; sx++ )
						{	
							int cur_state_write = states[sx];
							ntot[cur_state_write] += 1;	
							ntot_t[cur_state_write*max_time+cur_time] += 1;
						}
					}
					navc += 1;
				}

				if( !use_init )
					printf("%s\n", used_string );

				use_init = 1;

//				printf("hq: %le %le avc: %le\n", my_use_hq[0], my_use_hq[1], avc / navc );

				//printf("av_dz: %le av_dz_c: %le\n", av_dz/nav_dz, av_dz_c/nav_dz );
//				tot_av_c += avc / navc;
//				ntot += 1;
	
//				printf("frame %d %le %le\n",  frame, avc / navc, navc );

//				printf("blip %lf %lf.\n", hq[2*mode_select+0], hq[2*mode_select+1] );
				cur_time++;

			}
		
			if( track_mode )
			{
				if( time_between_frames > 0 )
				{
					printf("%le ", time_between_frames * frame );
				}

				if( bilayer_mode )
					printf("%lf %lf\n", hq[2*mode_select+0], hq[2*mode_select+1]);
				else
					printf("%lf %lf %lf %lf\n", hq[4*mode_select+0], hq[4*mode_select+1], hq[4*mode_select+2], hq[4*mode_select+3] );
			}	
			if( track_rho )
			{
				printf("%le %le %le %le %le (ic,is,oc,os)\n", time_between_frames * frame, rho_q_i_0, rho_q_i_1,rho_q_o_0,rho_q_o_1 );
			}	
			frame++;
		}
	
		if( do_curv )
		{
			for( int l = 0; l < nlipids; l++ )
			{
//				if( lipidInfo[l].navc > 0 )
//					printf("%s %s %s %le\n", lipidInfo[l].segid, lipidInfo[l].resName, lipidInfo[l].atomName, lipidInfo[l].avc/lipidInfo[l].navc );
			}
			
			if( do_gauss )
			{
				printf("<K>: %le <J>: %le (%d modes used).\n", avK/navK, avJ/navJ, g_nmodes_used);
			}
		}
			
		int nstates = max_nstates;

		if( mode_curv )
		{
			int state_lim = 0;
			
			for( int s = 0; s < nstates; s++ )
			for( int t = 0; t < max_time; t++ )
			{
				if( ntot_t_m[s*max_time*nmodes+t] >= 0  )
					state_lim = s+1;
			}

			if( do_config )
				printf("nmodes %d nconfigs %d num_per_state", nmodes, state_lim );
			else
				printf("nmodes %d nstates %d num_per_state", nmodes, state_lim );

			for( int s = 0; s < state_lim; s++ )
				printf(" %lf", ntot[s] );
			printf("\n");

			printf("mode q vals");
			for( int m = 0; m < nmodes; m++ )
				printf(" %le", sqrt(modes[2*m+0]*modes[2*m+0]+modes[2*m+1]*modes[2*m+1]) );
			printf("\n");

			for( int s = 0; s < state_lim; s++ )
			{
				int use = 0;
				for( int t = 0; t < max_time; t++ )
				{
					if( ntot_t_m[s*max_time*nmodes+t] > 0  )
						use = 1;
				}
				if( use  )
				{
					if( do_config )
					{
						char code = '0'; // 0, A-Z, a-z are the code sequences.
						if( s > 0 && s <= 26 )
							code = 'A' + s - 1;
						else if( s > 26 && s <= 52 )
							code = 'a' + s - 27;
						printf("state %c curv", code );
					}
					else if( do_multistate )
					{
						char code = 'X'; // 0, A-Z, a-z are the code sequences.
						if( s > 0 && s <= 9 )
							code = '0' + s - 1;
						printf("state %c curv", code );
					}
					else
						printf("state %d curv", s );
					
					for( int m = 0; m < nmodes; m++ )
						printf(" %le", tot_av_c_m[s*nmodes+m] / ntot_m[s*nmodes+m] );
	
					printf("\n");
				}
				else
				{
					printf("state %d curv", s );

					for( int m = 0; m < nmodes; m++ )
						printf(" %le", 0.0  );
	
					printf("\n");
					
				}	
			}

			for( int s = 0; s < state_lim; s++ )
			{
				int use = 0;
				for( int t = 0; t < max_time; t++ )
				{
					if( ntot_t_m[s*max_time*nmodes+t] >= 0  )
						use = 1;
				}
				if( use  )
				{
					if( do_config )
					{
						char code = '0'; // 0, A-Z, a-z are the code sequences.
						if( s > 0 && s <= 26 )
							code = 'A' + s - 1;
						else if( s > 26 && s <= 52 )
							code = 'a' + s - 27;
						printf("state %c ebar", code );
					}
					else if( do_multistate)
					{
						char code = 'X'; // 0, A-Z, a-z are the code sequences.
						if( s > 0 && s <= 9 )
							code = '0' + s - 1;
						printf("state %c ebar", code );
					}
					else
						printf("state %d ebar", s );
					for( int m = 0; m < nmodes; m++ )
					{
						int nebar = 7;
						
						double lavc = 0;
						double lavc2 = 0;
						double navc = 0;
		
						for( int t = 0; t < nebar; t++ )
						{
							double t_avc = 0;
							double t_navc = 0;
		
							for( int it = (t*cur_time/7); it < ((t+1)*cur_time/7) && it < cur_time; it++ )
							{
								t_avc += tot_av_c_t_m[s*max_time*nmodes+m*max_time+it];
								t_navc += ntot_t_m[s*max_time*nmodes+m*max_time+it]; 
							}
			
							if( t_navc > 50 )
							{
								t_avc /= t_navc;
		
								lavc += t_avc;
								lavc2 += t_avc*t_avc;
								navc += 1;
							}
						}
					
						double ebar = sqrt((lavc2 / navc) - (lavc/navc)*(lavc/navc));
						ebar /= sqrt(navc);
						printf(" %le", ebar );
					}	
					printf("\n");
				}
				else
				{
					printf("state %d ebar", s );
					for( int m = 0; m < nmodes; m++ )
						printf("%le ", 0.0 );
					printf("\n");
				}
			}
		}
		else if( do_curv  && ! do_gauss )	
		{

			double ebar[nstates];
			
			int nebar = 7;
			for( int s = 0; s < nstates; s++ )
			{
				double lavc = 0;
				double lavc2 = 0;
				double navc = 0;
				for( int t = 0; t < nebar; t++ )
				{
					double t_avc = 0;
					double t_navc = 0;
					for( int it = (t*cur_time/7); it < ((t+1)*cur_time/7) && it < cur_time; it++ )
					{
						t_avc += tot_av_c_t[s*max_time+it];
						t_navc += ntot_t[s*max_time+it]; 
					}

	
					if( t_navc > 50 )
					{
						t_avc /= t_navc;

						lavc += t_avc;
						lavc2 += t_avc*t_avc;
						navc += 1;
					}
				}

				ebar[s] = sqrt((lavc2 / navc) - (lavc/navc)*(lavc/navc));
				ebar[s] /= sqrt(navc);
			}

			if( nstates <= 1 )
				printf("<c>: %.14le nmodes %d\n", tot_av_c[0]/ntot[0], nmodes_used);
			else
			{
				for( int s = 0; s < nstates; s++ )
					printf("c_%c : %le ntot: %lf ebar %le nmodes %d\n", (s == 0 ? 'm' : '0' + s-1), tot_av_c[s] / ntot[s], ntot[s], ebar[s], nmodes_used  );
			}
		}
			//printf("<c>: %le %le\n", tot_av_c/ntot, av_dz_c / nav_dz  );
	}

	double qmode[nmodes];
	for( int i = 0; i < nmodes; i++ )
	{
		qmode[i] = sqrt(modes[2*i+0]*modes[2*i+0]+modes[2*i+1]*modes[2*i+1]);
	}
	int sorter[nmodes];
	for( int i = 0; i < nmodes; i++ )
		sorter[i]=i;
	int done = 0;
	while(!done)
	{
		done = 1;	
		for( int i = 0; i < nmodes-1;i++ )
		{
			if( qmode[sorter[i+1]] < qmode[sorter[i]] )
			{
				done = 0;
				int t = sorter[i];
				sorter[i] = sorter[i+1];
				sorter[i+1] =t;
			}				
		}
	}	

	if( for_kc )
	{
		for( int i = 0; i < nmodes; i++ )
			printf("%le %le\n", qmode[sorter[i]], mode_sq[sorter[i]]/(1e-14+n_mode_sq[sorter[i]]) ); 
	}
}

int matches_atomselect( const char *select, const char *target )
{
	int s = 0;
	char segment[256];
	
	while( s < strlen(select) )
	{
		int len = 0;
	
		
		while( s+len < strlen(select) && select[s+len] != '+' )
		{	
			if( len < 255 )
			{
				segment[len] = select[s+len];
				segment[len+1] = '\0';
			}
			len++;
		}

		if( !strcasecmp( segment, target ) )
		{
//			printf("%s with length %d matches %s.\n", select+s, len, target );
			return 1;
		}
		s += len;

		while( select[s] == '+' ) s++;
	}
			
//	printf("%s didn't match %s.\n", select, target );

	return 0;
}


