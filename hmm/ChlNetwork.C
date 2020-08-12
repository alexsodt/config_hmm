
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "HMM.h"

#define DISSOCIATION_STATE

extern int globalN_OBSERVABLES;


void MakeChlNetwork( HMMSystem &theSystem )
{
	char label[256];
	
	sprintf(label, "Lower" );
	HMMState *state1 = new HMMState('L', label, 1, 0 ); 
	
	sprintf(label, "MidLower" );
	HMMState *state2 = new HMMState('l', label, 1, 1 ); 

	sprintf(label, "MidUpper" );
	HMMState *state3 = new HMMState('u', label, 1, 2 ); 
	
	sprintf(label, "Upper" );
	HMMState *state4 = new HMMState('U', label, 1, 3 ); 
	

	state1->addTransition( state1 );
	state1->addTransition( state2 );

	state2->addTransition( state1 );
	state2->addTransition( state2 );
	state2->addTransition( state3 );
	
	state3->addTransition( state2 );
	state3->addTransition( state3 );
	state3->addTransition( state4 );
	
	state4->addTransition( state3 );
	state4->addTransition( state4 );

	theSystem.addState(state1);
	theSystem.addState(state2);
	theSystem.addState(state3);
	theSystem.addState(state4);
	
	theSystem.setup();	
}


