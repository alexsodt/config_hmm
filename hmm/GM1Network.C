
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "HMM.h"

#define DISSOCIATION_STATE

extern int globalN_OBSERVABLES;

/*
void MakeGM1Network( HMMSystem &theSystem )
{
	char label[256];
	
	sprintf(label, "State1" );
	HMMState *state1 = new HMMState('A', label, 1, 0 ); 

	sprintf(label, "State2" );
	HMMState *state2 = new HMMState('B', label, 1, 1 ); 
	

	state1->addTransition( state1 );
	state1->addTransition( state2 );

	state2->addTransition( state1 );
	state2->addTransition( state2 );

#ifdef DISSOCIATION_STATE
	sprintf(label, "StateX" );
	HMMState *state3 = new HMMState('x', label, 0, 2 ); 
	state1->addTransition( state3 );
	state2->addTransition( state3 );
	
	state3->addTransition( state1 );
	state3->addTransition( state2 );
	state3->addTransition( state3 );
#endif
	theSystem.addState(state1);
	theSystem.addState(state2);
#ifdef DISSOCIATION_STATE
	theSystem.addState(state3);
#endif
	theSystem.setup();	
}
*/

void MakeGM1Network( HMMSystem &theSystem, int nstates)
{
	char label[256];

	HMMState *states[nstates+1];

	for( int i = 1; i <= nstates; i++ )
	{
		sprintf(label, "State%d", 1 + i );
		states[i] = new HMMState('A'+i, label, 1, i ); 
	}

	sprintf(label, "StateX" );
	HMMState *stateX = new HMMState('x', label, 1, 0 ); 
	
	for( int i = 1; i <= nstates; i++ )
	{
		states[i]->addTransition( states[i] );
		states[i]->addTransition( stateX );
		stateX->addTransition(states[i] );
	}
		
	stateX->addTransition(stateX );
	
	for( int i = 1; i <= nstates; i++ )
	for( int j = 1; j <= nstates; j++ )
	{
		if( i == j ) continue;

		states[i]->addTransition( states[j] );
	}
	
	theSystem.addState(stateX);
	for( int i = 1; i <= nstates; i++ )
		theSystem.addState(states[i]);

	theSystem.setup();	
}
