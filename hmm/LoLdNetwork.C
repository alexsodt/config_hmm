
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "HMM.h"

extern int globalN_OBSERVABLES;

//#define PROTEIN_STATE // may be defined in HMM.h
//#define THREE_STATE
//#define SINGLE_STATE
//#define FOUR_STATE

void MakeLoLdNetwork( HMMSystem &theSystem, int ignore )
{
	char label[256];
	
	sprintf(label, "State1" );
	HMMState *state1 = new HMMState('o', label, 1, 0, 0 ); 

#ifdef SINGLE_STATE
	state1->addTransition(state1);
	theSystem.addState(state1);

#else
	sprintf(label, "State2" );
	HMMState *state2 = new HMMState('d', label, 1, 1, 1 ); 


#if defined(THREE_STATE) || defined(PROTEIN_STATE)

#if defined(THREE_STATE )
	sprintf(label, "Intermediate" );
	HMMState *state3 = new HMMState('i', label, 1, 2, 2 ); 
#else
	sprintf(label, "ProteinBorder" );
	HMMState *state3 = new HMMState('b', label, 1, 2, 2 ); 
#endif

	state1->addTransition( state3 );
	state2->addTransition( state3 );
	state3->addTransition( state3 );
	
	state1->addTransition( state2 );
	state2->addTransition( state2 );
	state3->addTransition( state2 );
	
	state1->addTransition( state1 );
	state2->addTransition( state1 );
	state3->addTransition( state1 );
	
	theSystem.addState(state3);
#elif defined(FOUR_STATE)
	sprintf(label, "Intermediate" );
	HMMState *state3 = new HMMState('i', label, 1, 2, 2 ); 
	sprintf(label, "Intermediate2" );
	HMMState *state4 = new HMMState('j', label, 1, 3, 3 );
 
	
	state1->addTransition( state1 );
	state2->addTransition( state1 );
	state3->addTransition( state1 );
	state4->addTransition( state1 );
	
	state1->addTransition( state2 );
	state2->addTransition( state2 );
	state3->addTransition( state2 );
	state4->addTransition( state2 );

	state1->addTransition( state3 );
	state2->addTransition( state3 );
	state3->addTransition( state3 );
	state4->addTransition( state3 );
	
	state1->addTransition( state4 );
	state2->addTransition( state4 );
	state3->addTransition( state4 );
	state4->addTransition( state4 );
	
	theSystem.addState(state3);
	theSystem.addState(state4);
#else
	state1->addTransition( state1 );
	state2->addTransition( state2 );

	state1->addTransition( state2 );
	state2->addTransition( state1 );
#endif
	theSystem.addState(state1);
	theSystem.addState(state2);
#endif

	theSystem.setup();	
}

