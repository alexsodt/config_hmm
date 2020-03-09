#ifndef __NR_LINMIN_H__
#define __NR_LINMIN_H__

void check_surface_consistency( HMMSystem &theSystem, HMMDirective &theDirective );
double surface_grad_opt( HMMSystem &theSystem, HMMDirective &theDirective );
double grad_opt( HMMSystem &theSystem, HMMDirective &theDirective );
double bwgrad_opt( HMMSystem &theSystem, HMMDirective &theDirective );
#endif
