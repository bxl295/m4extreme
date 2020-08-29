// NewtonRaphson.cpp: implementation of the NewtonRaphson class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include "./NewtonRaphson.h"

namespace Solver
{
NewtonRaphson::NewtonRaphson() {}

NewtonRaphson::~NewtonRaphson() {
            if(LSS!=0) delete LSS;
        }

NewtonRaphson::NewtonRaphson(
	Model::LocalState *LS_, 
	Model::Jet<1> *J_, 
	Solver::Linear::System<DOF> *S_, 
	set<Set::Manifold::Point *> *x_, UPDATE_TYPE type_) : 
	LS(LS_), J(J_), DE(0), DDE(0), S(S_), x(x_), 
	Print(false), Tol(0.0), NitMax(0), type(type_), LSS(0) {}

NewtonRaphson::NewtonRaphson(
	Model::LocalState *LS_, 
        Model::Energy<1> *rhs_DE, 
	Model::Jet<1> *J_, 
	Solver::Linear::System<DOF> *S_, 
	set<Set::Manifold::Point *> *x_, UPDATE_TYPE type_) : 
	LS(LS_), J(J_), DE(rhs_DE), DDE(0), S(S_), x(x_), 
	Print(false), Tol(0.0), NitMax(0), type(type_) 
        {
                LSS = new Solver::LineSearch::Secant(LS,DE,x);
                LSS->SetPrint() = false; 
		LSS->SetTol() = 1.0e-5; 
		LSS->SetNitMax() = 20;
		LSS->SetStep() = 1.0e-5; 
        }

NewtonRaphson::NewtonRaphson(
	Model::LocalState *LS_, 
	Model::Energy<1> *rhs_DE, 
	Model::Energy<2> *rhs_DDE, 
	Solver::Linear::System<DOF> *S_, 
	set<Set::Manifold::Point *> *x_, UPDATE_TYPE type_) : 
	LS(LS_), J(0), DE(rhs_DE), DDE(rhs_DDE), S(S_), x(x_), 
	Print(false), Tol(0.0), NitMax(0), type(type_) 
        {
                LSS = new Solver::LineSearch::Secant(LS,DE,x);
                LSS->SetPrint() = false; 
		LSS->SetTol() = 1.0e-5; 
		LSS->SetNitMax() = 20;
		LSS->SetStep() = 1.0e-5; 
        }

bool & 
NewtonRaphson::SetPrint()
{
	return Print;
}

unsigned int & 
NewtonRaphson::SetNitMax()
{
	return NitMax;
}

double & 
NewtonRaphson::SetTol()
{
	return Tol;
}

#if !defined(_M4EXTREME_NEWTONRAPHSON_FIXED_ITERATIONS_)
void 
NewtonRaphson::operator ++ ()
{
        typedef map<Set::Manifold::Point *, Set::Euclidean::Orthonormal::Point> map_type;
        typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;

	unsigned int k; double Error;
	if (Print) cout << "Newton-Raphson iteration:" << endl;
        if(J != 0){
            for (k=0; k<NitMax; k++)
            {
                    map_type xemb_old; ((Model::Static::LocalState *)LS)->Embed(*x,xemb_old); 
                    (*J)(*x,*S); 
		    try {
		      S->Solve(); 
		    }
		    catch (...) {
		      cerr << "Solver::NewtonRaphson Failed calling the SuperLU solver @ J!=0" << endl;
		      throw(0);
		    }
		    (*LS)(*S,*x); Error = S->Norm();
                    if (Error > Tol && LSS!=0 ){
                        map_type xemb; ((Model::Static::LocalState *)LS)->Embed(*x,xemb);
                        vector_type u;
                        map_type::iterator px1;
                        map_type::iterator px2;
                        for(px1 = xemb_old.begin(), px2 = xemb.begin(); px1 != xemb_old.end(); px1++, px2++){
                            Set::VectorSpace::Vector deltax = px2->second - px1->second;
                            u.insert(make_pair(px1->first,deltax));
                        }
                        vector_type uemb;
                        ((Model::Static::LocalState *)LS)->Submerge(u,uemb);
                        (*x) -= uemb;
                        (*LSS)(uemb);
                    }

                    if (Print) cout 
                            << "Iteration = " << k     << ", " 	
                            << "Error = "     << Error << endl;
                    if (Error < Tol) {if (type) ++(*LS); return;}
            }
        }
        else{
            for (k=0; k<NitMax; k++)
            {
                    map_type xemb_old; ((Model::Static::LocalState *)LS)->Embed(*x,xemb_old); 
                    (*DE)(*x,*S); (*DDE)(*x,*S);
		    try {
		      S->Solve();
		    }
		    catch (...) {
		      cerr << "Solver::NewtonRaphson Failed calling the SuperLU solver @ J==0" << endl;
		      throw(0);
		    }
		    (*LS)(*S,*x); Error = S->Norm();
                    if (Error > Tol && LSS!=0 ){
                        map_type xemb; ((Model::Static::LocalState *)LS)->Embed(*x,xemb);
                        vector_type u;
                        map_type::iterator px1;
                        map_type::iterator px2;
                        for(px1 = xemb_old.begin(), px2 = xemb.begin(); px1 != xemb_old.end(); px1++, px2++){
                            Set::VectorSpace::Vector deltax = px2->second - px1->second;
                            u.insert(make_pair(px1->first,deltax));
                        }
                        vector_type uemb;
                        ((Model::Static::LocalState *)LS)->Submerge(u,uemb);
                        (*x) -= uemb;
                        (*LSS)(uemb);
                    }
                    if (Print) cout 
                            << "Iteration = " << k     << ", " 	
                            << "Error = "     << Error << endl;
                    if (Error < Tol) {if (type) ++(*LS); return;}
            }
            
        }
        
        cerr << "Solver::NewtonRaphson Failed to converge @ Tol=" 
             << Tol << " and NitMax=" << NitMax << ": error=" << Error << endl;
	throw(0);
}

#else

void 
NewtonRaphson::operator ++ ()
{
        typedef map<Set::Manifold::Point *, Set::Euclidean::Orthonormal::Point> map_type;
        typedef map<Set::Manifold::Point *, Set::VectorSpace::Vector> vector_type;

        NitMax = 5;
        
	unsigned int k;
	if (Print) cout << "Newton-Raphson iteration:" << endl;
        if(J != 0){
            for (k=0; k<NitMax; k++)
            {
                    map_type xemb_old; ((Model::Static::LocalState *)LS)->Embed(*x,xemb_old); 
                    (*J)(*x,*S); 
		    try {
		      S->Solve(); 
		    }
		    catch (...) {
		      cerr << "Solver::NewtonRaphson Failed calling the SuperLU solver @ J!=0" << endl;
		      throw(0);
		    }
		    (*LS)(*S,*x); 
                    
		    if ( LSS != 0 ) {
		      map_type xemb; ((Model::Static::LocalState *)LS)->Embed(*x,xemb);
		      vector_type u;
		      map_type::iterator px1;
		      map_type::iterator px2;
		      for(px1 = xemb_old.begin(), px2 = xemb.begin(); px1 != xemb_old.end(); px1++, px2++){
			Set::VectorSpace::Vector deltax = px2->second - px1->second;
			u.insert(make_pair(px1->first,deltax));
		      }
		      vector_type uemb;
		      ((Model::Static::LocalState *)LS)->Submerge(u,uemb);
		      (*x) -= uemb;
		      (*LSS)(uemb);
                    }                 
            }	    
        }
        else{
            for (k=0; k<NitMax; k++)
            {
                    map_type xemb_old; ((Model::Static::LocalState *)LS)->Embed(*x,xemb_old); 
                    (*DE)(*x,*S); (*DDE)(*x,*S);
		    try {
		      S->Solve();
		    }
		    catch (...) {
		      cerr << "Solver::NewtonRaphson Failed calling the SuperLU solver @ J==0" << endl;
		      throw(0);
		    }
		    (*LS)(*S,*x);
                    if (LSS!=0 ){
                        map_type xemb; ((Model::Static::LocalState *)LS)->Embed(*x,xemb);
                        vector_type u;
                        map_type::iterator px1;
                        map_type::iterator px2;
                        for(px1 = xemb_old.begin(), px2 = xemb.begin(); px1 != xemb_old.end(); px1++, px2++){
                            Set::VectorSpace::Vector deltax = px2->second - px1->second;
                            u.insert(make_pair(px1->first,deltax));
                        }
                        vector_type uemb;
                        ((Model::Static::LocalState *)LS)->Submerge(u,uemb);
                        (*x) -= uemb;
                        (*LSS)(uemb);
                    }
             }
            
        }

	if (type) ++(*LS);
	return;
}
#endif

}
