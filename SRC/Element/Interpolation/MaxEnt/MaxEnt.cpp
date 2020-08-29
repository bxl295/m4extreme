// MaxEnt.cpp: implementation of the MaxEnt class.
// Copyright (c) 2006 by Bo Li - All rights Set::Manifold::Point *
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include <vector>
#include <algorithm>
#include <limits>
#include "./MaxEnt.h"
#include "floatpoint.h"

extern map<Set::Manifold::Point *, unsigned int> _gIndexSet;
//extern bool _gPrint;

namespace Element {
    namespace Interpolation {
        namespace MaxEnt {
            //////////////////////////////////////////////////////////////////////
            // Class Base
            //////////////////////////////////////////////////////////////////////

	    Base::Base() : _Beta(0.), _JInv(NULL), _Lambda(NULL), _x(NULL) {
            }

            Base::Base(const double &Beta, const nodeset_type & x)
            : _Beta(Beta), _x(&x) {

                const size_t n = _x->begin()->second.size();

                const size_t npp = n + 1;

                if (_x->size() == npp) {

                    Set::VectorSpace::Vector qp(n);

                    nodeset_type::const_iterator px;
                    for (px = _x->begin(); px != _x->end(); px++) {
                        qp += px->second / (double) npp;
                    }

                    vector<Set::VectorSpace::Vector> ya(npp, qp);
                    double ya2[npp];

                    px = _x->begin();
                    ya[0] -= px->second;
                    ya2[0] = ya[0](ya[0]);
                    px++;

                    Set::VectorSpace::Hom A(n);
                    Set::VectorSpace::Vector u(n);
                    size_t k = 1;

                    double minYa2 = ya2[0];
                    for (; px != _x->end(); px++, k++) {
                        ya[k] -= px->second;
                        ya2[k] = ya[k](ya[k]);
                        u[k - 1] = _Beta * (ya2[k] - ya2[0]);
                        A[k - 1] = ya[k] - ya[0];

                        if (minYa2 > ya2[k]) {
                            minYa2 = ya2[k];
                        }
                    }

                    _Lambda = new Set::VectorSpace::Vector(Inverse(Adjoint(A)) * u);

                    // check if beta has the right value                   
                    for (size_t i = 0; i < npp; i++) {
                        if (-_Beta * ya2[i] + ya[i](*_Lambda) < -60) {
                            (*_Lambda) /= _Beta;
                            _Beta = 1.0 / sqrt(minYa2);
                            (*_Lambda) *= _Beta;

                            break;
                        }
                    }

                } else {
                    _Lambda = new Set::VectorSpace::Vector(n);
                }

                _JInv = new Set::VectorSpace::Hom(n);

            }

            Base::~Base() {
                if (_JInv) {
                    delete _JInv;
                }

                if (_Lambda) {
                    delete _Lambda;
                }
            }

            const double & Base::GetBeta() const {
                return Base::_Beta;
            }

            void Base::SetBeta(const double & beta) {
                Base::_Beta = beta;
                return;
            }

            void Base::SetNodes(const nodeset_type * x) {
                _x = x;
            }

            const map<Set::Manifold::Point *, Set::VectorSpace::Vector> *
            Base::GetNodes() const {
                return _x;
            }

            const Set::VectorSpace::Vector & Base::GetLambda() const {
                return *_Lambda;
            }

            const vector<double> &
            Base::GetN() const {
                return _NNodes;
            }

            const Set::VectorSpace::Hom *
            Base::GetJInverse() const {
                return _JInv;
            }

            /* // original implementation in m4extreme
              void Base::operator() (const Set::VectorSpace::Vector & y) {

              assert(_JInv && _Lambda);
              _NNodes.clear();

              unsigned int n = y.size();	
              double Dx2 = 0.0;
              map<dof_type *, double> E; 
              map<dof_type *, Set::VectorSpace::Vector> Dx;
              map<dof_type *, Set::VectorSpace::Vector>::const_iterator px;
              for (px=_x->begin(); px!=_x->end(); px++)
              {
              Set::VectorSpace::Vector Dxa = y - px->second;
              double DxaDxa = Dxa(Dxa); 
              if (DxaDxa > Dx2) Dx2 = DxaDxa;
              E.insert( make_pair(px->first, _Beta*DxaDxa) );
              Dx.insert(make_pair(px->first,Dxa));
              }

              unsigned int k, Max=20;
              //double Tol=1.0e-12*sqrt(Dx2);
              double Tol = 1.0e-8;

              Set::VectorSpace::Vector r(n);	
              Set::VectorSpace::Vector DLambda(n);
              Set::VectorSpace::Hom J(n);
              for (k=0; k<Max; k++)
              {
              double Z=0.0; Null(r); Null(J);
              map<dof_type *, double>::iterator pE; 
              map<dof_type *, Set::VectorSpace::Vector>::iterator pDx;
              map<dof_type *, double> pa;
              for (pE=E.begin(), pDx=Dx.begin(); pE!=E.end(); pE++, pDx++)
              {
              double &Ea = pE->second;
              Set::VectorSpace::Vector &Dxa = pDx->second;
              double Za = exp( -Ea + (*_Lambda)(Dxa) );
              Set::VectorSpace::Vector Dy = Za*Dxa; 
              r += Dy; Z += Za; J += Dyadic(Dy,Dxa);

              pa.insert( make_pair(pE->first, Za) );
              }
              J /= Z; r /= Z; J -= Dyadic(r,r); 
              (*_JInv) = Inverse(J);
              DLambda = -1.0 * (*_JInv)(r); 
              (*_Lambda) += DLambda;

              if (Norm(DLambda) < Tol) {

              map<dof_type*, double>::iterator pPa;
              for ( pPa = pa.begin(); pPa != pa.end(); pPa++ ) {
              _NNodes.insert( make_pair(pPa->first, pPa->second / Z) );
              }

              return;
              }

              }
              throw(0);        
              }
             */

            void Base::operator() (const Set::VectorSpace::Vector &qp) {

                assert(_JInv && _Lambda);	  
		if ( !_NNodes.empty() ) _NNodes.clear();

#if !defined(_M4EXTREME_NEIGHBOR_FIXED_)
                int dim = _x->begin()->second.size();

#ifdef _M4EXTREME_DEBUG_III        
                std::cout << "dim = " << dim << "\t n = " << _x->begin()->second.size() << std::endl;
#endif

                switch (dim) {
                    case 1:
                        return _computeLambda1D(qp);
                    case 2:

#if defined(_M4EXTREME_VISCOUS_REGULARIZATION_)
                        return _computeLambda2D_Regularized(qp);
#else
			return _computeLambda2D(qp);
#endif

                    case 3:

#if defined(_M4EXTREME_VISCOUS_REGULARIZATION_)
                        return _computeLambda3D_Regularized(qp);
#else
			return _computeLambda3D(qp);
#endif
                    default:
                        assert(false);
                }
#else
                _computeExactLambda(qp);
#endif

            }

            void Base::_computeExactLambda(const Set::VectorSpace::Vector & qp) {

                if (_Lambda != NULL) {
                    delete _Lambda;
                }

                size_t n = _x->begin()->second.size();
                double pa = 1.0 / (double) (n + 1);

                nodeset_type::const_iterator px;
                px = _x->begin();
                Set::VectorSpace::Vector y0 = qp - px->second;
                px++;

                double norm_y0 = y0(y0);
                Set::VectorSpace::Hom A(n);
                Set::VectorSpace::Vector u(n);
                size_t k = 0;

                for (; px != _x->end(); px++, k++) {
                    Set::VectorSpace::Vector yloc = qp - px->second;
                    u[k] = _Beta * (yloc(yloc) - norm_y0);
                    A[k] = yloc - y0;
                }

                _Lambda = new Set::VectorSpace::Vector(Inverse(Adjoint(A)) * u);

                Set::VectorSpace::Hom J(n);
                for (px = _x->begin(); px != _x->end(); px++) {
		    _NNodes.push_back(pa);
                    Set::VectorSpace::Vector yloc = qp - px->second;
                    J += Dyadic(yloc, yloc) * pa;
                }

                (*_JInv) = Inverse(J);

                return;
            }

            void Base::_computeLambda1D(const Set::VectorSpace::Vector &qp) {

                vector<double> dx, dxdx;
                double qpx = qp[0];

                int numofNodes = _x->size();

                nodeset_type::const_iterator px;
                double dxa2 = 0.0, norm_dx = 0.0;
                for (px = _x->begin(); px != _x->end(); px++) {
                    double dxa = qpx - (px->second)[0];
                    dx.push_back(dxa);

                    dxa2 = dxa * dxa;
                    dxdx.push_back(dxa2);

                    if (norm_dx < dxa2) {
                        norm_dx = dxa2;
                    }
                }

                norm_dx = sqrt(norm_dx);
                //double detTol = 1.0e-8 * norm_dx;
                //double rTol   = 1.0e-4 * norm_dx;

                double detTol = 1.0e-8;
                double rTol = 1.0e-6;

                double & lambda = (*_Lambda)[0];

                // Iteration
                double err = 1.0;
                double r = 0.0;
                double tol_err = 1.0e-24;
                double tol_r = 1.0e-10;
                double dt = 0.1;
                double Z = 0.0;
                int k = 0;
                size_t gfcount = 0;
                size_t ResetCount = 100;
                double j = 0.0;

                vector<double> N(numofNodes, 0.0);

                do {
                    Z = 0.0;

                    for (unsigned int i = 0; i < numofNodes; i++) {

                        double pa = exp(-_Beta * dxdx[i] + lambda * dx[i]);

                        Z += pa;
                        N[i] = pa;
                    }

                    double shape_a;
                    j = 0.0;
                    r = 0.0;
                    for (unsigned int i = 0; i < numofNodes; i++) {
                        shape_a = N[i] / Z;
                        r += dx[i] * shape_a;
                        j += dxdx[i] * shape_a;
                    }

                    j -= r * r;

                    // Newton Raphson method fails here
                    // Using steepest gradient to find a close point near 
                    // the stationary point of the objective function log(Z);
                    if ((j > 0 ? j : -j) < detTol && fabs(r) > tol_r) {

#if !defined(_M4EXTREME_EXACT_LINE_SEARCH_)
                        double fold = log(Z);
#endif            

                        do {

#if defined(_M4EXTREME_EXACT_LINE_SEARCH_)
                            // compute the exact line search step length              
                            vector<double> cs;
                            for (int i = 0; i < numofNodes; i++) {
                                cs.push_back(-r * dx[i]);
                            }

                            double fl = 0.0;
                            for (int i = 0; i < numofNodes; i++) {
                                fl += N[i] * cs[i] / Z;
                            }

                            dt = _bisection(fl, cs, N, dxdx);
#endif	               

                            lambda -= r * dt;

                            Z = 0.0;
                            r = 0.0;

                            for (unsigned int i = 0; i < numofNodes; i++) {

                                const double & dxa = dx[i];

                                double pa = exp(-_Beta * dxdx[i] + lambda * dxa);

                                Z += pa;

                                r += dxa * pa;
                            }

                            r /= Z;

#if !defined(_M4EXTREME_EXACT_LINE_SEARCH_)              
                            double fnew = log(Z);
                            if (fnew > fold) {
                                dt *= 0.5;
                            } else {
                                dt *= 1.2;
                            }

                            fold = fnew;
#endif
                        } while (gfcount++ < ResetCount && fabs(r) > rTol);

                        if (gfcount > ResetCount) {
                            cout << "theta:=<" << qpx << ">;" << endl;
                            cout << "dx[1]:=<";
                            int i = 0;
                            for (; i < dx.size() - 1; i++) {
                                cout << dx[i] << ",";
                            }
                            cout << dx[i] << ">;" << endl;

                            //                return;
                            throw ("steepest descent failed!");
                        }

                        --k;
                    } else {

                        err = r / j; // This is the error!

#if defined(_M4EXTREME_MODIFIED_NEWTON_)
                        // compute the exact line search step length    

                        vector<double> cs;
                        for (int i = 0; i < numofNodes; i++) {
                            cs.push_back(-err * dx[i]);
                        }

                        double fl = 0.0;
                        for (int i = 0; i < numofNodes; i++) {
                            fl += N[i] * cs[i] / Z;
                        }

                        double sm = _bisection(fl, cs, N, dxdx);

                        // update lambda -- Modified Newton-Raphson method
                        lambda -= sm * err;
#else
                        // update lambda -- Newton-Raphson method
                        lambda -= err;
#endif	   

                    }

                } while (fabs(r) > tol_r && k++ < 20);

                if (k < 21) {
                    px = _x->begin();
                    for (unsigned int i = 0; i < numofNodes; i++, px++) {
		      _NNodes.push_back( N[i] / Z );
                    }

                    (*_JInv)[0][0] = 1.0 / j;

                    return;
                } else {
                    throw ("Newton Raphson failed!");
                }

            }

	  void Base::_computeLambda2D(const Set::VectorSpace::Vector &qp) {

	    vector<double> dx, dy, dxdx;
	    double qpx = qp[0];
	    double qpy = qp[1];

	    int numofNodes = _x->size();

	    double dxa2 = 0.0, norm_dx = 0.0;
	    double dxa, dya;
	    nodeset_type::const_iterator px;
	    for (px = _x->begin(); px != _x->end(); px++) {
	      dxa = qpx - (px->second)[0];
	      dya = qpy - (px->second)[1];
	      dx.push_back(dxa);
	      dy.push_back(dya);

	      dxa2 = dxa * dxa + dya * dya;
	      dxdx.push_back(dxa2);

	      if (norm_dx < dxa2) {
		norm_dx = dxa2;
	      }
	    }

	    norm_dx = sqrt(norm_dx);
	    double detTol = 1.0e-8;
	    double rTol = 1.0e-8;

	    double & lambda_1 = (*_Lambda)[0];
	    double & lambda_2 = (*_Lambda)[1];

	    // Iteration
	    m4extreme::floatpoint err(1.0);

	    double tol = 1.0e-24;
	    double r1 = 0.0, r2 = 0.0;
	    double tol_r = 1.0e-20;
	    double dt = 1.0;
	    double Z = 0.0, det;
	    int k = 0;
	    size_t gfcount = 0;
	    size_t ResetCount = 10000;
	    double i_j11, i_j12, i_j22;

	    vector<double> N(numofNodes, 0.0);

	    do {
	      Z = 0.0;

	      for (unsigned int i = 0; i < numofNodes; i++) {
		double pa = exp(-_Beta * dxdx[i] + lambda_1 * dx[i] + lambda_2 * dy[i]);

		Z += pa;
		N[i] = pa;
	      }

	      double j11 = 0.0, j12 = 0.0, j22 = 0.0;
	      double shape_a;

	      r1 = 0.0;
	      r2 = 0.0;
	      for (unsigned int i = 0; i < numofNodes; i++) {

		const double & dxa = dx[i];
		const double & dya = dy[i];

		shape_a = N[i] / Z;

		r1 += dxa * shape_a;
		r2 += dya * shape_a;

		j11 += dxa * dxa * shape_a;
		j12 += dxa * dya * shape_a;
		j22 += dya * dya * shape_a;
	      }

              err = r1 * r1 + r2 * r2;
              
	      j11 -= r1 * r1;
	      j12 -= r1 * r2;
	      j22 -= r2 * r2;

	      det = j11 * j22 - j12 * j12;

	      i_j11 = j22 / det;
	      i_j12 = -j12 / det;
	      i_j22 = j11 / det;

	      if (__builtin_expect(err < tol_r, 0) ) {
		    px = _x->begin();
		    for (unsigned int i = 0; i < numofNodes; i++, px++) {
		      double Na = N[i] / Z;
		      if (Na > 1.0) {
			cout << "Bug: Na = " << Na << endl;
			assert(false);
		      }
		      _NNodes.push_back(Na);
		    }

		    (*_JInv)[0][0] = i_j11;
		    (*_JInv)[1][1] = i_j22;
		    (*_JInv)[0][1] = (*_JInv)[1][0] = i_j12;

		    return;
		  } else {
		    // Newton Raphson method fails here, using
		    // steepest descent method to find a close point near
		    // the stationary point of the objective function log(Z);
		    if ((det > 0 ? det : -det) < detTol) {

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)
		      _loop_list L1, L2;
#endif

#if !defined(_M4EXTREME_EXACT_LINE_SEARCH_)
		      double fold = log(Z);
#endif            

		      do {

#if defined(_M4EXTREME_EXACT_LINE_SEARCH_)
			// compute the exact line search step length

			vector<double> cs;
			for (int i = 0; i < numofNodes; i++) {
			  cs.push_back(-r1 * dx[i] - r2 * dy[i]);
			}

			double fl = 0.0;
			for (int i = 0; i < numofNodes; i++) {
			  fl += N[i] * cs[i] / Z;
			}

			dt = _bisection(fl, cs, N, dxdx);
#endif	   

			lambda_1 -= r1 * dt;
			lambda_2 -= r2 * dt;

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)
			if (fabs(L1.head() - lambda_1) < tol &&
			    fabs(L2.head() - lambda_2) < tol) {
			  dt *= 1.2;
			  lambda_1 -= r1 * dt;
			  lambda_2 -= r2 * dt;
			}

			L1.head() = lambda_1;
			L2.head() = lambda_2;
                        ++L1; ++L2;
#endif

			Z = 0.0;
			r1 = 0.0;
			r2 = 0.0;

			for (unsigned int i = 0; i < numofNodes; i++) {

			  const double & dxa = dx[i];
			  const double & dya = dy[i];

			  double pa = exp(-_Beta * dxdx[i] + lambda_1 * dxa + lambda_2 * dya);

			  Z += pa;

			  r1 += dxa * pa;
			  r2 += dya * pa;
			}

			r1 /= Z;
			r2 /= Z;

			err = r1 * r1 + r2 * r2;

#if !defined(_M4EXTREME_EXACT_LINE_SEARCH_)              
			double fnew = log(Z);

			if (fnew > fold) {
			  dt *= 0.5;
			} else {
			  dt *= 1.2;
			}

			fold = fnew;
#endif
		      } while (gfcount++ < ResetCount && err > rTol);

		      if (gfcount > ResetCount || Z < 1.0e-24 || Z > 1.0e24) {

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)			        
			L1.reset(); L2.reset();
#endif

			tol_r = 1.0e-12;

			if (norm_dx < 1.0) {
			  _Beta *= norm_dx;
			} else {
			  _Beta *= 0.1;
			}

			gfcount = 0;
			ResetCount *= 1000;

			if (Z < 1.0e-24 || Z > 1.0e24) {
			  lambda_1 = 0.0;
			  lambda_2 = 0.0;
			}

			dt = 0.0;

			do {

			  Z = 0.0;
			  r1 = 0.0;
			  r2 = 0.0;

			  for (unsigned int i = 0; i < numofNodes; i++) {

			    const double & dxa = dx[i];
			    const double & dya = dy[i];

			    double pa = exp(-_Beta * dxdx[i] + lambda_1 * dxa + lambda_2 * dya);

			    Z += pa;
			    N[i] = pa;

			    r1 += dxa * pa;
			    r2 += dya * pa;
			  }

			  r1 /= Z;
			  r2 /= Z;

			  err = r1 * r1 + r2 * r2;

#ifdef _M4EXTREME_DEBUG_III                                    
			  cout << rTol << ":lambda_1 = " << lambda_1 << "\tlambda_2 = " << lambda_2 <<
			       << "\tr1 = " << r1 << "\tr2 = " << r2 << "\terr = " << err << endl;
#endif

			  if (err > tol_r) {
			    // compute the exact line search step length for next iteration
			    vector<double> cs;
			    for (int i = 0; i < numofNodes; i++) {
			      cs.push_back(-r1 * dx[i] - r2 * dy[i]);
			    }

			    double fl = 0.0;
			    for (int i = 0; i < numofNodes; i++) {
			      fl += N[i] * cs[i] / Z;
			    }

			    dt = _bisection(fl, cs, N, dxdx);

			    lambda_1 -= r1 * dt;
			    lambda_2 -= r2 * dt;

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)
			    if (fabs(L1.head() - lambda_1) < tol &&
				fabs(L2.head() - lambda_2) < tol) {
			      dt *= 1.2;
			      lambda_1 -= r1 * dt;
			      lambda_2 -= r2 * dt;
			    }

			    L1.head() = lambda_1;
			    L2.head() = lambda_2;
                            ++L1; ++L2;
#endif

			  } else {
#ifdef _M4EXTREME_DEBUG_II
			    cout << "Reset tolerance of the steepest descent solver" << endl;
#endif
			    break;
			  }

#ifdef _M4EXTREME_DEUBG_III
			  cout << rTol << ": dt = " << dt << endl;
#endif

			} while (gfcount++ < ResetCount);

			if (gfcount > ResetCount) {
			  cout << "theta:=<" << qpx << "," << qpy << ">;" << endl;
			  cout << "dx[1]:=<";
			  int i = 0;
			  for (; i < dx.size() - 1; i++) {
			    cout << dx[i] << ",";
			  }
			  cout << dx[i] << ">;" << endl;

			  cout << "dx[2]:=<";
			  i = 0;
			  for (; i < dy.size() - 1; i++) {
			    cout << dy[i] << ",";
			  }
			  cout << dy[i] << ">;" << endl;

			  cout << "lambda:=<" << lambda_1 << "," << lambda_2 << ">" << endl;

			  throw ("steepest descent failed!");
			}
		      }

		      --k;
		    } else {

		      double dl1 = i_j11 * r1 + i_j12 * r2;
		      double dl2 = i_j12 * r1 + i_j22 * r2;

#if !defined(_M4EXTREME_NEWTON_BACKTRACKING_LINE_SEARCH)
		      // compute the exact line search step length
		      vector<double> cs;
		      for (int i = 0; i < numofNodes; i++) {
			cs.push_back(-dl1 * dx[i] - dl2 * dy[i]);
		      }

		      double fl = 0.0;
		      for (int i = 0; i < numofNodes; i++) {
			fl += N[i] * cs[i] / Z;
		      }

		      double sm = _bisection(fl, cs, N, dxdx);

		      // update lambda -- ewton-Raphson method with exact line search
		      lambda_1 -= sm * dl1;
		      lambda_2 -= sm * dl2;
#else
		      // update lambda -- Newton-Raphson method with backtracking line search
		      double sm = 1.0, gamma = 0.5, alpha = 0.01;
		      double f0 = log(Z);
		      double df = alpha * (r1 * dl1 + r2 * dl2);
		      double fold, fnew;
		      int kb = 0;

		      do {
			fold = f0 - df * sm;

			lambda_1 -= sm * dl1;
			lambda_2 -= sm * dl2;

			fnew = 0.0;
			for (int i = 0; i < numofNodes; i++) {
			  fnew += exp(-_Beta * dxdx[i] + lambda_1 * dx[i] + lambda_2 * dy[i]);
			}

			fnew = log(fnew);
			sm *= gamma;

		      } while (fnew < fold && kb++ < 100);
#endif
		    }

		  }

	      } while (k++ < 20);

	      cout << "Newton-Raphson solver failed >>>>>>>>>>>>>>> " << endl
		   << "beta := " << _Beta << ";" << endl
		   << "theta:=<" << qpx << "," << qpy << ">;" << endl;

	      cout << "dx[1]:=<";
	      int i = 0;
	      for (; i < dx.size() - 1; i++) {
		cout << dx[i] << ",";
	      }
	      cout << dx[i] << ">;" << endl;

	      cout << "dx[2]:=<";
	      i = 0;
	      for (; i < dy.size() - 1; i++) {
		cout << dy[i] << ",";
	      }
	      cout << dy[i] << ">;" << endl;

	      cout << "lambda_1 = " << lambda_1 << "\tlambda_2 = " << lambda_2 << endl;
	      cout << "r1 = " << r1 << "\tr2 = " << r2
		   << "\terr=" << err << "\tdet = " << det << "\tZ = " << Z << endl;


	      cout << "i_j11 = " << i_j11 << "\ti_j12 = " << i_j12 << "i_j22 = " << i_j22 << endl;

	      cout << "Za : [\t";
	      for (unsigned int i = 0; i < numofNodes; i++) {
		cout << N[i] << "\t";
	      }
	      cout << "]" << endl;

	      throw ("Newton Raphson failed!");

            }
            
            void Base::_computeLambda3D(const Set::VectorSpace::Vector &qp) {

	        double beta = _Beta;

                vector<double> dx, dy, dz, dxdx;
                double qpx = qp[0];
                double qpy = qp[1];
                double qpz = qp[2];

                int numofNodes = _x->size();

                double dxa, dya, dza;
                double dxa2 = 0.0, norm_dx = 0.0;
                nodeset_type::const_iterator px;
                for (px = _x->begin(); px != _x->end(); px++) {
                    dxa = qpx - (px->second)[0];
                    dya = qpy - (px->second)[1];
                    dza = qpz - (px->second)[2];
                    dx.push_back(dxa);
                    dy.push_back(dya);
                    dz.push_back(dza);

                    dxa2 = dxa * dxa + dya * dya + dza * dza;
                    dxdx.push_back(dxa2);

                    if (norm_dx < dxa2) {
                        norm_dx = dxa2;
                    }
                }

                norm_dx = sqrt(norm_dx);
                //double detTol = 1.0e-8 * norm_dx;
		//double rTol   = 1.0e-8 * norm_dx;
                double detTol = 1.0e-8;		
                double rTol = 1.0e-8;

                double & lambda_1 = (*_Lambda)[0];
                double & lambda_2 = (*_Lambda)[1];
                double & lambda_3 = (*_Lambda)[2];

                // Iteration
		m4extreme::floatpoint err(1.0);

                double tol = 1.0e-24;
                double r1 = 0.0, r2 = 0.0, r3 = 0.0;
                double tol_r = 1.0e-20;
                double dt = 1.0;
                double Z = 0.0, det;
                int k = 0;
                size_t gfcount = 0;
                size_t ResetCount = 10000;
                double i_j11, i_j12, i_j13, i_j22, i_j23, i_j33;

                vector<double> N(numofNodes, 0.0);

                do {

                    Z = 0.0;
                    for (unsigned int i = 0; i < numofNodes; i++) {
                        double pa = exp(-beta * dxdx[i] + lambda_1 * dx[i] + lambda_2 * dy[i] + lambda_3 * dz[i]);

                        Z += pa;
                        N[i] = pa;
                    }

                    while (__builtin_expect(Z < 1.0e-32 || Z > 1.0e32, 0)) {
                        Z = 0.0;

                        if (norm_dx > 1.0) {
                            beta /= norm_dx;
                            lambda_1 /= norm_dx;
                            lambda_2 /= norm_dx;
                            lambda_3 /= norm_dx;
                        } else {
                            beta /= 2.0;
                            lambda_1 /= 2.0;
                            lambda_2 /= 2.0;
                            lambda_3 /= 2.0;
                        }

                        for (unsigned int i = 0; i < numofNodes; i++) {
                            double pa = exp(-beta * dxdx[i] + lambda_1 * dx[i] + lambda_2 * dy[i] + lambda_3 * dz[i]);

                            Z += pa;
                            N[i] = pa;
                        }
                    }

                    double j11 = 0.0, j12 = 0.0, j22 = 0.0;
                    double j13 = 0.0, j23 = 0.0, j33 = 0.0;
                    double shape_a;

                    r1 = 0.0;
                    r2 = 0.0;
                    r3 = 0.0;
                    for (unsigned int i = 0; i < numofNodes; i++) {

                        const double & dxa = dx[i];
                        const double & dya = dy[i];
                        const double & dza = dz[i];

                        shape_a = N[i] / Z;

                        r1 += dxa * shape_a;
                        r2 += dya * shape_a;
                        r3 += dza * shape_a;

                        j11 += dxa * dxa * shape_a;
                        j12 += dxa * dya * shape_a;
                        j22 += dya * dya * shape_a;
                        j13 += dxa * dza * shape_a;
                        j23 += dya * dza * shape_a;
                        j33 += dza * dza * shape_a;
                    }

                    err = r1 * r1 + r2 * r2 + r3 * r3; // this is the criterion

                    j11 -= r1 * r1;
                    j12 -= r1 * r2;
                    j22 -= r2 * r2;
                    j13 -= r1 * r3;
                    j23 -= r2 * r3;
                    j33 -= r3 * r3;

                    det = j11 * j22 * j33 - j11 * j23 * j23 + 2.0 * j12 * j23 * j13
                            - j12 * j12 * j33 - j22 * j13 * j13;

                    i_j11 = (j22 * j33 - j23 * j23) / det;
                    i_j12 = (j23 * j13 - j12 * j33) / det;
                    i_j13 = (j12 * j23 - j22 * j13) / det;
                    i_j22 = (-j13 * j13 + j11 * j33) / det;
                    i_j23 = (-j11 * j23 + j13 * j12) / det;
                    i_j33 = (-j12 * j12 + j11 * j22) / det;

                    if (__builtin_expect(err < tol_r, 0)) {
                        px = _x->begin();
                        for (unsigned int i = 0; i < numofNodes; i++, px++) {
			  double Na = N[i] / Z;
			  if ( Na > 1.0 ) {
			    cout << "Bug: Na = " << Na << endl;
			    assert(false);
			  }
			  _NNodes.push_back(Na);
                        }

                        (*_JInv)[0][0] = i_j11;
                        (*_JInv)[1][1] = i_j22;
                        (*_JInv)[0][1] = (*_JInv)[1][0] = i_j12;
                        (*_JInv)[2][2] = i_j33;
                        (*_JInv)[0][2] = (*_JInv)[2][0] = i_j13;
                        (*_JInv)[2][1] = (*_JInv)[1][2] = i_j23;			

                        return;
                    } else {

                        // Newton Raphson method fails here, using
                        // steepest descent method to find a closer point
                        // near the stationary point of the objective function log(Z)         
                        if ((det > 0 ? det : -det) < detTol) {

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)
			  _loop_list L1, L2, L3;
#endif

#if !defined(_M4EXTREME_EXACT_LINE_SEARCH_)
                            double fold = log(Z);
#endif

                            do {

#if defined(_M4EXTREME_EXACT_LINE_SEARCH_)
                                // compute the exact line search step length
                                vector<double> cs;
                                for (int i = 0; i < numofNodes; i++) {
                                    cs.push_back(-r1 * dx[i] - r2 * dy[i] - r3 * dz[i]);
                                }

                                double fl = 0.0;
                                for (int i = 0; i < numofNodes; i++) {
                                    fl += N[i] * cs[i] / Z;
                                }

                                dt = _bisection(fl, cs, N, dxdx);
#endif

                                lambda_1 -= r1 * dt;
                                lambda_2 -= r2 * dt;
                                lambda_3 -= r3 * dt;

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)
				if ( fabs(L1.head() - lambda_1) < tol &&
				     fabs(L2.head() - lambda_2) < tol &&
				     fabs(L3.head() - lambda_3) < tol ) {
				  dt *= 1.2;
				  lambda_1 -= r1 * dt;
				  lambda_2 -= r2 * dt;
				  lambda_3 -= r3 * dt;
				}

				L1.head() = lambda_1;
				L2.head() = lambda_2;
				L3.head() = lambda_3;
                                ++L1; ++L2; ++L3;
#endif

                                Z = 0.0;
                                r1 = 0.0;
                                r2 = 0.0;
                                r3 = 0.0;
                                for (unsigned int i = 0; i < numofNodes; i++) {

                                    const double & dxa = dx[i];
                                    const double & dya = dy[i];
                                    const double & dza = dz[i];

                                    double pa = exp(-beta * dxdx[i] + lambda_1 * dxa + lambda_2 * dya + lambda_3 * dza);

                                    Z += pa;

                                    r1 += dxa * pa;
                                    r2 += dya * pa;
                                    r3 += dza * pa;
                                }

                                r1 /= Z;
                                r2 /= Z;
                                r3 /= Z;

                                err = r1 * r1 + r2 * r2 + r3 * r3;

#if !defined(_M4EXTREME_EXACT_LINE_SEARCH_)
                                double fnew = log(Z);
                                if (fnew > fold) {
                                    dt *= 0.5;
                                } else {
                                    dt *= 1.2;
                                }

                                fold = fnew;
#endif              
                            } while (gfcount++ < ResetCount && err > rTol);

                            if (gfcount > ResetCount || Z < 1.0e-32 || Z > 1.0e32) {

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)			        
			        L1.reset(); L2.reset(); L3.reset();
#endif

                                tol_r = 1.0e-10;
                                
                                if (norm_dx < 1.0) {
                                    beta *= norm_dx;
                                } else {
                                    beta *= 0.1;
                                }                                                                

                                gfcount = 0;
                                ResetCount *= 1000;

                                dt = 0.0;

                                do {

                                    Z = 0.0;
                                    r1 = 0.0;
                                    r2 = 0.0;
                                    r3 = 0.0;
                                    for (unsigned int i = 0; i < numofNodes; i++) {

                                        const double & dxa = dx[i];
                                        const double & dya = dy[i];
                                        const double & dza = dz[i];

                                        double pa = exp(-beta * dxdx[i] + lambda_1 * dxa + lambda_2 * dya + lambda_3 * dza);

                                        Z += pa;
                                        N[i] = pa;

                                        r1 += dxa * pa;
                                        r2 += dya * pa;
                                        r3 += dza * pa;
                                    }

                                    r1 /= Z;
                                    r2 /= Z;
                                    r3 /= Z;

                                    err = r1 * r1 + r2 * r2 + r3 * r3;

#ifdef _M4EXTREME_DEBUG_III                                    
                                    cout << rTol << ":lambda_1 = " << lambda_1 << "\tlambda_2 = " << lambda_2 << "\tlambda_3 = " << lambda_3
                                            << "\tr1 = " << r1 << "\tr2 = " << r2 << "\tr3 = " << r3 << "\terr = " << err << endl;
#endif

                                    if (err > tol_r) {
                                        // compute the exact line search step length for next iteration
                                        vector<double> cs;
                                        for (int i = 0; i < numofNodes; i++) {
                                            cs.push_back(-r1 * dx[i] - r2 * dy[i] - r3 * dz[i]);
                                        }

                                        double fl = 0.0;
                                        for (int i = 0; i < numofNodes; i++) {
                                            fl += N[i] * cs[i] / Z;
                                        }

                                        dt = _bisection(fl, cs, N, dxdx);

                                        lambda_1 -= r1 * dt;
                                        lambda_2 -= r2 * dt;
                                        lambda_3 -= r3 * dt;

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)
					if ( fabs(L1.head() - lambda_1) < tol &&
					     fabs(L2.head() - lambda_2) < tol &&
					     fabs(L3.head() - lambda_3) < tol ) {
					  dt *= 1.2;
					  lambda_1 -= r1 * dt;
					  lambda_2 -= r2 * dt;
					  lambda_3 -= r3 * dt;
					}
					
					L1.head() = lambda_1;
					L2.head() = lambda_2;
					L3.head() = lambda_3;
                                        ++L1; ++L2; ++L3;
#endif

                                    } else { 
#ifdef _M4EXTREME_DEBUG_II
                                        cout << "Reset tolerance of the steepest descent solver" << endl;
#endif
                                        break;
                                    }

#ifdef _M4EXTREME_DEUBG_III
                                    cout << rTol << ": dt = " << dt << endl;
#endif

                                } while (gfcount++ < ResetCount); 
                                
                            if (gfcount > ResetCount) {
                                cout << "steepest descent solver failed >>>>>>>>>>>>>>> " << endl
                                        << "beta := " << beta << ";" << endl;

                                cout << "theta:=<" << qpx << "," << qpy << "," << qpz << ">;" << endl;

                                cout << "dx[1]:=<";
                                int i = 0;
                                for (; i < dx.size() - 1; i++) {
                                    cout << dx[i] << ",";
                                }
                                cout << dx[i] << ">;" << endl;

                                cout << "dx[2]:=<";
                                i = 0;
                                for (; i < dy.size() - 1; i++) {
                                    cout << dy[i] << ",";
                                }
                                cout << dy[i] << ">;" << endl;

                                cout << "dx[3]:=<";
                                i = 0;
                                for (; i < dz.size() - 1; i++) {
                                    cout << dz[i] << ",";
                                }
                                cout << dz[i] << ">;" << endl;

                                cout << "lambda_1 = " << lambda_1 << "\tlambda_2 = " << lambda_2 << "\tlambda_3 = " << lambda_3 << endl;
                                cout << "r1 = " << r1 << "\tr2 = " << r2 << "\tr3 = " << r3
                                        << "\terr=" << err << "\tdet = " << det << "\tZ = " << Z << endl;

                                cout << "j11 = " << j11 << "\tj12 = " << j12 << "\tj13 = " << j13 << endl
                                        << "j22 = " << j22 << "\tj23 = " << j23 << "\tj33 = " << j33 << endl;

                                cout << "i_j11 = " << i_j11 << "\ti_j12 = " << i_j12 << "\ti_j13 = " << i_j13 << endl
                                        << "i_j22 = " << i_j22 << "\ti_j23 = " << i_j23 << "\ti_j33 = " << i_j33 << endl;

                                cout << "Za : [\t";
                                for (unsigned int i = 0; i < numofNodes; i++) {
                                    cout << N[i] << "\t";
                                }
                                cout << "]" << endl;

                                throw ("steepest descent failed!");
                            }                            

                            }

                            --k;
                        } else {

                            double dl1 = (i_j11 * r1 + i_j12 * r2 + i_j13 * r3);
                            double dl2 = (i_j12 * r1 + i_j22 * r2 + i_j23 * r3);
                            double dl3 = (i_j13 * r1 + i_j23 * r2 + i_j33 * r3);

#if !defined(_M4EXTREME_NEWTON_BACKTRACKING_LINE_SEARCH)
                            // compute the exact line search step length
                            vector<double> cs;
                            for (int i = 0; i < numofNodes; i++) {
                                cs.push_back(-dl1 * dx[i] - dl2 * dy[i] - dl3 * dz[i]);
                            }

                            double fl = 0.0;
                            for (int i = 0; i < numofNodes; i++) {
                                fl += N[i] * cs[i] / Z;
                            }

                            double sm = _bisection(fl, cs, N, dxdx);

                            // update lambda -- ewton-Raphson method with exact line search
                            lambda_1 -= sm * dl1;
                            lambda_2 -= sm * dl2;
                            lambda_3 -= sm * dl3;
#else
                            // update lambda -- Newton-Raphson method with backtracking line search
                            double sm = 1.0, gamma = 0.5, alpha = 0.01;
                            double f0 = log(Z);
                            double df = alpha * (r1 * dl1 + r2 * dl2 + r3 * dl3);
                            double fold, fnew;
                            int kb = 0;

                            do {
                                fold = f0 - df * sm;

                                lambda_1 -= sm * dl1;
                                lambda_2 -= sm * dl2;
                                lambda_3 -= sm * dl3;

                                fnew = 0.0;
                                for (int i = 0; i < numofNodes; i++) {
                                    fnew += exp(-beta * dxdx[i] + lambda_1 * dx[i] + lambda_2 * dy[i] + lambda_3 * dz[i]);
                                }

                                fnew = log(fnew);

                                sm *= gamma;

                            } while (fnew < fold && kb++ < 100);
#endif
                        }

                    }

                } while (k++ < 20);

                cout << "Newton-Raphson solver failed >>>>>>>>>>>>>>> " << endl
                        << "beta := " << beta << ";" << endl
                        << "theta:=<" << qpx << "," << qpy << "," << qpz << ">;" << endl;

                cout << "dx[1]:=<";
                int i = 0;
                for (; i < dx.size() - 1; i++) {
                    cout << dx[i] << ",";
                }
                cout << dx[i] << ">;" << endl;

                cout << "dx[2]:=<";
                i = 0;
                for (; i < dy.size() - 1; i++) {
                    cout << dy[i] << ",";
                }
                cout << dy[i] << ">;" << endl;

                cout << "dx[3]:=<";
                i = 0;
                for (; i < dz.size() - 1; i++) {
                    cout << dz[i] << ",";
                }
                cout << dz[i] << ">;" << endl;

                cout << "lambda_1 = " << lambda_1 << "\tlambda_2 = " << lambda_2 << "\tlambda_3 = " << lambda_3 << endl;
                cout << "r1 = " << r1 << "\tr2 = " << r2 << "\tr3 = " << r3
                     << "\terr=" << err << "\tdet = " << det << "\tZ = " << Z << endl;


                cout << "i_j11 = " << i_j11 << "\ti_j12 = " << i_j12 << "\ti_j13 = " << i_j13 << endl
                        << "i_j22 = " << i_j22 << "\ti_j23 = " << i_j23 << "\ti_j33 = " << i_j33 << endl;

                cout << "Za : [\t";
                for (unsigned int i = 0; i < numofNodes; i++) {
                    cout << N[i] << "\t";
                }
                cout << "]" << endl;

                throw ("Newton Raphson failed!");

            }                      

            double Base::_bisection(double fl,
                    const vector<double> & cs,
                    const vector<double> & ps,
                    const vector<double> & dxdx) const {
                int BistMax = 100;
                int itBist = 0;

                double sl = 0.0, sr = 10.0, sm = 1.0;
                double fr = _computeLineSearchf(sr, cs, ps, dxdx);

                // correct the guess
                while (fl * fr > 0.0 && itBist++ < BistMax) {
                    sr *= 2.0;
                    fr = _computeLineSearchf(sr, cs, ps, dxdx);
                }

                // bisection method
                itBist = 0;
                while (sr - sl > 1.0e-8 && itBist++ < BistMax) {
                    sm = 0.5 * (sl + sr);
                    double fm = _computeLineSearchf(sm, cs, ps, dxdx);
                    if (fm * fl > 0.0) {
                        sl = sm;
                    } else {
                        sr = sm;
                    }
                }

                if (itBist < BistMax) {
                    return sm;
                } else {
                    return 1.0;
                }
            }

            double Base::_computeLineSearchf(double s, const vector<double> & cs,
                    const vector<double> & ps,
                    const vector<double> & dxdx) const {

                int numofNodes = ps.size();

                double z = 0.0;
                vector<double> N(numofNodes, 0.0);
                for (int i = 0; i < numofNodes; i++) {
                    N[i] = ps[i] * exp(s * cs[i]);
                    z += N[i];
                }

                double f = 0.0;
                for (int i = 0; i < numofNodes; i++) {
                    f += N[i] * cs[i] / z;
                }

                return f;
            }


#if defined (_M4EXTREME_VISCOUS_REGULARIZATION_)

            inline
            void Base::_get_psum(double** p, double * psum, int ndim) const
            {
		int i,j;
		double sum;               
                int mpts = ndim + 1;
		for (j=0;j<ndim;j++) {
			for (sum=0.0,i=0;i<mpts;i++)
				sum += p[i][j];
			psum[j]=sum;
		}
            }

	  /**
	   *
	   * 2D Solver
	   *
	   **/

            inline
            bool Base::_logZ(double & value,
                    double * Lambda,
                    double beta,
                    const vector<double> & dx,
                    const vector<double> & dy,
                    const vector<double> & dxdx) const {

                double Z = 0.0;
                size_t numofNodes = dx.size();

                for (size_t i = 0; i < numofNodes; ++i) {
                    Z += exp(-beta * dxdx[i] + Lambda[0] * dx[i] + Lambda[1] * dy[i]);
                }

                if ( Z > 0.0 && Z == Z ) {
                    value = log(Z);
                    return true;
                }
                else {
                    return false;
                }
            }

            double Base::_NMtry(double ** p, double * y, double * psum, int ihi, double fac,
                    double beta,
                    const vector<double> & dx,
                    const vector<double> & dy,
                    const vector<double> & dxdx) const {

                int j;
                double fac1, fac2, ytry;
                double TINY = 1.0e-12;

                int ndim = 2;
                double * ptry = new double[ndim];
                fac1 = (1.0 - fac) / ndim;
                fac2 = fac1 - fac;
                for (j = 0; j < ndim; j++)
                    ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;

                if ( _logZ(ytry, ptry, beta, dx, dy, dxdx) ) {
                    if (ytry < y[ihi]) {
                        y[ihi] = ytry;
                        for (j = 0; j < ndim; j++) {
                            psum[j] += ptry[j] - p[ihi][j];
                            p[ihi][j] = ptry[j];
                        }
                    }
                    
                    delete [] ptry;                    
                    return ytry;
                }
                else {
                    delete [] ptry;
                    return y[ihi] + TINY;
                }
            }


            void Base::_NelderMead(double & Z,
                    double beta,
                    double ftol,
                    const vector<double> & dx,
                    const vector<double> & dy,                    
                    const vector<double> & dxdx) {

                const int NMAX = 50;
                const double TINY = 1.0e-12;
                int i, ihi, ilo, inhi, j;
                double ysave, ytry;                

                if ( ftol < 1.0e-8 ) {
                    ftol = 1.0e-8;
                }
                
                int ndim = 2;
                int mpts = ndim + 1;

                double ** p = new double * [mpts];
                for ( size_t i = 0; i != mpts; ++i ) {
                    p[i] = new double[ndim];
                    for ( size_t j = 0; j != ndim; ++j ) p[i][j] = 0.0;
                }

                for ( size_t i = 1; i != mpts; ++i ) {
                    p[i][i-1] = (*_Lambda)[i-1];

                    if ( fabs(p[i][i-1]) > 1.0e32 || p[i][i-1] != p[i][i-1] ) { // check if lambda is inf or nan
                        p[i][i-1] = 1.0;
                    }
                }

                double * y = new double[mpts];
                if ( !_logZ(y[0], p[0], beta, dx, dy, dxdx) ) throw(1);

                for ( size_t i = 1; i != mpts; ++i ) {
                    int count = 0;
                    while ( !_logZ(y[i], p[i], beta, dx, dy, dxdx) && count++ < NMAX ) {
                        p[i][i-1] *= 0.5;
                    }

                    if ( count >= NMAX ) throw(1); // 2^100 ~ 10^30
                }

                double * psum = new double[ndim];
                int nfunk = 0;
                _get_psum(p, psum, ndim);
                for (;;) {
                    ilo = 0;
                    ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
                    for (i = 0; i < mpts; i++) {
                        if (y[i] <= y[ilo]) ilo = i;
                        if (y[i] > y[ihi]) {
                            inhi = ihi;
                            ihi = i;
                        } else if (y[i] > y[inhi] && i != ihi) inhi = i;
                    }
                    
                    if (fabs(y[ihi] - y[ilo]) < ftol || nfunk >= NMAX) {
                        Z = exp(y[ilo]);
                        for (i = 0; i < ndim; i++) (*_Lambda)[i] = p[ilo][i];
                        break;
                    }

                    nfunk += 2;
                    ytry = _NMtry(p, y, psum, ihi, -1.0, beta, dx, dy, dxdx);
                    if (ytry <= y[ilo])
                        ytry = _NMtry(p, y, psum, ihi, 2.0, beta, dx, dy, dxdx);
                    else if (ytry >= y[inhi]) {
                        ysave = y[ihi];
                        ytry = _NMtry(p, y, psum, ihi, 0.5, beta, dx, dy, dxdx);
                        if (ytry >= ysave) {
                            for (i = 0; i < mpts; i++) {
                                if (i != ilo) {
                                    for (j = 0; j < ndim; j++)
                                        p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
                                    _logZ(y[i], psum, beta, dx, dy, dxdx);
                                }
                            }
                            nfunk += ndim;
                            _get_psum(p, psum, ndim);
                        }
                    } else --nfunk;
                }

                delete [] psum;
                delete [] y;
                for (size_t i = 0; i != mpts; ++i ) {
                    delete [] p[i];
                }
                delete [] p;

                return;
            }

            void Base::_computeLambda2D_Regularized(const Set::VectorSpace::Vector &qp) {

                vector<double> dx, dy, dxdx;
                double qpx = qp[0];
                double qpy = qp[1];                

                int numofNodes = _x->size();

                double dxa, dya, dxa2;
		double max_dx = std::numeric_limits<double>::min();
                nodeset_type::const_iterator px;
                for (px = _x->begin(); px != _x->end(); px++) {
                    dxa = qpx - (px->second)[0];
                    dya = qpy - (px->second)[1];
                    dx.push_back(dxa);
                    dy.push_back(dya);

                    dxa2 = dxa * dxa + dya * dya;
                    dxdx.push_back(dxa2);
		    
		    if ( dxa2 > max_dx ) {
		      max_dx = dxa2;
		    }
                }

                double & lambda_1 = (*_Lambda)[0];
                double & lambda_2 = (*_Lambda)[1];

		double norm_dx = 4.0 * max_dx;
		max_dx = sqrt(max_dx);

                double beta = _Beta;
                double err;
                double r1 = 0.0, r2 = 0.0;
                double tol_r = 1.0e-8; // * max_dx;              
		double det;                
                double gamma = 4.0;
                double Z = 0.0;		
                vector<double> N(numofNodes, 0.0);
                double logZOld;

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)
                _loop_list L1, L2;
#endif

                int nDegenerater = 0;
                double ftol = 1.0e-4;
                while (nDegenerater++ < 5) {

                    int k = 0, nReset = 0;
                    double j11, j12, j22;

                    do {

                        Z = 0.0;
                        for (unsigned int i = 0; i < numofNodes; ++i) {
                            N[i] = exp(-beta * dxdx[i] + lambda_1 * dx[i] + lambda_2 * dy[i]);
                            Z += N[i];
                        }

                        nReset = 0;
                        double ftolloc = 1.0e-4;
                        while ( (Z < 1.0e-32 || Z > 1.0e32 || Z != Z) && nReset++ < 5 ) {                            
                            beta = gamma / norm_dx;
                            try {
                                _NelderMead(Z, beta, ftolloc, dx, dy, dxdx);
                            }
                            catch(...) {
                                //cerr << "NelderMead method failed to initialize" << endl;
                                goto EXCEPTION;
                            }
                            ftolloc *= 0.2;
                            gamma *= 0.5;
                        }

                        if (nReset >= 5) {
                            //cerr << "Z is zero or inf or nan @ beta = " << beta << endl;
                            goto EXCEPTION;
                        } else {
                            for (unsigned int i = 0; i < numofNodes; ++i) {
                                N[i] = exp(-beta * dxdx[i] + lambda_1 * dx[i] + lambda_2 * dy[i]);
                            }
                        }

                        j11 = 0.0;
                        j12 = 0.0;
                        j22 = 0.0;
                        double shape_a;

                        r1 = 0.0;
                        r2 = 0.0;
                        for (unsigned int i = 0; i < numofNodes; ++i) {

                            const double & dxa = dx[i];
                            const double & dya = dy[i];

                            shape_a = N[i] / Z;

                            r1 += dxa * shape_a;
                            r2 += dya * shape_a;

                            j11 += dxa * dxa * shape_a;
                            j12 += dxa * dya * shape_a;
                            j22 += dya * dya * shape_a;
                        }

                        j11 -= r1 * r1;
                        j12 -= r1 * r2;
                        j22 -= r2 * r2;

                        err = sqrt(r1 * r1 + r2 * r2);

                        if (err != err) {
                            //cout << "floating point is not a number" << endl;
                            goto EXCEPTION;
                        }

                        if (__builtin_expect(err < tol_r, 0)) {

                            for (unsigned int i = 0; i < numofNodes; i++) {
                                _NNodes.push_back(N[i] / Z);
                            }

                            det = j11 * j22 - j12 * j12;

                            (*_JInv)[0][0] = j22 / det;
                            (*_JInv)[1][1] = j11 / det;
                            (*_JInv)[0][1] = (*_JInv)[1][0] = -j12 / det;

                            return;
                        } else {

                            j11 += err;
                            j22 += err;                           
                            double det_h = j11 * j22 - j12 * j12;

                            double i_j11 = j22 / det_h;
                            double i_j12 =-j12 / det_h;
                            double i_j22 = j11 / det_h;

                            double dl1 = i_j11 * r1 + i_j12 * r2;
                            double dl2 = i_j12 * r1 + i_j22 * r2;

                            double dt = 1.0; //initial step size

                            lambda_1 -= dt * dl1;
                            lambda_2 -= dt * dl2;

                            logZOld = log(Z);

                            Z = 0.0;
                            for (unsigned int i = 0; i < numofNodes; ++i) {
                                Z += exp(-beta * dxdx[i] + lambda_1 * dx[i] + lambda_2 * dy[i]);
                            }

                            double logZ = log(Z);
                            if (logZ >= logZOld) {
                                // compute the exact line search step length
                                vector<double> cs;
                                for (int i = 0; i < numofNodes; ++i) {
                                    cs.push_back(-dl1 * dx[i] - dl2 * dy[i]);
                                }

                                double fl = 0.0;
                                for (int i = 0; i < numofNodes; ++i) {
                                    fl += N[i] * cs[i] / Z;
                                }

                                dt -= _bisection(fl, cs, N, dxdx);
                                lambda_1 += dt * dl1;
                                lambda_2 += dt * dl2;                                
                            }

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)
                            if (fabs(L1.head() - lambda_1) < tol_r &&
				fabs(L2.head() - lambda_2) < tol_r) {
                                dt *= 0.2;
                                lambda_1 -= dl1 * dt;
                                lambda_2 -= dl2 * dt;
                            }

                            L1.head() = lambda_1;
                            L2.head() = lambda_2;
                            ++L1;
                            ++L2;
#endif
                        }

                    } while (k++ < 100);

                    beta = gamma / norm_dx;
                    gamma *= 0.7;
                    ftol  *= 0.2;
                    try {
                        _NelderMead(Z, beta, ftol, dx, dy, dxdx);
                    } catch (...) {
                        //cerr << "NelderMead method failed to initialize" << endl;
                        goto EXCEPTION;
                    }

                    tol_r = 1.0e-4;
                }

                EXCEPTION:
#ifdef _M4EXTREME_DEBUG_II
                cout << "Newton-Raphson solver failed >>>>>>>>>>>>>>> " << endl
                        << "beta := " << beta << ";" << endl
                        << "ftol := " << ftol << ";" << endl
                        << "theta:=<" << qpx << "," << qpy << ">;" << endl;

                cout << "dx[1]:=<";
                int i = 0;
                for (; i < dx.size() - 1; i++) {
                    cout << dx[i] << ",";
                }
                cout << dx[i] << ">;" << endl;

                cout << "dx[2]:=<";
                i = 0;
                for (; i < dy.size() - 1; i++) {
                    cout << dy[i] << ",";
                }
                cout << dy[i] << ">;" << endl;

                cout << "lambda_1 = " << lambda_1 << "\tlambda_2 = " << lambda_2 << endl;
                cout << "r1 = " << r1 << "\tr2 = " << r2
                     << "\terr=" << err << "\tdet = " << det << "\tZ = " << Z << endl;

                cout << "Za : [\t";
                for (unsigned int i = 0; i < numofNodes; i++) {
                    cout << N[i] << "\t";
                }
                cout << "]" << endl;
#endif
                //cerr << "Newton Raphson failed@MaxEnt::Base" << endl;
                throw (1);

            }


	  /**
	   *
	   * 3D Solver
	   *
	   **/

            inline
            bool Base::_logZ(double & value,
                    double * Lambda,
                    double beta,
                    const vector<double> & dx,
                    const vector<double> & dy,
                    const vector<double> & dz,
                    const vector<double> & dxdx) const {

                double Z = 0.0;
                size_t numofNodes = dx.size();

                for (size_t i = 0; i < numofNodes; ++i) {
                    Z += exp(-beta * dxdx[i] + Lambda[0] * dx[i] + Lambda[1] * dy[i] + Lambda[2] * dz[i]);
                }

                if ( Z > 0.0 && Z == Z ) {
                    value = log(Z);
                    return true;
                }
                else {
                    return false;
                }
            }

            double Base::_NMtry(double ** p, double * y, double * psum, int ihi, double fac,
                    double beta,
                    const vector<double> & dx,
                    const vector<double> & dy,
                    const vector<double> & dz,
                    const vector<double> & dxdx) const {

                int j;
                double fac1, fac2, ytry;
                double TINY = 1.0e-12;

                int ndim = 3;
                double * ptry = new double[ndim];
                fac1 = (1.0 - fac) / ndim;
                fac2 = fac1 - fac;
                for (j = 0; j < ndim; j++)
                    ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;

                if ( _logZ(ytry, ptry, beta, dx, dy, dz, dxdx) ) {
                    if (ytry < y[ihi]) {
                        y[ihi] = ytry;
                        for (j = 0; j < ndim; j++) {
                            psum[j] += ptry[j] - p[ihi][j];
                            p[ihi][j] = ptry[j];
                        }
                    }
                    
                    delete [] ptry;                    
                    return ytry;
                }
                else {
                    delete [] ptry;
                    return y[ihi] + TINY;
                }
            }


            void Base::_NelderMead(double & Z,
                    double beta,
                    double ftol,
                    const vector<double> & dx,
                    const vector<double> & dy,
                    const vector<double> & dz,
                    const vector<double> & dxdx) {

                const int NMAX = 50;
                const double TINY = 1.0e-12;
                int i, ihi, ilo, inhi, j;
                double ysave, ytry;                

                if ( ftol < 1.0e-8 ) {
                    ftol = 1.0e-8;
                }
                
                int mpts = 4;
                int ndim = 3;

                double ** p = new double * [mpts];
                for ( size_t i = 0; i != mpts; ++i ) {
                    p[i] = new double[ndim];
                    for ( size_t j = 0; j != ndim; ++j ) p[i][j] = 0.0;
                }

                for ( size_t i = 1; i != mpts; ++i ) {
                    p[i][i-1] = (*_Lambda)[i-1];

                    if ( fabs(p[i][i-1]) > 1.0e32 || p[i][i-1] != p[i][i-1] ) { // check if lambda is inf or nan
                        p[i][i-1] = 1.0;
                    }
                }

                double * y = new double[mpts];
                if ( !_logZ(y[0], p[0], beta, dx, dy, dz, dxdx) ) throw(1);

                for ( size_t i = 1; i != mpts; ++i ) {
                    int count = 0;
                    while ( !_logZ(y[i], p[i], beta, dx, dy, dz, dxdx) && count++ < NMAX ) {
                        p[i][i-1] *= 0.5;
                    }

                    if ( count >= NMAX ) throw(1); // 2^100 ~ 10^30
                }

                double * psum = new double[ndim];
                int nfunk = 0;
                _get_psum(p, psum, ndim);
                for (;;) {
                    ilo = 0;
                    ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
                    for (i = 0; i < mpts; i++) {
                        if (y[i] <= y[ilo]) ilo = i;
                        if (y[i] > y[ihi]) {
                            inhi = ihi;
                            ihi = i;
                        } else if (y[i] > y[inhi] && i != ihi) inhi = i;
                    }
                    
                    if (fabs(y[ihi] - y[ilo]) < ftol || nfunk >= NMAX) {
                        Z = exp(y[ilo]);
                        for (i = 0; i < ndim; i++) (*_Lambda)[i] = p[ilo][i];
                        break;
                    }

                    nfunk += 2;
                    ytry = _NMtry(p, y, psum, ihi, -1.0, beta, dx, dy, dz, dxdx);
                    if (ytry <= y[ilo])
                        ytry = _NMtry(p, y, psum, ihi, 2.0, beta, dx, dy, dz, dxdx);
                    else if (ytry >= y[inhi]) {
                        ysave = y[ihi];
                        ytry = _NMtry(p, y, psum, ihi, 0.5, beta, dx, dy, dz, dxdx);
                        if (ytry >= ysave) {
                            for (i = 0; i < mpts; i++) {
                                if (i != ilo) {
                                    for (j = 0; j < ndim; j++)
                                        p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
                                    _logZ(y[i], psum, beta, dx, dy, dz, dxdx);
                                }
                            }
                            nfunk += ndim;
                            _get_psum(p, psum, ndim);
                        }
                    } else --nfunk;
                }

                delete [] psum;
                delete [] y;
                for (size_t i = 0; i != mpts; ++i ) {
                    delete [] p[i];
                }
                delete [] p;

                return;
            }

            void Base::_computeLambda3D_Regularized(const Set::VectorSpace::Vector &qp) {

                vector<double> dx, dy, dz, dxdx;
                double qpx = qp[0];
                double qpy = qp[1];
                double qpz = qp[2];

                int numofNodes = _x->size();

                double dxa, dya, dza;
                double dxa2 = 0.0, norm_dx = 0.0, min_dxa2 = std::numeric_limits<double>::max();
                size_t idx = 0, norm_idx = 0;
                nodeset_type::const_iterator px;
                for (px = _x->begin(); px != _x->end(); px++, idx++) {
                    dxa = qpx - (px->second)[0];
                    dya = qpy - (px->second)[1];
                    dza = qpz - (px->second)[2];
                    dx.push_back(dxa);
                    dy.push_back(dya);
                    dz.push_back(dza);

                    dxa2 = dxa * dxa + dya * dya + dza * dza;
                    dxdx.push_back(dxa2);

                    if (norm_dx < dxa2) {
                        norm_dx = dxa2;
                        norm_idx = idx;
                    }

                    if ( min_dxa2 > dxa2 ) {
                        min_dxa2 = dxa2;
                    }
                }

                norm_dx *= 4.0; // norm_dx = h^2 where h is the meshsize
                if ( norm_dx / min_dxa2 > 100.0 ) {
                    //To Do
                }

                if ( norm_dx > 1.0e12 ) {
                    //cerr << norm_dx << "is too big(>1.0e12)@MaxEnt::Base" << endl;
                    throw(1);
                }

                double & lambda_1 = (*_Lambda)[0];
                double & lambda_2 = (*_Lambda)[1];
                double & lambda_3 = (*_Lambda)[2];

                double beta = _Beta;
                double err;
                double r1 = 0.0, r2 = 0.0, r3 = 0.0;
                double tol_r = 1.0e-8;                
		double det;                
                double gamma = 4.0;
                double Z = 0.0;		
                vector<double> N(numofNodes, 0.0);
                double logZOld;

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)
                _loop_list L1, L2, L3;
#endif

                int nDegenerater = 0;
                double ftol = 1.0e-4;
                while (nDegenerater++ < 5) {

                    int k = 0, nReset = 0;
                    double j11, j12, j13, j22, j23, j33;

                    do {

                        Z = 0.0;
                        for (unsigned int i = 0; i < numofNodes; ++i) {
                            N[i] = exp(-beta * dxdx[i] + lambda_1 * dx[i] + lambda_2 * dy[i] + lambda_3 * dz[i]);
                            Z += N[i];
                        }

                        nReset = 0;
                        double ftolloc = 1.0e-4;
                        while ( (Z < 1.0e-32 || Z > 1.0e32 || Z != Z) && nReset++ < 5 ) {                            
                            beta = gamma / norm_dx;
                            try {
                                _NelderMead(Z, beta, ftolloc, dx, dy, dz, dxdx);
                            }
                            catch(...) {
                                //cerr << "NelderMead method failed to initialize" << endl;
                                goto EXCEPTION;
                            }
                            ftolloc *= 0.2;
                            gamma *= 0.5;
                        }

                        if (nReset >= 5) {
                            //cerr << "Z is zero or inf or nan @ beta = " << beta << endl;
                            goto EXCEPTION;
                        } else {
                            for (unsigned int i = 0; i < numofNodes; ++i) {
                                N[i] = exp(-beta * dxdx[i] + lambda_1 * dx[i] + lambda_2 * dy[i] + lambda_3 * dz[i]);
                            }
                        }

                        j11 = 0.0;
                        j12 = 0.0;
                        j22 = 0.0;
                        j13 = 0.0;
                        j23 = 0.0;
                        j33 = 0.0;
                        double shape_a;

                        r1 = 0.0;
                        r2 = 0.0;
                        r3 = 0.0;
                        for (unsigned int i = 0; i < numofNodes; ++i) {

                            const double & dxa = dx[i];
                            const double & dya = dy[i];
                            const double & dza = dz[i];

                            shape_a = N[i] / Z;

                            r1 += dxa * shape_a;
                            r2 += dya * shape_a;
                            r3 += dza * shape_a;

                            j11 += dxa * dxa * shape_a;
                            j12 += dxa * dya * shape_a;
                            j22 += dya * dya * shape_a;
                            j13 += dxa * dza * shape_a;
                            j23 += dya * dza * shape_a;
                            j33 += dza * dza * shape_a;
                        }

                        j11 -= r1 * r1;
                        j12 -= r1 * r2;
                        j22 -= r2 * r2;
                        j13 -= r1 * r3;
                        j23 -= r2 * r3;
                        j33 -= r3 * r3;

                        err = sqrt(r1 * r1 + r2 * r2 + r3 * r3);

                        if (err != err) {
                            //cout << "floating point is not a number" << endl;
                            goto EXCEPTION;
                        }

                        if (__builtin_expect(err < tol_r, 0)) {

                            for (unsigned int i = 0; i < numofNodes; i++) {
                                _NNodes.push_back(N[i] / Z);
                            }

                            det = j11 * j22 * j33 - j11 * j23 * j23 + 2.0 * j12 * j23 * j13
                                    - j12 * j12 * j33 - j22 * j13 * j13;

                            (*_JInv)[0][0] = (j22 * j33 - j23 * j23) / det;
                            (*_JInv)[1][1] = (-j13 * j13 + j11 * j33) / det;
                            (*_JInv)[0][1] = (*_JInv)[1][0] = (j23 * j13 - j12 * j33) / det;
                            (*_JInv)[2][2] = (-j12 * j12 + j11 * j22) / det;
                            (*_JInv)[0][2] = (*_JInv)[2][0] = (j12 * j23 - j22 * j13) / det;
                            (*_JInv)[2][1] = (*_JInv)[1][2] = (-j11 * j23 + j13 * j12) / det;

                            return;
                        } else {

                            j11 += err;
                            j22 += err;
                            j33 += err;
                            double det_h = j11 * j22 * j33 - j11 * j23 * j23 + 2.0 * j12 * j23 * j13
                                    - j12 * j12 * j33 - j22 * j13 * j13;

                            double i_j11 = (j22 * j33 - j23 * j23) / det_h;
                            double i_j12 = (j23 * j13 - j12 * j33) / det_h;
                            double i_j13 = (j12 * j23 - j22 * j13) / det_h;
                            double i_j22 = (-j13 * j13 + j11 * j33) / det_h;
                            double i_j23 = (-j11 * j23 + j13 * j12) / det_h;
                            double i_j33 = (-j12 * j12 + j11 * j22) / det_h;

                            double dl1 = i_j11 * r1 + i_j12 * r2 + i_j13 * r3;
                            double dl2 = i_j12 * r1 + i_j22 * r2 + i_j23 * r3;
                            double dl3 = i_j13 * r1 + i_j23 * r2 + i_j33 * r3;

                            double dt = 1.0; //initial step size

                            lambda_1 -= dt * dl1;
                            lambda_2 -= dt * dl2;
                            lambda_3 -= dt * dl3;

                            logZOld = log(Z);

                            Z = 0.0;
                            for (unsigned int i = 0; i < numofNodes; ++i) {
                                Z += exp(-beta * dxdx[i] + lambda_1 * dx[i] + lambda_2 * dy[i] + lambda_3 * dz[i]);
                            }

                            double logZ = log(Z);
                            if (logZ >= logZOld) {
                                // compute the exact line search step length
                                vector<double> cs;
                                for (int i = 0; i < numofNodes; ++i) {
                                    cs.push_back(-dl1 * dx[i] - dl2 * dy[i] - dl3 * dz[i]);
                                }

                                double fl = 0.0;
                                for (int i = 0; i < numofNodes; ++i) {
                                    fl += N[i] * cs[i] / Z;
                                }

                                dt -= _bisection(fl, cs, N, dxdx);
                                lambda_1 += dt * dl1;
                                lambda_2 += dt * dl2;
                                lambda_3 += dt * dl3;
                            }

#if defined(_M4EXTREME_MAXENT_SYMMETRIC_CASE_)
                            if (fabs(L1.head() - lambda_1) < tol_r &&
                                    fabs(L2.head() - lambda_2) < tol_r &&
                                    fabs(L3.head() - lambda_3) < tol_r) {
                                dt *= 0.2;
                                lambda_1 -= dl1 * dt;
                                lambda_2 -= dl2 * dt;
                                lambda_3 -= dl3 * dt;
                            }

                            L1.head() = lambda_1;
                            L2.head() = lambda_2;
                            L3.head() = lambda_3;
                            ++L1;
                            ++L2;
                            ++L3;
#endif
                        }

                    } while (k++ < 100);

                    beta = gamma / norm_dx;
                    gamma *= 0.7;
                    ftol  *= 0.2;
                    try {
                        _NelderMead(Z, beta, ftol, dx, dy, dz, dxdx);
                    } catch (...) {
                        //cerr << "NelderMead method failed to initialize" << endl;
                        goto EXCEPTION;
                    }

                    tol_r = 1.0e-4;
                }

                EXCEPTION:
#ifdef _M4EXTREME_DEBUG_II
                cout << "Newton-Raphson solver failed >>>>>>>>>>>>>>> " << endl
                        << "beta := " << beta << ";" << endl
                        << "ftol := " << ftol << ";" << endl
                        << "theta:=<" << qpx << "," << qpy << "," << qpz << ">;" << endl;

                cout << "dx[1]:=<";
                int i = 0;
                for (; i < dx.size() - 1; i++) {
                    cout << dx[i] << ",";
                }
                cout << dx[i] << ">;" << endl;

                cout << "dx[2]:=<";
                i = 0;
                for (; i < dy.size() - 1; i++) {
                    cout << dy[i] << ",";
                }
                cout << dy[i] << ">;" << endl;

                cout << "dx[3]:=<";
                i = 0;
                for (; i < dz.size() - 1; i++) {
                    cout << dz[i] << ",";
                }
                cout << dz[i] << ">;" << endl;

                cout << "lambda_1 = " << lambda_1 << "\tlambda_2 = " << lambda_2 << "\tlambda_3 = " << lambda_3 << endl;
                cout << "r1 = " << r1 << "\tr2 = " << r2 << "\tr3 = " << r3
                     << "\terr=" << err << "\tdet = " << det << "\tZ = " << Z << endl;

                cout << "Za : [\t";
                for (unsigned int i = 0; i < numofNodes; i++) {
                    cout << N[i] << "\t";
                }
                cout << "]" << endl;
#endif
                //cerr << "Newton Raphson failed@MaxEnt::Base" << endl;
                throw (1);

            }

#endif

      
            //////////////////////////////////////////////////////////////////////
            // Class Shape<0>
            //////////////////////////////////////////////////////////////////////

            Shape < 0 > ::Shape(Base *B_) : B(B_) {
            }

            Shape < 0 > ::~Shape() {
            }

            map<Set::Manifold::Point *, double>
            Shape < 0 > ::operator () (const Set::VectorSpace::Vector &y) const {
                try {
                    (*B)(y);
                } catch (int code) {
                    throw (code);
                }

		const map<Set::Manifold::Point *, Set::VectorSpace::Vector> * x = B->GetNodes();
		map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;
                const vector<double> & BN = B->GetN();

		map<Set::Manifold::Point *, double> shapefunction;
		int idx = 0;
                for (px = x->begin(); px != x->end(); px++) {
		  shapefunction.insert( make_pair(px->first, BN[idx++]) );
		}

		return shapefunction;
            }

            //////////////////////////////////////////////////////////////////////
            // Class Shape<1>
            //////////////////////////////////////////////////////////////////////

            Shape < 1 > ::Shape(Base *B_) : B(B_) {
            }

            Shape < 1 > ::Shape(const Shape < 0 > &rhs) : B(rhs.B) {
            }

            Shape < 1 > ::~Shape() {
            }

            map<Set::Manifold::Point *, Set::VectorSpace::Vector>
            Shape < 1 > ::operator () (const Set::VectorSpace::Vector &y) const {
                try {
                    (*B)(y);
                } catch (int code) {
                    throw (code);
                }

                const vector<double> & N = B->GetN();
                const Set::VectorSpace::Hom & JInv = *(B->GetJInverse());

                map<Set::Manifold::Point *, Set::VectorSpace::Vector> DN;

                const map<Set::Manifold::Point *, Set::VectorSpace::Vector> * x = B->GetNodes();
                map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;

		int idx = 0;
                for (px = x->begin(); px != x->end(); px++) {
                    Set::VectorSpace::Vector Dxa = px->second - y;
                    DN.insert(make_pair(px->first, N[idx++] * JInv(Dxa)));
                }

                return DN;
            }

            //////////////////////////////////////////////////////////////////////
            // Class Jet<0>
            //////////////////////////////////////////////////////////////////////

            Jet < 0 > ::Jet() {
            }

            Jet < 0 > ::Jet(Element::Interpolation::MaxEnt::Base * B) : _BNode(B) {
            }

            Jet < 0 > ::~Jet() {
            }

            void Jet < 0 > ::operator () (const domain_type & y, shape_type & N, shape_derivative_type & DN) const {

                try {
                    (*_BNode)(y);
                } catch (...) {
                    //cerr << "Caculating shape functions failed in Jet<0>" << endl;
		    throw(1);
                }

                const  vector<double> & BN = _BNode->GetN();
                const Set::VectorSpace::Hom & JInv = *(_BNode->GetJInverse());

                const map<Set::Manifold::Point *, Set::VectorSpace::Vector> * x = _BNode->GetNodes();
                map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator px;

                if ( !N.empty() ) N.clear();
                if ( !DN.empty() ) DN.clear();

		int idx = 0;
		int dim = y.size();
		Set::VectorSpace::Vector dshape_a(dim);
                for (px = x->begin(); px != x->end(); px++) {
                    double shape_a = BN[idx++];
                    N.insert(shape_type::value_type(px->first, shape_a));

                    Set::VectorSpace::Vector Dxa = px->second - y;
		    for ( int i = 0; i < dim; ++i ) {
		      double & dshape_ai = dshape_a[i];
		      dshape_ai = 0.0;
		      for ( int j = 0; j < dim; ++j ) {
			dshape_ai += JInv[i][j] * Dxa[j];
		      }
		      dshape_ai *= shape_a; 
		    }

                    DN.insert( shape_derivative_type::value_type(px->first, dshape_a) );
                }

#ifdef _M4EXTREME_DEBUG_III
		if ( /*_gPrint*/ true ) {
		  map<Set::Manifold::Point *, double>::const_iterator pN;
		  map<Set::Manifold::Point *, Set::VectorSpace::Vector>::const_iterator pDN;
		  cout << "New shape functions [";
		  for (pN = N.begin(); pN != N.end(); pN++) {
                    cout << "(" << _gIndexSet.find(pN->first)->second << ")" << pN->second << "\t";
		  }
		  cout << "] New lambda is [" << _BNode->GetLambda() << "]" << endl;
		  
		  cout << "New shape derivatives are\t";
		  for (pDN = DN.begin(); pDN != DN.end(); pDN++) {
                    cout << "(" << _gIndexSet.find(pDN->first)->second << ")" << "[" << pDN->second << "]\t";
		  }
		  cout << std::endl;
		}
#endif            

                return;

            }

        }

    }

}
