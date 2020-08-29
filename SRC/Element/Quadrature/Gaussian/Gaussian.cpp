// Gaussian.cpp: implementation of the Rule class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include "./Gaussian.h"

namespace Element
{
namespace Quadrature
{
namespace Gaussian
{
//////////////////////////////////////////////////////////////////////
// Class Rule
//////////////////////////////////////////////////////////////////////

Rule::Rule() {}

Rule::~Rule() {}

Rule::Rule(
	const unsigned int &rhs_n, 
	const unsigned int &rhs_N)
	: n(rhs_n), N(rhs_N) 
{
	assert(n > 0);
}

vector<Set::VectorSpace::Vector> 
Rule::GetQ() const
{
	vector<Set::VectorSpace::Vector> P1;

	switch (N)
	{
		case 0:  P1 = Q1(); break;
		case 1:  P1 = Q1(); break;
		case 2:  P1 = Q3(); break;
		case 3:  P1 = Q3(); break;
		case 4:  P1 = Q5(); break;
		case 5:  P1 = Q5(); break;
		case 6:  P1 = Q7(); break;
		case 7:  P1 = Q7(); break;
		case 8:  P1 = Q9(); break;
		case 9:  P1 = Q9(); break;
		case 10: P1 = Q11(); break;
		case 11: P1 = Q11(); break;
		case 12:  P1 = Q13(); break;
		case 13:  P1 = Q13(); break;
		case 14:  P1 = Q15(); break;
		case 15:  P1 = Q15(); break;
		case 16:  P1 = Q17(); break;
		case 17:  P1 = Q17(); break;
		case 18:  P1 = Q19(); break;
		case 19:  P1 = Q19(); break;
		case 20:  P1 = Q21(); break;
		case 21:  P1 = Q21(); break;
		case 22:  P1 = Q23(); break;
		case 23:  P1 = Q23(); break;
		case 24:  P1 = Q25(); break;
		case 25:  P1 = Q25(); break;
		case 26:  P1 = Q27(); break;
		case 27:  P1 = Q27(); break;
		case 28:  P1 = Q29(); break;
		case 29:  P1 = Q29(); break;
		case 30:  P1 = Q31(); break;
		case 31:  P1 = Q31(); break;

		default: throw(0);
	}

	if (n == 1) return P1;

	unsigned int i, j, k, l, lpp;
	vector<Set::VectorSpace::Vector> P=P1;
	vector<Set::VectorSpace::Vector> Q;
	for (l=1; l<n; l++)
	{
		lpp = l+1; Q.clear();
		for (j=0; j<P.size(); j++)
		{
			Set::VectorSpace::Vector A = P[j];
			for (i=0; i<P1.size(); i++)
			{
				Set::VectorSpace::Vector B(lpp); 
				for (k=0; k<l; k++) B[k] = A[k];
				B[l] = P1[i][0]; Q.push_back(B);
			}
		}
		P.clear(); P=Q;
	}

	return Q;
}

vector<double> 
Rule::GetW() const
{
	assert(n > 0);

	vector<double> P1;

	switch (N)
	{
		case 0:  P1 = W1(); break;
		case 1:  P1 = W1(); break;
		case 2:  P1 = W3(); break;
		case 3:  P1 = W3(); break;
		case 4:  P1 = W5(); break;
		case 5:  P1 = W5(); break;
		case 6:  P1 = W7(); break;
		case 7:  P1 = W7(); break;
		case 8:  P1 = W9(); break;
		case 9:  P1 = W9(); break;
		case 10: P1 = W11(); break;
		case 11: P1 = W11(); break;
		case 12:  P1 = W13(); break;
		case 13:  P1 = W13(); break;
		case 14:  P1 = W15(); break;
		case 15:  P1 = W15(); break;
		case 16:  P1 = W17(); break;
		case 17:  P1 = W17(); break;
		case 18:  P1 = W19(); break;
		case 19:  P1 = W19(); break;
		case 20:  P1 = W21(); break;
		case 21:  P1 = W21(); break;
		case 22:  P1 = W23(); break;
		case 23:  P1 = W23(); break;
		case 24:  P1 = W25(); break;
		case 25:  P1 = W25(); break;
		case 26:  P1 = W27(); break;
		case 27:  P1 = W27(); break;
		case 28:  P1 = W29(); break;
		case 29:  P1 = W29(); break;
		case 30: P1 = W31(); break;
		case 31: P1 = W31(); break;
		default: throw(0);
	}

	if (n == 1) return P1;

	unsigned int i, j, l, lpp;
	vector<double> P=P1;
	vector<double> Q;
	for (l=1; l<n; l++)
	{
		lpp = l+1; Q.clear();
		for (j=0; j<P.size(); j++)
		{
			double A = P[j];
			for (i=0; i<P1.size(); i++)
			{
				Q.push_back(P[j]*P1[i]);
			}
		}
		P.clear(); P=Q;
	}

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q1() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(1,Q0);

	Q[0][0] =   0.0e0;

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q3() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(2,Q0);

	Q[0][0] = - 0.577350269189626e0;
	Q[1][0] =   0.577350269189626e0;

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q5() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(3,Q0);

	Q[0][0] = - 0.774596669241483e0;
	Q[1][0] =   0.0e0;
	Q[2][0] =   0.774596669241483e0;

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q7() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(4,Q0);

	Q[0][0] = - 0.861136311594053e0;
	Q[1][0] = - 0.339981043584856e0;
	Q[2][0] =   0.339981043584856e0;
	Q[3][0] =   0.861136311594053e0;

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q9() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(5,Q0);

	Q[0][0] = - 0.906179845938664e0;
	Q[1][0] = - 0.538469310105683e0;
	Q[2][0] =   0.0e0;
	Q[3][0] =   0.538469310105683e0;
	Q[4][0] =   0.906179845938664e0;

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q11() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(6,Q0);

	Q[0][0] = - 0.932469514203152e0;
	Q[1][0] = - 0.661209386466265e0;
	Q[2][0] = - 0.238619186083197e0;
	Q[3][0] =   0.238619186083197e0;
	Q[4][0] =   0.661209386466265e0;
	Q[5][0] =   0.932469514203152e0;

	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q13() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(7,Q0);

        Q[0][0] = - 0.949107912342759e0; 
        Q[1][0] = - 0.741531185599394e0; 
        Q[2][0] = - 0.405845151377397e0; 
        Q[3][0] =   0.0e0;                
        Q[4][0] =   0.405845151377397e0;  
        Q[5][0] =   0.741531185599394e0;  
        Q[6][0] =   0.949107912342759e0;  
            
	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q15() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(8,Q0);

        Q[0][0] = - 0.960289856497536e0;
        Q[1][0] = - 0.796666477413627e0;
        Q[2][0] = - 0.525532409916329e0;
        Q[3][0] = - 0.183434642495650e0;
        Q[4][0] =   0.183434642495650e0; 
        Q[5][0] =   0.525532409916329e0; 
        Q[6][0] =   0.796666477413627e0; 
        Q[7][0] =   0.960289856497536e0; 
            
	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q17() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(9,Q0);

        Q[0][0] = - 0.968160239507626e0; 
        Q[1][0] = - 0.836031107326636e0; 
        Q[2][0] = - 0.613371432700590e0; 
        Q[3][0] = - 0.324253423403809e0; 
        Q[4][0] =   0.0e0;                
        Q[5][0] =   0.324253423403809e0;  
        Q[6][0] =   0.613371432700590e0;  
        Q[7][0] =   0.836031107326636e0;  
        Q[8][0] =   0.968160239507626e0;  
            
	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q19() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(10,Q0);

        Q[0][0] = - 0.973906528517172e0; 
        Q[1][0] = - 0.865063366688985e0; 
        Q[2][0] = - 0.679409568299024e0; 
        Q[3][0] = - 0.433395394129247e0; 
        Q[4][0] = - 0.148874338981631e0; 
        Q[5][0] =   0.148874338981631e0;  
        Q[6][0] =   0.433395394129247e0;  
        Q[7][0] =   0.679409568299024e0;  
        Q[8][0] =   0.865063366688985e0;  
        Q[9][0] =   0.973906528517172e0;  
            
	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q21() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(11,Q0);

        Q[0][0] = - 0.978228658146057e0; 
        Q[1][0] = - 0.887062599768095e0; 
        Q[2][0] = - 0.730152005574049e0; 
        Q[3][0] = - 0.519096129206812e0; 
        Q[4][0] = - 0.269543155952345e0; 
        Q[5][0] =   0.0e0;                
        Q[6][0] =   0.269543155952345e0;  
        Q[7][0] =   0.519096129206812e0; 
        Q[8][0] =   0.730152005574049e0; 
        Q[9][0] =   0.887062599768095e0; 
        Q[10][0] =  0.978228658146057e0; 
            
	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q23() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(12,Q0);

        Q[0][0] = - 0.981560634246719e0;
        Q[1][0] = - 0.904117256370475e0;
        Q[2][0] = - 0.769902674194305e0;
        Q[3][0] = - 0.587317954286617e0;
        Q[4][0] = - 0.367831498918180e0;
        Q[5][0] = - 0.125233408511469e0;
        Q[6][0] =   0.125233408511469e0; 
        Q[7][0] =   0.367831498918180e0; 
        Q[8][0] =   0.587317954286617e0; 
        Q[9][0] =   0.769902674194305e0; 
        Q[10][0] =  0.904117256370475e0; 
        Q[11][0] =  0.981560634246719e0; 
            
	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q25() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(13,Q0);

        Q[0][0] = - 0.984183054718588e0; 
        Q[1][0] = - 0.917598399222978e0; 
        Q[2][0] = - 0.801578090733310e0; 
        Q[3][0] = - 0.642349339440340e0; 
        Q[4][0] = - 0.448492751036447e0; 
        Q[5][0] = - 0.230458315955135e0; 
        Q[6][0] =   0.0e0;                
        Q[7][0] =   0.230458315955135e0;  
        Q[8][0] =   0.448492751036447e0;  
        Q[9][0] =   0.642349339440340e0;  
        Q[10][0] =  0.801578090733310e0;  
        Q[11][0] =  0.917598399222978e0;  
        Q[12][0] =  0.984183054718588e0;  
            
	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q27() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(14,Q0);

        Q[0][0] = - 0.986283808696812e0;  
        Q[1][0] = - 0.928434883663574e0;  
        Q[2][0] = - 0.827201315069765e0;  
        Q[3][0] = - 0.687292904811685e0;  
        Q[4][0] = - 0.515248636358154e0;  
        Q[5][0] = - 0.319112368927890e0;  
        Q[6][0] = - 0.108054948707344e0;  
        Q[7][0] =   0.108054948707344e0;   
        Q[8][0] =   0.319112368927890e0;   
        Q[9][0] =   0.515248636358154e0;           
        Q[10][0] =  0.687292904811685e0;           
        Q[11][0] =  0.827201315069765e0;   
        Q[12][0] =  0.928434883663574e0;   
        Q[13][0] =  0.986283808696812e0;     
            
	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q29() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(15,Q0);

        Q[0][0] = - 0.987992518020485e0;  
        Q[1][0] = - 0.937273392400706e0;  
        Q[2][0] = - 0.848206583410427e0;  
        Q[3][0] = - 0.724417731360170e0;  
        Q[4][0] = - 0.570972172608539e0;  
        Q[5][0] = - 0.394151347077563e0;  
        Q[6][0] = - 0.201194093997435e0;  
        Q[7][0] =   0.0e0;                 
        Q[8][0] =   0.201194093997435e0;   
        Q[9][0] =   0.394151347077563e0;           
        Q[10][0] =  0.570972172608539e0;           
        Q[11][0] =  0.724417731360170e0;   
        Q[12][0] =  0.848206583410427e0;   
        Q[13][0] =  0.937273392400706e0;   
        Q[14][0] =  0.987992518020485e0;       
            
	return Q;
}

const vector<Set::VectorSpace::Vector>
Rule::Q31() const
{
	Set::VectorSpace::Vector Q0(1);
	vector<Set::VectorSpace::Vector> Q(16,Q0);

        Q[0][0] = - 0.989400934991650e0;  
        Q[1][0] = - 0.944575023073233e0;  
        Q[2][0] = - 0.865631202387832e0;  
        Q[3][0] = - 0.755404408355003e0;  
        Q[4][0] = - 0.617876244402644e0;  
        Q[5][0] = - 0.458016777657227e0;  
        Q[6][0] = - 0.281603550779259e0;  
        Q[7][0] = - 0.095012509837637e0;  
        Q[8][0] =   0.095012509837637e0;   
        Q[9][0] =   0.281603550779259e0;           
        Q[10][0] =  0.458016777657227e0;           
        Q[11][0] =  0.617876244402644e0;   
        Q[12][0] =  0.755404408355003e0;   
        Q[13][0] =  0.865631202387832e0;   
        Q[14][0] =  0.944575023073233e0;   
        Q[15][0] =  0.989400934991650e0;   
            
	return Q;
}

const vector<double>
Rule::W1() const
{
	vector<double> W(1);

	W[0]    =   2.0e0;

	return W;
}

const vector<double>
Rule::W3() const
{
	vector<double> W(2);

	W[0]    =   1.0e0;
	W[1]    =   1.0e0;

	return W;
}

const vector<double>
Rule::W5() const
{
	vector<double> W(3);

	W[0]    =   0.555555555555556e0;
	W[1]    =   0.888888888888889e0;
	W[2]    =   0.555555555555556e0;

	return W;
}

const vector<double>
Rule::W7() const
{
	vector<double> W(4);

	W[0]    =   0.347854845137454e0;
	W[1]    =   0.652145154862546e0;
	W[2]    =   0.652145154862546e0;
	W[3]    =   0.347854845137454e0;

	return W;
}

const vector<double>
Rule::W9() const
{
	vector<double> W(5);

	W[0]    =   0.236926885056189e0;
	W[1]    =   0.478628670499366e0;
	W[2]    =   0.568888888888889e0;
	W[3]    =   0.478628670499366e0;
	W[4]    =   0.236926885056189e0;

	return W;
}

const vector<double>
Rule::W11() const
{
	vector<double> W(6);

	W[0]    =   0.171324492379170e0;
	W[1]    =   0.360761573048139e0;
	W[2]    =   0.467913934572691e0;
	W[3]    =   0.467913934572691e0;
	W[4]    =   0.360761573048139e0;
	W[5]    =   0.171324492379170e0;

	return W;
}

const vector<double>
Rule::W13() const
{
	vector<double> W(7);

        W[0]    =   0.129484966168870e0;
        W[1]    =   0.279705391489277e0;
        W[2]    =   0.381830050505119e0;
        W[3]    =   0.417959183673469e0;
        W[4]    =   0.381830050505119e0;
        W[5]    =   0.279705391489277e0;
        W[6]    =   0.129484966168870e0;

	return W;
}

const vector<double>
Rule::W15() const
{
	vector<double> W(8);

        W[0]    =   0.101228536290376e0;
        W[1]    =   0.222381034453374e0;
        W[2]    =   0.313706645877887e0;
        W[3]    =   0.362683783378362e0;
        W[4]    =   0.362683783378362e0;
        W[5]    =   0.313706645877887e0;
        W[6]    =   0.222381034453374e0;
        W[7]    =   0.101228536290376e0;

	return W;
}

const vector<double>
Rule::W17() const
{
	vector<double> W(9);

        W[0]    =   0.081274388361574e0;
        W[1]    =   0.180648160694857e0;
        W[2]    =   0.260610696402935e0;
        W[3]    =   0.312347077040003e0;
        W[4]    =   0.330239355001260e0;
        W[5]    =   0.312347077040003e0;
        W[6]    =   0.260610696402935e0;
        W[7]    =   0.180648160694857e0;
        W[8]    =   0.081274388361574e0;

	return W;
}

const vector<double>
Rule::W19() const
{
	vector<double> W(10);

        W[0]    =   0.066671344308688e0;
        W[1]    =   0.149451349150581e0;
        W[2]    =   0.219086362515982e0;
        W[3]    =   0.269266719309996e0;
        W[4]    =   0.295524224714753e0;
        W[5]    =   0.295524224714753e0;
        W[6]    =   0.269266719309996e0;
        W[7]    =   0.219086362515982e0;
        W[8]    =   0.149451349150581e0;
        W[9]    =   0.066671344308688e0;

	return W;
}

const vector<double>
Rule::W21() const
{
	vector<double> W(11);

        W[0]    =   0.055668567116174e0;
        W[1]    =   0.125580369464905e0;
        W[2]    =   0.186290210927734e0;
        W[3]    =   0.233193764591990e0;
        W[4]    =   0.262804544510247e0;
        W[5]    =   0.272925086777901e0;
        W[6]    =   0.262804544510247e0;
        W[7]    =   0.233193764591990e0;
        W[8]    =   0.186290210927734e0;
        W[9]    =   0.125580369464905e0;
        W[10]   =   0.055668567116174e0;

	return W;
}

const vector<double>
Rule::W23() const
{
	vector<double> W(12);

        W[0]    =   0.047175336386512e0;
        W[1]    =   0.106939325995318e0;
        W[2]    =   0.160078328543346e0;
        W[3]    =   0.203167426723066e0;
        W[4]    =   0.233492536538355e0;
        W[5]    =   0.249147045813403e0;
        W[6]    =   0.249147045813403e0;
        W[7]    =   0.233492536538355e0;
        W[8]    =   0.203167426723066e0;
        W[9]    =   0.160078328543346e0;
        W[10]   =   0.106939325995318e0;
        W[11]   =   0.047175336386512e0;

	return W;
}

const vector<double>
Rule::W25() const
{
	vector<double> W(13);

        W[0]    =   0.040484004765316e0;
        W[1]    =   0.092121499837728e0;
        W[2]    =   0.138873510219787e0;
        W[3]    =   0.178145980761946e0;
        W[4]    =   0.207816047536889e0;
        W[5]    =   0.226283180262897e0;
        W[6]    =   0.232551553230874e0;
        W[7]    =   0.226283180262897e0;
        W[8]    =   0.207816047536889e0;
        W[9]    =   0.178145980761946e0;
        W[10]   =   0.138873510219787e0;
        W[11]   =   0.092121499837728e0;
        W[12]   =   0.040484004765316e0;

	return W;
}

const vector<double>
Rule::W27() const
{
	vector<double> W(14);

        W[0]    =   0.035119460331752e0;
        W[1]    =   0.080158087159760e0;
        W[2]    =   0.121518570687903e0;
        W[3]    =   0.157203167158194e0;
        W[4]    =   0.185538397477938e0;
        W[5]    =   0.205198463721296e0;
        W[6]    =   0.215263853463158e0;
        W[7]    =   0.215263853463158e0;
        W[8]    =   0.205198463721296e0;
        W[9]    =   0.185538397477938e0;
        W[10]   =   0.157203167158194e0;
        W[11]   =   0.121518570687903e0;
        W[12]   =   0.080158087159760e0;
        W[13]   =   0.035119460331752e0;

	return W;
}

const vector<double>
Rule::W29() const
{
	vector<double> W(15);

        W[0]    =   0.030753241996117e0;
        W[1]    =   0.070366047488108e0;
        W[2]    =   0.107159220467172e0;
        W[3]    =   0.139570677926154e0;
        W[4]    =   0.166269205816994e0;
        W[5]    =   0.186161000015562e0;
        W[6]    =   0.198431485327111e0;
        W[7]    =   0.202578241925561e0;
        W[8]    =   0.198431485327111e0;
        W[9]    =   0.186161000015562e0;
        W[10]   =   0.166269205816994e0;
        W[11]   =   0.139570677926154e0;
        W[12]   =   0.107159220467172e0;
        W[13]   =   0.070366047488108e0;
        W[14]   =   0.030753241996117e0;

	return W;
}

const vector<double>
Rule::W31() const
{
	vector<double> W(16);

        W[0]    =   0.027152459411754e0;
        W[1]    =   0.062253523938648e0;
        W[2]    =   0.095158511682493e0;
        W[3]    =   0.124628971255534e0;
        W[4]    =   0.149595988816577e0;
        W[5]    =   0.169156519395003e0;
        W[6]    =   0.182603415044924e0;
        W[7]    =   0.189450610455069e0;
        W[8]    =   0.189450610455069e0;
        W[9]    =   0.182603415044924e0;
        W[10]   =   0.169156519395003e0;
        W[11]   =   0.149595988816577e0;
        W[12]   =   0.124628971255534e0;
        W[13]   =   0.095158511682493e0;
        W[14]   =   0.062253523938648e0;
        W[15]   =   0.027152459411754e0;

	return W;
}

}

}

}
