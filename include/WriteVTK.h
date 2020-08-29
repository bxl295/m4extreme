//
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#ifndef _M4EXTREME_WRITEVTKFILE_H_
#define _M4EXTREME_WRITEVTKFILE_H_

#include "Factory/ModelBuilder.h"

namespace m4extreme {

  enum ElementType
  {
    // Linear cell types defined in VTK
    VTK_VERTEX = 1,
    VTK_POLY_VERTEX = 2,
    VTK_LINE = 3,
    VTK_POLY_LINE = 4,
    VTK_TRIANGLE = 5,
    VTK_TRIANGLE_STRIP = 6,
    VTK_POLYGON = 7,
    VTK_PIXEL = 8,
    VTK_QUAD = 9,
    VTK_TETRA = 10,
    VTK_VOXEL = 11,
    VTK_HEXAHEDRON = 12,
    VTK_WEDGE = 13,
    VTK_PYRAMID = 14,

    // Non-linear cell types defined in VTK
    VTK_QUADRATIC_EDGE = 21,
    VTK_QUADRATIC_TRIANGLE = 22,
    VTK_QUADRATIC_QUAD = 23,
    VTK_QUADRATIC_TETRA = 24,
    VTK_QUADRATIC_HEXAHEDRON = 25,
  };

  bool isabadnumber(double value) {
    if ( value != value || 
	 value > std::numeric_limits<double>::max() ||
	 value < -std::numeric_limits<double>::max() ) {
      return true;
    }
    else {
      return false;
    }
  }

  bool isabadvector(const Set::VectorSpace::Vector & value) {
    bool result = false;
    for ( int i = 0; i < value.size(); ++i ) {
      if ( value[i] != value[i] || 
	   value[i] > std::numeric_limits<double>::max() ||
	   value[i] < -std::numeric_limits<double>::max() ) {
	result = true;
	break;
      }      
    }

    return result;
  }

  //-----------------------------------------------------------------------
  //
  // Output Header information
  //
  //------------------------------------------------------------------------
  void _write_header(ofstream & ofs) {
    ofs << "# vtk DataFile Version 2.0" << std::endl
	<< "VTK file" << std::endl
	<< "ASCII" << std::endl
	<< "DATASET UNSTRUCTURED_GRID" << std::endl;
    return;
  }

  //-----------------------------------------------------------------------
  //
  // Output geometric information
  //
  //------------------------------------------------------------------------  
  void _write_geometry(ofstream & ofs, int elementType,
		       const map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point> & P,
		       Geometry::CellComplex* solid) {
    int dim = P.begin()->second.size();

    // Point information
    int numofNodes = P.size();
    ofs << "POINTS " << numofNodes << " double" << endl;

    map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point>::const_iterator pP;
    map<Geometry::Cell *, int> Nodesidmap;
    int i = 0;
    for( pP = P.begin(); pP != P.end(); pP++ ) {
      const Set::Euclidean::Orthonormal::Point & pLocal = pP->second;
      ofs << pLocal[0] << " " << pLocal[1];
      if (dim > 2) {
	ofs << " " << pLocal[2];
      }
      else {
	ofs << " " << 0.0;
      }
      ofs << std::endl;
      Nodesidmap.insert(make_pair(pP->first,i));
      i++;
    }

    ofs << std::endl;

 
    // Cell information
    int cell_dim = solid->size() - 1;
    const set<Geometry::Cell*> & emloc = (*solid)[cell_dim];
    int numofEle = emloc.size();
    int size = 0;
    int tag = 0;
    switch ( elementType ) {
    case VTK_TRIANGLE:
      {
	size = numofEle * 4;
	tag = 3;
      }
      break;
    case VTK_QUAD:
      {
	size = numofEle * 5;
	tag = 4;
      }
      break;
    case VTK_TETRA:
      {
	size = numofEle * 5;
	tag = 4;
      }
      break;
    }
  
    ofs << "CELLS " << numofEle << " " << size << std::endl;

    //Output Cells
    for(set<Geometry::Cell*>::iterator pEN = emloc.begin(); pEN != emloc.end(); pEN++) {
      set<Geometry::Cell *> v_set;
      m4extreme::getVertices(*pEN, v_set);
      ofs << tag;
      for(set<Geometry::Cell *>::iterator pV = v_set.begin(); pV != v_set.end(); pV++) {
	ofs << " " << Nodesidmap.find(*pV)->second;
      }
      ofs << std::endl;
    }
  
    ofs << std::endl;

    ofs << "CELL_TYPES " << numofEle << endl;
    for (int count = 0; count < numofEle; ++count) {
      ofs << elementType << std::endl;
    }

    ofs << std::endl;

    return;
  }

  //-----------------------------------------------------------------------
  //
  // Output point scalar information
  //
  //------------------------------------------------------------------------  
  void _write_data(ofstream & ofs,
		   const map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point> & P,
		   const map< std::string, map<Geometry::Cell*, double>* > & nodal_scalars) {
   
    if ( !nodal_scalars.empty() ) {
      map<std::string, map<Geometry::Cell*, double>* >::const_iterator it_field;
      map<Geometry::Cell*, double>::const_iterator it_data;

      for ( it_field = nodal_scalars.begin(); it_field != nodal_scalars.end(); ++it_field) {
	ofs << "SCALARS " << it_field->first.c_str() << " double 1" << std::endl;
	ofs << "LOOKUP_TABLE default" << std::endl;
	const map<Geometry::Cell*, double>* pData = it_field->second;

	map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point>::const_iterator pP;
	for( pP = P.begin(); pP != P.end(); pP++ ) {
	  if ( (it_data=pData->find(pP->first)) != pData->end() ) {
	    double value = it_data->second;
	    if ( isabadnumber(value) ) {
	      value = 0.0;
	    }

	    ofs << value << std::endl;
	  }
	  else {
	    ofs << 0.0 << std::endl;
	  }
	}
      }
    }

    return;
  }

  //-----------------------------------------------------------------------
  //
  // Output point vector information
  //
  //------------------------------------------------------------------------  
  void _write_data(ofstream & ofs,
		   const map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point> & P,
		   const map< std::string, map<Geometry::Cell*, Set::VectorSpace::Vector>* > & nodal_vectors) {

    if ( !nodal_vectors.empty() ) {
      map<std::string, map<Geometry::Cell*, Set::VectorSpace::Vector>* >::const_iterator it_vectorfield;
      map<Geometry::Cell*, Set::VectorSpace::Vector>::const_iterator it_vector;
      
      for ( it_vectorfield = nodal_vectors.begin(); it_vectorfield != nodal_vectors.end(); ++it_vectorfield) {
	ofs << "VECTORS " << it_vectorfield->first.c_str() << " double" << std::endl;
	const map<Geometry::Cell*, Set::VectorSpace::Vector>* pVector = it_vectorfield->second;
	
	int dim = P.begin()->second.size();
	map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point>::const_iterator pP;
	for( pP = P.begin(); pP != P.end(); pP++ ) {
	  if ( (it_vector=pVector->find(pP->first)) != pVector->end() ) {
	    const Set::VectorSpace::Vector & vloc = it_vector->second;
	    if ( isabadvector(vloc) ) {
	      ofs << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
	    }
	    else {
	      if ( dim == 3 ) {
		ofs << vloc[0] << " " << vloc[1] << " " << vloc[2] << std::endl;	      
	      }
	      else {
		ofs << vloc[0] << " " << vloc[1] << " " << 0.0 << std::endl;
	      }
	    }
	  }
	  else {
	    ofs << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
	  }
	}
      }
    }    
    
    return;
  }

  //-----------------------------------------------------------------------
  //
  // Output cell scalar information
  //
  //------------------------------------------------------------------------  
  void _write_data(ofstream & ofs,
		   Geometry::CellComplex* solid,
		   const map< std::string, map<Geometry::Cell*, double>* > & mpt_scalars) {
    
    if ( !mpt_scalars.empty() ) {
      map<std::string, map<Geometry::Cell*, double>* >::const_iterator it_field;
      map<Geometry::Cell*, double>::const_iterator it_data; 

      int dim = solid->size() - 1;
      const set<Geometry::Cell*> & emloc = (*solid)[dim];

      for ( it_field = mpt_scalars.begin(); it_field != mpt_scalars.end(); ++it_field) {
	ofs << "SCALARS " << it_field->first.c_str() << " double 1" << std::endl;
	ofs << "LOOKUP_TABLE default" << std::endl;
	const map<Geometry::Cell*, double>* pData = it_field->second;
	
	for(set<Geometry::Cell*>::iterator pEN = emloc.begin(); pEN != emloc.end(); pEN++) {
	  if ( (it_data=pData->find(*pEN)) != pData->end() ) {
	    double value = it_data->second;
	    if ( isabadnumber(value) ) {
	      ofs << 0.0 << std::endl;
	    }
	    else {
	      ofs << value << std::endl;
	    }
	  }
	  else {
	    ofs << 0.0 << std::endl;
	  }
	}
      }
    }

    return;
  }

  //-----------------------------------------------------------------------
  //
  // Output cell vector information
  //
  //------------------------------------------------------------------------  
  void _write_data(ofstream & ofs,
		   Geometry::CellComplex* solid,
		   const map< std::string, map<Geometry::Cell*, Set::VectorSpace::Vector>* > & mpt_vectors) {

    if ( !mpt_vectors.empty() ) {
      map<std::string, map<Geometry::Cell*, Set::VectorSpace::Vector>* >::const_iterator it_vectorfield;
      map<Geometry::Cell*, Set::VectorSpace::Vector>::const_iterator it_vector;

      int dim = solid->size() - 1;
      const set<Geometry::Cell*> & emloc = (*solid)[dim];
      
      for ( it_vectorfield = mpt_vectors.begin(); it_vectorfield != mpt_vectors.end(); ++it_vectorfield) {
	ofs << "VECTORS " << it_vectorfield->first.c_str() << " double" << std::endl;
	const map<Geometry::Cell*, Set::VectorSpace::Vector>* pVector = it_vectorfield->second;
	
	for(set<Geometry::Cell*>::iterator pEN = emloc.begin(); pEN != emloc.end(); pEN++) {
	  if ( (it_vector=pVector->find(*pEN)) != pVector->end() ) {
	    const Set::VectorSpace::Vector & vloc = it_vector->second;
	    
	    if ( isabadvector(vloc) ) {
	      ofs << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
	    }
	    else {
	      if ( dim == 3 ) {
		ofs << vloc[0] << " " << vloc[1] << " " << vloc[2] << std::endl;
	      }
	      else {
		ofs << vloc[0] << " " << vloc[1] << " " << 0.0 << std::endl;
	      }
	    }
	  }
	  else {
	    ofs << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
	  }
	}
      }
    }    
    
    return;
  }

  //-----------------------------------------------------------------------
  //
  // get data field from pModel
  //
  //------------------------------------------------------------------------  
  void _get_data(int index,
		 m4extreme::MEMPModelBuilder * pModel,
		 map<string, map<Geometry::Cell*, double>* > & nodal_field, 
		 map<string, map<Geometry::Cell*, double>* > & mpt_field,
		 map<string, map<Geometry::Cell*, Set::VectorSpace::Vector>* > & nodal_vectorfield) {

    int dim = pModel->_DIM;
    assert(dim == 2 || dim == 3);

    /* nodal_field.insert( make_pair("pressure",  new map<Geometry::Cell*, double>) ); */
    /* nodal_field.insert( make_pair("mises",     new map<Geometry::Cell*, double>) ); */
    nodal_field.insert( make_pair("status",    new map<Geometry::Cell*, double>) );
    /* nodal_field.insert( make_pair("s11",       new map<Geometry::Cell*, double>) ); */
    /* nodal_field.insert( make_pair("s22",       new map<Geometry::Cell*, double>) ); */
    /* nodal_field.insert( make_pair("s12",       new map<Geometry::Cell*, double>) ); */
    nodal_vectorfield.insert( make_pair("velocity", new map<Geometry::Cell*, Set::VectorSpace::Vector>) );

    /* if ( dim > 2 ) { */
    /*   nodal_field.insert( make_pair("s13", new map<Geometry::Cell*, double>) ); */
    /*   nodal_field.insert( make_pair("s23", new map<Geometry::Cell*, double>) ); */
    /*   nodal_field.insert( make_pair("s33", new map<Geometry::Cell*, double>) );    */
    /* } */

    pModel->GetNodeData(*nodal_vectorfield["velocity"], m4extreme::MEMPModelBuilder::VELOCITY);
    pModel->GetNodeData(index, *nodal_field["status"], m4extreme::MEMPModelBuilder::STATUS);
      
    /* map<Geometry::Cell*, Set::VectorSpace::Hom> cauchy; */
    /* pModel->GetNodeData(index, cauchy, m4extreme::MEMPModelBuilder::CAUCHY_STRESS); */

    /* map<Geometry::Cell*, Set::VectorSpace::Hom>::iterator pC; */
    /* for (pC = cauchy.begin(); pC != cauchy.end(); pC++) { */
    /*   Set::VectorSpace::Hom & sigma = pC->second; */
    /*   nodal_field["s11"]->insert( make_pair(pC->first, sigma[0][0]) ); */
    /*   nodal_field["s22"]->insert( make_pair(pC->first, sigma[1][1]) ); */
    /*   nodal_field["s12"]->insert( make_pair(pC->first, sigma[0][1]) ); */

    /*   if ( dim > 2 ) { */
    /* 	nodal_field["s33"]->insert( make_pair(pC->first, sigma[2][2]) ); */
    /* 	nodal_field["s13"]->insert( make_pair(pC->first, sigma[0][2]) ); */
    /* 	nodal_field["s23"]->insert( make_pair(pC->first, sigma[1][2]) ); */
	  
    /* 	double mean = ( sigma[0][0] + sigma[1][1] + sigma[2][2] ) / 3.0; */
    /* 	nodal_field["pressure"]->insert( make_pair(pC->first, mean) ); */
	  
    /* 	double sigv = 0.0; */

    /* 	double tempsig = sigma[0][0] - sigma[1][1]; */
    /* 	tempsig *= tempsig; */
    /* 	sigv += tempsig; */

    /* 	tempsig = sigma[1][1] - sigma[2][2]; */
    /* 	tempsig *= tempsig; */
    /* 	sigv += tempsig; */

    /* 	tempsig = sigma[2][2] - sigma[0][0]; */
    /* 	tempsig *= tempsig; */
    /* 	sigv += tempsig; */

    /* 	tempsig  = sigma[1][2] * sigma[1][2]; */
    /* 	tempsig += sigma[0][2] * sigma[0][2]; */
    /* 	tempsig += sigma[0][1] * sigma[0][1]; */
    /* 	tempsig *= 6.0; */
    /* 	sigv += tempsig; */
    /* 	sigv *= 0.5; */

    /* 	nodal_field["mises"]->insert( make_pair(pC->first, sqrt(sigv)) ); */
    /*   } */
    /*   else { */
    /* 	double mean = ( sigma[0][0] + sigma[1][1] ) / 2.0; */
    /* 	nodal_field["pressure"]->insert( make_pair(pC->first, mean) ); */

    /* 	double sum = sigma[0][0]*sigma[0][0]; */
    /* 	sum -= sigma[0][0] * sigma[1][1]; */
    /* 	sum += 3.0 * sigma[0][1]*sigma[0][1]; */
    /* 	sum += sigma[1][1]*sigma[1][1]; */
    /* 	nodal_field["mises"]->insert( make_pair(pC->first, sqrt(sum)) ); */
    /*   } */
    /* } */

    // mpt field
    map<Geometry::Cell*, Set::VectorSpace::Hom> mpt_cauchy;
    map<Geometry::Cell*, Set::VectorSpace::Hom>::iterator mpt_pC;

    mpt_field.insert( make_pair("temperature",  new map<Geometry::Cell*, double>) );
    mpt_field.insert( make_pair("pressure",  new map<Geometry::Cell*, double>) );
    mpt_field.insert( make_pair("mises",     new map<Geometry::Cell*, double>) );
    mpt_field.insert( make_pair("weight",    new map<Geometry::Cell*, double>) );
    mpt_field.insert( make_pair("sigma_max", new map<Geometry::Cell*, double>) );
    mpt_field.insert( make_pair("sigma_min", new map<Geometry::Cell*, double>) );

    if ( dim > 2 ) {
      mpt_field.insert( make_pair("sigma_median", new map<Geometry::Cell*, double>) );
    }

    /* mpt_field.insert( make_pair("s11", new map<Geometry::Cell*, double>) ); */
    /* mpt_field.insert( make_pair("s22", new map<Geometry::Cell*, double>) ); */
    /* mpt_field.insert( make_pair("s12", new map<Geometry::Cell*, double>) ); */

    /* if ( dim > 2 ) { */
    /*   mpt_field.insert( make_pair("s33", new map<Geometry::Cell*, double>) );   */
    /*   mpt_field.insert( make_pair("s13", new map<Geometry::Cell*, double>) ); */
    /*   mpt_field.insert( make_pair("s23", new map<Geometry::Cell*, double>) ); */
    /* } */

    pModel->GetQPData(index, *mpt_field["temperature"], m4extreme::MEMPModelBuilder::TEMPERATURE);  	
    pModel->GetQPData(index, *mpt_field["weight"], m4extreme::MEMPModelBuilder::WEIGHT);
 
#if defined(_M4EXTREME_EIGEN_FRACTURE_)
    mpt_field.insert( make_pair("status", new map<Geometry::Cell*, double>) );
    pModel->GetQPData(index, *mpt_field["status"], m4extreme::MEMPModelBuilder::STATUS);
#endif
    
    pModel->GetQPData(index, mpt_cauchy, m4extreme::MEMPModelBuilder::CAUCHY_STRESS);      
    for (mpt_pC = mpt_cauchy.begin(); mpt_pC != mpt_cauchy.end(); mpt_pC++) {
      Set::VectorSpace::Hom & sigma = mpt_pC->second;
      /* mpt_field["s11"]->insert( make_pair(mpt_pC->first, sigma[0][0]) ); */
      /* mpt_field["s22"]->insert( make_pair(mpt_pC->first, sigma[1][1]) ); */
      /* mpt_field["s12"]->insert( make_pair(mpt_pC->first, sigma[0][1]) ); */

      if ( dim > 2 ) {
	/* mpt_field["s33"]->insert( make_pair(mpt_pC->first, sigma[2][2]) ); */
	/* mpt_field["s13"]->insert( make_pair(mpt_pC->first, sigma[0][2]) ); */
	/* mpt_field["s23"]->insert( make_pair(mpt_pC->first, sigma[1][2]) ); */
	double mpt_mean = ( sigma[0][0] + sigma[1][1] + sigma[2][2] ) / 3.0;
	mpt_field["pressure"]->insert( make_pair(mpt_pC->first, mpt_mean) );
	  
	/* double sigv = 0.0; */

	/* double tempsig = sigma[0][0] - sigma[1][1]; */
	/* tempsig *= tempsig; */
	/* sigv += tempsig; */

	/* tempsig = sigma[1][1] - sigma[2][2]; */
	/* tempsig *= tempsig; */
	/* sigv += tempsig; */

	/* tempsig = sigma[2][2] - sigma[0][0]; */
	/* tempsig *= tempsig; */
	/* sigv += tempsig; */

	/* tempsig  = sigma[1][2] * sigma[1][2]; */
	/* tempsig += sigma[0][2] * sigma[0][2]; */
	/* tempsig += sigma[0][1] * sigma[0][1]; */
	/* tempsig *= 6.0; */
	/* sigv += tempsig; */
	/* sigv *= 0.5; */

	/* mpt_field["mises"]->insert( make_pair(mpt_pC->first, sqrt(sigv)) ); */
	double s1, s2, s3;
	double Tol = 1.0e-3;
	if (fabs(sigma[0][1]) > Tol || fabs(sigma[0][2]) > Tol || fabs(sigma[1][2]) > Tol ) {
	  double I1 = 3.0 * mpt_mean;
	  double I2 = sigma[0][0]*sigma[1][1] + sigma[2][2]*sigma[1][1]
	    + sigma[0][0]*sigma[2][2] - sigma[0][1]*sigma[0][1]
	    - sigma[1][2]*sigma[1][2] - sigma[0][2]*sigma[0][2];
	  double I3 = sigma[0][0]*sigma[1][1]*sigma[2][2]
	    - sigma[0][0]*sigma[1][2]*sigma[1][2]
	    - sigma[1][1]*sigma[0][2]*sigma[0][2]
	    - sigma[2][2]*sigma[0][1]*sigma[0][1]
	    + 2.0*sigma[0][1]*sigma[1][2]*sigma[0][2];

	  double phi;
	  double R = sqrt( fabs(I1 * I1 - 3.0 * I2) );

	  if ( R < 1.0e-12 ) {
	    s1 = 0.0;
	    s2 = 0.0;
	    s3 = 0.0;
	  }
	  else {
	    double R3 = R;
	    R3 *= R;
	    R3 *= R;
	    double q = I1*I1*I1 - 4.5*I1*I2 + 13.5*I3;
	    q /= R3;

	    double PI = 3.1415926;
	    if ( q > 1.0 ) {
	      phi = 0.0;
	    }
	    else if ( q < -1.0 ) {
	      phi = PI;
	    }
	    else {
	      phi = acos(q);
	      phi /= 3.0;
	    }	  

	    double I1_third = I1/3.0;
	    double vloc = 2.0 * R / 3.0;
	    s1 = I1_third + vloc * cos(phi);
	    s2 = I1_third + vloc * cos(phi - 2.0 * PI / 3.0);
	    s3 = I1_third + vloc * cos(phi - 4.0 * PI / 3.0);
	  }
	}
	else {
	  s1 = sigma[0][0];
	  s2 = sigma[1][1];
	  s3 = sigma[2][2];
	}

	double smax = 0.0;
	double smin = 0.0;
	double smedian = s3;

	if ( s1 > s2 ) {
	  smax = s1;
	  smin = s2;
	}
	else {
	  smax = s2;
	  smin = s1;
	}
	
	if ( smax < s3 ) {
	  smedian = smax;
	  smax = s3;	  
	}
	
	if ( smin > s3 ) {
	  smedian = smin;
	  smin = s3;
	}
	
	mpt_field["sigma_max"]->insert( make_pair(mpt_pC->first, smax) ); 	
	mpt_field["sigma_median"]->insert( make_pair(mpt_pC->first, smedian) ); 	
	mpt_field["sigma_min"]->insert( make_pair(mpt_pC->first, smin) ); 

	double sigv = 0.0;
	double tmp = s1 - s2;
	tmp *= tmp;
	sigv += tmp;

	tmp = s1 - s3;
	tmp *= tmp;
	sigv += tmp;

	tmp = s2 - s3;
	tmp *= tmp;
	sigv += tmp;
	sigv *= 0.5;
	mpt_field["mises"]->insert( make_pair(mpt_pC->first, sqrt(sigv)) );
      }
      else {
	double mpt_mean = ( sigma[0][0] + sigma[1][1] ) / 2.0;
	mpt_field["pressure"]->insert( make_pair(mpt_pC->first, mpt_mean) );

	double s2 = sigma[0][0] - sigma[1][1];
	s2 *= 0.5;
	s2 *= s2;
	s2 += sigma[0][1] * sigma[0][1];
	s2 = sqrt(s2);

	double smax = mpt_mean + s2;
	double smin = mpt_mean - s2;
	mpt_field["sigma_max"]->insert( make_pair(mpt_pC->first, smax) ); 	
	mpt_field["sigma_min"]->insert( make_pair(mpt_pC->first, smin) ); 


	mpt_field["mises"]->insert( make_pair(mpt_pC->first, sqrt(smax*smax - smax*smin + smin*smin)) );

	/* double sum = sigma[0][0]*sigma[0][0]; */
	/* sum -= sigma[0][0] * sigma[1][1]; */
	/* sum += 3.0 * sigma[0][1]*sigma[0][1]; */
	/* sum += sigma[1][1]*sigma[1][1]; */
	/* mpt_field["mises"]->insert( make_pair(mpt_pC->first, sqrt(sum)) ); */
      }
    }

    return;
  }

  /**************************************************************************************/
  // No data field
  /*************************************************************************************/

  void WriteVTKFile(const map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point> & P,
		    Geometry::CellComplex* solid,
		    const char * name,
		    int elementType) {
    std::ofstream ofs(name);
    if (ofs.fail()) {
      std::cerr << "could not open the file" << name << ". Abort writing!" << std::endl;
      assert(false);
    }

    _write_header(ofs);
    _write_geometry(ofs, elementType, P, solid);

    ofs.close();
    return;
  }

  /**************************************************************************************/
  // with data field input 
  /*************************************************************************************/

  void WriteVTKFile(const map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point> & P,
		    Geometry::CellComplex* solid,
		    const char * name,
		    int elementType,
		    const map< std::string, map<Geometry::Cell*, double>* > & scalars,
		    bool isNodalData) {
    std::ofstream ofs(name);
    if (ofs.fail()) {
      std::cerr << "could not open the file" << name << ". Abort writing!" << std::endl;
      assert(false);
    }

    _write_header(ofs);
    _write_geometry(ofs, elementType, P, solid);

    if ( isNodalData ) {
      int numofNodes = P.size();
      ofs << "POINT_DATA " << numofNodes << std::endl;
      _write_data(ofs, P, scalars);
    }
    else {
      int dim = solid->size() - 1;
      int numofEle = (*solid)[dim].size();
      ofs << "CELL_DATA " << numofEle << std::endl;
      _write_data(ofs, solid, scalars);
    }

    ofs.close();
    return;
  }

  /**************************************************************************************/
  // with data field input 
  /*************************************************************************************/

  void WriteVTKFile(const map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point> & P,
		    Geometry::CellComplex* solid,
		    const char * name,
		    int elementType,
		    const map< std::string, map<Geometry::Cell*, double>* > & nodal_scalars,
		    const map< std::string, map<Geometry::Cell*, double>* > & mpt_scalars) {
    std::ofstream ofs(name);
    if (ofs.fail()) {
      std::cerr << "could not open the file" << name << ". Abort writing!" << std::endl;
      assert(false);
    }

    _write_header(ofs);
    _write_geometry(ofs, elementType, P, solid);

    int numofNodes = P.size();
    ofs << "POINT_DATA " << numofNodes << std::endl;
    _write_data(ofs, P, nodal_scalars);

    int dim = solid->size() - 1;
    int numofEle = (*solid)[dim].size();
    ofs << "CELL_DATA " << numofEle << std::endl;
    _write_data(ofs, solid, mpt_scalars);

    ofs.close();
    return;
  }

  /**************************************************************************************/
  // with data field input 
  /*************************************************************************************/

  void WriteVTKFile(const map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point> & P,
		    Geometry::CellComplex* solid,
		    const char * name,
		    int elementType,
		    const map< std::string, map<Geometry::Cell*, double>* > & scalars,
		    const map< std::string, map<Geometry::Cell*, Set::VectorSpace::Vector>* > & vectors,
		    bool isNodalData) {
    std::ofstream ofs(name);
    if (ofs.fail()) {
      std::cerr << "could not open the file" << name << ". Abort writing!" << std::endl;
      assert(false);
    }

    _write_header(ofs);
    _write_geometry(ofs, elementType, P, solid);

    if ( isNodalData ) {
      int numofNodes = P.size();
      ofs << "POINT_DATA " << numofNodes << std::endl;
      _write_data(ofs, P, scalars);
      _write_data(ofs, P, vectors);
    }
    else {
      int dim = solid->size() - 1;
      int numofEle = (*solid)[dim].size();
      ofs << "CELL_DATA " << numofEle << std::endl;
      _write_data(ofs, solid, scalars);
      _write_data(ofs, solid, vectors);
    }

    ofs.close();
    return;
  }

  /**************************************************************************************/
  // with data field input 
  /*************************************************************************************/

  void WriteVTKFile(const map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point> & P,
		    Geometry::CellComplex* solid,
		    const char * name,
		    int elementType,
		    const map< std::string, map<Geometry::Cell*, double>* > & nodal_scalars,
		    const map< std::string, map<Geometry::Cell*, Set::VectorSpace::Vector>* > & nodal_vectors,
		    const map< std::string, map<Geometry::Cell*, double>* > & mpt_scalars,
		    const map< std::string, map<Geometry::Cell*, Set::VectorSpace::Vector>* > & mpt_vectors) {
    std::ofstream ofs(name);
    if (ofs.fail()) {
      std::cerr << "could not open the file" << name << ". Abort writing!" << std::endl;
      assert(false);
    }

    _write_header(ofs);
    _write_geometry(ofs, elementType, P, solid);

    int numofNodes = P.size();
    ofs << "POINT_DATA " << numofNodes << std::endl;
    _write_data(ofs, P, nodal_scalars);
    _write_data(ofs, P, nodal_vectors);

    int dim = solid->size() - 1;
    int numofEle = (*solid)[dim].size();
    ofs << "CELL_DATA " << numofEle << std::endl;
    _write_data(ofs, solid, mpt_scalars);
    _write_data(ofs, solid, mpt_vectors);

    ofs.close();
    return;
  }

  /**************************************************************************************/
  // data field obtained from pModel
  /*************************************************************************************/

  void WriteVTKFile(int index,
		    const map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point> & P,
		    Geometry::CellComplex* solid,
		    const char * name,
		    m4extreme::MEMPModelBuilder * pModel,
		    int elementType) {

    // check filename
    std::ofstream ofs(name);
    if (ofs.fail()) {
      std::cerr << "could not open the file" << name << ". Abort writing!" << std::endl;
      assert(false);
    } 

    // output information
    _write_header(ofs);
    _write_geometry(ofs, elementType, P, solid);


    // get data fields
    map<string, map<Geometry::Cell*, double>* > nodal_field, mpt_field;
    map<string, map<Geometry::Cell*, Set::VectorSpace::Vector>* > nodal_vectorfield;
    _get_data(index, pModel, nodal_field, mpt_field, nodal_vectorfield);

    int numofNodes = P.size();
    ofs << "POINT_DATA " << numofNodes << std::endl;
    _write_data(ofs, P, nodal_field);
    _write_data(ofs, P, nodal_vectorfield);

    int dim = solid->size() - 1;
    int numofEle = (*solid)[dim].size();
    ofs << "CELL_DATA " << numofEle << std::endl;
    _write_data(ofs, solid, mpt_field);

    ofs.close();

    for (map<string, map<Geometry::Cell*, double>* >::iterator pF = nodal_field.begin();
	 pF != nodal_field.end(); pF++) {
      delete pF->second;
    }

    for (map<string, map<Geometry::Cell*,  Set::VectorSpace::Vector>* >::iterator pF = nodal_vectorfield.begin();
	 pF != nodal_vectorfield.end(); pF++) {
      delete pF->second;
    }

    for (map<string, map<Geometry::Cell*, double>* >::iterator pF = mpt_field.begin();
	 pF != mpt_field.end(); pF++) {
      delete pF->second;
    }

    return;
  }

  /**************************************************************************************/
  // data field obtained from pModel + extra inputs
  /*************************************************************************************/

  void WriteVTKFile(int index,
		    const map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point> & P,
		    Geometry::CellComplex* solid,
		    const char * name,
		    m4extreme::MEMPModelBuilder * pModel,
		    int elementType,
		    const map< std::string, map<Geometry::Cell*, double>* > & scalarData,
		    bool isNodalData) {

    // check filename
    std::ofstream ofs(name);
    if (ofs.fail()) {
      std::cerr << "could not open the file" << name << ". Abort writing!" << std::endl;
      assert(false);
    } 

    // output information
    _write_header(ofs);
    _write_geometry(ofs, elementType, P, solid);


    // get data fields
    map<string, map<Geometry::Cell*, double>* > nodal_field, mpt_field;
    map<string, map<Geometry::Cell*, Set::VectorSpace::Vector>* > nodal_vectorfield;
    _get_data(index, pModel, nodal_field, mpt_field, nodal_vectorfield);

    int numofNodes = P.size();
    ofs << "POINT_DATA " << numofNodes << std::endl;
    _write_data(ofs, P, nodal_field);
    if ( isNodalData ) {
      _write_data(ofs, P, scalarData);
    }
    _write_data(ofs, P, nodal_vectorfield);

    int dim = solid->size() - 1;
    int numofEle = (*solid)[dim].size();
    ofs << "CELL_DATA " << numofEle << std::endl;
    _write_data(ofs, solid, mpt_field);
    if ( !isNodalData ) {
      _write_data(ofs, solid, scalarData);
    }

    ofs.close();

    for (map<string, map<Geometry::Cell*, double>* >::iterator pF = nodal_field.begin();
	 pF != nodal_field.end(); pF++) {
      delete pF->second;
    }

    for (map<string, map<Geometry::Cell*,  Set::VectorSpace::Vector>* >::iterator pF = nodal_vectorfield.begin();
	 pF != nodal_vectorfield.end(); pF++) {
      delete pF->second;
    }

    for (map<string, map<Geometry::Cell*, double>* >::iterator pF = mpt_field.begin();
	 pF != mpt_field.end(); pF++) {
      delete pF->second;
    }

    return;
  }

  /**************************************************************************************/
  // data field obtained from pModel + extra inputs
  /*************************************************************************************/

  void WriteVTKFile(int index,
		    const map<Geometry::Cell *, Set::Euclidean::Orthonormal::Point> & P,
		    Geometry::CellComplex* solid,
		    const char * name,
		    m4extreme::MEMPModelBuilder * pModel,
		    int elementType,
		    const map< std::string, map<Geometry::Cell*, double>* > & nodeData,
		    const map< std::string, map<Geometry::Cell*, double>* > & mptData) {

    // check filename
    std::ofstream ofs(name);
    if (ofs.fail()) {
      std::cerr << "could not open the file" << name << ". Abort writing!" << std::endl;
      assert(false);
    } 

    // output information
    _write_header(ofs);
    _write_geometry(ofs, elementType, P, solid);


    // get data fields
    map<string, map<Geometry::Cell*, double>* > nodal_field, mpt_field;
    map<string, map<Geometry::Cell*, Set::VectorSpace::Vector>* > nodal_vectorfield;
    _get_data(index, pModel, nodal_field, mpt_field, nodal_vectorfield);

    int numofNodes = P.size();
    ofs << "POINT_DATA " << numofNodes << std::endl;
    _write_data(ofs, P, nodal_field);
    _write_data(ofs, P, nodeData);
    _write_data(ofs, P, nodal_vectorfield);

    int dim = solid->size() - 1;
    int numofEle = (*solid)[dim].size();
    ofs << "CELL_DATA " << numofEle << std::endl;
    _write_data(ofs, solid, mpt_field);
    _write_data(ofs, solid, mptData);

    ofs.close();

    for (map<string, map<Geometry::Cell*, double>* >::iterator pF = nodal_field.begin();
	 pF != nodal_field.end(); pF++) {
      delete pF->second;
    }

    for (map<string, map<Geometry::Cell*,  Set::VectorSpace::Vector>* >::iterator pF = nodal_vectorfield.begin();
	 pF != nodal_vectorfield.end(); pF++) {
      delete pF->second;
    }

    for (map<string, map<Geometry::Cell*, double>* >::iterator pF = mpt_field.begin();
	 pF != mpt_field.end(); pF++) {
      delete pF->second;
    }

    return;
  }

}

#endif
