// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <fstream>
#include "./GeneralGeometry.h"
#include "Geometry/Utils/Simplex/Simplex.h"

#if defined(_M4EXTREME_VTK_)
#include "Geometry/Meshing/vtkDelaunay2m4extreme/vtkDelaunay2m4extreme.h"
#endif

#if defined(_M4EXTREME_MPI_)
#include "mpi/MPI_Core.h"
#endif

namespace Geometry {
    namespace Solid {
        namespace Mesh {

            GeneralGeometry::GeneralGeometry()
#if defined(_M4EXTREME_VTK_)
	      :_DEGEN_TOL(1.0e-4)
#endif
	    { return;}

            GeneralGeometry::~GeneralGeometry() {
	      for ( int i = 0; i < _SPatches.size(); ++i ) {
		delete _SPatches[i];
	      }
            }

            /**
             * Read mesh from a Femap neutral file (sequential version) or m4extreme data file (mpi version)
             */
	    GeneralGeometry::GeneralGeometry(int dim, 
					     const char* filename, 
					     const FILE_TYPE & type,
					     int flag):
#if defined(_M4EXTREME_VTK_)
	      _DEGEN_TOL(1.0e-4),
#endif
	      _type(type), _flag(flag) {

                assert(dim > 1 && dim < 4);

                bool isSuccessful = true;
                if ( _type == M4EXTREME_NEUTRAL ||
		     _type == M4EXTREME_NEUTRAL_MINI ) {                    
		  isSuccessful = this->_readm4extremeNeutralFormat(dim, filename);
		}
		else if ( _type == M4EXTREME_NEUTRAL_FULL ){
		  isSuccessful = this->_readm4extremeNeutralFull(dim, filename);
		}
		else {
		  assert(false);
                }

                assert(isSuccessful);
            }

            const CellComplex &
            GeneralGeometry::GetCellComplex() const {
                return _SC;
            }

            CellComplex &
            GeneralGeometry::GetCellComplex() {
                return _SC;
            }

            const map<Cell *, Cell *> &
            GeneralGeometry::GetBinding() const {
                return _h;
            }

            map<Cell *, Cell *> &
            GeneralGeometry::GetBinding() {
                return _h;
            }

            const map<Cell *, Set::Euclidean::Orthonormal::Point> &
            GeneralGeometry::GetPointSet() const {
                return _XNodes;
            }

            map<Cell *, Set::Euclidean::Orthonormal::Point> &
            GeneralGeometry::GetPointSet() {
                return _XNodes;
            }

	    const map<int, Cell *> & 
	    GeneralGeometry::GetIdMap() const {		  
	       return _node_idmap;
	    }
	  
	    const map<Cell *, int> & 
	    GeneralGeometry::GetInverseIdMap() const {		  
	       return _node_inverse_idmap;
	    }
	  
	    const map<int, Cell *> & 
	    GeneralGeometry::GetElementIds() const {		  
	       return _elm_idmap;
	    }

	    const map<Cell *, int> & 
	    GeneralGeometry::GetInverseElementIds() const {		  
	       return _elm_inverse_idmap;
	    }
	  
  	    bool GeneralGeometry::FindElementGroup(const string & name) const {
	      map< string, set<Cell*> >::const_iterator pc = _solid_groups.find(name);
	      if (pc == _solid_groups.end()) {
		cerr << "ERROR: couldnot find element group " << name.c_str() << endl;
		return false;
	      }
	      else {
		return true;
	      }	      
	    }

            const set<Cell*> & GeneralGeometry::GetElementGroup(const string & name) const {
                map< string, set<Cell*> >::const_iterator pc = _solid_groups.find(name);
                if (pc == _solid_groups.end()) {
		    cerr << "ERROR: couldnot find element group " << name.c_str() << endl;
                    assert(false);
                }
                return pc->second;
            }

  	    bool GeneralGeometry::FindNodeGroup(const string & name) const {
	      map< string, set<Cell*> >::const_iterator pc = _node_groups.find(name);
	      if (pc == _node_groups.end()) {
		cerr << "ERROR: couldnot find node group " << name.c_str() << endl;
		return false;
	      }
	      else {
		return true;
	      }	      
	    }

            const set<Cell*> & GeneralGeometry::GetNodeGroup(const string & name) const {
                map< string, set<Cell*> >::const_iterator pc = _node_groups.find(name);
                if (pc == _node_groups.end()) {
		    cerr << "ERROR: couldnot find node group " << name.c_str() << endl;
                    assert(false);
                }
                return pc->second;
            }

  	    int GeneralGeometry::GetElementGroupSize() const {
	        return _solid_groups.size();
	    }

	    int GeneralGeometry::GetNodeGroupSize() const {
	        return _node_groups.size();
	    }

	    const map<string, set<Cell*> > & 
	    GeneralGeometry::GetAllElmentGroups() const {
	      return _solid_groups;
	    }

	    const map<string, set<Cell*> > & 
	    GeneralGeometry::GetAllNodeGroups() const {
	      return _node_groups;
	    }

	    bool GeneralGeometry::_error_info(const string & str, ifstream & fs) {
	      std::cerr << "invalid input format \'" << str.c_str() << "\' . Abort!" << std::endl;
	      fs.close();
	      return false;
	    }

            /**
             * Reading m4extreme Neutral File Format, a complete m4extreme Neutral file format is shown below,
	     * but it can be simplified, as long as the basic information is presented, i.e.,
	     * the coordinates of nodes and the connectivity table.
             * $
             * $ dimension element_type  number_of_nodes number_of_element
             * $ id x y z groupid //nodal id, nodal coordinates(for instance 3 doubles for 3-dimensional nodes), group id
	     * $                  // (id must be non-negative integer)
             * $ ...
             * $ id id0 id1 id2 id3 groupid //element id, connectivity table(for instance 4 node tetrahedron), group id
             * $ ...
             */
            bool GeneralGeometry::_readm4extremeNeutralFormat(int dim, const char * filename) {

                std::ifstream ifs(filename);

                if (ifs.is_open()) {
                    vector<point_type> v;
                    vector< vector<unsigned int> > conn;
                    vector<unsigned int> ids;		    

                    //printf("Succeed in opening file %s\n\n", filename);
                    std::string str;
                    std::getline(ifs, str);

                    //printf("core %d reading string: %s\n", _rank, str.c_str());

                    int dimension, type, numofNodes, numofMpts;

                    if (sscanf(str.c_str(), "%d %d %d %d", &dimension, &type, &numofNodes, &numofMpts) != 4) {
		      str += " at reading header";
		      return _error_info(str,ifs);
                    }

                    assert(dimension == dim);

                    //printf("core %d reading data: %d %d\n", _rank, numofNodes, numofMpts);
		    std::getline(ifs, str);

		    /**
		     * node information
		     */

                    int count = 0;
                    int idloc;
                    point_type ploc(dim);
                    map<int, int> bind;

                    if (dim == 3) {

		        switch(_type) {
			case M4EXTREME_NEUTRAL: //format: id x y z
			    while (!ifs.eof() && count < numofNodes) {
			      //printf("core %d reading string %d: %s\n", _rank, count, str.c_str());
			      if ( !str.empty() ) {
				if (sscanf(str.c_str(), "%d %lf %lf %lf", &idloc, &ploc[0], &ploc[1], &ploc[2]) != 4) {
				  str += " at reading nodal information";
				  return _error_info(str,ifs);
				}
				//printf("core %d reading line %d: %d %lf %lf %lf\n", _rank, count, idloc, ploc[0], ploc[1], ploc[2]);
			      
				ids.push_back(idloc);
				v.push_back(ploc);
				bind.insert(make_pair(idloc, count));
			      }

			      std::getline(ifs, str);
			      ++count;
			    }
			    break;
			case M4EXTREME_NEUTRAL_MINI: //format: x y z
			    while (!ifs.eof() && count < numofNodes) {
			      //printf("core %d reading string %d: %s\n", _rank, count, str.c_str());
			      if ( !str.empty() ) {
				if (sscanf(str.c_str(), "%lf %lf %lf", &ploc[0], &ploc[1], &ploc[2]) != 3) {
				  str += " at reading nodal information";
				  return _error_info(str,ifs);
				}
				//printf("core %d reading line %d: %d %lf %lf %lf\n", _rank, count, idloc, ploc[0], ploc[1], ploc[2]);
				
				idloc = count;
				
				ids.push_back(idloc);
				v.push_back(ploc);
				bind.insert(make_pair(idloc, count));
			      }

			      std::getline(ifs, str);
			      ++count;
			    }
			    break;
			default:
			  assert(false);
			}


			/**
			 * element information
			 */
                        const int mptType = type;
                        int idsloc[mptType];
                        std::vector<unsigned int> mpt(mptType);

                        count = 0;
                        while (!ifs.eof() && count < numofMpts) {

			  if ( !str.empty() ) {
			    if ( std::sscanf(str.c_str(), "%d %d %d %d",
					     idsloc, idsloc + 1, idsloc + 2, idsloc + 3) == mptType) {
			      mpt[0] = bind.find(idsloc[0])->second;
			      mpt[1] = bind.find(idsloc[1])->second;
			      mpt[2] = bind.find(idsloc[2])->second;
			      mpt[3] = bind.find(idsloc[3])->second;
			    }
			    else if ( std::sscanf(str.c_str(), "%d %d %d %d %d %d %d %d",
						  idsloc, idsloc + 1, idsloc + 2, idsloc + 3, 
						  idsloc + 4, idsloc + 5, idsloc + 6, idsloc + 7 ) == mptType) {
			      mpt[0] = bind.find(idsloc[0])->second;
			      mpt[1] = bind.find(idsloc[1])->second;
			      mpt[2] = bind.find(idsloc[2])->second;
			      mpt[3] = bind.find(idsloc[3])->second;
			      mpt[4] = bind.find(idsloc[4])->second;
			      mpt[5] = bind.find(idsloc[5])->second;
			      mpt[6] = bind.find(idsloc[6])->second;
			      mpt[7] = bind.find(idsloc[7])->second;
			    }
			    else {
			      str += " at reading element information";
			      return _error_info(str,ifs);
			    }			    


                            conn.push_back(mpt);
			  }

			  std::getline(ifs, str);
			  ++count;
                        }
                    } else if (dim == 2) {
		      
		      switch(_type) {
			case M4EXTREME_NEUTRAL: //format: id x y
			    while (!ifs.eof() && count < numofNodes) {
			      //printf("core %d reading string %d: %s\n", _rank, count, str.c_str());
			      if ( !str.empty() ) {
				if (sscanf(str.c_str(), "%d %lf %lf", &idloc, &ploc[0], &ploc[1]) != 3) {
				  str += " at reading nodal information";
				  return _error_info(str,ifs);
				}
				//printf("core %d reading line %d: %d %lf %lf \n", _rank, count, idloc, ploc[0], ploc[1]);
			      
				ids.push_back(idloc);
				v.push_back(ploc);
				bind.insert(make_pair(idloc, count));
			      }

			      std::getline(ifs, str);
			      ++count;
			    }
			    break;
			case M4EXTREME_NEUTRAL_MINI: //format: x y
			    while (!ifs.eof() && count < numofNodes) {
			      //printf("core %d reading string %d: %s\n", _rank, count, str.c_str());
			      if ( !str.empty() ) {
				if (sscanf(str.c_str(), "%lf %lf", &ploc[0], &ploc[1]) != 2) {
				  str += " at reading nodal information";
				  return _error_info(str,ifs);
				}
				//printf("core %d reading line %d: %d %lf %lf\n", _rank, count, idloc, ploc[0], ploc[1]);
				
				idloc = count;
				
				ids.push_back(idloc);
				v.push_back(ploc);
				bind.insert(make_pair(idloc, count));
			      }

			      std::getline(ifs, str);
			      ++count;
			    }
			    break;
			default:
			  assert(false);
			}


		      /**
		       * element information
		       */
		      const int mptType = type;
		      int idsloc[mptType];
		      std::vector<unsigned int> mpt(mptType);
		      
		      count = 0;
		      while (!ifs.eof() && count < numofMpts) {
			
			if ( !str.empty() ) {
			  if (std::sscanf(str.c_str(), "%d %d %d",
					  idsloc, idsloc + 1, idsloc + 2) == mptType) {
			    mpt[0] = bind.find(idsloc[0])->second;
			    mpt[1] = bind.find(idsloc[1])->second;
			    mpt[2] = bind.find(idsloc[2])->second;
			  }
			  else if (std::sscanf(str.c_str(), "%d %d %d",
					  idsloc, idsloc + 1, idsloc + 2, idsloc + 3) == mptType) {
			    mpt[0] = bind.find(idsloc[0])->second;
			    mpt[1] = bind.find(idsloc[1])->second;
			    mpt[2] = bind.find(idsloc[2])->second;
			    mpt[3] = bind.find(idsloc[3])->second;
			  } 
			  else {
			    str += " at reading element information";
			    return _error_info(str,ifs);
			  }			    
			  
			  conn.push_back(mpt);
			}
			
			std::getline(ifs, str);
			++count;
		      }
                        
                    } else {
		        cerr << "dimension " << dim << " is not supported!" << endl;
                        assert(false);
                    }

                    ifs.close();

                    if (!v.empty()) {
                        vector<Geometry::Cell *> V = _SC.ReadTable(conn);

                        for (unsigned int i = 0; i < v.size(); i++) {
                            _XNodes.insert(nodeset_type::value_type(V[i], v[i]));
                            _node_idmap.insert(make_pair(ids[i], V[i]));
                        }

                    } else {
                        for (int i = 0; i <= dim; ++i) {
                            _SC.push_back(set<Cell*>());
                        }
                    }

                    return true;
                } else {
                    printf("couldnot open file %s", filename);
                    return false;
                }
            }

            bool GeneralGeometry::_readm4extremeNeutralFull(int dim, const char * filename) {

                std::ifstream ifs(filename);

                if (ifs.is_open()) {
		    map<unsigned int, point_type> v;
                    map<unsigned int, vector<unsigned int> > conn;
                    vector<unsigned int> ids;		    

                    //printf("Succeed in opening file %s\n\n", filename);
                    std::string str;
                    std::getline(ifs, str);

                    //printf("core %d reading string: %s\n", _rank, str.c_str());

                    int dimension, type, numofNodes, numofMpts;

                    if (sscanf(str.c_str(), "%d %d %d %d", &dimension, &type, &numofNodes, &numofMpts) != 4) {
		      str += " at reading header";
		      return _error_info(str,ifs);
                    }

                    assert(dimension == dim);

                    //printf("core %d reading data: %d %d\n", _rank, numofNodes, numofMpts);
		    std::getline(ifs, str);

                    int count = 0;
                    int idloc;		    
		    point_type ploc(dim);

                    if (dim == 3) {

		        /**
			 * node information
			 */
		        while (!ifs.eof() && count < numofNodes) {		    
			  
			  if ( !str.empty() ) {
			    //printf("core %d reading string %d: %s\n", _rank, count, str.c_str());
			    if (sscanf(str.c_str(), "%d %lf %lf %lf", &idloc, &ploc[0], &ploc[1], &ploc[2]) != 4) {
			      str += " at reading nodal information";
			      return _error_info(str,ifs);						      
			    }
			    //printf("core %d reading line %d: %d %lf %lf %lf\n", _rank, count, idloc, ploc[0], ploc[1], ploc[2]);
			    			    
			    v.insert( make_pair(idloc, ploc) );			    
			  }
			  std::getline(ifs, str);
			  ++count;
			}

			/**
			 * element information
			 */
                        const int mptType = type;
                        int idsloc[mptType];
                        std::vector<unsigned int> mpt(mptType);
			int numInfo = mptType + 1; //add additional information id and groupid

                        count = 0;
                        while (!ifs.eof() && count < numofMpts) {
			  if ( !str.empty() ) {
			    // four nodes per element--tetrahedron
			    if (std::sscanf(str.c_str(), "%d %d %d %d %d",
					    &idloc, idsloc, idsloc + 1, idsloc + 2, idsloc + 3) == numInfo) {
			      mpt[0] = idsloc[0];
			      mpt[1] = idsloc[1];
			      mpt[2] = idsloc[2];
			      mpt[3] = idsloc[3];
			    }
			    // eight nodes per element--hexagon
			    else if (std::sscanf(str.c_str(), "%d %d %d %d %d %d %d %d %d",
						 &idloc, idsloc, idsloc + 1, idsloc + 2, idsloc + 3,
						 idsloc + 4, idsloc + 5, idsloc + 6, idsloc + 7) == numInfo) {
			      mpt[0] = idsloc[0];
			      mpt[1] = idsloc[1];
			      mpt[2] = idsloc[2];
			      mpt[3] = idsloc[3];
			      mpt[4] = idsloc[4];
			      mpt[5] = idsloc[5];
			      mpt[6] = idsloc[6];
			      mpt[7] = idsloc[7];
			    }
			    else {
			      str += " at reading element information";
			      return _error_info(str,ifs);
			    }

                            conn.insert( make_pair(idloc, mpt) );
			  }

			  std::getline(ifs, str);
			  ++count;
                        }
                    } else if (dim == 2) {
                        while (!ifs.eof() && count < numofNodes) {
                            //printf("core %d reading string %d: %s\n", _rank, count, str.c_str());
			  if ( !str.empty() ) {
			    if (sscanf(str.c_str(), "%d %lf %lf", &idloc, &ploc[0], &ploc[1]) != 3) {
			      str += " at reading nodal information";
			      return _error_info(str,ifs);
                            }
                            //printf("core %d reading line %d: %d %lf %lf %lf\n", _rank, count, idloc, ploc[0], ploc[1], ploc[2]);
			    v.insert( make_pair(idloc, ploc) );			    
			  }
			  std::getline(ifs, str);
			  ++count;                            
                        }

                        const int mptType = type;
                        int idsloc[mptType];
                        std::vector<unsigned int> mpt(mptType);
			int numInfo = mptType + 1;

                        count = 0;
                        while (!ifs.eof() && count < numofMpts) {

			  if ( !str.empty() ) {
			    // three nodes per element -- linear triangle
                            if (std::sscanf(str.c_str(), "%d %d %d %d",
					    &idloc, idsloc, idsloc + 1, idsloc + 2) == numInfo) {
			      mpt[0] = idsloc[0];
			      mpt[1] = idsloc[1];
			      mpt[2] = idsloc[2];
			    }
			    // four nodes per element -- bilinear rectangle
			    else if (std::sscanf(str.c_str(), "%d %d %d %d %d",
						 &idloc, idsloc, idsloc + 1, idsloc + 2, idsloc + 3) == numInfo) {
			      mpt[0] = idsloc[0];
			      mpt[1] = idsloc[1];
			      mpt[2] = idsloc[2];
			      mpt[3] = idsloc[3];
			    }
			    else {
			      str += " at reading element information";
			      return _error_info(str,ifs);
                            }

			    conn.insert( make_pair(idloc, mpt) );
			  }
			  std::getline(ifs, str);
			  ++count;
                        }
                    } else {
                        cerr << "dimension " << dim << " is not supported!" << endl;
                        assert(false);
                    }


                    if (!v.empty()) {
			_node_idmap.clear();
			_elm_idmap.clear();
			_node_inverse_idmap.clear();
			_elm_inverse_idmap.clear();

			if ( _flag == 0 ) {
			  _readTable(v, conn, _node_idmap, _elm_idmap, _node_inverse_idmap, _elm_inverse_idmap);
			}
			else {
			  _readTable_general(v, conn, _node_idmap, _elm_idmap, _node_inverse_idmap, _elm_inverse_idmap);
			}
			
			/**
			 * group information
			 */
			int grouptype;		    
			char substr[128];
			int numofEntity;
			set<Cell*> entityloc;
			map<int, Cell*>::iterator pC;
			
			while ( !ifs.eof() ) {
			  entityloc.clear();			  
			  //while (std::sscanf(str.c_str(), "%d %[^\r] %d", &grouptype, groupname, &numofEntity) == 3) {
			  //while (std::sscanf(str.c_str(), "%d %s %d", &grouptype, groupname, &numofEntity) == 3) {
			  while (std::sscanf(str.c_str(), "%d %[^\n]", &grouptype, substr) == 2) {
			    string groupname(substr);
			    int found = groupname.rfind(' ');
			    int last  = groupname.length() - 1;
			    if ( found != std::string::npos ) {
			      numofEntity = atoi(groupname.substr(found+1, last).c_str());
			    }
			    else {
			      break;
			    }
				
			    if ( groupname.find('\r') != std::string::npos ) {			    
			      groupname.erase(found-1, last);
			    }
			    else {
			      groupname.erase(found, last);
			    }
						    
			    if ( grouptype == 7 ) {
			        for ( int i = 0; i < numofEntity; ++i ) {
				    std::getline(ifs, str);
				    idloc = atoi(str.c_str());
			      
				    if ( (pC = _node_idmap.find(idloc)) != _node_idmap.end() ) {
				      entityloc.insert(pC->second);
				    }
				}
		      
				//if ( !entityloc.empty() ) {
				_node_groups.insert( make_pair(groupname, entityloc) );
				//}
			    }
			    else if ( grouptype == 8 ) {
			        for ( int i = 0; i < numofEntity; ++i ) {
				    std::getline(ifs, str);
				    idloc = atoi(str.c_str());
			      
				    if ( (pC = _elm_idmap.find(idloc)) != _elm_idmap.end() ) {
				      entityloc.insert(pC->second);
				    }
				}
		      
				//if ( !entityloc.empty() ) {
				_solid_groups.insert( make_pair(groupname, entityloc) );
				//}
			    }
			    else {
			        break;
			    }
			  }

			  std::getline(ifs, str);
			}	    
			
                    } else {
                        for (int i = 0; i <= dim; ++i) {
                            _SC.push_back(set<Cell*>());
                        }

			/**
			 * group information
			 */
			int grouptype;
			char substr[128];
			int numofEntity;

			while ( !ifs.eof() ) {
			  //while (std::sscanf(str.c_str(), "%d %[^\r] %d", &grouptype, groupname, &numofEntity) == 3) {    
			  //while (std::sscanf(str.c_str(), "%d %s %d", &grouptype, groupname, &numofEntity) == 3) {
			  while (std::sscanf(str.c_str(), "%d %[^\n]", &grouptype, substr) == 2) {
			    string groupname(substr);
			    int found = groupname.rfind(' ');
			    int last  = groupname.length() - 1;
			    if ( found != std::string::npos ) {
			      numofEntity = atoi(groupname.substr(found+1, last).c_str());
			    }
			    else {
			      break;
			    }
				
			    if ( groupname.find('\r') != std::string::npos ) {			    
			      groupname.erase(found-1, last);
			    }
			    else {
			      groupname.erase(found, last);
			    }

			    if ( grouptype == 7 ) {
			      for ( int i = 0; i < numofEntity; ++i ) {
				std::getline(ifs, str);
			      }
			      _node_groups.insert( make_pair(groupname, set<Cell*>()) );
			    }
			    else if ( grouptype == 8 ) {
			      for ( int i = 0; i < numofEntity; ++i ) {
				std::getline(ifs, str);
			      }
			      _solid_groups.insert( make_pair(groupname, set<Cell*>()) );
			    }
			    else {
			        break;
			    }
			  }

			  std::getline(ifs, str);
			}	    
                    }

                    ifs.close();

                    return true;
                } else {
                    printf("couldnot open file %s", filename);
                    return false;
                }
            }


            void GeneralGeometry::_readTable(const map<unsigned int, point_type> & v,
					     const map<unsigned int, vector<unsigned int> > & conn,
					     map<int, Cell *> & node_idmap,
					     map<int, Cell *> & elm_idmap,
					     map<Cell *, int> & node_inverse_idmap,
					     map<Cell *, int> & elm_inverse_idmap) {

                assert(conn.size() > 0);
                unsigned int i, j;

                map< unsigned int, point_type > ::const_iterator pv;
                for (pv = v.begin(); pv != v.end(); pv++) {
                    Cell * e = new Cell(0);
                    node_idmap.insert( make_pair(pv->first, e));
		    node_inverse_idmap.insert( make_pair(e, pv->first) );
                    _XNodes.insert(make_pair(e, pv->second));
                }

                set<set<Cell *> > CT;
                map< set<Cell *>, unsigned int > ids;
                map< unsigned int, vector<unsigned int> > ::const_iterator pc;
                for (pc = conn.begin(); pc != conn.end(); pc++) {
                    set<Cell *> k;
                    const vector<unsigned int> & CA = pc->second;
                    for (j = 0; j < CA.size(); j++) {
                        k.insert(node_idmap.find(CA[j])->second);
                    }

                    CT.insert(k);
                    ids.insert(make_pair(k, pc->first));
                }

                if (CT.size() == 0) return;

		set<set<Cell *> >::const_iterator pCT;
                unsigned int n = (unsigned int) CT.begin()->size();
                for (pCT = CT.begin(); pCT != CT.end(); pCT++) {
                    if (pCT->size() > n) n = pCT->size();
                }

                unsigned int dim = n - 1;
                if (dim == 0) return;

                unsigned int p, q;
                if (_SC.size() == 0) {
                    for (p = 0; p <= dim; p++)
                        _SC.push_back(set<Cell *>());
                } else {
                    assert(n == _SC.size());
                    set<Cell *>::iterator pE;
                    for (pE = _SC[0].begin(); pE != _SC[0].end(); pE++)
                        (*pE)->Coboundary().clear();
                    for (p = 1; p <= dim; p++) {
                        set<Cell *>::iterator pE;
                        for (pE = _SC[p].begin(); pE != _SC[p].end(); pE++) delete *pE;
                        _SC[p].clear();
                    }
                }

                vector<set<set<Cell *> > > K;
                for (p = 0; p <= dim; p++)
                    K.push_back(set<set<Cell *> >());
                for (pCT = CT.begin(); pCT != CT.end(); pCT++)
		  K[pCT->size() - 1].insert(*pCT);		  

                vector<map<set<Cell *>, Cell *> > e_map;
                for (p = 0; p <= dim; p++)
                    e_map.push_back(map<set<Cell *>, Cell *>());

                set<set<Cell *> >::iterator pK;
                set<Cell *>::iterator pk;

                if (dim > 1) {
                    for (p = dim; p > 1; p--) {
                        q = p - 1;
                        for (pK = K[p].begin(); pK != K[p].end(); pK++) {
                            set<Cell *> k = *pK;
                            for (pk = k.begin(); pk != k.end(); pk++) {
                                set<Cell *> m = k;
                                m.erase(*pk);
                                K[q].insert(m);
                            }
                        }
                    }
                }

                if (_SC[0].size() == 0) {
                    for (pK = K[dim].begin(); pK != K[dim].end(); pK++) {
                        set<Cell *> k = *pK;
                        for (pk = k.begin(); pk != k.end(); pk++)
                            _SC[0].insert(*pk);
                    }
                }

                for (p = 1; p <= dim; p++) {
                    for (pK = K[p].begin(); pK != K[p].end(); pK++) {
                        Cell *e = new Cell(p);
                        e_map[p][*pK] = e;
                        _SC[p].insert(e);
                    }
                }

                for (pK = K[dim].begin(); pK != K[dim].end(); pK++) {
                    elm_idmap.insert(make_pair(ids.find(*pK)->second, e_map[dim].find(*pK)->second));
		    elm_inverse_idmap.insert( make_pair(e_map[dim].find(*pK)->second, ids.find(*pK)->second) );
                }

                int sign;

                for (pK = K[1].begin(); pK != K[1].end(); pK++) {
                    set<Cell *> k = *pK;
                    set<Cell *>::iterator pk;
                    Cell *e = e_map[1][k];
                    Chain &e_bo = e->Boundary();
                    for (pk = k.begin(), sign = 1; pk != k.end(); pk++, sign = -sign) {
                        set<Cell *> m = k;
                        m.erase(*pk);
                        Cell *f = *m.begin();
                        e_bo[f] = sign;
                        Cochain &f_co = f->Coboundary();
                        f_co[e] = sign;
                    }
                }

                if (dim > 1) {
                    for (p = 2; p <= dim; p++) {
                        q = p - 1;
                        for (pK = K[p].begin(); pK != K[p].end(); pK++) {
                            set<Cell *> k = *pK;
                            set<Cell *>::iterator pk;
                            Cell *e = e_map[p][k];
                            Chain &e_bo = e->Boundary();
                            for (pk = k.begin(), sign = 1; pk != k.end(); pk++, sign = -sign) {
                                set<Cell *> m = k;
                                m.erase(*pk);
                                Cell *f = e_map[q][m];
                                e_bo[f] = sign;
                                Cochain &f_co = f->Coboundary();
                                f_co[e] = sign;
                            }
                        }
                    }
                }

                return;
            }

            void GeneralGeometry::_readTable_general(const map<unsigned int, point_type> & v,
						     const map<unsigned int, vector<unsigned int> > & conn,
						     map<int, Cell *> & node_idmap,
						     map<int, Cell *> & elm_idmap,
						     map<Cell *, int> & node_inverse_idmap,
						     map<Cell *, int> & elm_inverse_idmap) {

	        assert(conn.size() > 0);
                unsigned int i, j;
		unsigned int dim = v.begin()->second.size();
		
                map< unsigned int, point_type > ::const_iterator pv;
                for (pv = v.begin(); pv != v.end(); pv++) {
                    Cell * e = new Cell(0);
                    node_idmap.insert(make_pair(pv->first, e));
		    node_inverse_idmap.insert( make_pair(e, pv->first) );
                    _XNodes.insert(make_pair(e, pv->second));
                }

                set<set<Cell *> > CT;
                map< set<Cell *>, unsigned int > ids;
                map< unsigned int, vector<unsigned int> > ::const_iterator pc;
                for (pc = conn.begin(); pc != conn.end(); pc++) {
                    set<Cell *> k;
                    const vector<unsigned int> & CA = pc->second;
                    for (j = 0; j < CA.size(); j++) {
                        k.insert(node_idmap.find(CA[j])->second);
                    }

                    CT.insert(k);
                    ids.insert(make_pair(k, pc->first));
                }

                if (CT.size() == 0) return;

		set<set<Cell *> >::const_iterator pCT;
                unsigned int p, q;
                if (_SC.size() == 0) {
                    for (p = 0; p <= dim; p++)
                        _SC.push_back(set<Cell *>());
                } else {
                    assert(dim+1 == _SC.size());
                    set<Cell *>::iterator pE;
                    for (pE = _SC[0].begin(); pE != _SC[0].end(); pE++)
                        (*pE)->Coboundary().clear();
                    for (p = 1; p <= dim; p++) {
                        set<Cell *>::iterator pE;
                        for (pE = _SC[p].begin(); pE != _SC[p].end(); pE++) delete *pE;
                        _SC[p].clear();
                    }
                }

                vector<set<set<Cell *> > > K;
                for (p = 0; p <= dim; p++)
                    K.push_back(set<set<Cell *> >());
                for (pCT = CT.begin(); pCT != CT.end(); pCT++)
		    K[dim].insert(*pCT);		  

                vector<map<set<Cell *>, Cell *> > e_map;
                for (p = 0; p <= dim; p++)
                    e_map.push_back(map<set<Cell *>, Cell *>());

                set<set<Cell *> >::iterator pK;
                set<Cell *>::iterator pk;

                if (dim > 1) {
		  set<Cell *>::iterator pk1, pk2;
		  for (pK = K[dim].begin(); pK != K[dim].end(); pK++){
		    set<Cell*> kd = *pK;                		
		    for (pk1 = kd.begin(); pk1 != kd.end(); ) {
		      for (pk2 = kd.begin(); pk2 != kd.end(); pk2++) {
			if ( *pk2 != *pk1 ) {
			  set<Cell*> m;
			  m.insert(*pk1);
			  m.insert(*pk2);
			  K[1].insert(m);
			}
		      }

		      kd.erase(*pk1);
		      pk1++;
		    }
		  }
                }

		if ( dim > 2 ) {
		  assert(false);
                }

                if (_SC[0].size() == 0) {
                    for (pK = K[dim].begin(); pK != K[dim].end(); pK++) {
                        set<Cell *> k = *pK;
                        for (pk = k.begin(); pk != k.end(); pk++)
                            _SC[0].insert(*pk);
                    }
                }

                for (p = 1; p <= dim; p++) {
                    for (pK = K[p].begin(); pK != K[p].end(); pK++) {
                        Cell *e = new Cell(p);
                        e_map[p][*pK] = e;
                        _SC[p].insert(e);
                    }
                }

                for (pK = K[dim].begin(); pK != K[dim].end(); pK++) {
                    elm_idmap.insert(make_pair(ids.find(*pK)->second, e_map[dim].find(*pK)->second));
		    elm_inverse_idmap.insert( make_pair(e_map[dim].find(*pK)->second, ids.find(*pK)->second) );
                }

                int sign;

                for (pK = K[1].begin(); pK != K[1].end(); pK++) {
                    set<Cell *> k = *pK;
                    set<Cell *>::iterator pk;
                    Cell *e = e_map[1][k];
                    Chain &e_bo = e->Boundary();
                    for (pk = k.begin(), sign = 1; pk != k.end(); pk++, sign = -sign) {
                        set<Cell *> m = k;
                        m.erase(*pk);
                        Cell *f = *m.begin();
                        e_bo[f] = sign;
                        Cochain &f_co = f->Coboundary();
                        f_co[e] = sign;
                    }
                }

                if (dim == 2) {
		  set<Cell *>::iterator pk1, pk2;
		  for (pK = K[dim].begin(); pK != K[dim].end(); pK++){
		    set<Cell*> kd = *pK;
		    Cell *e = e_map[dim][kd];
		    Chain &e_bo = e->Boundary();
		    for (pk1 = kd.begin(), sign=1; pk1 != kd.end(); sign=-sign) {
		      for (pk2 = kd.begin(); pk2 != kd.end(); pk2++) {
			if ( *pk2 != *pk1 ) {
			  set<Cell*> m;
			  m.insert(*pk1);
			  m.insert(*pk2);
			  Cell *f = e_map[1][m];
			  e_bo[f] = sign;
			  Cochain &f_co = f->Coboundary();
			  f_co[e] = sign;
			}
		      }

		      kd.erase(*pk1);
		      pk1++;
		    }
		  }		  
		}
		
                return;
            }
	  
	  void GeneralGeometry::_getVertices(Geometry::Cell* e, 
					     set<Geometry::Cell*> & vs) { // note make sure vs is empty
	    if (e->dim() == 0) {
	      vs.insert(e);
	    } else {
	      Geometry::Chain & bo = e->Boundary();
	      for (Geometry::Chain::iterator pe = bo.begin();
		   pe != bo.end(); pe++) {
                _getVertices(pe->first, vs);
	      }
	    }
	    
	    return;
	  }

	  //
	  // Patch the geometry
	  //
#if defined(_M4EXTREME_VTK_)

	  void GeneralGeometry::SetDegenTol(double tol) {
	    _DEGEN_TOL = tol;
	  }

	  // initialize the mesh of the patch geometry;
	  // cloned patch_nodes are translated by vector N;
	  // the vtkDelaunay mesh generator is called to generate 
	  // a delaunay mesh for the convex hull of the nodal set 
	  // {original patch_nodes, support_nodes, translated patch_nodes}

#if defined(_M4EXTREME_MPI_)
	  //!!!bug: need to add support for boundary conditions for patch nodes
	  void GeneralGeometry::InitializePatch(const set<Cell*> & patch_nodes,
						const set<Cell*> & support_nodes,
						const Set::VectorSpace::Vector & N,
						int id_base,
						int numofCores) {
	    int dim = N.size();

	    /**
	     * patch nodes
	     */
	    vector<double> coords;
	    vector<int> ids;
	    vector<point_type> PatchNodes;
	    vector<Cell*> PatchCells, newCells;

	    // get the patch nodes from all the partitions
	    for ( set<Cell*>::const_iterator pC = patch_nodes.begin(); 
		  pC != patch_nodes.end(); pC++) {
	      const point_type & xloc = _XNodes.find(*pC)->second;	      
	      ids.push_back(_node_inverse_idmap.find(*pC)->second);
	      for ( int i = 0; i < dim; ++i ) {
		coords.push_back(xloc[i]);
	      }
	    }

	    int sendcount = coords.size();
	    int recvcount[numofCores];
	    MPI_Allgather(&sendcount, 1, MPI_INT, recvcount, 1, MPI_INT, MPI_COMM_WORLD);

	    int displs[numofCores];
	    int totalrecv= 0;
	    for (int i = 0; i < numofCores; i++) {
	      displs[i] = totalrecv;
	      totalrecv += recvcount[i];
	    }

	    double * recv = (double*)malloc(totalrecv * sizeof(double));

	    int ERR = MPI_SUCCESS;
	    if ( (ERR = MPI_Allgatherv(&coords.front(), sendcount, MPI_DOUBLE, recv, recvcount, displs, 
				       MPI_DOUBLE, MPI_COMM_WORLD)) != MPI_SUCCESS ){
	      MPI_Abort(MPI_COMM_WORLD, ERR);
	    }

	    // get the ids of the patch nodes from all the partitions
	    sendcount /= dim;
	    totalrecv = 0;
	    for (int i = 0; i < numofCores; i++) {
	      recvcount[i] /= dim;
	      displs[i] = totalrecv;
	      totalrecv += recvcount[i];
	    }

	    int * ids_recv = (int*)malloc(totalrecv * sizeof(int));
	    if ( (ERR = MPI_Allgatherv(&ids.front(), sendcount, MPI_INT, ids_recv, recvcount, displs, 
				       MPI_INT, MPI_COMM_WORLD)) != MPI_SUCCESS ){
	      MPI_Abort(MPI_COMM_WORLD, ERR);
	    }

	    // construct the id pool for new patch nodes and 
	    // merge the patch nodes shared by different partitions    
	    map<int, Cell*>::iterator pID;
	    Cell * etemp = NULL;	   
	    map<int, int> idpool;
	    map<int, point_type> all_patch_nodes;
	    int count = 0;
	    
	    point_type temp(dim);
	    for ( int i = 0; i < totalrecv; i++ ){
	      if ( idpool.find(ids_recv[i]) == idpool.end() ) {
		idpool.insert( make_pair(ids_recv[i], id_base+count) );

		for ( int j = 0; j < dim; j++ ){
		  temp[j] =  recv[i*dim + j];
		}
		all_patch_nodes.insert( make_pair(ids_recv[i], temp) );

		++count;
	      }
	    }

	    free(ids_recv);
	    free(recv);

	    // insert new patch nodes (patch back)
	    // note: idpool and all_patch_nodes store the elements
	    //       in the same order based on the value of ids
	    map<int, point_type>::iterator pN = all_patch_nodes.begin();
	    for ( map<int, int>::iterator pI = idpool.begin();
		  pI != idpool.end(); pI++, pN++ ) {
	      // original patch nodes (patch front)
	      if ( (pID = _node_idmap.find(pI->first)) == _node_idmap.end() ) {
		etemp = new Cell(0);		
		PatchCells.push_back(etemp);
		newCells.push_back(etemp);
	      }
	      else {
		PatchCells.push_back(pID->second);
	      }

	      PatchNodes.push_back(pN->second);

	      // new patch nodes by a translation in N
	      etemp = new Cell(0);
	      PatchCells.push_back(etemp);
 
	      point_type xnew = pN->second + N;
	      PatchNodes.push_back(xnew);

	      // only save the new patch nodes corresponding to the 
	      // original patch nodes on the current partition
	      if ( pID != _node_idmap.end() ) {
		_node_idmap.insert( make_pair(idpool.find(pID->first)->second, etemp) );
		_XNodes.insert( make_pair(etemp, xnew) );
		_XPatch.push_back(xnew);
		_CPatch_front.push_back(pID->second);
		_CPatch_back.push_back(etemp);
	      }
	      else {
		newCells.push_back(etemp);
	      }
	    }

	    all_patch_nodes.clear();
	    idpool.clear();

	    /**
	     * support nodes
	     */
	    coords.clear();
	    ids.clear();
	    for ( set<Cell*>::const_iterator pC = support_nodes.begin(); 
		  pC != support_nodes.end(); pC++) {
	      const point_type & xloc = _XNodes.find(*pC)->second;
	      ids.push_back(_node_inverse_idmap.find(*pC)->second);	      
	      for ( int i = 0; i < dim; ++i ) {
		coords.push_back(xloc[i]);
	      }
	    }

	    sendcount = coords.size();	    
	    MPI_Allgather(&sendcount, 1, MPI_INT, recvcount, 1, MPI_INT, MPI_COMM_WORLD);

	    totalrecv= 0;
	    for (int i = 0; i < numofCores; i++) {
	      displs[i] = totalrecv;
	      totalrecv += recvcount[i];
	    }

	    recv = (double*)malloc(totalrecv * sizeof(double));
	    if ( (ERR = MPI_Allgatherv(&coords.front(), sendcount, MPI_DOUBLE, recv, recvcount, displs, 
				       MPI_DOUBLE, MPI_COMM_WORLD)) != MPI_SUCCESS ){
	      MPI_Abort(MPI_COMM_WORLD, ERR);
	    }
	 
	    sendcount /= dim;
	    totalrecv = 0;
	    for (int i = 0; i < numofCores; i++) {
	      recvcount[i] /= dim;
	      displs[i] = totalrecv;
	      totalrecv += recvcount[i];
	    }

	    ids_recv = (int*)malloc(totalrecv * sizeof(int));
	    if ( (ERR = MPI_Allgatherv(&ids.front(), sendcount, MPI_INT, ids_recv, recvcount, displs, 
				       MPI_INT, MPI_COMM_WORLD)) != MPI_SUCCESS ){
	      MPI_Abort(MPI_COMM_WORLD, ERR);
	    }

	    set<int> all_ids;
	    for ( int i = 0; i < totalrecv; i++ ) {
	      if ( all_ids.find(ids_recv[i]) == all_ids.end() ) {
		all_ids.insert(ids_recv[i]);

		for ( int j = 0; j < dim; j++ ){
		  temp[j] =  recv[i*dim+j];
		}
		PatchNodes.push_back(temp);

		if ( (pID = _node_idmap.find(ids_recv[i])) == _node_idmap.end() ) {
		  etemp = new Cell(0);
		  PatchCells.push_back(etemp);
		  newCells.push_back(etemp);
		}
		else {
		  PatchCells.push_back(pID->second);
		}
	      }
	    }

	    all_ids.clear();
	    free(ids_recv);
	    free(recv);
  
	    /**
	     * Delaunay triangulation
	     */
	    vtkDelaunay2m4extreme PatchTriangulation;
	    vector< vector<unsigned int> > PatchConn = PatchTriangulation(dim, PatchNodes);

	    _DEGEN_TOL *= Norm(N);
	    for (int i=0; i<PatchConn.size(); i++)
	    {
	      set<Cell *> k;
	      vector<Set::VectorSpace::Vector> PSloc;
	      for (int j=0; j<PatchConn[i].size(); j++)
	      { 		
		etemp = PatchCells[PatchConn[i][j]];
		if ( patch_nodes.find(etemp) == patch_nodes.end() &&
		     support_nodes.find(etemp) == support_nodes.end() &&
		     std::find(_CPatch_back.begin(), _CPatch_back.end(), etemp) == _CPatch_back.end() ) {
		  k.clear();
		  break;
		}
		else {
		  k.insert(etemp);
		  PSloc.push_back(PatchNodes[PatchConn[i][j]]);
		}
	      }

	      if ( !k.empty() ) {
		// check if it is a degenerate case
		double volume = Geometry::Utils::Simplex(PSloc).Volume();
	      
		if (  volume > _DEGEN_TOL ) {
		  _CTPatch.insert(k);
		}
		else {
		  cout << "ignore patch element with volume = " << volume << endl;
		}
	      }
	    }

	    _SPatches.push_back(new SimplicialComplex(dim));		
	    _SPatches.back()->ReadTable(_CTPatch); 

	    // clean up
	    for ( int i = 0; i < newCells.size(); i++ ) {
	      delete newCells[i];
	    }
	    
	    return;
	  }

#else

	  void GeneralGeometry::InitializePatch(const map<Cell*, Set::VectorSpace::Hom*> & patch_nodes,
						const set<Cell*> & support_nodes,
						const Set::VectorSpace::Vector & N,
						int id_base) {

	    if ( patch_nodes.empty() ) return;
	    
	    double normN = Norm(N);
	    if ( normN < 1.0e-8 ) return;

	    int dim = _XNodes.begin()->second.size();
	    int count = 0;
	    _DEGEN_TOL *= normN;	    
	    
	    vector<point_type> PatchNodes;
	    vector<Cell*> PatchCells;

	    for ( map<Cell*, Set::VectorSpace::Hom*>::const_iterator pC = patch_nodes.begin(); 
		  pC != patch_nodes.end(); pC++) {

	      const point_type & xloc = _XNodes.find(pC->first)->second;
	      PatchNodes.push_back(xloc);
	      PatchCells.push_back(pC->first);

	      Cell * e = new Cell(0);	      
	      point_type xnew = xloc + N;
	      
	      _node_idmap.insert( make_pair(id_base + count, e) );
	      _XNodes.insert( make_pair(e, xnew) );
	      _XPatch.push_back(xnew);
	      _CPatch_front.push_back(pC->first);
	      _CPatch_back.push_back(e);
	      _CPatch_back_BCs.push_back(pC->second);

	      PatchNodes.push_back(xnew);
	      PatchCells.push_back(e);

	      count++;
	    }
	    
	    for ( set<Cell*>::const_iterator pC = support_nodes.begin();
		  pC != support_nodes.end(); pC++) {
	      const point_type & xloc = _XNodes.find(*pC)->second;
	      PatchNodes.push_back(xloc);
	      PatchCells.push_back(*pC);
	    }

	    vtkDelaunay2m4extreme PatchTriangulation;
	    vector< vector<unsigned int> > PatchConn = PatchTriangulation(dim, PatchNodes);
	   
	    for (int i=0; i<PatchConn.size(); i++)
	    {
	      set<Cell *> k;
	      vector<Set::VectorSpace::Vector> PSloc;
	      for (int j=0; j<PatchConn[i].size(); j++)
	      { 		
		k.insert(PatchCells[PatchConn[i][j]]);
		PSloc.push_back(PatchNodes[PatchConn[i][j]]);
	      }

	      // check if it is a degenerate case
	      double volume = Geometry::Utils::Simplex(PSloc).Volume();
	      
	      if (  volume > _DEGEN_TOL ) {
		_CTPatch.insert(k);
	      }
	      else {
		cout << "ignore patch element with volume = " << volume << endl;
	      }
	    }
		
	    _SPatches.push_back(new SimplicialComplex(dim));		
	    _SPatches.back()->ReadTable(_CTPatch); 

	    return;
	  }	  

#endif
	  //!!! need to add support for boundary conditions
	  void GeneralGeometry::InitializePatch(const set<Cell*> & front_nodes,
						const set<Cell*> & back_nodes,
						set<Cell*> & SD) {

	    int dim = _XNodes.begin()->second.size();

	    for ( set<Cell*>::const_iterator pC = front_nodes.begin();
		  pC != front_nodes.end(); pC++){
	      _CPatch_front.push_back(*pC);
	    }

	    for ( set<Cell*>::const_iterator pC = back_nodes.begin();
		  pC != back_nodes.end(); pC++){
	      _CPatch_back.push_back(*pC);
	      _XPatch.push_back(_XNodes.find(*pC)->second);
	    }

	    for ( set<Cell*>::iterator pC = SD.begin();
		  pC != SD.end(); pC++ ) {
	      set<Cell*> k;
	      _getVertices(*pC, k);
	      _CTPatch.insert(k);
	    }

	    _SPatches.push_back(new SimplicialComplex(dim));		
	    _SPatches.back()->ReadTable(_CTPatch); 

	    return;
	  }
	  

	  // insert a new patch with new nodes located at the coordinates
	  // obtrained in the reference configuration
	  void GeneralGeometry::InsertPatch(
#if defined(_M4EXTREME_MPI_)
				 int numofCores,
#endif
				 int id_base) {
 
	    int count = 0;

#if defined(_M4EXTREME_MPI_)
	    map<int, int> idpool;
	    vector<int> ids;
	    for ( int i = 0; i < _CPatch_back.size(); ++i ) {
	      ids.push_back(_node_inverse_idmap.find(_CPatch_back[i])->second);
	    }

	    int sendcount = ids.size();
	    int recvcount[numofCores];
	    MPI_Allgather(&sendcount, 1, MPI_INT, recvcount, 1, MPI_INT, MPI_COMM_WORLD);
	    
	    int displs[numofCores];
	    int totalrecv= 0;
	    for (int i = 0; i < numofCores; i++) {
	      displs[i] = totalrecv;
	      totalrecv += recvcount[i];
	    }
	    
	    int * ids_recv = (int*)malloc(totalrecv * sizeof(int));

	    int ERR = MPI_SUCCESS;
	    if ( (ERR = MPI_Allgatherv(&ids.front(), sendcount, MPI_INT, ids_recv, recvcount, displs, 
				       MPI_INT, MPI_COMM_WORLD)) != MPI_SUCCESS ){
	      MPI_Abort(MPI_COMM_WORLD, ERR);
	    }    
	    
	    for ( int i = 0; i < totalrecv; ++i ) {
	      if ( idpool.find(ids_recv[i]) == idpool.end() ) {
		idpool.insert( make_pair(ids_recv[i], id_base+count) );
		++count;
	      }
	    }
#endif

	    if ( _XPatch.empty() ) return;

	    int dim = _XPatch[0].size();
	    map<Cell*, Cell*> cellmap;	   
	    
	    for ( int i = 0; i < _XPatch.size(); ++i ) {
	      Cell * e = new Cell(0);

#if defined(_M4EXTREME_MPI_)
	      _node_idmap.insert( make_pair(idpool.find(ids[i])->second, e) );
#else
	      _node_idmap.insert( make_pair(id_base + count, e) );
#endif	  

	      _XNodes.insert( make_pair(e, _XPatch[i]) );

	      cellmap.insert( make_pair(_CPatch_front[i], e) );
	      _CPatch_front[i] = _CPatch_back[i];
	      _CPatch_back[i]  = e;

	      count++;
	    }

	    map<Cell*, Cell*>::iterator pPC;
	    for ( set<set<Cell*> >::iterator pS = _CTPatch.begin(); 
		  pS != _CTPatch.end(); ) {	      

	      set<Cell*> vertices;
	      for ( set<Cell*>::iterator pC = pS->begin(); pC != pS->end(); pC++) {
		if ( (pPC = cellmap.find(*pC)) != cellmap.end() ) {
		  vertices.insert(pPC->second);
		}
		else {
		  vertices.insert(*pC);
		}
	      }
	      
	      _CTPatch.erase(pS++);
	      _CTPatch.insert(vertices);
	    }		

	    _SPatches.push_back(new SimplicialComplex(dim));		
	    _SPatches.back()->ReadTable(_CTPatch);    

	    return;
	  }

	  // insert a new patch with nodes whose coordinates are calculated
	  // as xpatch[i] = _CPatch_back[i]+translation
	  void GeneralGeometry::InsertPatch(
#if defined(_M4EXTREME_MPI_)
					    int numofCores,
#endif
					    int id_base, 
					    const Set::VectorSpace::Vector & N) {

	    int count = 0;

#if defined(_M4EXTREME_MPI_)
	    map<int, int> idpool;
	    vector<int> ids;
	    for ( int i = 0; i < _CPatch_back.size(); ++i ) {
	      ids.push_back(_node_inverse_idmap.find(_CPatch_back[i])->second);
	    }

	    int sendcount = ids.size();
	    int recvcount[numofCores];
	    MPI_Allgather(&sendcount, 1, MPI_INT, recvcount, 1, MPI_INT, MPI_COMM_WORLD);
	    
	    int displs[numofCores];
	    int totalrecv= 0;
	    for (int i = 0; i < numofCores; i++) {
	      displs[i] = totalrecv;
	      totalrecv += recvcount[i];
	    }
	    
	    int * ids_recv = (int*)malloc(totalrecv * sizeof(int));

	    int ERR = MPI_SUCCESS;
	    if ( (ERR = MPI_Allgatherv(&ids.front(), sendcount, MPI_INT, ids_recv, recvcount, displs, 
				       MPI_INT, MPI_COMM_WORLD)) != MPI_SUCCESS ){
	      MPI_Abort(MPI_COMM_WORLD, ERR);
	    }    
	    
	    for ( int i = 0; i < totalrecv; ++i ) {
	      if ( idpool.find(ids_recv[i]) == idpool.end() ) {
		idpool.insert( make_pair(ids_recv[i], id_base+count) );
		++count;
	      }
	    }
#endif

	    if ( _CPatch_back.empty() ) return;

	    int dim = _XPatch[0].size();
	    map<Cell*, Cell*> cellmap;	    

	    _XPatch.clear();
	    
	    for ( int i = 0; i < _CPatch_back.size(); ++i ) {
	      Cell * e = new Cell(0);
	      point_type xnew = _XNodes.find(_CPatch_back[i])->second + N;

#if defined(_M4EXTREME_MPI_)
	      _node_idmap.insert( make_pair(idpool.find(ids[i])->second, e) );
#else
	      _node_idmap.insert( make_pair(id_base + count, e) );
#endif	    

	      _XPatch.push_back(xnew);
	      _XNodes.insert( make_pair(e, xnew) );

	      cellmap.insert( make_pair(_CPatch_front[i], e) );// !!!bug: assume xpatch and CPatch_front in the same order
	      _CPatch_front[i] = _CPatch_back[i];
	      _CPatch_back[i]  = e;
	      count++;
	    }

	    map<Cell*, Cell*>::iterator pPC;
	    for ( set<set<Cell*> >::iterator pS = _CTPatch.begin(); 
		  pS != _CTPatch.end(); ) {	      

	      set<Cell*> vertices;
	      for ( set<Cell*>::iterator pC = pS->begin(); pC != pS->end(); pC++) {
		if ( (pPC = cellmap.find(*pC)) != cellmap.end() ) {
		  vertices.insert(pPC->second);
		}
		else {
		  vertices.insert(*pC);
		}
	      }
	      
	      _CTPatch.erase(pS++);
	      _CTPatch.insert(vertices);
	    }		

	    _SPatches.push_back(new SimplicialComplex(dim));		
	    _SPatches.back()->ReadTable(_CTPatch); 

	    return;
	  }
	  
	  CellComplex * GeneralGeometry::GetPatchCellComplex() const {
	    if ( !_SPatches.empty() ) {
	      return _SPatches.back();
	    }
	    else {
	      return NULL;
	    }
	  }

	  const vector<Cell*> & GeneralGeometry::GetPatchBackNodes() const {
	    return _CPatch_back;
	  }

	  const vector<Set::VectorSpace::Hom*> & GeneralGeometry::GetPatchBackNodesBC() const {
	    return _CPatch_back_BCs;
	  }

	  const vector<Cell*> & GeneralGeometry::GetPatchFrontNodes() const {
	    return _CPatch_front;
	  }

	  void GeneralGeometry::GetPatchNodes(map<Cell*, point_type> & patch_nodes) const {
	    patch_nodes.clear();

	    for ( int i = 0; i < _XPatch.size(); ++i ) {
	      patch_nodes.insert( make_pair(_CPatch_back[i], _XPatch[i]) );
	    }
	    
	    return;
	  }

#if defined(_M4EXTREME_MPI_)
	  int GeneralGeometry::GetMaxId(int numofCores) {
	    vector<int> ids(_node_idmap.size());
	    for ( map<int, Geometry::Cell*>::const_iterator pI = _node_idmap.begin();
		  pI != _node_idmap.end(); pI++ ) {
	      ids.push_back(pI->first);
	    }
	    
	    int sendcount = ids.size();
	    int recvcount[numofCores];
	    MPI_Allgather(&sendcount, 1, MPI_INT, recvcount, 1, MPI_INT, MPI_COMM_WORLD);
	    int displs[numofCores];
	    int totalrecv= 0;
	    for (int i = 0; i < numofCores; i++) {
	      displs[i] = totalrecv;
	      totalrecv += recvcount[i];
	    }
	    
	    int * ids_recv = (int*)malloc(totalrecv * sizeof(int));

	    int ERR = MPI_SUCCESS;
	    if ( (ERR = MPI_Allgatherv(&ids.front(), sendcount, MPI_INT, ids_recv, recvcount, displs, 
				       MPI_INT, MPI_COMM_WORLD)) != MPI_SUCCESS ){
	      MPI_Abort(MPI_COMM_WORLD, ERR);
	    }

	    return *std::max_element(ids_recv, ids_recv+totalrecv);
	  }
#else
	  int GeneralGeometry::GetMaxId() {
	    int maxId = 0;
	    for ( map<int, Cell*>::const_iterator pI = _node_idmap.begin();
		  pI != _node_idmap.end(); pI++ ) {
	      if( maxId < pI->first) {
		maxId = pI->first;
	      }
	    }
	    
	    return maxId;
	  }
#endif

#endif //_M4EXTREME_VTK_

        } // end_of_namespace_Mesh
    } // end_of_namespace_Solid
} // end_of_namespace_Geometry
