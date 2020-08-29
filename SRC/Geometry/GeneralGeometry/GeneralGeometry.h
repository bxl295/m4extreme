// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
////////////////////////////////////////////////////////////////////////////

#if !defined(GEOMETRY_MESH_GENERALGEOMETRY_H__INCLUDED_)
#define GEOMETRY_MESH_GENERALGEOMETRY_H__INCLUDED_

#pragma once

#include "Geometry/Solid/Mesh/Mesh.h"
#include "Geometry/Algebraic/SimplicialComplex/SimplicialComplex.h"

using namespace std;

namespace Geometry {
    namespace Solid {
        namespace Mesh {

            class GeneralGeometry : public Geometry::Mesh {
            public:
                typedef Set::Euclidean::Orthonormal::Point point_type;
                typedef map<Cell *, point_type> nodeset_type;

                enum FILE_TYPE {
		  M4EXTREME_NEUTRAL = 0, M4EXTREME_NEUTRAL_MINI = 1, M4EXTREME_NEUTRAL_FULL = 2, M4EXTREME_NEUTRAL_FULL_II = 3,
		  FEMAP_NEUTRAL = 4, FEMAP_NEUTRAL_FULL = 5
                };

            private:
		FILE_TYPE _type;
		int _flag; //Simplicial complex: _flag=0; Cell complex: _flag != 0
                SimplicialComplex _SC;
                map<Cell *, Cell *> _h;
                map<Cell *, point_type> _XNodes;

                map<string, set<Cell*> > _node_groups, _solid_groups;
                map<int, Cell *> _node_idmap;
		map<Cell*, int>  _node_inverse_idmap;
		map<int, Cell *> _elm_idmap;
		map<Cell*, int>  _elm_inverse_idmap;
		
		// patches
		vector<point_type> _XPatch;
		vector<Cell*> _CPatch_front, _CPatch_back;
		vector<Set::VectorSpace::Hom*> _CPatch_back_BCs;
		set< set<Cell*> > _CTPatch;
		vector<SimplicialComplex*> _SPatches;	

            public:
                GeneralGeometry();
                virtual ~GeneralGeometry();

                GeneralGeometry(int, const char*, const FILE_TYPE &, int flag = 0);

                const CellComplex & GetCellComplex() const;
                CellComplex & GetCellComplex();

                const map<Cell *, Cell *> & GetBinding() const;
                map<Cell *, Cell *> & GetBinding();
                
                const map<Cell *, point_type> & GetPointSet() const;
                map<Cell *, point_type> & GetPointSet();

		bool FindElementGroup(const string & name) const;
                const set<Cell*> & GetElementGroup(const string & name) const;

		bool FindNodeGroup(const string & name) const;
                const set<Cell*> & GetNodeGroup(const string & name) const;
   
		int GetElementGroupSize() const;
		int GetNodeGroupSize() const;
                const map<string, set<Cell*> > & GetAllElmentGroups() const;
		const map<string, set<Cell*> > & GetAllNodeGroups() const;

                const map<int, Cell *> & GetIdMap() const;
		const map<Cell *, int> & GetInverseIdMap() const;
                const map<int, Cell *> & GetElementIds() const;
		const map<Cell *, int> & GetInverseElementIds() const;		

		//
		// patch the geometry
		//

#if defined(_M4EXTREME_VTK_)
#if defined(_M4EXTREME_MPI_)
		void InitializePatch(const set<Cell*>  & patch_nodes,
				     const set<Cell*> & support_nodes,
				     const Set::VectorSpace::Vector & N,
				     int id_base,
				     int numofCores);
#else
		void InitializePatch(const map<Cell*, Set::VectorSpace::Hom*> & patch_nodes,
				     const set<Cell*> & support_nodes,
				     const Set::VectorSpace::Vector & N,
				     int id_base);
#endif
		void InitializePatch(const set<Cell*> & front_nodes,
				     const set<Cell*> & back_nodes,
				     set<Cell*> & SD);

		void InsertPatch(
#if defined(_M4EXTREME_MPI_)
				 int numofCores,
#endif				 
				 int id_base);
		void InsertPatch(
#if defined(_M4EXTREME_MPI_)
				 int numofCores,
#endif
				 int id_base, 
				 const Set::VectorSpace::Vector & N);

		CellComplex * GetPatchCellComplex() const;
		void GetPatchNodes(map<Cell*, point_type> & patch_nodes) const;
		const vector<Cell*> & GetPatchFrontNodes() const;
		const vector<Cell*> & GetPatchBackNodes() const;
		const vector<Set::VectorSpace::Hom*> & GetPatchBackNodesBC() const;
		
		void SetDegenTol(double);

#if defined(_M4EXTREME_MPI_)
		int GetMaxId(int numofCores);
#else
		int GetMaxId();
#endif

	  private:
		double _DEGEN_TOL;

#endif

	  private:
		bool _error_info(const string &, ifstream &);

                bool _readm4extremeNeutralFormat(int, const char *);
		bool _readm4extremeNeutralFull(int, const char *);

                void _readTable(const map<unsigned int, point_type> &,
				const map<unsigned int, vector<unsigned int> > &,
				map<int, Cell *> &,
				map<int, Cell *> &,
				map<Cell *, int> &,
				map<Cell *, int> &);
		void _readTable_general(const map<unsigned int, point_type> &,
					const map<unsigned int, vector<unsigned int> > &,
					map<int, Cell *> &,
					map<int, Cell *> &,
					map<Cell *, int> &,
					map<Cell *, int> &);

		void _getVertices(Geometry::Cell* e, 
				  set<Geometry::Cell*> & vs);

                GeneralGeometry(const GeneralGeometry &); // Not implemented
                GeneralGeometry & operator =(const GeneralGeometry &); // Not implemented
            };

        }

    }

}

#endif // !defined(GEOMETRY_MESH_GENERALGEOMETRY_H__INCLUDED_)
