// -*- C++ -*-

/*!
  \file geom/orq/Octree.h
  \brief A class for an octree in 3-D.
*/

#if !defined(__geom_orq_Octree_h__)
#define __geom_orq_Octree_h__

#include "ORQ.h"

#include "../../ads/functor/Dereference.h"

#include <string>

namespace geom {

//
//----------------------------OctreeNode--------------------------------
//

//! A node in an octree.
template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
class OctreeNode {
public:

   //
   // Public types.
   //

   //! A pointer to the record type.
   typedef _Record Record;
   //! The multi-key type.
   typedef _MultiKey MultiKey;
   //! The key type.
   typedef _Key Key;
   //! The multi-key accessor functor.
   typedef _MultiKeyAccessor MultiKeyAccessor;
   //! The output iterator for records.
   typedef _RecordOutputIterator RecordOutputIterator;

   //! A Cartesian point.
   typedef std::tr1::array<Key, 3> Point;
   //! Bounding box.
   typedef geom::BBox<Key, 3> BBox;

protected:

   //! The domain of this node.
   BBox _domain;

protected:

   //! The multi-key accessor.
   static MultiKeyAccessor _multiKeyAccessor;

public:

   //--------------------------------------------------------------------------
   //! \name Destructor.
   //@{

   //! Trivial destructor.
   virtual
   ~OctreeNode() {}

   //@}
   //--------------------------------------------------------------------------
   //! \name Accesors.
   //@{

   //! Return the domain.
   const BBox&
   getDomain() const {
      return _domain;
   }

   //@}
   //--------------------------------------------------------------------------
   //! \name Insert records.
   //@{

   //! Add a record. Return the OctreeNode containing the record.
   virtual
   OctreeNode*
   insert(Record record, std::size_t leafSize) = 0;

   //@}
   //--------------------------------------------------------------------------
   //! \name Window queries.
   //@{

   //! Get the records in the node and children.  Return the # of records.
   virtual
   std::size_t
   report(RecordOutputIterator iter) const = 0;

   //! Get the records in the window.  Return the # of records inside.
   virtual
   std::size_t
   computeWindowQuery(RecordOutputIterator iter, const BBox& window)
   const = 0;

   //! Get the records in the window.  Return the # of records inside.
   virtual
   std::size_t
   computeWindowQueryCheckDomain(RecordOutputIterator iter,
                                 const BBox& window) const = 0;

   //@}
   //--------------------------------------------------------------------------
   //! \name File I/O.
   //@{

   //! Print the octree node
   virtual
   void
   put(std::ostream& out) const = 0;

   //! Print the octree node as a tree.
   virtual
   void
   print(std::ostream& out, std::string tabbing) const = 0;

   //@}
   //--------------------------------------------------------------------------
   //! \name Memory usage.
   //@{

   //! Return the memory usage of this node and its children.
   virtual
   std::size_t
   getMemoryUsage() const = 0;

   //@}
   //--------------------------------------------------------------------------
   //! \name Validity.
   //@{

   // Check the validity of the node.
   virtual
   bool
   isValid() const = 0;

   //@}
};


//! Static member variable.
template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
_MultiKeyAccessor OctreeNode < _Record, _MultiKey, _Key, _MultiKeyAccessor,
                  _RecordOutputIterator >::_multiKeyAccessor;


//! Write to a file stream.
/*! \relates OctreeNode */
template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
std::ostream&
operator<<(std::ostream& out,
           const OctreeNode < _Record, _MultiKey, _Key, _MultiKeyAccessor,
           _RecordOutputIterator > & node) {
   node.put(out);
   return out;
}



//
//-------------------------OctreeBranch------------------------------
//

//
// Forward declarations
//

template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
class OctreeLeaf;

//! A branch in an octree.
template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
class OctreeBranch :
   public OctreeNode < _Record, _MultiKey, _Key, _MultiKeyAccessor,
      _RecordOutputIterator > {
private:

   typedef OctreeNode < _Record, _MultiKey, _Key, _MultiKeyAccessor,
           _RecordOutputIterator > Base;
   typedef OctreeNode < _Record, _MultiKey, _Key, _MultiKeyAccessor,
           _RecordOutputIterator > Node;
   typedef OctreeLeaf < _Record, _MultiKey, _Key, _MultiKeyAccessor,
           _RecordOutputIterator > Leaf;

public:

   //
   // Public types.
   //

   //! A pointer to the record type.
   typedef _Record Record;
   //! The multi-key type.
   typedef _MultiKey MultiKey;
   //! The key type.
   typedef _Key Key;
   //! The multi-key accessor functor.
   typedef _MultiKeyAccessor MultiKeyAccessor;
   //! The output iterator for records.
   typedef _RecordOutputIterator RecordOutputIterator;

   //! A Cartesian point.
   typedef typename Base::Point Point;
   //! Bounding box.
   typedef typename Base::BBox BBox;

private:

   //! The midpoint of the domain.
   Point _midpoint;
   //! The eight octants.
   Node* _octant[8];

private:

   //
   // Not implemented
   //

   //! Default constructor not implemented.
   OctreeBranch();

   //! Copy constructor not implemented
   OctreeBranch(const OctreeBranch&);

   //! Assignment operator not implemented
   OctreeBranch&
   operator=(const OctreeBranch&);

public:

   //--------------------------------------------------------------------------
   //! \name Constructors and destructor.
   //@{

   //! Construct from the Cartesian domain.
   OctreeBranch(const BBox& domain);

   //! Destructor.
   virtual
   ~OctreeBranch();

   //@}
   //--------------------------------------------------------------------------
   //! \name Insert records.
   //@{

   //! Add a record.
   Node*
   insert(Record record, std::size_t leafSize);

   //@}
   //--------------------------------------------------------------------------
   //! \name Window queries.
   //@{

   //! Get the records in the node and children.  Return the # of records.
   std::size_t
   report(RecordOutputIterator iter) const;

   //! Get the records in the window.  Return the # of records inside.
   std::size_t
   computeWindowQuery(RecordOutputIterator iter, const BBox& window) const;

   //! Get the records in the window.  Return the # of records inside.
   std::size_t
   computeWindowQueryCheckDomain(RecordOutputIterator iter,
                                 const BBox& window)
   const;

   //@}
   //--------------------------------------------------------------------------
   //! \name File I/O.
   //@{

   //! Print the octree.
   void
   put(std::ostream& out) const;

   //! Print the octree as a tree.
   void
   print(std::ostream& out, std::string tabbing) const;

   //@}
   //--------------------------------------------------------------------------
   //! \name Memory usage.
   //@{

   //! Return the memory usage of this branch and its children.
   std::size_t
   getMemoryUsage() const;

   //@}
   //--------------------------------------------------------------------------
   //! \name Validity check.
   //@{

   bool
   isValid() const;

   //@}

private:

   //! Return the octant index containing the record.
   std::size_t
   getOctantIndex(Record record) const;

   //! Return the domain for the given octant.
   BBox
   getOctantDomain(std::size_t index) const;

};


//
//----------------------------OctreeLeaf--------------------------------
//

//! A leaf in an octree.
template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
class OctreeLeaf :
   public OctreeNode < _Record, _MultiKey, _Key, _MultiKeyAccessor,
      _RecordOutputIterator > {
private:

   typedef OctreeNode < _Record, _MultiKey, _Key, _MultiKeyAccessor,
           _RecordOutputIterator > Base;
   typedef OctreeNode < _Record, _MultiKey, _Key, _MultiKeyAccessor,
           _RecordOutputIterator > Node;
   typedef OctreeBranch < _Record, _MultiKey, _Key, _MultiKeyAccessor,
           _RecordOutputIterator > Branch;

public:

   //
   // Public types.
   //

   //! A pointer to the record type.
   typedef _Record Record;
   //! The multi-key type.
   typedef _MultiKey MultiKey;
   //! The key type.
   typedef _Key Key;
   //! The multi-key accessor functor.
   typedef _MultiKeyAccessor MultiKeyAccessor;
   //! The output iterator for records.
   typedef _RecordOutputIterator RecordOutputIterator;

   //! A Cartesian point.
   typedef typename Base::Point Point;
   //! Bounding box.
   typedef typename Base::BBox BBox;


private:

   //
   // Private types.
   //

   typedef std::vector<Record> Container;
   typedef typename Container::iterator Iterator;
   typedef typename Container::const_iterator ConstIterator;

private:

   //! The records
   Container _records;

private:

   //
   // Not implemented
   //

   //! Default constructor not implemented.
   OctreeLeaf();

   //! Copy constructor not implemented
   OctreeLeaf(const OctreeLeaf&);

   //! Assignment operator not implemented
   OctreeLeaf&
   operator=(const OctreeLeaf&);

public:

   //--------------------------------------------------------------------------
   //! \name Constructors and destructor.
   //@{

   //! Construct from a Cartesian domain.
   OctreeLeaf(const BBox& domain) :
      _records() {
      Base::_domain = domain;
   }

   //! Trivual destructor.
   virtual
   ~OctreeLeaf() {}

   //@}
   //--------------------------------------------------------------------------
   //! \name Insert records.
   //@{

   //! Add a record.
   Node*
   insert(Record record, std::size_t leafSize);

   //@}
   //--------------------------------------------------------------------------
   //! \name Window queries.
   //@{

   //! Get the records in the node and children.  Return the # of records.
   std::size_t
   report(RecordOutputIterator iter) const;

   //! Get the records in the window.  Return the # of records inside.
   std::size_t
   computeWindowQuery(RecordOutputIterator iter, const BBox& window) const;

   //! Get the records in the window.  Return the # of records inside.
   std::size_t
   computeWindowQueryCheckDomain(RecordOutputIterator iter,
                                 const BBox& window) const;

   //@}
   //--------------------------------------------------------------------------
   //! \name File I/O.
   //@{

   //! Print the octree.
   void
   put(std::ostream& out) const;

   //! Print the octree in tree form.
   void
   print(std::ostream& out, std::string tabbing) const;

   //@}
   //--------------------------------------------------------------------------
   //! \name Memory usage.
   //@{

   //! Return the memory usage of this leaf.
   std::size_t
   getMemoryUsage() const {
      return (sizeof(OctreeLeaf) + _records.size() * sizeof(Record));
   }

   //@}
   //--------------------------------------------------------------------------
   //! \name Validity.
   //@{

   bool
   isValid() const {
      for (ConstIterator i = _records.begin(); i != _records.end(); ++i) {
         if (! Base::_domain.isIn(_multiKeyAccessor(*i))) {
            return false;
         }
      }
      return true;
   }

   //@}
};




//
//---------------------------Octree class---------------------------------
//

//! An octree in 3-D.
/*!
  An octree in 3-D.

  RecordOutputIterator is the typename that you will use in the
  computeWindowQuery() function.  You have to specify it as a template
  parameter because of the way KDTree is implemented.  C++ does
  not allow templated virtual functions.

  Use computeWindowQuery() if the number of records returned by a typical
  window query is small. Use computeWindowQueryCheckDomain only if the
  number of records returned is much larger than the leaf size.

  This implementation of Octree stores the domain information at each
  branch and leaf.  This choice increases the memory usage but avoids
  the computational cost of having to compute the domain as the
  window query progresses.
*/
template < typename _Record,
         typename _MultiKey =
         typename std::iterator_traits<_Record>::value_type,
         typename _Key = typename _MultiKey::value_type,
         typename _MultiKeyAccessor = ads::Dereference<_Record>,
         typename _RecordOutputIterator = _Record* >
class Octree :
   public ORQ<3, _Record, _MultiKey, _Key, _MultiKeyAccessor> {
private:

   //
   // Private types.
   //

   typedef ORQ<3, _Record, _MultiKey, _Key, _MultiKeyAccessor> Base;

public:

   //
   // Public types.
   //

   //! A pointer to the record type.
   typedef _Record Record;
   //! The multi-key type.
   typedef _MultiKey MultiKey;
   //! The key type.
   typedef _Key Key;
   //! The record output iterator.
   typedef _RecordOutputIterator RecordOutputIterator;

   //! A Cartesian point.
   typedef typename Base::Point Point;
   //! Bounding box.
   typedef typename Base::BBox BBox;

private:

   //
   // private typedefs
   //

   typedef OctreeNode < _Record, _MultiKey, _Key, _MultiKeyAccessor,
           _RecordOutputIterator > Node;
   typedef OctreeBranch < _Record, _MultiKey, _Key, _MultiKeyAccessor,
           _RecordOutputIterator > Branch;
   typedef OctreeLeaf < _Record, _MultiKey, _Key, _MultiKeyAccessor,
           _RecordOutputIterator > Leaf;

   //
   // Member data
   //

   //! The root of the octree.
   Node* _root;

   //! The leaf size.
   std::size_t _leafSize;

protected:

   using Base::_multiKeyAccessor;

private:

   //
   // Not implemented
   //

   //! Default constructor not implemented.
   Octree();

   //! Copy constructor not implemented
   Octree(const Octree&);

   //! Assignment operator not implemented
   Octree&
   operator=(const Octree&);

public:

   //--------------------------------------------------------------------------
   //! \name Constructors and destructor.
   //@{

   //! Construct from a Cartesian domain.
   /*!
     leafSize is the maximum number of records that are stored in a leaf.
     Choose this value to be about the number records that you expect a
     window query to return.  The default value is 8.
   */
   Octree(const BBox& domain, const std::size_t leafSize = 8) :
      Base(),
      _leafSize(leafSize) {
      _root = new Leaf(domain);
   }

   //! Construct from a Cartesian domain and a range of records.
   /*!
     leafSize is the maximum number of records that are stored in a leaf.
     Choose this value to be about the number records that you expect a
     window query to return.  The default value is 8.
   */
   template<class InputIterator>
   Octree(const BBox& domain,
          InputIterator first, InputIterator last,
          const std::size_t leafSize = 8) :
      Base(),
      _leafSize(leafSize) {
      _root = new Leaf(domain);
      insert(first, last);
   }

   //! Destructor.
   ~Octree() {
      delete _root;
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Accesors.
   // @{

   //
   // Inherited.
   //

   // Return the number of records.
   using Base::getSize;

   // Return true if the grid is empty.
   using Base::isEmpty;

   //
   // New.
   //

   //! Return the domain spanned by the records.
   const BBox&
   getDomain() const {
      return _root->getDomain();
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Insert records.
   // @{

   //! Add a single record.
   void
   insert(const Record record) {
      _root = _root->insert(record, _leafSize);
      incrementSize();
   }

   //! Add a range of records.
   template<typename InputIterator>
   void
   insert(InputIterator first, InputIterator last) {
      while (first != last) {
         insert(first);
         ++ first;
      }
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Window queries.
   // @{

   //! Get the records in the node and children.  Return the # of records.
   std::size_t
   report(RecordOutputIterator iter) const {
      return _root->report(iter);
   }

   //! Get the records in the window.  Return the # of records inside.
   std::size_t
   computeWindowQuery(RecordOutputIterator iter, const BBox& window) const {
      return _root->computeWindowQuery(iter, window);
   }

   //! Get the records in the window.  Return the # of records inside.
   /*!
     Use this version of window query only if the number of records
     returned is much larger than the leaf size.
   */
   std::size_t
   computeWindowQueryCheckDomain(RecordOutputIterator iter,
                                 const BBox& window) const {
      return _root->computeWindowQueryCheckDomain(iter, window);
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name File I/O.
   // @{

   //! Print the octree.
   void
   put(std::ostream& out) const;

   //! Print the structure of the octree.
   void
   print(std::ostream& out) const;

   // @}
   //--------------------------------------------------------------------------
   //! \name Memory usage.
   // @{

   //! Return the memory usage of the octree.
   std::size_t
   getMemoryUsage() const {
      return (sizeof(Octree) + _root->getMemoryUsage());
   }

   // @}
   //--------------------------------------------------------------------------
   //! \name Validity.
   // @{

   //! Check the validity of the octree.
   bool
   isValid() const {
      return _root->isValid();
   }

   // @}

protected:

   //! Increment the number of records.
   using Base::incrementSize;
};

//! Write to a file stream.
/*! \relates Octree */
template < typename _Record,
         typename _MultiKey,
         typename _Key,
         typename _MultiKeyAccessor,
         typename _RecordOutputIterator >
inline
std::ostream&
operator<<(std::ostream& out,
           const Octree<_Record, _MultiKey, _Key, _MultiKeyAccessor>& x) {
   x.put(out);
   return out;
}


} // namespace geom

#define __geom_orq_Octree_ipp__
#include "Octree.ipp"
#undef __geom_orq_Octree_ipp__

#endif
