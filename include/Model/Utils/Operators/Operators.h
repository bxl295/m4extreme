// Operators.h: interface for the Operators class.
// Copyright (c) 2017-2018 Extreme Computation Technology and Solutions, LLC 
// All rights reserved
// see file License.txt for license details
//////////////////////////////////////////////////////////////////////

#if !defined(MODEL_UTILS_OPERATORS_H__INCLUDED_)
#define MODEL_UTILS_OPERATORS_H__INCLUDED_

#pragma once

#include "../../Model.h"

void operator += (
	set<Set::Manifold::Point *> &,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &);
void operator -= (
	set<Set::Manifold::Point *> &,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &);
void operator *= (
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
	const double &);
void operator /= (
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
	const double &);
map<Set::Manifold::Point *, Set::VectorSpace::Vector> operator + (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &);
map<Set::Manifold::Point *, Set::VectorSpace::Vector> operator - (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &);
map<Set::Manifold::Point *, Set::VectorSpace::Vector> operator - (
	const map<Set::Manifold::Point *, Set::Euclidean::Orthonormal::Point> &,
	const map<Set::Manifold::Point *, Set::Euclidean::Orthonormal::Point> &);
map<Set::Manifold::Point *, Set::VectorSpace::Vector> operator - (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &);
map<Set::Manifold::Point *, Set::VectorSpace::Vector> operator * (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
	const double &);
map<Set::Manifold::Point *, Set::VectorSpace::Vector> operator * (
	const double &,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &);
map<Set::Manifold::Point *, Set::VectorSpace::Vector> operator / (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
	const double &);
double operator * (
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &,
	const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &);
void Null(
	const set<Set::Manifold::Point *> &,
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> &);
void Random(
	const set<Set::Manifold::Point *> &,
	map<Set::Manifold::Point *, Set::VectorSpace::Vector> &);
double Norm(const map<Set::Manifold::Point *, Set::VectorSpace::Vector> &);

#endif // !defined(MODEL_UTILS_OPERATORS_H__INCLUDED_)
