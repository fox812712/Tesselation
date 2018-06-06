/*
DelFEM (Finite Element Analysis)
Copyright (C) 2009  Nobuyuki Umetani    n.umetani@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*! @file
@brief ２次元メッシュのUtility関数群
@author Nobuyuki Umetani
*/


#if !defined(MESH_KERNEL_2D_H)
#define MESH_KERNEL_2D_H

#if defined(__VISUALC__)
	#pragma warning( disable : 4786 )
#endif

#include <vector>

#include "vector2d.h"

namespace Msh
{

/*! 
@addtogroup Msh2D
*/
//! @{

////////////////////////////////////////////////
// ここから本当は共有させたい

const double MIN_TRI_AREA = 1.0e-10;

const unsigned int nNoEd  = 2;	//!< 線分にいくら点があるか

const unsigned int nNoTri = 3;	//!< ３角形にいくら頂点があるか
const unsigned int nEdTri = 3;	//!< ３角形にいくら辺があるか
//! ３角形の各辺の頂点番号
const unsigned int noelTriEdge[nEdTri][nNoEd] = {	
	{ 1, 2 },
	{ 2, 0 },
	{ 0, 1 },
};
//! ３角形の隣接関係
const unsigned int relTriTri[3][3] = {
	{ 0, 2, 1 }, //  0
	{ 2, 1, 0 }, //  1 
	{ 1, 0, 2 }, //  2
};

const unsigned int nNoQuad = 4;	//!< ４角形にいくら頂点があるか
const unsigned int nEdQuad = 4;	//!< ４角形にいくら辺があるか
//! ４角形の各辺の頂点番号
const unsigned int noelQuadEdge[nEdQuad][nNoEd] = {	
	{ 0, 1 },
	{ 1, 2 },
	{ 2, 3 },
	{ 3, 0 }
};
//! ４角形の隣接関係
static const unsigned int relQuadQuad[nNoQuad][nNoQuad] = {
	{ 0, 3, 2, 1 }, //  0
	{ 1, 0, 3, 2 }, //  1
	{ 2, 1, 0, 3 }, //  2
	{ 3, 2, 1, 0 }, //  3
};


struct SBar{
	unsigned int v[2];	//!< 頂点番号
	////////////////
    // 左側が０右側が1
	unsigned int s2[2];
	unsigned int r2[2];
};


//! ２次元３角形要素構造体
struct STri2D{
	unsigned int v[3];	//!< 頂点Index
    int g2[3];			//!< 隣接する要素配列ID(-1:隣接要素なし、-2:自分の要素配列に隣接)
	unsigned int s2[3];	//!< 隣接要素Index
    unsigned int r2[3];	//!< 隣接関係
};


class CPoint2D{
public:
	CPoint2D(){}
	CPoint2D( const CPoint2D& rhs )
		: p(rhs.p), e(rhs.e), d(rhs.d){}
	CPoint2D(double x, double y, int ielem, unsigned int idir)
		: p(x,y), e(ielem), d(idir){}
public:
    Com::CVector2D p;   //!< 点の座標
    int e;              //!< 点を囲む要素のうち一つの番号(孤立している点なら-1が入る)
    unsigned int d;     //!< 点を囲む要素eの要素内節点番号
};

////////////////////////////////////////////////////////////////

bool MakePointSurTri(const std::vector<Msh::STri2D>& aTri, const unsigned int npoin, 
		unsigned int* const elsup_ind, unsigned int& nelsup, unsigned int*& elsup );


bool MakePointSurBar( const std::vector<Msh::SBar>& aBar, const unsigned int npoin, 
		unsigned int* const elsup_ind, unsigned int& nelsup, unsigned int*& elsup );

bool MakeInnerRelationTri( std::vector<Msh::STri2D>& aTri, const unsigned int npoin, 
	const unsigned int* elsup_ind, const unsigned int nelsup, const unsigned int* elsup);

bool MakeOuterBoundTri( const std::vector<Msh::STri2D>& aTri, std::vector<Msh::SBar>& aBar );

bool CheckTri( const std::vector<Msh::STri2D>& aTri );

bool CheckTri( const std::vector<Msh::CPoint2D>& po, const std::vector<Msh::STri2D>& tri );

bool InsertPoint_Elem( const unsigned int ipo_ins, const unsigned int itri_ins,
					  std::vector<Msh::CPoint2D>& po, std::vector<Msh::STri2D>& tri );

bool InsertPoint_ElemEdge( const unsigned int ipo_ins,
    const unsigned int itri_ins, const unsigned int ied_ins,
	std::vector<Msh::CPoint2D>& po, std::vector<Msh::STri2D>& tri );


bool DelaunayAroundPoint( unsigned int ipo0,
						 std::vector<Msh::CPoint2D>& po, std::vector<Msh::STri2D>& tri );

bool DelaunayAroundPoint( unsigned int itri, unsigned int inotri, 
						 std::vector<Com::CVector2D>& aVec, std::vector<Msh::STri2D>& tri,
						 unsigned int& num_flip);

bool FlipEdge( unsigned int itri0, unsigned int ied0,  std::vector<Msh::STri2D>& aTri );

bool FlipEdge( unsigned int itri0, unsigned int ied0,
			  std::vector<Msh::CPoint2D>& aPo, std::vector<Msh::STri2D>& aTri );

bool FindEdge( const unsigned int& ipo0, const unsigned int& ipo1,
	unsigned int& itri0, unsigned int& inotri0, unsigned int& inotri1,
	const std::vector<Msh::CPoint2D>& po, const std::vector<Msh::STri2D>& tri );

bool FindEdgePoint_AcrossEdge( const unsigned int& ipo0, const unsigned int& ipo1,
	unsigned int& itri0, unsigned int& inotri0, unsigned int& inotri1, double& ratio,
	std::vector<Msh::CPoint2D>& po, std::vector<Msh::STri2D>& tri );


bool Tesselate_GiftWrapping(
	unsigned int ipo0, unsigned int ipo1,
	std::vector< std::pair<unsigned int,unsigned int> >& aTriEd,
	const std::vector<Com::CVector2D>& aPo, std::vector<Msh::STri2D>& aTri);


bool DeletePointFromMesh(unsigned int ipoin, std::vector<Msh::CPoint2D>& aPo, std::vector<Msh::STri2D>& aTri);

bool DeleteTri(unsigned int itri0, std::vector<Msh::CPoint2D>& aPo, std::vector<Msh::STri2D>& aTri);


void LaplacianSmoothing( std::vector<Msh::CPoint2D>& aPo, const std::vector<Msh::STri2D>& aTri,
	const std::vector<unsigned int>& aflg_isnt_move );

void LaplaceDelaunaySmoothing( std::vector<Msh::CPoint2D>& aPo, std::vector<Msh::STri2D>& aTri,
	const std::vector<unsigned int>& aflg_isnt_move );

void PliantBossenHeckbertSmoothing( double elen, std::vector<Msh::CPoint2D>& aPo, std::vector<Msh::STri2D>& aTri );

void ColorCodeBarAry( const std::vector<SBar>& aBar, const std::vector<Com::CVector2D>& aVec, 
					 std::vector< std::vector<unsigned int> >& aIndBarAry );

} 



#endif // MSH_H
