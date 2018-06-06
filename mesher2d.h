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
@brief interface of 2D mesher class (Msh::CMesher2D)
@author Nobuyuki Umetani
*/

#if !defined(MSHER_2D_H)
#define MSHER_2D_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <vector>
#include <set>
#include <map>

#include "commonstruct.h"
#include "vector2d.h"
#include "cad2d_interface.h"
#include "meshkernel2d.h"

////////////////////////////////////////////////

namespace Msh{

	enum MSH_TYPE{
		VERTEX,	//!< 頂点要素
		BAR,	//!< 線分要素
		TRI,	//!< 三角形要素
		QUAD,	//!< ４角形要素
		TET,	//!< 四面体要素
		HEX		//!< 六面体要素
	};

struct SVertex{
public:
    SVertex() : id(0), id_v_cad(0), ilayer(0){}
public:
    unsigned int id;	//!< ID
    unsigned int id_v_cad;	//!< vertex id in CAD（0 if not related to CAD）
    int ilayer;
    unsigned int v;	//!< index of node
};

//! array of line element
class CBarAry{
public:
    CBarAry() : id(0), id_e_cad(0), ilayer(0){}
public:
    unsigned int id;	//!< ID
    unsigned int id_e_cad;	//!< CADの辺ID（CADに関連されてなければ０）
    unsigned int id_se[2];
    unsigned int id_lr[2];
    int ilayer;
    std::vector<SBar> m_aBar;	//!< array of line element
};

class CTriAry2D{
public:
    CTriAry2D() : id(0), id_l_cad(0), ilayer(0){}
public:
    unsigned int id;	//!< ID
    unsigned int id_l_cad;	//!< CADの面ID（CADに関連されてなければ０）
    int ilayer;
    std::vector<STri2D> m_aTri;	//!< array of 2d triangle element
};


class CMesher2D
{
public:
	CMesher2D();

	virtual ~CMesher2D();

	virtual void AddIdLCad_CutMesh(unsigned int id_l_cad);

	virtual bool IsIdLCad_CutMesh(unsigned int id_l_cad) const;

	virtual void RemoveIdLCad_CutMesh(unsigned int id_l_cad);

	virtual std::vector<unsigned int> GetIdLCad_CutMesh();

	void Get_TesslationResult(std::vector<D2Point>& points, std::vector<Triangle>& triangles);

	void set_Elementlength(const double& edgelength, const double& arealength);

	virtual bool Meshing(const Cad::ICad2D_Msh& cad_2d);

    virtual void Clear();
    unsigned int GetElemID_FromCadID(unsigned int id_cad, Cad::CAD_ELEM_TYPE type_cad) const;
    bool IsID(unsigned int id) const;


    const std::vector<CTriAry2D>& GetTriArySet() const { return m_aTriAry; }
    const std::vector<CBarAry>& GetBarArySet() const { return m_aBarAry; }
    const std::vector<SVertex>& GetVertexAry() const { return m_aVertex; }
    const std::vector<Com::CVector2D>& GetVectorAry() const { return aVec2D; }

protected:
	void ClearMeshData();

    bool Meshing_ElemLength(const Cad::ICad2D_Msh& cad_2d, const std::vector<unsigned int>& aIdLoop);
    unsigned int FindMaxID() const;
    int CheckMesh();	// 異常がなければ０を返す
    unsigned int GetFreeObjID();
    bool MakeMesh_Edge(const Cad::ICad2D_Msh& cad_2d, unsigned int id_e, const double len);
    bool MakeMesh_Loop(const Cad::ICad2D_Msh& cad_2d, unsigned int id_l, const double len);
    bool Tesselate_Loop(const Cad::ICad2D_Msh& cad_2d, const unsigned int id_l);
    bool FindElemLocType_CadIDType(
        unsigned int& iloc, unsigned int& itype,
        unsigned int id_cad_part, Cad::CAD_ELEM_TYPE itype_cad);

private:
    std::set<unsigned int> setIdLCad_CutMesh;
protected:
    std::vector<int> m_ElemType;	// vertex(0) bar(1) tri(2) quad(3)	always valid (return -1 if no corresponding ID)
    std::vector<int> m_ElemLoc;		// index of elem_ary : always valid (return -1 if no corresponding ID)
    std::vector<SVertex> m_aVertex;		// type(0)
    std::vector<CBarAry> m_aBarAry;		// type(1)
    std::vector<CTriAry2D> m_aTriAry;	// type(2)
    std::vector<Com::CVector2D> aVec2D;
	double m_pedgelength;
	double m_parealength;
};

// @}
}

#endif
