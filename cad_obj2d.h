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
@brief interface of 2D cad class (Cad::CCadObj2Dm)
@author Nobuyuki Umetani
*/

#if !defined(CAD_OBJ_2D_H)
#define CAD_OBJ_2D_H

#include <vector>

#include "vector2d.h"
#include "cad2d_interface.h"
#include "objset.h"
#include "brep2d.h"
#include "cad_elem2d.h"

namespace Cad{

class CVertex2D;
class CLoop2D;
class CEdge2D;
class CTopology;


class CCadObj2D : public Cad::ICad2D_Msh
{
public:
    class CResAddVertex{
    public:
        CResAddVertex(){
            id_v_add = 0;
            id_e_add = 0;
        }
    public:
        unsigned int id_v_add;
        unsigned int id_e_add;
    };
    class CResAddPolygon{
    public:
        CResAddPolygon(){
            id_l_add = 0;
        }
        CResAddPolygon(const CResAddPolygon& lhs){
            this->id_l_add = lhs.id_l_add;
            this->aIdV = lhs.aIdV;
            this->aIdE = lhs.aIdE;
        }
    public:
        unsigned int id_l_add;
        std::vector<unsigned int> aIdV;
        std::vector<unsigned int> aIdE;
    };


    CCadObj2D();
    virtual ~CCadObj2D();
    void Clear();

    virtual std::auto_ptr<Cad::IItrLoop> GetPtrItrLoop(unsigned int id_l) const{
        return std::auto_ptr<Cad::IItrLoop>(new CBRepSurface::CItrLoop(m_BRep, id_l));	// instance
    }
    virtual bool IsElemID(Cad::CAD_ELEM_TYPE, unsigned int id) const;
    virtual const std::vector<unsigned int> GetAryElemID(Cad::CAD_ELEM_TYPE itype) const;
    virtual bool GetIdVertex_Edge(unsigned int &id_v_s, unsigned int& id_v_e, unsigned int id_e) const;
    unsigned int GetIdVertex_Edge(unsigned int id_e, bool is_s) const;
    virtual bool GetIdLoop_Edge(unsigned int &id_l_l, unsigned int& id_l_r, unsigned int id_e) const;
    CBRepSurface::CItrVertex GetItrVertex(unsigned int id_v) const { return CBRepSurface::CItrVertex(m_BRep, id_v); }
    CBRepSurface::CItrLoop GetItrLoop(unsigned int id_l) const { return CBRepSurface::CItrLoop(m_BRep, id_l); }

    virtual int GetLayer(Cad::CAD_ELEM_TYPE, unsigned int id) const;

    double GetMinClearance() const { return min_clearance; }

    double SignedDistPointLoop(unsigned int id_l1, const Com::CVector2D& point, unsigned int id_v_ignore = 0) const;
    virtual bool SetColor_Loop(unsigned int id_l, const double color[3]);

    const CEdge2D& GetEdge(unsigned int id_e) const;


    virtual int GetEdgeCurveType(const unsigned int& id_e) const;

	virtual void Set_Edge_InsertIndex(const unsigned int& id_e, const unsigned& id_edge_s, const unsigned& id_edge_e);

    virtual bool GetCurveAsPolyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo, double elen = -1) const;

    virtual bool GetCurveAsPolyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo,
        unsigned int ndiv, const Com::CVector2D& ps, const Com::CVector2D& pe) const;

    Com::CVector2D GetNearestPoint(unsigned int id_e, const Com::CVector2D& po_in) const;


    virtual Com::CVector2D GetVertexCoord(unsigned int id_v) const;

    CResAddPolygon AddPolygon(const std::vector<Com::CVector2D>& vec_ary, unsigned int id_l = 0);

    CResAddVertex AddVertex(Cad::CAD_ELEM_TYPE itype, unsigned int id_elem, const Com::CVector2D& vec);

    bool RemoveElement(Cad::CAD_ELEM_TYPE itype, unsigned int id_elem);

    CBRepSurface::CResConnectVertex ConnectVertex(CEdge2D edge);
    CBRepSurface::CResConnectVertex ConnectVertex_Line(unsigned int id_v1, unsigned int id_v2){
        Cad::CEdge2D e(id_v1, id_v2);
        return this->ConnectVertex(e);
    }
	bool DrawSpline(unsigned int id_e, const std::vector<double>& ctrls);
protected:

    CEdge2D& GetEdgeRef(unsigned int id_e);
    int AssertValid() const;
    bool CheckIsPointInside_ItrLoop(CBRepSurface::CItrLoop& itrl, const Com::CVector2D& p1) const;
    double DistPointItrLoop(CBRepSurface::CItrLoop& itrl, const Com::CVector2D& point) const;

    unsigned int CheckInOut_ItrLoopPoint_ItrLoop(CBRepSurface::CItrLoop& itrl1, CBRepSurface::CItrLoop& itrl2) const;

    int CheckLoop(unsigned int id_l) const;
    CBRepSurface::CItrVertex FindCorner_HalfLine(unsigned int id_v, const Com::CVector2D& point) const;
    bool CheckIntersection_Loop(unsigned int id_l = 0) const;
    bool CheckIntersection_EdgeAgainstLoop(const CEdge2D& edge, unsigned int id_l = 0) const;
    double GetArea_ItrLoop(CBRepSurface::CItrLoop& itrl) const;
protected:
    ////////////////
    Com::CObjSet<CLoop2D>   m_LoopSet;
    Com::CObjSet<CEdge2D>   m_EdgeSet;
    Com::CObjSet<CVertex2D> m_VertexSet;
    ////////////////
    CBRepSurface m_BRep;
    double min_clearance;
};

}

#endif
