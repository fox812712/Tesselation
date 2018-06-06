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
@brief interface of the class (Cad::CBRep) wich represents topology with B-Rep data strcutre
@author Nobuyuki Umetani
*/

#if !defined(B_REP_H)
#define B_REP_H

#ifdef __VISUALC__
#pragma warning( disable : 4786 )
#endif

#include "objset_cad.h"

namespace Cad{

class CUseLoop{
public:
    CUseLoop(const CUseLoop& rhs)
        : id(rhs.id), id_l(rhs.id_l), id_he(rhs.id_he), id_ul_c(rhs.id_ul_c), id_ul_p(rhs.id_ul_p){}
    CUseLoop(const unsigned int id,
        const unsigned int id_he, const unsigned int id_ul_c, const unsigned int id_ul_p)
        : id(id), id_l(0), id_he(id_he), id_ul_c(id_ul_c), id_ul_p(id_ul_p){}
public:
	unsigned int id;
	unsigned int id_l;
    unsigned int id_he; 
    unsigned int id_ul_c;	
    unsigned int id_ul_p;	
};

class CHalfEdge{
public:
    CHalfEdge(const CHalfEdge& rhs)
        : id(rhs.id),
        id_uv(rhs.id_uv),
        id_he_f(rhs.id_he_f), id_he_b(rhs.id_he_b), id_he_o(rhs.id_he_o),
        id_ul(rhs.id_ul),
        id_e(rhs.id_e), is_same_dir(rhs.is_same_dir){}
    CHalfEdge(const unsigned int id,
        const unsigned int id_uv,
        const unsigned int id_he_f, const unsigned int id_he_b, const unsigned int id_he_o,
        const unsigned int id_ul)
        : id(id),
        id_uv(id_uv),
        id_he_f(id_he_f), id_he_b(id_he_b), id_he_o(id_he_o),
        id_ul(id_ul),
        id_e(0), is_same_dir(true){}
public:
    unsigned int id;        
    unsigned int id_uv;    
    unsigned int id_he_f;   
    unsigned int id_he_b;   
    unsigned int id_he_o;	  
    unsigned int id_ul;     

    unsigned int id_e;     
    bool is_same_dir;       
};


class CUseVertex{
public:
    CUseVertex(const unsigned int id, const unsigned int id_he)
        : id(id), id_he(id_he), id_v(0){}
    CUseVertex(const CUseVertex& rhs)
        : id(rhs.id), id_he(rhs.id_he), id_v(rhs.id_v){}
public:
    unsigned int id;    
    unsigned int id_he; 
    unsigned int id_v;  
};

class CBRep{
public:
    void Clear();

    bool IsID_UseLoop(unsigned int id_ul) const { return m_UseLoopSet.IsObjID(id_ul); }
    bool IsID_HalfEdge(unsigned int id_he) const { return m_HalfEdgeSet.IsObjID(id_he); }
    bool IsID_UseVertex(unsigned int id_uv) const { return m_UseVertexSet.IsObjID(id_uv); }
    std::vector<unsigned int> GetAry_UseVertexID() const { return m_UseVertexSet.GetAry_ObjID(); }

    const CUseLoop& GetUseLoop(unsigned int id_ul) const;
    const CUseVertex& GetUseVertex(unsigned int id_uv) const;
    const CHalfEdge& GetHalfEdge(unsigned int id_he) const;

    bool SetLoopIDtoUseLoop(unsigned int id_ul, unsigned int id_l);
    bool SetVertexIDtoUseVertex(unsigned int id_uv, unsigned int id_v);
    bool SetEdgeIDtoHalfEdge(unsigned int id_he, unsigned int id_e, bool is_same_dir);

    int AssertValid_Use() const;
    std::vector<unsigned int> FindHalfEdge_Edge(const unsigned int& id_e) const;
    std::vector<unsigned int> FindHalfEdge_Vertex(const unsigned int& id_v) const;

    bool MEVVL(unsigned int& id_he_add1, unsigned int& id_he_add2,
        unsigned int& id_uv_add1, unsigned int& id_uv_add2, unsigned int& id_ul_add);
    bool MEL(unsigned int& id_he_add1, unsigned int& id_he_add2, unsigned int& id_ul_add,
        const unsigned int id_he1, const unsigned int id_he2);
    bool KEL(const unsigned int id_he_rem1);
    bool MEV(unsigned int& id_he_add1, unsigned int& id_he_add2, unsigned int& id_uv_add,
        const unsigned int id_he);
    bool KVE(unsigned int& id_he_rem1);
    bool MVE(unsigned int& id_he_add1, unsigned int& id_he_add2, unsigned int& id_uv_add,
        const unsigned int id_he);
    bool MEKL(unsigned int& id_he_add1, unsigned int& id_he_add2,
        const unsigned int id_he1, const unsigned int id_he2);
    bool KEML(unsigned int& id_ul_add1,
        const unsigned int& id_he1);
    bool MEKL_OneFloatingVertex(unsigned int& id_he_add1,
        const unsigned int id_he1, const unsigned int id_he2);
    bool MEKL_TwoFloatingVertex(const unsigned int id_he1, const unsigned int id_he2);
    bool KEML_OneFloatingVertex(unsigned int& id_ul_add,
        const unsigned int id_he1);
    bool KEML_TwoFloatingVertex(unsigned int& id_ul_add,
        const unsigned int id_he1);
    bool MVEL(unsigned int& id_uv_add, unsigned int& id_he_add, unsigned int& id_ul_add,
        const unsigned int id_ul_p);
    bool KVEL(const unsigned int id_uv_rem);
    bool SwapUseLoop(unsigned int id_ul1, unsigned int id_ul2);
    bool MoveUseLoop(unsigned int id_ul1, unsigned int id_ul2);
private:
    bool SwapUseLoop_CP_SameLoop(unsigned int id_ul1, unsigned int id_ul2);
    bool SwapUseLoop_CP_DifferentLoop(unsigned int id_ul1, unsigned int id_ul2);
public:
    CObjSetCad<CUseLoop>   m_UseLoopSet;  
    CObjSetCad<CHalfEdge>  m_HalfEdgeSet;  
    CObjSetCad<CUseVertex> m_UseVertexSet; 
};
}

#endif

