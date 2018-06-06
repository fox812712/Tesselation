#pragma once
#pragma execution_character_set("utf-8")

#include "Tesselation.h"

Tesselation::Tesselation(){}

Tesselation::~Tesselation(){}

void Tesselation::LoadData(const TesLoop& outerloop, const std::vector<TesLoop>& innerloop,
	const double& edgelength,const double& arealength,
	std::vector<D2Point>& d2points, std::vector<Triangle>& triangles)
{
	cad_2d.Clear();
	msh_2d.Clear();

	msh_2d.set_Elementlength(edgelength, arealength);

	unsigned int id_loop = AddOutlineloop(cad_2d, msh_2d, outerloop);

	for (auto var : innerloop){
		AddInnerloop(cad_2d, msh_2d, id_loop, var);
	}

	msh_2d.Get_TesslationResult(d2points, triangles);
}


void Tesselation::LoadEdgeIndex(unsigned int& id_edge, unsigned int& id_e_s, unsigned int& id_e_e){
	const std::vector<CBarAry>& barsh = msh_2d.GetBarArySet();
	CBarAry b;
	for (auto var : barsh){
		if (var.id_e_cad == id_edge){ b = var; break; }
	}
	id_e_s = b.m_aBar.front().v[1];
	id_e_e = b.m_aBar.back().v[0];
	int tr = 0;
	
}

unsigned int Tesselation::AddOutlineloop(Cad::CCadObj2D& cad2d, Msh::CMesher2D& msh, const TesLoop& loop){
	std::vector<Com::CVector2D> vec_ary;

	for (auto var : loop.edges){
		vec_ary.push_back(Com::CVector2D(var.pos.x,var.pos.y));
		for (auto t : var.ctrlpts){
			vec_ary.push_back(Com::CVector2D(t.x,t.y));
		}
	}
	Cad::CCadObj2D::CResAddPolygon res = cad2d.AddPolygon(vec_ary);
	
	if (res.aIdE.size() == 0 || res.aIdV.size() == 0)return 0;

	msh.AddIdLCad_CutMesh(res.id_l_add);

	msh.Meshing(cad2d);

	return res.id_l_add;
}

void Tesselation::AddInnerloop(Cad::CCadObj2D& cad2d, Msh::CMesher2D& msh, const unsigned int& id_loop, const TesLoop& loop)
{

	std::vector<Com::CVector2D> vec_ary;

	for (auto var : loop.edges){
		vec_ary.push_back(Com::CVector2D(var.pos.x, var.pos.y));
	}

	cad2d.AddPolygon(vec_ary, id_loop);
	msh.Meshing(cad2d);
}