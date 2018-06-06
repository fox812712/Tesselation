#pragma once
#pragma execution_character_set("utf-8")

#include <vector>


#include "commonstruct.h"
#include "mesher2d.h"
#include "cad2d_interface.h"
#include "cad_obj2d.h"
//#include "cad_obj2d_move.h"


using namespace Cad;
using namespace Com;

using namespace Msh;

class Tesselation{
public:

	Tesselation();
	~Tesselation();
	void LoadData(const TesLoop& outerloop, const std::vector<TesLoop>& innerloop,
		               const double& edgelength,const double& arealength,
		               std::vector<D2Point>& d2points, std::vector<Triangle>& triangles);

	void LoadEdgeIndex(unsigned int& id_edge, unsigned int& id_e_s, unsigned int& id_e_e);


protected:
	unsigned int AddOutlineloop(Cad::CCadObj2D& cad2d, Msh::CMesher2D& msh, const TesLoop& loop);
	void AddInnerloop(Cad::CCadObj2D& cad2d, Msh::CMesher2D& msh, const unsigned int& id_loop, const TesLoop& loop);

private:
	Cad::CCadObj2D cad_2d;
	Msh::CMesher2D msh_2d;
};
