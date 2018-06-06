#pragma once
#pragma execution_character_set("utf-8")

#include <osgViewer/Viewer>
#include <osg/Node>
#include <osg/Geometry>
#include <osg/Drawable>
#include <osgViewer/ViewerEventHandlers>
#include <osgGA/StandardManipulator>
#include <osgGA/StateSetManipulator>

#ifdef _DEBUG

#pragma  comment(lib,"osgViewerd.lib")
#pragma comment(lib,"osgGAd.lib")
#pragma  comment(lib,"osgd.lib")

#else
#pragma  comment(lib,"osgViewer.lib")
#pragma  comment(lib,"osg.lib")
#endif // DEBUG

#include "Tesselation.h"

int main(){


	osg::ref_ptr<osgViewer::Viewer> viewer = new osgViewer::Viewer;

	viewer->addEventHandler(new osgViewer::StatsHandler);

	viewer->addEventHandler(new osgGA::StateSetManipulator(viewer->getCamera()->getOrCreateStateSet()));

	osg::ref_ptr<osg::Group> gp = new osg::Group;

	osg::ref_ptr<osg::Geode> geode = new osg::Geode;

	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;

	viewer->setSceneData(gp);
	gp->addChild(geode);
	geode->addDrawable(geom);

	Tesselation tel;

	std::vector<D2Point> pts;
	std::vector<D2Point> calcpts;
	std::vector<Triangle> triangles;
	
	TesLoop outloop;
	std::vector<TesLoop> innerloops;
	{
		TesEdge a; a.pos = (D2Point(-0.96589, 1.78663));
		{
			//a.ctrlpts.push_back(D2Point(-0.83285, 1.78663));
			//a.ctrlpts.push_back(D2Point(-0.69981, 1.78663));
			//a.ctrlpts.push_back(D2Point(-0.56677, 1.78663));
			//a.ctrlpts.push_back(D2Point(-0.43373, 1.78663));
			//a.ctrlpts.push_back(D2Point(-0.30069, 1.78663));
			//a.ctrlpts.push_back(D2Point(-0.16765, 1.78663));
			//a.ctrlpts.push_back(D2Point(-0.03461, 1.78663));
			//a.ctrlpts.push_back(D2Point(0.09842, 1.78663));
		}

		TesEdge b; b.pos=(D2Point(0.23146, 1.78663));
		TesEdge c; c.pos=(D2Point(0.23146, 0.90254));
		TesEdge d; d.pos=(D2Point(-0.96589, 0.90254));

		outloop.edges.push_back(a);
		outloop.edges.push_back(b);
		outloop.edges.push_back(c);
		outloop.edges.push_back(d);
	}

	{
		TesLoop innerloop1;
		TesEdge a; a.pos=(D2Point(-0.83362, 1.65088));
		TesEdge b; b.pos=(D2Point(-0.62826, 1.65088));
		TesEdge c; c.pos=(D2Point(-0.62826, 1.44552));
		TesEdge d; d.pos=(D2Point(-0.83362, 1.44552));

		innerloop1.edges.push_back(a);
		innerloop1.edges.push_back(b);
		innerloop1.edges.push_back(c);
		innerloop1.edges.push_back(d);

		innerloops.push_back(innerloop1);
	}
	{
		TesLoop innerloop2;
		TesEdge a; a.pos=(D2Point(-0.29411, 1.62652));
		TesEdge b; b.pos=(D2Point(0.04002, 1.62652));
		TesEdge c; c.pos=(D2Point(0.04002, 1.36547));
		TesEdge d; d.pos=(D2Point(-0.29411, 1.36547));

		innerloop2.edges.push_back(a);
		innerloop2.edges.push_back(b);
		innerloop2.edges.push_back(c);
		innerloop2.edges.push_back(d);

		innerloops.push_back(innerloop2);

	}
	
	tel.LoadData(outloop, innerloops,0.03,0.03, calcpts, triangles);

	osg::ref_ptr<osg::Vec3dArray> arry = new osg::Vec3dArray;

	int size = calcpts.size();

	for (auto var : calcpts){
		arry->push_back(osg::Vec3d(var.x, var.y, 0));
	}

	osg::ref_ptr<osg::Vec4dArray> acolor_ary = new osg::Vec4dArray;

	unsigned int begin, end;
	unsigned int edge_id = 3;
	tel.LoadEdgeIndex(edge_id, begin, end);
	for (int i = 0; i < size; i++){
		if (i <= 3)
			acolor_ary->push_back(osg::Vec4d(1.0, 0.0, 0.0, 1.0));
		else if (i > 3 && i <= 7)
			acolor_ary->push_back(osg::Vec4d(1.0, 0.0, 1.0, 1.0));
		else if (i > 7 && i <= 11)
			acolor_ary->push_back(osg::Vec4d(1.0, 1.0, 0.0, 1.0));
		else if (i >= begin && i <= end)
			acolor_ary->push_back(osg::Vec4d(0.5,1.0,0.3,1.0));
		else
			acolor_ary->push_back(osg::Vec4d(1.0,1.0,1.0,1.0));
	}
	

	geom->setVertexArray(arry);

	geom->setColorArray(acolor_ary);

	geom->setColorBinding(osg::Geometry::BIND_PER_VERTEX);

	osg::ref_ptr<osg::DrawElements> ps = NULL;
	size_t  vertsize = arry->size();
	if (vertsize < 0xff){
		ps = new osg::DrawElementsUByte;
	}
	else if (vertsize < 0xffff){
		ps = new osg::DrawElementsUShort;
	}
	else{
		ps = new osg::DrawElementsUInt;
	}

	for (auto var : triangles){
		ps->addElement(var.a);
		ps->addElement(var.b);
		ps->addElement(var.c);
	}

	ps->setMode(GL_TRIANGLES);
	geom->addPrimitiveSet(ps);


	
	return viewer->run();
}