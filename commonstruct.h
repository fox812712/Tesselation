#pragma once
#pragma execution_character_set("utf-8")


#ifndef COMMONSTRUCT_H
#define COMMONSTRUCT_H



typedef struct _D2Point{

	double x;
	double y;
	_D2Point(){ x = 0.0; y = 0.0; }
	_D2Point(double _x, double _y){ x = _x; y = _y; }

	_D2Point& operator =(const _D2Point& rhs){
		x = rhs.x;
		y = rhs.y;
		return *this;
	}

}D2Point;

typedef struct _TesEdge{
	D2Point pos;
	std::vector<D2Point> ctrlpts;
	_TesEdge(){}
}TesEdge;

typedef struct _TesLoop{

	std::vector<TesEdge> edges;

}TesLoop;

typedef struct _Triangle{
	int a;
	int b;
	int c;

	_Triangle(int _a, int _b, int _c){ a = _a; b = _b; c = _c; }

}Triangle;

#endif
