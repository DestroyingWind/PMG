#pragma once

using namespace std;

struct point;
struct face;
struct region;

class PMG
{
public:
	PMG();
	~PMG();

	bool readOM();
	bool divideM();
	bool writeNM();

	string infilename="box.poly";

	int nx=1,ny=1,nz=1;
private:
	vector<point> PMG_point;
	vector<face> PMG_face;
	vector<point> PMG_hole;
	vector<region> PMG_region;
};

