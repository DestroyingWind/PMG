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

	string infilename="example.poly";

	int dx=1,dy=1,dz=1;
private:
	vector<point> PMG_point;
	vector<face> PMG_face;
	vector<point> PMG_hole;
	vector<region> PMG_region;
};

