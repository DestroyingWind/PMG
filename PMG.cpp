#include<vector>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<string>
#include "PMG.h"

using namespace std;

struct point
{
	double cood[3];
	int tag;
};
struct face
{
	vector<int> vertex;
	int vn;
	vector<int> hole;
	vector<int> holen;
	int hn = 0;
	vector<point> hv;
	bool boundmark;
};
struct region
{
	point rp;
	int rs;
	double ra;
};


PMG::PMG()
{
}

PMG::~PMG()
{
}



bool PMG::readOM()
{
	ifstream file(infilename);
	cout << infilename << endl;
	if (!file)
	{
		cout << "fail to open the file" << endl;
		return false;
	}
	string line;
	int pn = -1, fn = -1, hn = -1,rn=-1, i = 0;
	point tp;
	face tf;
	region tr;
	while (1)
	{
		file >> line;
		if (file.eof())
			break;
		if (line[0] == '#')
		{
			getline(file, line);
			continue;
		}
		if (pn == -1)
		{
			pn = stoi(line);
			getline(file, line);
			continue;
		}
		if (i < pn)
		{
			i++;
			file >> line;
			tp.cood[0] = stof(line);
			file >> line;
			tp.cood[1] = stof(line);
			file >> line;
			tp.cood[2] = stof(line);
			file >> line;
			tp.tag = stoi(line);
			PMG_point.push_back(tp);
			getline(file, line);
			continue;
		}
		if (fn == -1)
		{
			fn = stoi(line);
			getline(file, line);
			continue;
		}
		if (i < pn + fn)
		{
			i++;
			file >> line;
			tf.hn = stoi(line);
			file >> line;
			tf.boundmark = (stoi(line) == 1 ? true : false);
			getline(file, line);
			file >> line;
			tf.vn = stoi(line);
			for (int iter = 0; iter < tf.vn; iter++)
			{
				file >> line;
				tf.vertex.push_back(stoi(line));
			}
			getline(file, line);
			for (int iter = 0; iter < tf.hn; iter++)
			{
				file >> line;
				tf.holen.push_back(stoi(line));
				for (int iter2 = 0; iter2 < tf.holen[iter]; iter2++)
				{
					file >> line;
					tf.hole.push_back(stoi(line));
				}
				getline(file, line);
				file >> line;
				file >> line;
				tp.cood[0] = stof(line);
				file >> line;
				tp.cood[1] = stof(line);
				file >> line;
				tp.cood[2] = stof(line);
				tp.tag = -1;
				tf.hv.push_back(tp);
				getline(file, line);
			}
			PMG_face.push_back(tf);
			tf.vertex.clear();
			tf.hole.clear();
			tf.holen.clear();
			continue;
		}
		if (hn == -1)
		{
			hn = stoi(line);
			getline(file, line);
			continue;
		}
		if (i < pn + fn + hn)
		{
			i++;
			file >> line;
			tp.cood[0] = stof(line);
			file >> line;
			tp.cood[1] = stof(line);
			file >> line;
			tp.cood[2] = stof(line);
			tp.tag = -1;
			PMG_hole.push_back(tp);
			getline(file, line);
			continue;
		}
		if (rn == -1)
		{
			rn = stoi(line);
			getline(file, line);
			continue;
		}
		if (i < pn + fn + hn + rn)
		{
			i++;
			file >> line;
			tr.rp.cood[0] = stof(line);
			file >> line;
			tr.rp.cood[1] = stof(line);
			file >> line;
			tr.rp.cood[2] = stof(line);
			file >> line;
			tr.rs = stoi(line); 
			file >> line;
			tr.ra = stof(line);
			PMG_region.push_back(tr);
			getline(file, line);
			continue;
		}
	}
	return true;
}

template <int shift,bool mm>
bool compare(point left, point right)
{
	if (left.cood[shift] < right.cood[shift] )
		return mm;
	else if (left.cood[shift] == right.cood[shift])
	{
		if (left.cood[(shift + 1) % 3] < right.cood[(shift + 1) % 3])
			return mm;
		else if (left.cood[(shift + 1) % 3] == right.cood[(shift + 1) % 3])
		{
			if (left.cood[(shift + 2) % 3] < right.cood[(shift + 2) % 3])
				return mm;
			else
				return !mm;
		}
		else
			return !mm;
	}
	else
		return !mm;
}

bool PMG::divideM()
{
	double box[3][2];
	vector<point>::iterator tpi;
	tpi = max_element(PMG_point.begin(), PMG_point.end(), compare<0, true>);
	box[0][1] = (*tpi).cood[0];
	tpi = max_element(PMG_point.begin(), PMG_point.end(), compare<1, true>);
	box[1][1] = (*tpi).cood[1];
	tpi = max_element(PMG_point.begin(), PMG_point.end(), compare<2, true>);
	box[2][1] = (*tpi).cood[2];
	tpi = max_element(PMG_point.begin(), PMG_point.end(), compare<0, false>);
	box[0][0] = (*tpi).cood[0];
	tpi = max_element(PMG_point.begin(), PMG_point.end(), compare<1, false>);
	box[1][0] = (*tpi).cood[1];
	tpi = max_element(PMG_point.begin(), PMG_point.end(), compare<2, false>);
	box[2][0] = (*tpi).cood[2];

	
	
	return true;
}