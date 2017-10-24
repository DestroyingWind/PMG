#include<vector>
#include<iostream>
#include<fstream>
#include<algorithm>
#include<string>
#include<assert.h>
#include<map>

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
				tf.vertex.push_back(stoi(line)-1);
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

	vector<point> inner_point;
	vector<int> bound_point;
	vector<int> inner_sn;
	int position = 0;
	double inner_box[3][2];

	for (tpi = PMG_point.begin(); tpi < PMG_point.end(); tpi++)
	{
		bool xyz[3];
		for (int i = 0; i < 3; i++)
			xyz[i] = ((*tpi).cood[i] == box[i][0]) || ((*tpi).cood[i] == box[i][1]);
		if (xyz[0] || xyz[1] || xyz[2])
		{
			bound_point.push_back(position);
			position++;
			continue;
		}
		else
		{
			inner_point.push_back(*tpi);
			inner_sn.push_back(position);
		}
		position++;
	}

	tpi = max_element(inner_point.begin(), inner_point.end(), compare<0, true>);
	inner_box[0][1] = (*tpi).cood[0];
	tpi = max_element(inner_point.begin(), inner_point.end(), compare<1, true>);
	inner_box[1][1] = (*tpi).cood[1];
	tpi = max_element(inner_point.begin(), inner_point.end(), compare<2, true>);
	inner_box[2][1] = (*tpi).cood[2];
	tpi = max_element(inner_point.begin(), inner_point.end(), compare<0, false>);
	inner_box[0][0] = (*tpi).cood[0];
	tpi = max_element(inner_point.begin(), inner_point.end(), compare<1, false>);
	inner_box[1][0] = (*tpi).cood[1];
	tpi = max_element(inner_point.begin(), inner_point.end(), compare<2, false>);
	inner_box[2][0] = (*tpi).cood[2];

	for (int i = 0; i < 3; i++)
	{
		double temp[2] = { inner_box[i][0],inner_box[i][1] };
		inner_box[i][0] = temp[0] * 1.0 - temp[1] * 0.0;
		inner_box[i][1] = temp[1] * 1.0 - temp[0] * 0.0;
	}

	for (int i = 0; i < 3; i++)
	{
		assert(inner_box[i][1]<box[i][1] && inner_box[i][0]>box[i][0]);
	}

	for (int x = -1; x < 2; x++)
	{
		for (int y = -1; y < 2; y++)
		{
			for (int z = -1; z < 2; z++)
			{
				char buf[1024];
				sprintf_s(buf, "temp_file%d,%d,%d.poly", x, y, z);
				ofstream f(buf);
				if (x*x + y*y + z*z == 0)
				{
					vector<int>::iterator op,fn;
					vector<face>::iterator fi;
					for (fi = PMG_face.begin(); fi < PMG_face.end();)
					{
						for (fn = (*fi).vertex.begin(); fn < (*fi).vertex.end(); fn++)
						{
							for (op = bound_point.begin(); op < bound_point.end(); op++)
							{
								if ((*op) == (*fn))
								{
									fi = PMG_face.erase(fi);
									break;
								}
							}
							if (op < bound_point.end())
								break;
						}
						if(op==bound_point.end())
							fi++;
					}
					map<int, int> dic;
					int count = 1;
					for (op=inner_sn.begin();op<inner_sn.end();op++)
					{
						if ((dic.insert(map<int, int>::value_type((*op), count))).second)
							count++;
					}
					f << 8 + inner_sn.size() << "  3  0  1\n";
					for (op = inner_sn.begin(); op<inner_sn.end(); op++)
					{
						f << dic[*op] << "  " << PMG_point[*op].cood[0] << "  " << PMG_point[*op].cood[1] << "  " << PMG_point[*op].cood[2] << endl;
					}
					for (int i = 0; i < 2; i++)
					{
						for (int j = 0; j < 2; j++)
						{
							for (int k = 0; k < 2; k++)
							{
								f << inner_sn.size() + 1 + i * 4 + j * 2 + k << "  " << inner_box[0][i] << "  " << inner_box[0][j] << "  " << inner_box[0][k]<<endl;
							}
						}
					}
					f << 6 + PMG_face.size() << "  1\n";
					for (fi = PMG_face.begin(); fi < PMG_face.end(); fi++)
					{
						f << "1  0  1\n" << fi->vn<<"  ";
						for (int i = 0; i < fi->vn; i++)
							f << "  " << dic[fi->vertex[i]];
						f << endl;
					}
					f << "1  0  1\n" << "4    " << inner_sn.size() + 1 << "  " << inner_sn.size() + 2 << "  " << inner_sn.size() + 4 << "  " << inner_sn.size() + 3 << endl;
					f << "1  0  1\n" << "4    " << inner_sn.size() + 5 << "  " << inner_sn.size() + 6 << "  " << inner_sn.size() + 8 << "  " << inner_sn.size() + 7 << endl;
					f << "1  0  1\n" << "4    " << inner_sn.size() + 1 << "  " << inner_sn.size() + 2 << "  " << inner_sn.size() + 6 << "  " << inner_sn.size() + 5 << endl;
					f << "1  0  1\n" << "4    " << inner_sn.size() + 3 << "  " << inner_sn.size() + 4 << "  " << inner_sn.size() + 8 << "  " << inner_sn.size() + 7 << endl;
					f << "1  0  1\n" << "4    " << inner_sn.size() + 1 << "  " << inner_sn.size() + 3 << "  " << inner_sn.size() + 7 << "  " << inner_sn.size() + 5 << endl;
					f << "1  0  1\n" << "4    " << inner_sn.size() + 2 << "  " << inner_sn.size() + 4 << "  " << inner_sn.size() + 8 << "  " << inner_sn.size() + 6 << endl;
					f << PMG_hole.size()<<endl;
					for (int i = 0; i < PMG_hole.size();i++)
						f << i + 1 << "  " << PMG_hole[i].cood[0] << "  " << PMG_hole[i].cood[1] << "  " << PMG_hole[i].cood[2]<<endl;
					f << 0;
				}
				if (x*x + y*y + z*z == 1)
				{
					f << "8  3  0  1\n";
					if (x != 0)
					{
						f << 1 << "  " << box[0][(x + 1) / 2] << "  " << inner_box[1][1] << "  " << inner_box[2][1] << endl;
						f << 2 << "  " << box[0][(x + 1) / 2] << "  " << inner_box[1][1] << "  " << inner_box[2][0] << endl;
						f << 3 << "  " << box[0][(x + 1) / 2] << "  " << inner_box[1][0] << "  " << inner_box[2][1] << endl;
						f << 4 << "  " << box[0][(x + 1) / 2] << "  " << inner_box[1][0] << "  " << inner_box[2][0] << endl;
						f << 5 << "  " << inner_box[0][(x + 1) / 2] << "  " << inner_box[1][1] << "  " << inner_box[2][1] << endl;
						f << 6 << "  " << inner_box[0][(x + 1) / 2] << "  " << inner_box[1][1] << "  " << inner_box[2][0] << endl;
						f << 7 << "  " << inner_box[0][(x + 1) / 2] << "  " << inner_box[1][0] << "  " << inner_box[2][1] << endl;
						f << 8 << "  " << inner_box[0][(x + 1) / 2] << "  " << inner_box[1][0] << "  " << inner_box[2][0] << endl;
					}
					if (y != 0)
					{
						f << 1 << "  " << inner_box[0][1] << "  " << box[1][(y + 1) / 2] << "  " << inner_box[2][1] << endl;
						f << 2 << "  " << inner_box[0][1] << "  " << box[1][(y + 1) / 2] << "  " << inner_box[2][0] << endl;
						f << 3 << "  " << inner_box[0][1] << "  " << inner_box[1][(y + 1) / 2] << "  " << inner_box[2][1] << endl;
						f << 4 << "  " << inner_box[0][1] << "  " << inner_box[1][(y + 1) / 2] << "  " << inner_box[2][0] << endl;
						f << 5 << "  " << inner_box[0][0] << "  " << box[1][(y + 1) / 2] << "  " << inner_box[2][1] << endl;
						f << 6 << "  " << inner_box[0][0] << "  " << box[1][(y + 1) / 2] << "  " << inner_box[2][0] << endl;
						f << 7 << "  " << inner_box[0][0] << "  " << inner_box[1][(y + 1) / 2] << "  " << inner_box[2][1] << endl;
						f << 8 << "  " << inner_box[0][0] << "  " << inner_box[1][(y + 1) / 2] << "  " << inner_box[2][0] << endl;
					}
					if (z != 0)
					{
						f << 1 << "  " << inner_box[0][1] << "  " << inner_box[1][1] << "  " << box[2][(z + 1) / 2] << endl;
						f << 2 << "  " << inner_box[0][1] << "  " << inner_box[1][1] << "  " << inner_box[2][(z + 1) / 2] << endl;
						f << 3 << "  " << inner_box[0][1] << "  " << inner_box[1][0] << "  " << box[2][(z + 1) / 2] << endl;
						f << 4 << "  " << inner_box[0][1] << "  " << inner_box[1][0] << "  " << inner_box[2][(z + 1) / 2] << endl;
						f << 5 << "  " << inner_box[0][0] << "  " << inner_box[1][1] << "  " << box[2][(z + 1) / 2] << endl;
						f << 6 << "  " << inner_box[0][0] << "  " << inner_box[1][1] << "  " << inner_box[2][(z + 1) / 2] << endl;
						f << 7 << "  " << inner_box[0][0] << "  " << inner_box[1][0] << "  " << box[2][(z + 1) / 2] << endl;
						f << 8 << "  " << inner_box[0][0] << "  " << inner_box[1][0] << "  " << inner_box[2][(z + 1) / 2] << endl;
					}
					f << "6  1\n";
					f << "1  0  1\n" << "4    1  2  4  3\n";
					f << "1  0  1\n" << "4    5  6  8  7\n";
					f << "1  0  1\n" << "4    1  2  6  5\n";
					f << "1  0  1\n" << "4    3  4  8  7\n";
					f << "1  0  1\n" << "4    1  3  7  5\n";
					f << "1  0  1\n" << "4    2  4  8  6\n";
					f << "0\n0";
				}
				if (x*x + y*y + z*z == 2)
				{
					f << "8  3  0  1\n";
					if(x==0)
					{
						f << 1 << "  " << inner_box[0][1] << "  " << box[1][(y + 1) / 2] << "  " << box[2][(z + 1) / 2] << endl;
						f << 2 << "  " << inner_box[0][1] << "  " << box[1][(y + 1) / 2] << "  " << inner_box[2][(z + 1) / 2] << endl;
						f << 3 << "  " << inner_box[0][1] << "  " << inner_box[1][(y + 1) / 2] << "  " << box[2][(z + 1) / 2] << endl;
						f << 4 << "  " << inner_box[0][1] << "  " << inner_box[1][(y + 1) / 2] << "  " << inner_box[2][(z + 1) / 2] << endl;
						f << 5 << "  " << inner_box[0][0] << "  " << box[1][(y + 1) / 2] << "  " << box[2][(z + 1) / 2] << endl;
						f << 6 << "  " << inner_box[0][0] << "  " << box[1][(y + 1) / 2] << "  " << inner_box[2][(z + 1) / 2] << endl;
						f << 7 << "  " << inner_box[0][0] << "  " << inner_box[1][(y + 1) / 2] << "  " << box[2][(z + 1) / 2] << endl;
						f << 8 << "  " << inner_box[0][0] << "  " << inner_box[1][(y + 1) / 2] << "  " << inner_box[2][(z + 1) / 2] << endl;
					}
					if (y == 0)
					{
						f << 1 << "  " << box[0][(x + 1) / 2] << "  " << inner_box[1][1] << "  " << box[2][(z + 1) / 2] << endl;
						f << 2 << "  " << box[0][(x + 1) / 2] << "  " << inner_box[1][1] << "  " << inner_box[2][(z + 1) / 2] << endl;
						f << 3 << "  " << box[0][(x + 1) / 2] << "  " << inner_box[1][0] << "  " << box[2][(z + 1) / 2] << endl;
						f << 4 << "  " << box[0][(x + 1) / 2] << "  " << inner_box[1][0] << "  " << inner_box[2][(z + 1) / 2] << endl;
						f << 5 << "  " << inner_box[0][(x + 1) / 2] << "  " << inner_box[1][1] << "  " << box[2][(z + 1) / 2] << endl;
						f << 6 << "  " << inner_box[0][(x + 1) / 2] << "  " << inner_box[1][1] << "  " << inner_box[2][(z + 1) / 2] << endl;
						f << 7 << "  " << inner_box[0][(x + 1) / 2] << "  " << inner_box[1][0] << "  " << box[2][(z + 1) / 2] << endl;
						f << 8 << "  " << inner_box[0][(x + 1) / 2] << "  " << inner_box[1][0] << "  " << inner_box[2][(z + 1) / 2] << endl;
					}
					if (z == 0)
					{
						f << 1 << "  " << box[0][(x + 1) / 2] << "  " << box[1][(y + 1) / 2] << "  " << inner_box[2][1] << endl;
						f << 2 << "  " << box[0][(x + 1) / 2] << "  " << box[1][(y + 1) / 2] << "  " << inner_box[2][0] << endl;
						f << 3 << "  " << box[0][(x + 1) / 2] << "  " << inner_box[1][(y + 1) / 2] << "  " << inner_box[2][1] << endl;
						f << 4 << "  " << box[0][(x + 1) / 2] << "  " << inner_box[1][(y + 1) / 2] << "  " << inner_box[2][0] << endl;
						f << 5 << "  " << inner_box[0][(x + 1) / 2] << "  " << box[1][(y + 1) / 2] << "  " << inner_box[2][1] << endl;
						f << 6 << "  " << inner_box[0][(x + 1) / 2] << "  " << box[1][(y + 1) / 2] << "  " << inner_box[2][0] << endl;
						f << 7 << "  " << inner_box[0][(x + 1) / 2] << "  " << inner_box[1][(y + 1) / 2] << "  " << inner_box[2][1] << endl;
						f << 8 << "  " << inner_box[0][(x + 1) / 2] << "  " << inner_box[1][(y + 1) / 2] << "  " << inner_box[2][0] << endl;
					}
					f << "6  1\n";
					f << "1  0  1\n" << "4    1  2  4  3\n";
					f << "1  0  1\n" << "4    5  6  8  7\n";
					f << "1  0  1\n" << "4    1  2  6  5\n";
					f << "1  0  1\n" << "4    3  4  8  7\n";
					f << "1  0  1\n" << "4    1  3  7  5\n";
					f << "1  0  1\n" << "4    2  4  8  6\n";
					f << "0\n0";
				}
				if (x*x + y*y + z*z == 3)
				{
					f << "8  3  0  1\n";
					f << 1 << "  " << box[0][(x + 1) / 2] << "  " << box[1][(y + 1) / 2] << "  " << box[2][(z + 1) / 2] << endl;
					f << 2 << "  " << box[0][(x + 1) / 2] << "  " << box[1][(y + 1) / 2] << "  " << inner_box[2][(z + 1) / 2] << endl;
					f << 3 << "  " << box[0][(x + 1) / 2] << "  " << inner_box[1][(y + 1) / 2] << "  " << box[2][(z + 1) / 2] << endl;
					f << 4 << "  " << box[0][(x + 1) / 2] << "  " << inner_box[1][(y + 1) / 2] << "  " << inner_box[2][(z + 1) / 2] << endl;
					f << 5 << "  " << inner_box[0][(x + 1) / 2] << "  " << box[1][(y + 1) / 2] << "  " << box[2][(z + 1) / 2] << endl;
					f << 6 << "  " << inner_box[0][(x + 1) / 2] << "  " << box[1][(y + 1) / 2] << "  " << inner_box[2][(z + 1) / 2] << endl;
					f << 7 << "  " << inner_box[0][(x + 1) / 2] << "  " << inner_box[1][(y + 1) / 2] << "  " << box[2][(z + 1) / 2] << endl;
					f << 8 << "  " << inner_box[0][(x + 1) / 2] << "  " << inner_box[1][(y + 1) / 2] << "  " << inner_box[2][(z + 1) / 2] << endl;
					f << "6  1\n";
					f << "1  0  1\n" << "4    1  2  4  3\n";
					f << "1  0  1\n" << "4    5  6  8  7\n";
					f << "1  0  1\n" << "4    1  2  6  5\n";
					f << "1  0  1\n" << "4    3  4  8  7\n";
					f << "1  0  1\n" << "4    1  3  7  5\n";
					f << "1  0  1\n" << "4    2  4  8  6\n";
					f << "0\n0";
				}
				f.close();
			}
		}
	}

	return true;
}