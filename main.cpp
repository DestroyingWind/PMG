#include<fstream>
#include<iostream>
#include<vector>
#include<math.h>
#include<time.h>
#include<string>
#include<vector>

#include "PMG.h"

using namespace std;


int main(int argc,char* argv[])
{
	PMG test;
	test.readOM();
	test.divideM();
	test.writeNM();
	return 0;
}