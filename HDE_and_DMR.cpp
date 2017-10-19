#include<stdio.h>
#include<assert.h>
#include<math.h>
#include<random>
#include<time.h>
#include<vector>
#include<algorithm>

#include "detri2.h"
#include "pred2d.h"

using namespace detri2;


//ADD for high dimensional embedding is all under this line.
void Triangulation::lift_vrt()
{
	double nor;
	nor = (xmax - xmin) > (ymax - ymin) ? (xmax - xmin) : (ymax - ymin);
	double size = nor / (vrts->objects) / 1000;
	int j = 0;
	for (int i = 0; j<vrts->objects - unused_vertices; i++)
	{
		if (vrts->is_deleted((Vertex*)(vrts->lookup(i))) || ((Vertex*)(*vrts)[i])->get_vrttype() == 0)
			continue;
		((Vertex*)(vrts->lookup(i)))->highdi[0] = s*embed_fun(((Vertex*)(vrts->lookup(i)))->crd[0], ((Vertex*)(vrts->lookup(i)))->crd[1]);
		((Vertex*)(vrts->lookup(i)))->highdi[1] = s*(embed_fun(((Vertex*)(vrts->lookup(i)))->crd[0] + size, ((Vertex*)(vrts->lookup(i)))->crd[1]) - embed_fun(((Vertex*)(vrts->lookup(i)))->crd[0] - size, ((Vertex*)(vrts->lookup(i)))->crd[1])) / 2 / size;
		((Vertex*)(vrts->lookup(i)))->highdi[2] = s*(embed_fun(((Vertex*)(vrts->lookup(i)))->crd[0], ((Vertex*)(vrts->lookup(i)))->crd[1] + size) - embed_fun(((Vertex*)(vrts->lookup(i)))->crd[0], ((Vertex*)(vrts->lookup(i)))->crd[1] - size)) / 2 / size;
		j++;
		if (j == 1)
		{
			fmax = fmin = ((Vertex*)(*vrts)[i])->highdi[0];
			gxmax = gxmin = ((Vertex*)(*vrts)[i])->highdi[1];
			gymax = gymin = ((Vertex*)(*vrts)[i])->highdi[2];
		}
		if (fmax < ((Vertex*)(*vrts)[i])->highdi[0])
			fmax = ((Vertex*)(*vrts)[i])->highdi[0];
		if (fmin >((Vertex*)(*vrts)[i])->highdi[0])
			fmin = ((Vertex*)(*vrts)[i])->highdi[0];
		if (gxmax < ((Vertex*)(*vrts)[i])->highdi[1])
			gxmax = ((Vertex*)(*vrts)[i])->highdi[1];
		if (gxmin >((Vertex*)(*vrts)[i])->highdi[1])
			gxmin = ((Vertex*)(*vrts)[i])->highdi[1];
		if (gymax < ((Vertex*)(*vrts)[i])->highdi[2])
			gymax = ((Vertex*)(*vrts)[i])->highdi[2];
		if (gymin >((Vertex*)(*vrts)[i])->highdi[2])
			gymin = ((Vertex*)(*vrts)[i])->highdi[2];
	}
}

void Triangulation::edge_flip()
{
	chktris = new arraypool(sizeof(EncEle), 8);
	chktris2 = new arraypool(sizeof(EncEle), 8);

	TriEdge E;
	int count = 0;
	for (int i = 0; count<tris->objects; i++)
	{
		E.tri = tris->get(i);
		if (!(tris->is_deleted(E.tri)))
		{
			enq_triangle(chktris, E);
			count++;
		}
	}

	int loc = 1;

	while (chktris->objects>0)
	{
		arraypool *swapary = chktris;
		chktris = chktris2;
		chktris2 = swapary;
		for (int i = 0; i<chktris2->objects; i++)
		{
			EncEle *paryEle = (EncEle*)(*chktris2)[i];
			if (!(paryEle->is_changed()))
			{
				for (int j = 0; j < 3; j++)
				{
					if (check_spaceangle(paryEle))
					{
						//printf("check pass\n");
						flip_space22(paryEle);
						break;
					}
					else
					{
						//printf("failed\n");
						paryEle->Ele.enextself();
					}
				}
			}
		}
		chktris2->restart();
	}
	delete chktris;
	delete chktris2;
}

bool Triangulation::check_spaceangle(EncEle *paryEle)
{
	if (paryEle->Ele.is_segment())
		return 0; //paryEle->Ele should not be a boundary edge

	if (paryEle->Ele.org()->idx == -1 || paryEle->Ele.dest()->idx == -1)
		return 0; //paryEle->Ele should not be a TriEdge of two virtue triangle

	TriEdge tt[4];
	tt[0] = paryEle->Ele.enext2();
	tt[3] = paryEle->Ele.enext();
	tt[1] = paryEle->Ele.sym().enext();
	tt[2] = paryEle->Ele.sym().enext2();
	

	Vertex *vertem[4];
	vertem[0] = paryEle->Ele.apex();
	vertem[1] = paryEle->Ele.org();
	vertem[2] = paryEle->Ele.sym().apex();
	vertem[3] = paryEle->Ele.dest();

	double plangle[2] = { 0.0, 0.0 };
	double vectem[6][2];//vec of the two triangles,with ccw
	vectem[0][0] = vertem[1]->crd[0] - vertem[0]->crd[0]; vectem[0][1] = vertem[1]->crd[1] - vertem[0]->crd[1];
	vectem[1][0] = vertem[3]->crd[0] - vertem[1]->crd[0]; vectem[1][1] = vertem[3]->crd[1] - vertem[1]->crd[1];
	vectem[2][0] = vertem[0]->crd[0] - vertem[3]->crd[0]; vectem[2][1] = vertem[0]->crd[1] - vertem[3]->crd[1];
	vectem[3][0] = vertem[3]->crd[0] - vertem[2]->crd[0]; vectem[3][1] = vertem[3]->crd[1] - vertem[2]->crd[1];
	vectem[4][0] = vertem[1]->crd[0] - vertem[3]->crd[0]; vectem[4][1] = vertem[1]->crd[1] - vertem[3]->crd[1];
	vectem[5][0] = vertem[2]->crd[0] - vertem[1]->crd[0]; vectem[5][1] = vertem[2]->crd[1] - vertem[1]->crd[1];

	plangle[0] += acos((-1.0)*(vectem[0][0] * vectem[1][0] + vectem[0][1] * vectem[1][1]) / sqrt((vectem[0][0] * vectem[0][0] + vectem[0][1] * vectem[0][1])*(vectem[1][0] * vectem[1][0] + vectem[1][1] * vectem[1][1])));
	plangle[1] += acos((-1.0)*(vectem[2][0] * vectem[1][0] + vectem[2][1] * vectem[1][1]) / sqrt((vectem[2][0] * vectem[2][0] + vectem[2][1] * vectem[2][1])*(vectem[1][0] * vectem[1][0] + vectem[1][1] * vectem[1][1])));
	plangle[1] += acos((-1.0)*(vectem[3][0] * vectem[4][0] + vectem[3][1] * vectem[4][1]) / sqrt((vectem[3][0] * vectem[3][0] + vectem[3][1] * vectem[3][1])*(vectem[4][0] * vectem[4][0] + vectem[4][1] * vectem[4][1])));
	plangle[0] += acos((-1.0)*(vectem[5][0] * vectem[4][0] + vectem[5][1] * vectem[4][1]) / sqrt((vectem[5][0] * vectem[5][0] + vectem[5][1] * vectem[5][1])*(vectem[4][0] * vectem[4][0] + vectem[4][1] * vectem[4][1])));

	if (plangle[0]>PI*0.9 || plangle[1]>PI*0.9)
		return 0; //no angle bigger than 180 degree


	
	double distan[2];
	distan[0] = hi_distance(vertem[0], vertem[2]);
	distan[1] = hi_distance(vertem[1], vertem[3]);
	if (distan[1] <= distan[0])
		return 0;

	return 1;
}

double Triangulation::cal_angle(TriEdge &tt1, TriEdge &tt2)
{
	//tt1.org()->crd[0],tt1.org()->crd[1],tt1.org()->highdi[0],tt1.org()->highdi[1],tt1.org()->highdi[2]
	double vec[2][5];
	vec[0][0] = (tt1.org()->crd[0] - tt1.dest()->crd[0]) / (xmax - xmin);
	vec[0][1] = (tt1.org()->crd[1] - tt1.dest()->crd[1]) / (ymax - ymin);
	vec[0][2] = (tt1.org()->highdi[0] - tt1.dest()->highdi[0]) / (fmax - fmin);
	vec[0][3] = (tt1.org()->highdi[1] - tt1.dest()->highdi[1]) / (gxmax - gxmin);
	vec[0][4] = (tt1.org()->highdi[2] - tt1.dest()->highdi[2]) / (gymax - gymin);
	vec[1][0] = (tt2.dest()->crd[0] - tt2.org()->crd[0]) / (xmax - xmin);
	vec[1][1] = (tt2.dest()->crd[1] - tt2.org()->crd[1]) / (ymax - ymin);
	vec[1][2] = (tt2.dest()->highdi[0] - tt2.org()->highdi[0]) / (fmax - fmin);
	vec[1][3] = (tt2.dest()->highdi[1] - tt2.org()->highdi[1]) / (gxmax - gxmin);
	vec[1][4] = (tt2.dest()->highdi[2] - tt2.org()->highdi[2]) / (gymax - gymin);

	double product = 0.0;
	double vecpro[2] = { 0.0, 0.0 };
	for (int i = 0; i<5; i++)
	{
		product += vec[0][i] * vec[1][i];
		vecpro[0] += vec[0][i] * vec[0][i];
		vecpro[1] += vec[1][i] * vec[1][i];
	}
	double angle = acos(product / sqrt(vecpro[0] * vecpro[1]));//printf("angle is :%f\n",angle);
	return angle;
}

void Triangulation::flip_space22(EncEle *paryEle)
{
	TriEdge tt[2];
	tt[0] = paryEle->Ele;
	tt[1] = paryEle->Ele.sym();

	int hullflag = (tt[0].tri->is_hulltri() || tt[1].tri->is_hulltri());
	Vertex *pa = tt[0].org();
	Vertex *pb = tt[0].dest();
	Vertex *pc = tt[0].apex();
	Vertex *pd = tt[1].apex();

	TriEdge nn[4];
	nn[0] = (tt[0].enext()).sym();
	nn[1] = (tt[0].enext2()).sym();
	nn[2] = (tt[1].enext()).sym();
	nn[3] = (tt[1].enext2()).sym();

	tris->dealloc_tri(tt[0].tri);
	tris->dealloc_tri(tt[1].tri);

	tt[0].tri = tris->alloc_tri();
	tt[1].tri = tris->alloc_tri();

	tt[0].set_vertices(pc, pd, pb);
	tt[1].set_vertices(pd, pc, pa);

	if (hullflag)
	{
		if (pa == infvrt) {
			tt[1].tri->set_hullflag();
			hullsize -= 1; // Remove 2 hull edges, add 1. (Convex Assumption)
		}
		else if (pb == infvrt) {
			tt[0].tri->set_hullflag();
			hullsize -= 1;
		}
		else if (pc == infvrt) {
			tt[0].tri->set_hullflag();
			tt[1].tri->set_hullflag();
			hullsize += 1; // Remove 1 hull edges, add 2. (Convex Assumption)
		}
		else if (pd == infvrt) {
			tt[0].tri->set_hullflag();
			tt[1].tri->set_hullflag();
			hullsize += 1;
		}
		else {
			assert(0); // Not possible.
		}
	}

	tt[0].bond2(tt[1]);
	nn[0].bond2(tt[0].enext2()); // [b,c]
	nn[1].bond2(tt[1].enext());  // [c,a]
	nn[2].bond2(tt[1].enext2()); // [a,d]
	nn[3].bond2(tt[0].enext());  // [d,b]
	for (int i = 0; i < 4; i++)
	{
		if (nn[i].is_segment())
		{
			nn[i].sym().set_segment();
		}
	}

	pa->adj = tt[1].tri;
	pb->adj = tt[0].tri;
	pc->adj = tt[0].tri;
	pd->adj = tt[0].tri;

	enq_triangle(chktris, tt[0]);
	enq_triangle(chktris, tt[1]);
}

void Triangulation::edge_contraction()
{
	chktris = new arraypool(sizeof(EncEle), 8);
	chktris2 = new arraypool(sizeof(EncEle), 8);

	TriEdge E;
	for (int i = 0, j = 0; j < tris->objects; i++)
	{
		E.tri = tris->get(i);
		if (!tris->is_deleted(E.tri))
		{
			enq_triangle(chktris, E);
			j++;
		}
	}
	if (lawsonflag)
		lawsonflag = 0;
	int contracase;
	int loc = 0;
	while (chktris->objects>0)
	{
		//save_to_ucd(loc, 0);
		//loc++;
		//printf("chktris->objects is %d\n", chktris->objects);
		arraypool *sawapary = chktris;
		chktris = chktris2;
		chktris2 = sawapary;
		for (int i = 0; i < chktris2->objects; i++)
		{
			EncEle *paryEle = (EncEle*)(*chktris2)[i];
			if (!(paryEle->is_changed()))
			{
				for (int j = 0; j < 3; j++)
				{
					if (contracase = chk_contra(paryEle))
					{
						//printf("contracase is %d  ",contracase);
						contract_edge(paryEle, contracase);
						//save_to_ucd(100 + loc, 0);
						//loc++;
						break;
					}
					else
					{
						paryEle->Ele.enextself();
					}
				}
			}
		}
		chktris2->restart();
	}
	delete chktris;
	delete chktris2;
}

int Triangulation::chk_contra(EncEle *paryEle)
{
	TriEdge E = paryEle->Ele;
	int orgflag = 0, destflag = 0;
	if (hi_distance(E) > 0.5*L)
		return 0;
	if (E.is_segment())//TO BE improve
	{
		int c = 3;
		if (E.org()->get_vrttype() != 2)
			c = c + 1;
		if (E.dest()->get_vrttype() != 2)
			c = c + 2;
		if (c != 3)
			return c;
		else
			return 0;
	}
	else
	{
		if (E.org()->idx == -1 || E.dest()->idx == -1)
			return 0;
		TriEdge N = E;
		Vertex* org = E.org();
		Vertex* dest = E.dest();
		N.enext2self();
		N.symself();
		while (1)
		{
			if (N.dest() == dest)
			{
				break;
			}
			if (N.is_segment())
			{
				orgflag = 1;
				break;
			}
			N.enext2self();
			N.symself();
		}
		N.enextself();
		N.symself();
		while (1)
		{
			if (N.org() == org)
			{
				break;
			}
			if (N.is_segment())
			{
				destflag = 1;
				break;
			}
			N.enextself();
			N.symself();
		}
		return 3 - (2 * orgflag + destflag);
	}
}

void Triangulation::contract_edge(EncEle *paryEle, int contracase)
{
	TriEdge tt[4];
	Vertex *org, *dest;;
	if (contracase == 1)
	{
		dest = paryEle->Ele.dest();
		arraypool bound(sizeof(EncEle), 8);
		TriEdge E = paryEle->Ele;
		int count = 0;
		while (1)
		{
			((EncEle*)(bound.alloc()))->Ele = E.enext2().sym();
			count++;
			E.enextself();
			E.symself();
			if (E.org() == paryEle->Ele.org())
			{
				break;
			}
		}
		remove_point(dest, tt);
		for (int i = 0; i < count; i++)
		{
			enq_triangle(chktris, ((EncEle*)bound[i])->Ele.sym());
		}
	}
	if (contracase == 2 || contracase == 3)
	{
		org = paryEle->Ele.org();
		arraypool bound(sizeof(EncEle), 8);
		TriEdge E = paryEle->Ele;
		int count = 0;
		while (1)
		{
			((EncEle*)(bound.alloc()))->Ele = E.enext().sym();
			count++;
			E.enext2self();
			E.symself();
			if (E.dest() == paryEle->Ele.dest())
			{
				break;
			}
		}
		remove_point(org, tt);
		for (int i = 0; i < count; i++)
		{
			enq_triangle(chktris, ((EncEle*)bound[i])->Ele.sym());
		}
	}
	if (contracase == 4)
	{
		org = paryEle->Ele.org();
		arraypool bound(sizeof(EncEle), 8);
		TriEdge E = paryEle->Ele;
		int count = 0;
		while (1)
		{
			((EncEle*)(bound.alloc()))->Ele = E.enext().sym();
			count++;
			E.enext2self();
			E.symself();
			if (E.dest() == paryEle->Ele.dest())
			{
				break;
			}
		}
		remove_point(org, tt);
		for (int i = 0; i < count; i++)
		{
			enq_triangle(chktris, ((EncEle*)bound[i])->Ele.sym());
		}
	}
	if (contracase == 5 || contracase == 6)
	{
		dest = paryEle->Ele.dest();
		arraypool bound(sizeof(EncEle), 8);
		TriEdge E = paryEle->Ele;
		int count = 0;
		while (1)
		{
			((EncEle*)(bound.alloc()))->Ele = E.enext2().sym();
			count++;
			E.enextself();
			E.symself();
			if (E.org() == paryEle->Ele.org())
			{
				break;
			}
		}
		remove_point(dest, tt);
		for (int i = 0; i < count; i++)
		{
			enq_triangle(chktris, ((EncEle*)bound[i])->Ele.sym());
		}
	}
}

void Triangulation::edge_splitting()
{
	chktris = new arraypool(sizeof(EncEle), 8);
	chktris2 = new arraypool(sizeof(EncEle), 8);

	TriEdge E;
	int count = 0;
	for (int i = 0; count<tris->objects; i++)
	{
		E.tri = tris->get(i);
		if (!(tris->is_deleted(E.tri)))
		{
			enq_triangle(chktris, E);
			count++;
		}
	}
	if (lawsonflag)
		lawsonflag = 0;
	int loc = 0;
	while (chktris->objects > 0)
	{
		arraypool *swapary = chktris;
		chktris = chktris2;
		chktris2 = swapary;
		for (int i = 0; i < chktris2->objects; i++)
		{
			EncEle *paryEle = (EncEle*)(*chktris2)[i];
			if (!(paryEle->is_changed()))
			{
				for (int j = 0; j < 3; j++)
				{
					if (chk_split(paryEle))
					{
						/*Vertex* a = paryEle->Ele.org();
						Vertex* b = paryEle->Ele.dest();
						Vertex* c = paryEle->Ele.apex();

						printf("%f,%f\t%f,%f\t%f,%f\n%d\t%d\t%d\n%d\t%d\t%d\n",a->crd[0],a->crd[1],b->crd[0],b->crd[1],c->crd[0],c->crd[1],a->idx,b->idx,c->idx,a->tag,b->tag,c->tag);
						save_to_ucd(loc, 0);*/
						split_edge(paryEle->Ele);//split the edge and change the arraypool of triangals
						//loc++;
						break;
					}
					else
						paryEle->Ele.enextself();
				}
			}
		}
		chktris2->restart();
	}
	delete chktris;
	delete chktris2;
}

bool Triangulation::chk_split(EncEle *paryEle)
{
	if (paryEle->Ele.org()->idx == -1 || paryEle->Ele.dest()->idx == -1)
		return false;
	double distance = hi_distance(paryEle->Ele);
	if (distance>1.5*L)
		return true;
	return false;
}

int Triangulation::split_edge(TriEdge &tt)
{
	TriEdge tt2 = tt.sym();
	if ((Orient2d(tt.org()->crd, tt.dest()->crd, tt.apex()->crd) == 0 || Orient2d(tt2.org()->crd, tt2.dest()->crd, tt2.apex()->crd) == 0) && (!tt.is_segment()) && !(tt.apex()->idx == -1 || tt2.apex()->idx == -1))
	{
		TriEdge nn[4];
		nn[0] = tt;
		nn[1] = tt.sym();
		flip22(nn);
		return 0;
	}
	Vertex *mdvr = vrts->alloc_vrt();
	mdvr->idx = firstindex + (vrts->objects - 1);
	mdvr->crd[0] = (tt.org()->crd[0] + tt.dest()->crd[0]) / 2;
	mdvr->crd[1] = (tt.org()->crd[1] + tt.dest()->crd[1]) / 2;
	mdvr->crd[2] = 0;
	//printf("tt.org:%d,%f,%f;tt.dest:%d,%f,%f\n",tt.org()->idx,tt.org()->crd[0],tt.org()->crd[1],tt.dest()->idx,tt.dest()->crd[0],tt.dest()->crd[1]);
	TriEdge nn[4], mm[4];
	nn[0] = tt;

	mm[0] = tt.enext().sym();
	mm[1] = tt.enext2().sym();
	mm[2] = tt.sym().enext().sym();
	mm[3] = tt.sym().enext2().sym();
	insert_point(mdvr, nn, 2);
	unused_vertices++;
	mdvr->highdi[0] = s*embed_fun(mdvr->crd[0], mdvr->crd[1]);
	mdvr->highdi[1] = s*(embed_fun(mdvr->crd[0] + 1e-6, mdvr->crd[1]) - embed_fun(mdvr->crd[0] - 1e-6, mdvr->crd[1])) / 2 / 1e-6;
	mdvr->highdi[2] = s*(embed_fun(mdvr->crd[0], mdvr->crd[1] + 1e-6) - embed_fun(mdvr->crd[0], mdvr->crd[1] - 1e-6)) / 2 / 1e-6;
	for (int i = 0; i < 4; i++)
	{
		mm[i].symself();
		enq_triangle(chktris, mm[i]);
	}
	return 1;
}

void Triangulation::node_smoothing(double alph)
{
	arraypool *smoothing_direction = new arraypool(2 * sizeof(double), 8);
	TriEdge *tt = new TriEdge;
	Vertex *vrt;
	double dis = 0.0;
	double *direct;
	for (int i = 0, j = 0; j < vrts->objects - unused_vertices; i++)
	{
		if (((Vertex*)(*vrts)[i])->get_vrttype() == 0 || vrts->is_deleted((Vertex*)(*vrts)[i]))
			continue;
		if (((Vertex*)(*vrts)[i])->get_vrttype() == 6 || ((Vertex*)(*vrts)[i])->get_vrttype() == 2)
		{
			j++;
			continue;
		}
		vrt = (Vertex*)(*vrts)[i];
		direct = (double*)smoothing_direction->alloc();
		direct[0] = 0.0;
		direct[1] = 0.0;
		tt->tri = vrt->adj;
		for (int k = 0; k < 3; k++)
		{
			tt->ver = k;
			if (tt->org() == vrt)
				break;
		}
		while (1)
		{
			if (tt->dest() == infvrt)
			{
				printf("tt->org's vrttype is %d\nidx is %d\ncrd is (%f,%f)\n,i=%d,j=%d\n", tt->org()->get_vrttype(), tt->org()->idx, tt->org()->crd[0], tt->org()->crd[1], i, j);
				assert(0);
			}
			//dis = hi_distance(*tt);
			dis = sqrt(pow(tt->dest()->crd[0] - tt->org()->crd[0], 2) + pow(tt->dest()->crd[1] - tt->org()->crd[1], 2));
			direct[0] += (tt->dest()->crd[0] - tt->org()->crd[0]) / dis*(1 - pow(dis, 4)*exp(0.0 - dis));//here the smoothing method
			direct[1] += (tt->dest()->crd[1] - tt->org()->crd[1]) / dis*(1 - pow(dis, 4)*exp(0.0 - dis));
			tt->symself();
			tt->enextself();
			assert(tt->org() == vrt);
			if (tt->tri == vrt->adj)
				break;
		}
		j++;
		/*tt->enextself();
		double temppo[2];
		while (1)
		{
		assert(tt->apex() == vrt);
		temppo[0] = vrt->crd[0] + alph*direct[0];
		temppo[1] = vrt->crd[1] + alph*direct[1];
		if (Orient2d(tt->org()->crd, tt->dest()->crd, temppo)!=1);
		{
		direct[0] = 0.0;
		direct[1] = 0.0;
		break;
		}
		tt->enextself();
		tt->symself();
		tt->enext2self();
		if (tt->tri == vrt->adj)
		break;
		}*/
	}
	for (int i = 0, j = 0; j < smoothing_direction->objects; i++)
	{
		if (((Vertex*)(*vrts)[i])->get_vrttype() == 0 || vrts->is_deleted((Vertex*)(*vrts)[i]))
			continue;
		if (((Vertex*)(*vrts)[i])->get_vrttype() == 6 || ((Vertex*)(*vrts)[i])->get_vrttype() == 2)
		{
			continue;
		}
		direct = (double*)(*smoothing_direction)[j];
		((Vertex*)(*vrts)[i])->crd[0] = ((Vertex*)(*vrts)[i])->crd[0] + alph*direct[0];
		((Vertex*)(*vrts)[i])->crd[1] = ((Vertex*)(*vrts)[i])->crd[1] + alph*direct[1];
		int tagg = 0;
		while (1)
		{
			if (tagg == 0)
			{
				tt->tri = ((Vertex*)(*vrts)[i])->adj;
				for (int k = 0; k < 3; k++)
				{
					tt->ver = k;
					if (tt->apex() == ((Vertex*)(*vrts)[i]))
						break;
				}
			}
			assert(tt->apex() == (Vertex*)(*vrts)[i]);
			if (tagg != 0 && tt->tri == ((Vertex*)(*vrts)[i])->adj)
				break;
			tagg += 1;
			if (Orient2d(tt->org()->crd, tt->dest()->crd, ((Vertex*)(*vrts)[i])->crd) <= 0)
			{
				((Vertex*)(*vrts)[i])->crd[0] = ((Vertex*)(*vrts)[i])->crd[0] - alph*direct[0];
				((Vertex*)(*vrts)[i])->crd[1] = ((Vertex*)(*vrts)[i])->crd[1] - alph*direct[1];
				break;
			}
			else
			{
				tt->enextself();
				tt->symself();
				tt->enextself();
			}
		}
		((Vertex*)(*vrts)[i])->highdi[0] = s*embed_fun(((Vertex*)(*vrts)[i])->crd[0], ((Vertex*)(*vrts)[i])->crd[1]);
		((Vertex*)(*vrts)[i])->highdi[1] = s*(embed_fun(((Vertex*)(*vrts)[i])->crd[0] + ((xmax - xmin) / (vrts->objects) / 1000), ((Vertex*)(*vrts)[i])->crd[1]) - embed_fun(((Vertex*)(*vrts)[i])->crd[0] - ((xmax - xmin) / (vrts->objects) / 1000), ((Vertex*)(*vrts)[i])->crd[1])) / 2 / ((xmax - xmin) / (vrts->objects) / 1000);
		((Vertex*)(*vrts)[i])->highdi[2] = s*(embed_fun(((Vertex*)(*vrts)[i])->crd[0], ((Vertex*)(*vrts)[i])->crd[1] + ((xmax - xmin) / (vrts->objects) / 1000)) - embed_fun(((Vertex*)(*vrts)[i])->crd[0], ((Vertex*)(*vrts)[i])->crd[1] - ((xmax - xmin) / (vrts->objects) / 1000))) / 2 / ((xmax - xmin) / (vrts->objects) / 1000);
		j++;
	}
	delete smoothing_direction;
	delete tt;
}

double Triangulation::interpolation_error_sum()
{
	double sum = 0;
	for (int i = 0, j = 0; i < tris->objects + j; i++)
	{
		if (!(((Triangle*)(*tris)[i])->is_hulltri()) && !(tris->is_deleted((Triangle*)(*tris)[i])))
			sum = sum + interpolation_error((Triangle*)(*tris)[i]);
		else if ((((Triangle*)(*tris)[i])->is_hulltri()))
			continue;
		else
		{
			j++;
			continue;
		}
	}
	return sum;
}

double Triangulation::interpolation_error(Triangle *tri)
{
	double error = 0.0;
	int m = 3, n = 3;
	double coef[3];//the coefficient of the linear interpolation function.Inter-fun(x,y)=coef[0]*x+coef[1]*y+coef[2]
	double x[2], y[2];//x[0]=xmin x[1]=xmax y si similar
	Vertex* temp[3];
	for (int i = 0; i < 3; i++)
	{
		temp[i] = tri->vrt[i];
	}
	formula_solver(coef, temp);
	if (coef[0] == 0.0 && coef[1] == 0.0 && coef[2] == 0.0)
		return 0.0;

	//have a test for the 4 formulas
	x[0] = x[1] = temp[0]->crd[0];
	y[0] = y[1] = temp[0]->crd[1];
	for (int i = 1; i < 3; i++)
	{
		if (temp[i]->crd[0] < x[0])
			x[0] = temp[i]->crd[0];
		if (temp[i]->crd[0]>x[1])
			x[1] = temp[i]->crd[0];
		if (temp[i]->crd[1] < y[0])
			y[0] = temp[i]->crd[1];
		if (temp[i]->crd[1]>y[1])
			y[1] = temp[i]->crd[1];
	}
	double h = (x[1] - x[0]) / (double)m;
	double g = (y[1] - y[0]) / (double)n;


	//ti xing gong shi
	//error = (x[1] - x[0])*(y[1] - y[0]) / 4 * (extra_func(temp, coef, x[0], y[0]) + extra_func(temp, coef, x[0], y[1]) + extra_func(temp, coef, x[1], y[0]) + extra_func(temp, coef, x[1], y[1]));


	//fu hua ti xing gong  shi
	/*for (int i = 0; i <= m; i++)
	{
	for (int j = 0; j <= n; j++)
	{
	if ((i == 0 || i == m) && (j == 0 || j == n))
	{
	error += extra_func(temp, coef, x[0] + i*h, y[0] + j*g);
	}
	else if (i == 0 || i == m || j == 0 || j == n)
	{
	error += 2 * extra_func(temp, coef, x[0] + i*h, y[0] + j*g);
	}
	else
	{
	error += 4 * extra_func(temp, coef, x[0] + i*h, y[0] + j*g);
	}
	}
	}
	error = error*h*g / 4.0;*/


	//xin pu sen gong shi
	/*double far = extra_func(temp, coef, x[0], y[0]) + extra_func(temp, coef, x[0], y[1]) + extra_func(temp, coef, x[1], y[0]) + extra_func(temp, coef, x[1], y[1]);
	double near = extra_func(temp, coef, x[0], (y[0] + y[1]) / 2) + extra_func(temp, coef, x[1], (y[0] + y[1]) / 2) + extra_func(temp, coef, (x[0] + x[1])/2, y[0]) + extra_func(temp, coef, (x[0] + x[1])/2, y[1]);
	double mid = extra_func(temp, coef, (x[0] + x[1]) / 2, (y[0] + y[1]) / 2);
	error = (x[1] - x[0])*(y[1] - y[0]) / 36 * (far + 4 * near + 16 * mid);*/


	//fu hua xin pu sen gong shi
	/*m = 2 * m;
	n = 2 * n;
	h = h / 2;
	g = g / 2;
	for (int i = 0; i <= m; i++)
	{
	for (int j = 0; j <= n; j++)
	{
	if ((i == 0 || i == m) && (j == 0 || j == n))
	error += extra_func(temp, coef, x[0] + i*h, y[0] + j*g);
	else if (((i == 0 || i == m) && j%2==1)||(i%2==1 && (j==0 || j==n)))
	error += 4 * extra_func(temp, coef, x[0] + i*h, y[0] + j*g);
	else if (i==0 || i==m || j==0 || j==n)
	error += 2 * extra_func(temp, coef, x[0] + i*h, y[0] + j*g);
	else if (i%2==0 && j%2==0)
	error += 4 * extra_func(temp, coef, x[0] + i*h, y[0] + j*g);
	else if ((i%2==0 && j%2==1)||(i%2==1 && j%2==0))
	error += 8 * extra_func(temp, coef, x[0] + i*h, y[0] + j*g);
	else
	error += 16 * extra_func(temp, coef, x[0] + i*h, y[0] + j*g);
	}
	}*/
	//4 kind fi
	//gauss inter point
	double area = (temp[2]->crd[0] - temp[0]->crd[0])*(temp[1]->crd[1] - temp[0]->crd[1]) - (temp[1]->crd[0] - temp[0]->crd[0])*(temp[2]->crd[1] - temp[0]->crd[1]);
	area = area / -2.0;
	double* gqpc;
	double area_crd[2];
	//1 D
	area_crd[0] = 1.0 / 3.0;
	area_crd[1] = 1.0 / 3.0;
	gqpc = area_to_axis(area_crd, *tri);
	error = area*extra_func(temp, coef, gqpc[0], gqpc[1]);
	return error;
}

int Triangulation::formula_solver(double *solve, Vertex* *vrt)
{
	double mat[3][3];
	double right[3];
	double temp[3][3];
	for (int i = 0; i < 3; i++)
	{
		mat[i][0] = vrt[i]->crd[0];
		mat[i][1] = vrt[i]->crd[1];
		mat[i][2] = 1.0;
		right[i] = embed_fun(vrt[i]->crd[0], vrt[i]->crd[1]);
	}

	double det = detmination(mat);
	if (det < 1e-12)
	{
		for (int i = 0; i < 3; i++)
		{
			printf("cord of one point is:%f %f %f\n", mat[i][0], mat[i][1], right[i]);
		}
		for (int i = 0; i < 3; i++)
		{
			solve[i] = 0.0;
		}
		return 0;
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				if (k == i)
					temp[j][k] = right[k];
				else
					temp[j][k] = mat[j][k];
			}
		}
		solve[i] = detmination(temp) / det;
	}
	return 0;
}

double inline Triangulation::detmination(double mat[3][3])
{
	return mat[0][0] * mat[1][1] * mat[2][2] + mat[2][0] * mat[0][1] * mat[1][2] + mat[1][0] * mat[2][1] * mat[0][2] - mat[2][0] * mat[1][1] * mat[0][2] - mat[1][0] * mat[0][1] * mat[2][2] - mat[0][0] * mat[2][1] * mat[1][2];
}

double Triangulation::extra_func(Vertex* temp[3], double coef[3], double x, double y)
{
	double error;
	double pointscrd[3][2];
	for (int i = 0; i < 3; i++)
	{
		pointscrd[i][0] = temp[i]->crd[0];
		pointscrd[i][1] = temp[i]->crd[1];
	}
	double crd[2] = { x, y };
	assert(Orient2d(pointscrd[0], pointscrd[1], crd) >= 0 && Orient2d(pointscrd[1], pointscrd[2], crd) >= 0 && Orient2d(pointscrd[2], pointscrd[0], crd) >= 0);
	error = (coef[0] * x + coef[1] * y + coef[2] - embed_fun(x, y));
	return error*error;
}

double* Triangulation::area_to_axis(double area[2], Triangle tr)
{
	double crd[3][2];
	for (int i = 0; i < 3; i++)
	{
		crd[i][0] = tr.vrt[i]->crd[0];
		crd[i][1] = tr.vrt[i]->crd[1];
	}
	double answer[2] = { area[0], area[1] };
	answer[0] = answer[0] * (crd[0][0] * (crd[1][1] - crd[2][1]) + crd[0][1] * (crd[2][0] - crd[1][0]) + crd[1][0] * crd[2][1] - crd[2][0] * crd[1][1]) - crd[1][0] * crd[2][1] + crd[2][0] * crd[1][1];
	answer[1] = answer[1] * (crd[1][0] * (crd[2][1] - crd[0][1]) + crd[1][1] * (crd[0][0] - crd[2][0]) + crd[2][0] * crd[0][1] - crd[0][0] * crd[2][1]) - crd[2][0] * crd[0][1] + crd[0][0] * crd[2][1];
	double static temp[2];
	temp[0] = (answer[0] * (crd[0][0] - crd[2][0]) - answer[1] * (crd[2][0] - crd[1][0])) / ((crd[1][1] - crd[2][1])*(crd[0][0] - crd[2][0]) - (crd[2][1] - crd[0][1])*(crd[2][0] - crd[1][0]));
	temp[1] = (answer[0] * (crd[2][1] - crd[0][1]) - answer[1] * (crd[1][1] - crd[2][1])) / ((crd[2][0] - crd[1][0])*(crd[2][1] - crd[0][1]) - (crd[0][0] - crd[2][0])*(crd[1][1] - crd[2][1]));
	return temp;
}





double Triangulation::embed_fun(double x, double y)
{
	//return tanh(2 * (x+0.6 - y) - 1.0);
	//return tanh(2*(x+para[0]-(y+para[1]))-1.0);
	return 1.0 / ((x + para[0])*(x + para[0]) + (y + para[1])*(y + para[1]) + 1.0);
}

void Triangulation::high_dimen_embedding()
{
	printf("\ninitial interpolation error is: %f\n", interpolation_error_sum());
	s = 1;
	L = 0.1;
	lift_vrt();
	for (int i = 0; i < 5; i++)
	{
		edge_flip();
		if (maxarea == 0.0)
		{
			save_to_ucd(i * 5 + 0);
			edge_contraction();
			save_to_ucd(i * 5 + 1);
			edge_flip();
			save_to_ucd(i * 5 + 2);
			if (idx2seglist != NULL)
			{
				delete[] idx2seglist;
				delete[] segperverlist;
				idx2seglist = NULL;
				segperverlist = NULL;
			}
			make_vertex_to_segment_map();
			edge_splitting();
			save_to_ucd(i * 5 + 3);
			edge_flip();
			save_to_ucd(i * 5 + 4);
		}
		printf("total interpolation error is: %f\n", interpolation_error_sum());
		//node_smoothing(0.2);
		//edge_flip();
	}
	getchar();
	//printf("total interpolation error is: %f\n", interpolation_error_sum());
	save_to_ucd(totalflipcount + count + 2, 1);
	printf("hd_embed fi\n");
}

double Triangulation::hi_distance(TriEdge E)
{
	double length = 0;
	double a[5];
	a[0] = (E.org()->crd[0] - E.dest()->crd[0]) / (xmax - xmin);
	a[1] = (E.org()->crd[1] - E.dest()->crd[1]) / (ymax - ymin);
	a[2] = (E.org()->highdi[0] - E.dest()->highdi[0]) / (fmax - fmin);
	a[3] = (E.org()->highdi[1] - E.dest()->highdi[1]) / (gxmax - gxmin);
	a[4] = (E.org()->highdi[2] - E.dest()->highdi[2]) / (gymax - gymin);
	for (int i = 0; i < 5; i++)
	{
		length += a[i] * a[i];
	}
	return sqrt(length);
}
double Triangulation::hi_distance(Vertex* pt1, Vertex* pt2)
{
	double length = 0.0;
	double a[5];
	a[0] = (pt1->crd[0] - pt2->crd[0]) / (xmax - xmin);
	a[1] = (pt1->crd[1] - pt2->crd[1]) / (ymax - ymin);
	a[2] = (pt1->highdi[0] - pt2->highdi[0]) / (fmax - fmin);
	a[3] = (pt1->highdi[1] - pt2->highdi[1]) / (gxmax - gxmin);
	a[4] = (pt1->highdi[2] - pt2->highdi[2]) / (gymax - gymin);
	for (int i = 0; i < 5; i++)
	{
		length += a[i] * a[i];
	}
	return sqrt(length);
}



















































//ADD for locally remesh
void Triangulation::load_move_vertices(double *a, double *b, int *c, int k)
{
	MovVertex *mvv;
	int i , j = 0;
	for (i = 0; j < vrts->objects - unused_vertices; i++)
	{
		if (vrts->is_deleted((Vertex*)(*vrts)[i]) || ((Vertex*)(*vrts)[i])->get_vrttype() == 0)
			continue;
		if (((Vertex*)(*vrts)[i])->tag != 0)
		{
			for (int l = 0; l < k; l++)
			{
				if (((Vertex*)(*vrts)[i])->tag == c[l])
				{
					mvv = (MovVertex*)(moving_vrts->alloc());
					mvv->mov_vrt = (Vertex*)((*vrts)[i]);
					mvv->vec[0] = a[l];
					mvv->vec[1] = b[l];
				}
			}
		}
		j++;
	}
}

void Triangulation::findandremove_segpoint()
{
	TriEdge tt[4];
	for (int i = 0, j = 0; j < vrts->objects - unused_vertices; i++)
	{
		if (vrts->is_deleted((Vertex*)(vrts->lookup(i))) || ((Vertex*)(*vrts)[i])->get_vrttype() == 0)
		{
			continue;
		}
		if (((Vertex*)(*vrts)[i])->get_vrttype() == 6)
		{
			Vertex* vr = (Vertex*)(*vrts)[i];
			TriEdge *E = new TriEdge;
			E->tri = vr->adj;
			for (int k = 0; k < 3; k++)
			{
				E->ver = k;
				if (E->org() == vr)
				{
					break;
				}
			}
			while (1)
			{
				if (E->is_segment())
				{
					break;
				}
				else
				{
					E->symself();
					E->enextself();
				}
			}
			Triangle *tr = get_segment(*E);
			for (int k = 0; k < moving_vrts->objects; k++)
			{
				if (tr->vrt[0] == ((MovVertex*)(*moving_vrts)[k])->mov_vrt || tr->vrt[1] == ((MovVertex*)(*moving_vrts)[k])->mov_vrt)
				{

					remove_point((Vertex*)(*vrts)[i], tt);
					j--;
					break;
				}
			}
			delete E;
		}
		j++;
	}
}

void Triangulation::delete_collision_points()
{
	int searching_list[30000];
	int searching_amount;
	TriEdge *E = new TriEdge;
	double dir1[2], dir2[2];
	double pt1[2], pt2[2], pt3[2], pt4[2];
	int panju = 0;
	double EL12[2], EL23[2], EL34[2], EL41[2];
	TriEdge tt[4];
	double xb[2], yb[2];
	double xt, yt,xd,yd;
	//get the bounding box
	for (int i = 0; i < moving_vrts->objects; i++)
	{
		xt = ((MovVertex*)(*moving_vrts)[i])->mov_vrt->crd[0];
		yt = ((MovVertex*)(*moving_vrts)[i])->mov_vrt->crd[1];
		xd = ((MovVertex*)(*moving_vrts)[i])->vec[0] + xt;
		yd = ((MovVertex*)(*moving_vrts)[i])->vec[1] + yt;
		if (i == 0)
		{
			xb[0] = xb[1] = xt;
			yb[0] = yb[1] = yt;
		}
		if (xt> xb[1])
			xb[1]=xt;
		if (xt < xb[0])
			xb[0] = xt;
		if (yt > yb[1])
			yb[1] = yt;
		if (yt < yb[0])
			yb[0] = yt;
		if (xd> xb[1])
			xb[1] = xd;
		if (xd < xb[0])
			xb[0] = xd;
		if (yd > yb[1])
			yb[1] = yd;
		if (yd < yb[0])
			yb[0] = yd;
	}
	xd = xb[0]; xt = xb[1];
	yd = yb[0]; yt = yb[1];
	//xb[0] = 1.25*xd - 0.25*xt;
	//xb[1] = 1.25*xt - 0.25*xd;
	//yb[0] = 1.25*yd - 0.25*yt;
	//yb[1] = 1.25*yt - 0.25*yd;
	xb[0] -= 0.5;
	xb[1] += 0.5;
	yb[0] -= 0.5;
	yb[1] += 0.5;
	//gotten
	//get the searching list
	searching_amount = 0;
	double jiji = (double)ji;
	for (int k = 0; k < vrts->objects ; k++)
	{
		if (vrts->is_deleted((Vertex*)(*vrts)[k]) || ((Vertex*)(*vrts)[k])->get_vrttype() == 0)
			continue;
		if (((Vertex*)(*vrts)[k])->get_vrttype() == 6)
			continue;
		if (((Vertex*)(*vrts)[k])->tag != 0)
			continue;
		xt = ((Vertex*)(*vrts)[k])->crd[0];
		yt = ((Vertex*)(*vrts)[k])->crd[1];
		if (xt<xb[0] || xt>xb[1] || yt<yb[0] || yt>yb[1])
			continue;
		else
		{
			searching_list[searching_amount] = k;
			searching_amount++;
		}
		/*if ((xt > 2.2 - jiji*0.05 && xt<4.5 - jiji*0.05 && yt>-0.4 && yt < 1.0) || (xt<-2.2 + jiji*0.05 && xt>-4.5 + jiji*0.05 && yt<0.4 && yt>-1.0))
		{
			searching_list[searching_amount] = k;
			searching_amount++;
		}*/
	}
	ji++;
	//gotten
	for (int i = 0; i < searching_amount; i++)
	{
		remove_point((Vertex*)(*vrts)[searching_list[i]], tt);
	}





	//for (int i = 0; i < moving_vrts->objects; i++)
	//{
	//	E->tri = ((MovVertex*)(*moving_vrts)[i])->mov_vrt->adj;
	//	for (int j = 0; j < 3; j++)
	//	{
	//		E->ver = j;
	//		if (E->org() == ((MovVertex*)(*moving_vrts)[i])->mov_vrt)
	//		{
	//			break;
	//		}
	//	}
	//	Vertex* temp = E->dest();
	//	while (1)
	//	{
	//		if (E->is_segment())
	//		{
	//			if(E->org()->tag < E->dest()->tag || E->dest()->tag==0)
	//			{
	//				jilu[0][ji] = E->org();
	//				jilu[1][ji] = E->dest();
	//				ji++;
	//			}
	//		}
	//		E->symself();
	//		E->enextself();
	//		if (E->dest() == temp)
	//		{
	//			break;
	//		}
	//	}
	//}
	//for (int j = 0; j < ji; j++)
	//{
	//	dir1[0] = dir1[1] = dir2[0] = dir2[1] = 0;
	//	get_edge(jilu[0][j], jilu[1][j], *E);
	//	for (int k = 0; k < moving_vrts->objects; k++)
	//	{
	//		if (E->org() == ((MovVertex*)(*moving_vrts)[k])->mov_vrt)
	//		{
	//			dir1[0] = ((MovVertex*)(*moving_vrts)[k])->vec[0];
	//			dir1[1] = ((MovVertex*)(*moving_vrts)[k])->vec[1];
	//		}
	//		if (E->dest() == ((MovVertex*)(*moving_vrts)[k])->mov_vrt)
	//		{
	//			dir2[0] = ((MovVertex*)(*moving_vrts)[k])->vec[0];
	//			dir2[1] = ((MovVertex*)(*moving_vrts)[k])->vec[1];
	//		}
	//	}

	//	pt1[0] = E->org()->crd[0]; pt1[1] = E->org()->crd[1];
	//	pt2[0] = E->dest()->crd[0]; pt2[1] = E->dest()->crd[1];
	//	pt3[0] = pt2[0] + dir2[0]; pt3[1] = pt2[1] + dir2[1];
	//	pt4[0] = pt1[0] + dir1[0]; pt4[1] = pt1[1] + dir1[1];

	//	EL12[0] = pt1[0] - pt2[0]; EL12[1] = pt1[1] - pt2[1];
	//	EL23[0] = pt2[0] - pt3[0]; EL23[1] = pt2[1] - pt3[1];
	//	EL34[0] = pt3[0] - pt4[0]; EL34[1] = pt3[1] - pt4[1];
	//	EL41[0] = pt4[0] - pt1[0]; EL41[1] = pt4[1] - pt1[1];

	//	xt = EL12[0] * EL12[0] + EL12[1] * EL12[1];
	//	xt = sqrt(xt);
	//	yt = (EL23[0] * EL23[0] + EL41[0] * EL41[0] + EL23[1] * EL23[1] + EL41[1] * EL41[1])/2.0;
	//	yt = sqrt(yt);
	//	if (yt == 0.0)
	//		continue;
	//	xd = xt / yt/2.0;
	//	double Len[4]={EL12[0]*EL12[0]+EL12[1]*EL12[1],EL23[0]*EL23[0]+EL23[1]*EL23[1],EL34[0]*EL34[0]+EL34[1]*EL34[1],EL41[0]*EL41[0]+EL41[1]*EL41[1]};

	//	int ori[2]={Orient2d(pt1,pt2,pt3),Orient2d(pt1,pt2,pt4)};
	//	if(ori[0]==ori[1])
	//	{
	//		if (ori[0] == 0)//remove nothing
	//			continue;
	//		else
	//		{
	//			pt1[0] = pt1[0] + 0.5*EL12[0]; pt1[1] = pt1[1] + 0.5*EL12[1];
	//			pt2[0] = pt2[0] - 0.5*EL12[0]; pt2[1] = pt2[1] - 0.5*EL12[1];
	//			pt3[0] = pt3[0] - xd*EL23[0] + 0.5*EL34[0]; pt3[1] = pt3[1] - xd*EL23[1] + 0.5*EL34[1];
	//			pt4[0] = pt4[0] - 0.5*EL34[0] + xd*EL41[0]; pt4[1] = pt4[1] - 0.5*EL34[1] + xd*EL41[1];
	//			for (int k = 0; k < searching_amount; k++)
	//			{
	//				int sn = searching_list[k];
	//				panju = Orient2d(pt1, pt2, ((Vertex*)(*vrts)[sn])->crd) + Orient2d(pt2, pt3, ((Vertex*)(*vrts)[sn])->crd) + Orient2d(pt3, pt4, ((Vertex*)(*vrts)[sn])->crd) + Orient2d(pt4, pt1, ((Vertex*)(*vrts)[sn])->crd);
	//				if (panju == 4 || panju == -4)
	//				{
	//					remove_point((Vertex*)(*vrts)[sn], tt);
	//					check_mesh();
	//				}
	//			}
	//			//for (int k = 0, l = 0; l < vrts->objects - unused_vertices; k++)
	//			//{
	//			//	if (vrts->is_deleted((Vertex*)(*vrts)[k]) || ((Vertex*)(*vrts)[k])->get_vrttype() == 0)
	//			//	{
	//			//		continue;
	//			//	}
	//			//	if ((Vertex*)(*vrts)[k] == E->org() || (Vertex*)(*vrts)[k] == E->dest() || ((Vertex*)(*vrts)[k])->get_vrttype() == 6)
	//			//	{
	//			//		l++;
	//			//		continue;
	//			//	}
	//			//	if (((Vertex*)(*vrts)[k])->tag != 0)
	//			//	{
	//			//		l++;
	//			//		continue;
	//			//	}
	//			//	int x = 0;
	//			//	/*for (int m = 0; m < moving_vrts->objects; m++)
	//			//	{
	//			//		MovVertex* vt = ((MovVertex*)(*moving_vrts)[m]);
	//			//		double dis1 = pow((vt->mov_vrt->crd[0] - ((Vertex*)(*vrts)[k])->crd[0]), 2) + pow((vt->mov_vrt->crd[1] - ((Vertex*)(*vrts)[k])->crd[1]), 2);
	//			//		double dis2 = pow(vt->vec[0], 2) + pow(vt->vec[1], 2);
	//			//		if (dis1 < dis2)
	//			//		{
	//			//			remove_point((Vertex*)(*vrts)[k], tt);
	//			//			x = 1;
	//			//			break;
	//			//		}
	//			//	}
	//			//	if (x == 1)
	//			//		continue;*/

	//			//	panju = Orient2d(pt1, pt2, ((Vertex*)(*vrts)[k])->crd) + Orient2d(pt2, pt3, ((Vertex*)(*vrts)[k])->crd) + Orient2d(pt3, pt4, ((Vertex*)(*vrts)[k])->crd) + Orient2d(pt4, pt1, ((Vertex*)(*vrts)[k])->crd);
	//			//	if (panju == 4 || panju == -4)
	//			//	{
	//			//		remove_point((Vertex*)(*vrts)[k], tt);
	//			//		continue;
	//			//		//j--;
	//			//		//break;
	//			//	}
	//			//	l++;
	//			//}
	//		}
	//	}
	//	else
	//	{
	//		if(ori[0]==0 || ori[1]==0)
	//		{
	//			pt3[0] = pt3[0] - 1.0*EL23[0] + 0.5*EL34[0]; pt3[1] = pt3[1] - 1.0*EL23[1] + 0.5*EL34[1];
	//			pt4[0] = pt4[0] - 0.5*EL34[0] + 1.0*EL41[0]; pt4[1] = pt4[1] - 0.5*EL34[1] + 1.0*EL41[1];
	//			for (int k = 0, l = 0; l < vrts->objects - unused_vertices; k++)
	//			{
	//				if (vrts->is_deleted((Vertex*)(*vrts)[k]) || ((Vertex*)(*vrts)[k])->get_vrttype() == 0)
	//				{
	//					continue;
	//				}
	//				if ((Vertex*)(*vrts)[k] == E->org() || (Vertex*)(*vrts)[k] == E->dest() || ((Vertex*)(*vrts)[k])->get_vrttype() == 6)
	//				{
	//					l++;
	//					continue;
	//				}
	//				if (((Vertex*)(*vrts)[k])->tag != 0)
	//				{
	//					l++;
	//					continue;
	//				}

	//				panju = Orient2d(pt1, pt2, ((Vertex*)(*vrts)[k])->crd) + Orient2d(pt2, pt3, ((Vertex*)(*vrts)[k])->crd) + Orient2d(pt3, pt4, ((Vertex*)(*vrts)[k])->crd) + Orient2d(pt4, pt1, ((Vertex*)(*vrts)[k])->crd);
	//				if (panju == 4 || panju == -4)
	//				{
	//					remove_point((Vertex*)(*vrts)[k], tt);
	//					//j--;
	//					//break;
	//					l--;
	//					continue;
	//				}
	//				l++;
	//			}
	//		}
	//		else
	//		{
	//			pt1[0] = pt1[0] + 0.5*EL12[0]; pt1[1] = pt1[1] + 0.5*EL12[1];
	//			pt2[0] = pt2[0] - 0.5*EL12[0]; pt2[1] = pt2[1] - 0.5*EL12[1];
	//			pt3[0] = pt3[0] - 1.0*EL23[0] + 0.5*EL34[0]; pt3[1] = pt3[1] - 1.0*EL23[1] + 0.5*EL34[1];
	//			pt4[0] = pt4[0] - 0.5*EL34[0] + 1.0*EL41[0]; pt4[1] = pt4[1] - 0.5*EL34[1] + 1.0*EL41[1];

	//			for (int k = 0, l = 0; l < vrts->objects - unused_vertices; k++)
	//			{
	//				if (vrts->is_deleted((Vertex*)(*vrts)[k]) || ((Vertex*)(*vrts)[k])->get_vrttype() == 0)
	//				{
	//					continue;
	//				}
	//				if ((Vertex*)(*vrts)[k] == E->org() || (Vertex*)(*vrts)[k] == E->dest() || ((Vertex*)(*vrts)[k])->get_vrttype() == 6)
	//				{
	//					l++;
	//					continue;
	//				}
	//				if (((Vertex*)(*vrts)[k])->tag != 0)
	//				{
	//					l++;
	//					continue;
	//				}

	//				int coun[2]={Orient2d(pt1, pt2, ((Vertex*)(*vrts)[k])->crd) + Orient2d(pt2, pt3, ((Vertex*)(*vrts)[k])->crd) + Orient2d(pt3, pt4, ((Vertex*)(*vrts)[k])->crd),Orient2d(pt3, pt4, ((Vertex*)(*vrts)[k])->crd) + Orient2d(pt4, pt1, ((Vertex*)(*vrts)[k])->crd) + Orient2d(pt1, pt2, ((Vertex*)(*vrts)[k])->crd)};
	//				if ( coun[0]== 3 || coun[0] == -3 || coun[1]==3 || coun[1]==-3)
	//				{
	//					remove_point((Vertex*)(*vrts)[k], tt);
	//					//j--;
	//					//break;
	//					l--;
	//					continue;
	//				}
	//				l++;
	//			}
	//		}
	//	}		
	//}
	delete E;
}



void Triangulation::change_vrts()
{
	Vertex* pttest = new Vertex;
	TriEdge* edtest = new TriEdge[4];
	int flag;
	int ori[2];
	for (int i = 0; i < moving_vrts->objects; i++) 
	{
		pttest->crd[0] = ((MovVertex*)(*moving_vrts)[i])->mov_vrt->crd[0] + ((MovVertex*)(*moving_vrts)[i])->vec[0];
		pttest->crd[1] = ((MovVertex*)(*moving_vrts)[i])->mov_vrt->crd[1] + ((MovVertex*)(*moving_vrts)[i])->vec[1];
		flag = -1;
		while (true)
		{
			if (flag == -1)
			{
				edtest[0].tri = ((MovVertex*)(*moving_vrts)[i])->mov_vrt->adj;
				for (int j = 0; j < 3; j++)
				{
					edtest[0].ver = j;
					if (edtest[0].apex() == ((MovVertex*)(*moving_vrts)[i])->mov_vrt)
						break;
				}
			}
			if (edtest[0].tri == ((MovVertex*)(*moving_vrts)[i])->mov_vrt->adj && flag != -1)
				break;
			if (edtest[0].org() == infvrt || edtest[0].dest() == infvrt)
			{
				edtest[0] = edtest[0].enext().sym().enext();
				flag = 0;
				continue;
			}
			edtest[1] = edtest[0].enext();
			edtest[2] = edtest[0].enext2();
			ori[0] = Orient2d(edtest[1].org()->crd, edtest[1].dest()->crd, pttest->crd);
			ori[1] = Orient2d(edtest[2].org()->crd, edtest[2].dest()->crd, pttest->crd);
			if (ori[0] == 1 && ori[1] == 1)
			{
				ori[0] = Orient2d(edtest[0].org()->crd, edtest[0].dest()->crd, pttest->crd);
				if (ori[0] > 0)
					flag = 1;
				else
				{
					//enforce_flip(edtest[0]);
					edtest[1] = edtest[0].sym();
					flip22(edtest);
					flag = -1;
					continue;
				}
			}
			else if (ori[0] == 0 || ori[1] == 0)
			{
				if (Orient2d(edtest[0].org()->crd, edtest[0].dest()->crd, pttest->crd) > 0)
				{
					flag = 1;
				}
				else if (Orient2d(edtest[0].org()->crd, edtest[0].dest()->crd, pttest->crd) == 0)
				{
					//coincidence
					if (ori[0] == 0)
						remove_point(edtest[1].org(), edtest);
					if (ori[1] == 0)
						remove_point(edtest[2].dest(), edtest);
					flag = -1;
					continue;
				}
				else
				{
					if (ori[0] == 0)
						remove_point(edtest[1].org(), edtest);
					if (ori[1] == 0)
						remove_point(edtest[2].dest(), edtest);
					flag = -1;
					continue;
				}
			}
			else
			{
				edtest[0] = edtest[0].enext().sym().enext();
				flag = 0;
				continue;
			}
			if (flag = 1)
				break;
		}
	}
	//crossing edge removed.
	for (int i = 0; i <moving_vrts->objects; i++)
	{
		((MovVertex*)(*moving_vrts)[i])->mov_vrt->crd[0] = ((MovVertex*)(*moving_vrts)[i])->mov_vrt->crd[0] + ((MovVertex*)(*moving_vrts)[i])->vec[0];
		((MovVertex*)(*moving_vrts)[i])->mov_vrt->crd[1] = ((MovVertex*)(*moving_vrts)[i])->mov_vrt->crd[1] + ((MovVertex*)(*moving_vrts)[i])->vec[1];
	}
}


void Triangulation::locally_improvement()
{

	TriEdge tt[4];
	double ccent[3];
	int loc, i;
	chktris = new arraypool(sizeof(EncEle), 8);
	chktris2 = new arraypool(sizeof(EncEle), 8);
	for (int i = 0; i < moving_vrts->objects; i++)
	{
		tt[0].tri = ((MovVertex*)(*moving_vrts)[i])->mov_vrt->adj;
		for (int j = 0; j < 3; j++)
		{
			tt[0].ver = j;
			if (tt[0].org() == ((MovVertex*)(*moving_vrts)[i])->mov_vrt)
				break;
		}
		tt[1].tri = tt[0].tri;
		tt[1].ver = tt[0].ver;
		while (1)
		{
			if (tt[0].dest() != infvrt && tt[0].apex() != infvrt)
			{
				enq_triangle(chktris, tt[0]);
			}
			tt[0].symself();
			tt[0].enextself();
			if (tt[0].tri == tt[1].tri)
				break;
		}
	}
	/*for (int i = 0; i < ji; i++)
	{
		get_edge(jilu[0][i], jilu[1][i], tt[0]);
		assert(tt[0].is_segment());
		if (tt[0].apex() == infvrt)
			tt[0].symself();
		enq_triangle(chktris, tt[0]);
	}*/

	if (verbose) {
		printf("Adding Steiner points to enforce quality.\n");
	}
	make_vertex_to_segment_map();
	encsegs = new arraypool(sizeof(EncEle), 4);
	encsegs2 = new arraypool(sizeof(EncEle), 4);

	lawsonflag = 1;
	encsegflag = 1;


	if (quality) { // -q
		

		// Calculate the squre of the maximum radius-edge ratio.
		// r/d = 1/(2*sin(theta_min));
		maxratio2 = 0.5 / sin(PI * minangle / 180.0);
		maxratio2 *= maxratio2;

		enqtriflag = 1;

		repair_triangles();

		enqtriflag = 0;

		delete chktris;
		delete chktris2;
		chktris = NULL;
		chktris2 = NULL;
	}

	lawsonflag = 0;
	encsegflag = 0;

	delete encsegs;
	delete encsegs2;
	encsegs = NULL;
	encsegs2 = NULL;

	delete[] idx2seglist;
	delete[] segperverlist;
	idx2seglist = NULL;
	segperverlist = NULL;
}




































void Triangulation::get_flip_queue()
{
	//get star and push it to the flip queue
	TriEdge TE, TF, *TT;

	for (int i = 0; i < (moving_vrts->objects); i++)
	{
		if (((MovVertex*)(*moving_vrts)[i])->cla != 1)
		{
			TE.tri = ((MovVertex*)(*moving_vrts)[i])->mov_vrt->adj;
			for (int j = 0; j < 3; j++)
			{
				if ((TE.tri->vrt[j]) == ((MovVertex*)(*moving_vrts)[i])->mov_vrt)
				{
					TE.ver = (j + 2) % 3;
				}
			}
			TF = TE;
			do
			{
				for (int k = 0; k < 3; k++)
				{
					TT = (TriEdge*)flip_queue->alloc();
					TT->tri = TF.tri;
					TT->ver = k;

					//TT->clear_segment();
				}
				TF.symself();
				TF.enextself();
			} while (TE.tri != TF.tri);
		}
	}
	//get the trace of the moving points and push the TriEdge into the flip queue
	for (int i = 0; i < (moving_vrts->objects); i++)
	{
		if (((MovVertex*)(*moving_vrts)[i])->cla != 1)
		{
			Vertex* pt;
			pt = new Vertex;
			TE.tri = (Triangle*)((MovVertex*)(*moving_vrts)[i])->mov_vrt->adj;
			TE.ver = 0;
			pt->crd[0] = ((MovVertex*)(*moving_vrts)[i])->mov_vrt->crd[0] + ((MovVertex*)(*moving_vrts)[i])->vec[0];
			pt->crd[1] = ((MovVertex*)(*moving_vrts)[i])->mov_vrt->crd[1] + ((MovVertex*)(*moving_vrts)[i])->vec[1];
			int j;
			if (TE.tri->is_hulltri()) {
				// Get a non-hull triangle.
				for (j = 0; j < 3; j++) {
					if (is_hulledge(TE)) break;
					TE.enextself();
				}
				assert(j < 3);
				TE.symself();
				assert(!TE.tri->is_hulltri());
			}

			// Select an edge such that pt lies to CCW of it.
			for (j = 0; j < 3; j++) {
				if (Orient2d((TE.org())->crd, (TE.dest())->crd, pt->crd) > 0) break;
				TE.enextself();
			}
			for (int k = 0; k < 3; k++)
			{
				TT = (TriEdge*)flip_queue->alloc();
				TT->tri = TE.tri;
				TT->ver = k;
			}

			// Let E = [a,b,c] and p lies to the CCW of [a->b].
			int ori1, ori2;

			do {

				ori1 = Orient2d((TE.dest())->crd, (TE.apex())->crd, pt->crd);
				ori2 = Orient2d((TE.apex())->crd, (TE.org())->crd, pt->crd);

				if (ori1 > 0) {
					if (ori2 > 0) {
						break; // Found.
					}
					else if (ori2 < 0) {
						TE.enext2self();
					}
					else { // ori2 == 0
						TE.enext2self(); // p lies on edge [c,a]
						for (int k = 0; k < 3; k++)
						{
							TT = (TriEdge*)flip_queue->alloc();
							TT->tri = TE.tri;
							TT->ver = k;
						}
					}
				}
				else if (ori1 < 0) {
					if (ori2 > 0) {
						TE.enextself();
					}
					else if (ori2 < 0) {
						// Randomly choose one.
						if (rand() % 2) { // flipping a coin.
							TE.enextself();
						}
						else {
							TE.enext2self();
						}
					}
					else { // ori2 == 0
						TE.enextself();
					}
				}
				else { // ori1 == 0
					if (ori2 > 0) {
						TE.enextself(); // p lies on edge [b,c].
						for (int k = 0; k < 3; k++)
						{
							TT = (TriEdge*)flip_queue->alloc();
							TT->tri = TE.tri;
							TT->ver = k;
						}
					}
					else if (ori2 < 0) {
						TE.enext2self();
					}
					else { // ori2 == 0
						TE.enext2self(); // p is coincident with apex.
						for (int k = 0; k < 3; k++)
						{
							TT = (TriEdge*)flip_queue->alloc();
							TT->tri = TE.tri;
							TT->ver = k;
						}
					}
				}

				if (encsegflag) {
					if (TE.is_segment()) {
						for (int k = 0; k < 3; k++)
						{
							TT = (TriEdge*)flip_queue->alloc();
							TT->tri = TE.tri;
							TT->ver = k;
						}
					}
				}

				TE.symself();
			} while (!TE.tri->is_hulltri());
		}
	}
}

struct anglecell
{
	double angle;
	double coordinate[3][2];
};

bool comp(anglecell a, anglecell b)
{
	return a.angle < b.angle;
}


void Triangulation::angle_statistic()
{
	Triangle* tr;
	std::vector<anglecell> smallangles;
	double vt[3][2], vd[2][2];
	double inner, length,at;
	anglecell ac;
	for (int i = 0; i < tris->objects; i++)
	{
		at = 1.0;
		tr = (Triangle*)(*tris)[i];
		if (tris->is_deleted(tr))
			continue;
		if (tr->is_hulltri())
			continue;
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 2; k++)
				vt[j][k] = tr->vrt[j]->crd[k];
		}
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				for (int l = 0; l < 2; l++)
					vd[k][l] = vt[(j + 2 - k)%3][l] - vt[j][l];
			}
			inner=vd[0][0] * vd[1][1] - vd[0][1] * vd[1][0];
			length = sqrt(vd[0][0] * vd[0][0] + vd[0][1] * vd[0][1])*sqrt(vd[1][0] * vd[1][0] + vd[1][1] * vd[1][1]);
			inner = fabs(inner)/length;
			inner < at ? at = inner : 1;
		}
		ac.angle = at;
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 2; k++)
				ac.coordinate[j][k] = vt[j][k];
		}
		smallangles.push_back(ac);
	}
	std::sort(smallangles.begin(), smallangles.end(),comp);
	printf("smallest angle is arcsin(%f) at %f,%f  %f,%f  %f,%f\n", smallangles[0].angle, smallangles[0].coordinate[0][0], smallangles[0].coordinate[0][1], smallangles[0].coordinate[1][0], smallangles[0].coordinate[1][1], smallangles[0].coordinate[2][0], smallangles[0].coordinate[2][1]);
	printf("biggest angle is arcsin(%f) at %f,%f  %f,%f  %f,%f\n", smallangles[smallangles.size()-1].angle, smallangles[smallangles.size() - 1].coordinate[0][0], smallangles[smallangles.size() - 1].coordinate[0][1], smallangles[smallangles.size() - 1].coordinate[1][0], smallangles[smallangles.size() - 1].coordinate[1][1], smallangles[smallangles.size() - 1].coordinate[2][0], smallangles[smallangles.size() - 1].coordinate[2][1]);
}















void Triangulation::remesh_Triangulation(double *x, double *y, int *c, int k)
{
	printf("start time %d\n", clock());
	if (moving_vrts == NULL)
		moving_vrts = (new arraypool(sizeof(MovVertex), 4));
	load_move_vertices(x, y, c, k);
	if (gabriel || quality || maxarea > 0.0 || minedgelen)
	{
		findandremove_segpoint();
	}
	//printf("remove seg time %d\n", clock());
	lawsonflag = 1;
	delete_collision_points();
	//printf("delete point time %d\n", clock());
	change_vrts();
	//lawson_flip(NULL);
	if ((gabriel || quality) && (!hdflag))
	{
		locally_improvement();
	}
	//printf("refine tri time %d\n", clock());
	lawsonflag = 0;
	if (dump_to_ucd)
	{
		save_to_ucd(totalflipcount + count + 2, 1);
	}
	count++;
	moving_vrts->restart();
}