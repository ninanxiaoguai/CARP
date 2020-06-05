
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include<bits/stdc++.h>
//#include <algorithm.h>
#include "functions.h"
#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void rand_selection(individual *p1, individual *p2, individual *pop)
/* pop is sorted increasingly already */
{
	int k1, k2;
	int candi[MAX_POPSIZE+1];
	candi[0] = popsize;
	for (int i = 1; i <= popsize; i++)
	{
		candi[i] = i-1;
	}

	k1 = rand_choose(candi[0]);
	delete_element(candi, k1);
	k2 = rand_choose(candi[0]);
	//printf("k1 = %d, k2 = %d, popsize = %d\n", k1, k2, popsize);
	indi_copy(p1, &pop[candi[k1]]);
	indi_copy(p2, &pop[candi[k2]]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void tour_selection(individual *p1, individual *p2, individual *pop)
/* pop is sorted increasingly already */
{
	int k1, k2;
	int candi1[MAX_POPSIZE+1], candi2[MAX_POPSIZE+1];
	candi1[0] = popsize;
	for (int i = 1; i <= popsize; i++)
	{
		candi1[i] = i-1;
	}
	memcpy(candi2, candi1, sizeof(candi1));

	k1 = rand_choose(candi1[0]);
	delete_element(candi1, k1);
	k2 = rand_choose(candi1[0]);
	if (k1 < k2)
	{
		indi_copy(p1, &pop[candi1[k1]]);
		delete_element(candi2, k1);
	}
	else
	{
		indi_copy(p1, &pop[candi1[k2]]);
		delete_element(candi2, k2);
	}

	k1 = rand_choose(candi2[0]);
	delete_element(candi2, k1);
	k2 = rand_choose(candi2[0]);
	if (k1 < k2)
	{
		indi_copy(p2, &pop[candi2[k1]]);
	}
	else
	{
		indi_copy(p2, &pop[candi2[k2]]);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SBX(individual *xed_child, individual *p1, individual *p2, const task *inst_tasks)
{
	int i, j, k, l, m, n, w, k1, k2, NO_LeftTasks, IVLoad;

	int S1[MAX_TASK_SEQ_LENGTH], S2[MAX_TASK_SEQ_LENGTH], SubPath1[MAX_TASK_SEQ_LENGTH], SubPath2[MAX_TASK_SEQ_LENGTH];
	int Routes1[MAX_SEG_TAG_LENGTH][MAX_TASK_SEQ_LENGTH], Routes2[MAX_SEG_TAG_LENGTH][MAX_TASK_SEQ_LENGTH];
	int XCLds[MAX_SEG_TAG_LENGTH], LeftTasks[MAX_TASK_SEQ_LENGTH], Positions[MAX_SEG_TAG_LENGTH];
	
	struct CandSelection
	{
		int RouteID;
		int Pos;
	};
	
	struct CandSelection CandSLCTList1[MAX_TASK_SEQ_LENGTH], CandSLCTList2[MAX_TASK_SEQ_LENGTH];
	int LLength1 = 0, LLength2 = 0;
	
	memcpy(S1, p1->sequence, sizeof(p1->sequence));
	memcpy(S2, p2->sequence, sizeof(p2->sequence));
		
	find_ele_positions(Positions, S1, 0);
	Routes1[0][0] = Positions[0]-1;
	for (i = 1; i < Positions[0]; i++)
	{
		copy_sub_array(Routes1[i], S1, Positions[i], Positions[i+1]);
	}
	
	find_ele_positions(Positions, S2, 0);
	Routes2[0][0] = Positions[0]-1;
	for (i = 1; i < Positions[0]; i++)
	{
		copy_sub_array(Routes2[i], S2, Positions[i], Positions[i+1]);
	}
	
	memcpy(XCLds, p1->route_seg_load, sizeof(p1->route_seg_load));
	
	for (i = 1; i <= Routes1[0][0]; i++)
	{
		for (j = 2; j < Routes1[i][0]; j++)
		{
			LLength1 ++;
			CandSLCTList1[LLength1].RouteID = i;
			CandSLCTList1[LLength1].Pos = j;
		}
	}
	
	for (i = 1; i <= Routes2[0][0]; i++)
	{
		for (j = 2; j < Routes2[i][0]; j++)
		{
			LLength2 ++;
			CandSLCTList2[LLength2].RouteID = i;
			CandSLCTList2[LLength2].Pos = j;
		}
	}
	
	k1 = rand_choose(LLength1);
	k2 = rand_choose(LLength2);
	
	copy_sub_array(SubPath1, Routes1[CandSLCTList1[k1].RouteID], 1, CandSLCTList1[k1].Pos);
	copy_sub_array(SubPath2, Routes2[CandSLCTList2[k2].RouteID], CandSLCTList2[k2].Pos, Routes2[CandSLCTList2[k2].RouteID][0]);	
	copy_sub_array(LeftTasks, Routes1[CandSLCTList1[k1].RouteID], CandSLCTList1[k1].Pos+1, Routes1[CandSLCTList1[k1].RouteID][0]-1);
	
	int Checked[MAX_TASK_SEQ_LENGTH];
	memset(Checked, 0, sizeof(Checked));
	
	for (i = 1; i < SubPath2[0]; i++)
	{
		if (Checked[i])
			continue;
		
		for (j = SubPath1[0]; j > 1; j--)
		{
			if (SubPath1[j] == SubPath2[i] || SubPath1[j] == inst_tasks[SubPath2[i]].inverse)
			{
				delete_element(SubPath1, j);
				Checked[i] = 1;
				break;
			}
		}
	}
	
	for (i = 1; i < SubPath2[0]; i++)
	{
		if (Checked[i])
			continue;
		
		for (j = LeftTasks[0]; j > 0; j--)
		{
			if (LeftTasks[j] == SubPath2[i] || LeftTasks[j] == inst_tasks[SubPath2[i]].inverse)
			{
				delete_element(LeftTasks, j);
				Checked[i] = 1;
				break;
			}
		}
	}
	
	for (i = 1; i < SubPath2[0]; i++)
	{
		if (Checked[i])
			continue;
		
		for (j = 1; j <= Routes1[0][0]; j++)
		{
			if (j == CandSLCTList1[k1].RouteID)
			//if (j == ChangedRouteID)
				continue;
			
			for (k = Routes1[j][0]; k > 1; k--)
			{
				if (Routes1[j][k] == SubPath2[i] || Routes1[j][k] == inst_tasks[SubPath2[i]].inverse)
				{
					delete_element(Routes1[j], k);
					XCLds[j] -= inst_tasks[SubPath2[i]].demand;
					Checked[i] = 1;
					break;
				}
			}
			
			if (Checked[i])
				break;
		}
	}
	
	link_array(SubPath1, SubPath2);
	memcpy(Routes1[CandSLCTList1[k1].RouteID], SubPath1, sizeof(SubPath1));
	XCLds[CandSLCTList1[k1].RouteID] = 0;
	for (i = 2; i < Routes1[CandSLCTList1[k1].RouteID][0]; i++)
	{
		XCLds[CandSLCTList1[k1].RouteID] += inst_tasks[Routes1[CandSLCTList1[k1].RouteID][i]].demand;
	}
	
	NO_LeftTasks = LeftTasks[0];
	
	// insert left tasks
	
	struct Insert
	{
		int InsertedTask;
		int InsertRouteID;
		int InsertPos;
		int InsertCost;
		int InsertVioLoad;
	};
	
	struct Insert CandInsertions[6000];
	int NO_CandInsertions;
	
	struct Insert ParetoSetInsertions[6000];
	int ParetoSetSize;
	int Out[6000], Add;
	
	struct Insert BestInsertion;

	int NO_CurrRoutes = 0;
	for (int j = 1; j <= Routes1[0][0]; j++)
	{
		if (Routes1[j][0] == 2)
			continue;

		NO_CurrRoutes ++;
	}
	
	for (n = 1; n <= NO_LeftTasks; n++)
	{
		NO_CandInsertions = 0;
		ParetoSetSize = 0;
	  
		for (j = 1; j <= Routes1[0][0]; j++)
		{
			if (Routes1[j][0] == 2)
				continue;

			if (XCLds[j] > capacity)
			{
				IVLoad = inst_tasks[LeftTasks[n]].demand;
			}
			else if (XCLds[j] > capacity-inst_tasks[LeftTasks[n]].demand)
			{
				IVLoad = XCLds[j]+inst_tasks[LeftTasks[n]].demand-capacity;
			}
			else 
			{
				IVLoad = 0;
			}
   		
			for (k = 2; k <= Routes1[j][0]; k++)
			{
				NO_CandInsertions ++;
				CandInsertions[NO_CandInsertions].InsertedTask = LeftTasks[n];
				CandInsertions[NO_CandInsertions].InsertRouteID = j;
				CandInsertions[NO_CandInsertions].InsertPos = k;
				CandInsertions[NO_CandInsertions].InsertCost = min_cost[inst_tasks[Routes1[j][k-1]].tail_node][inst_tasks[LeftTasks[n]].head_node]+
					min_cost[inst_tasks[LeftTasks[n]].tail_node][inst_tasks[Routes1[j][k]].head_node]-
					min_cost[inst_tasks[Routes1[j][k-1]].tail_node][inst_tasks[Routes1[j][k]].head_node];
				CandInsertions[NO_CandInsertions].InsertVioLoad = IVLoad;
   			
				Out[0] = 0;
				Add = 1;
   			
				for (m = 1; m <= ParetoSetSize; m++)
				{
					if (CandInsertions[NO_CandInsertions].InsertCost > ParetoSetInsertions[m].InsertCost
						&& CandInsertions[NO_CandInsertions].InsertVioLoad > ParetoSetInsertions[m].InsertVioLoad)
					{
						Add = 0;
						break;
					}
					else if (CandInsertions[NO_CandInsertions].InsertCost < ParetoSetInsertions[m].InsertCost
						&& CandInsertions[NO_CandInsertions].InsertVioLoad < ParetoSetInsertions[m].InsertVioLoad)
					{
						Out[0] ++;
						Out[Out[0]] = m;
					}
				}
   			
				if (Add)
				{
					for (m = Out[0]; m > 0; m--)
					{
						for (l = Out[m]; l < ParetoSetSize; l++)
						{
							ParetoSetInsertions[l] = ParetoSetInsertions[l+1];
						}
						ParetoSetSize --;
					}
     			
					ParetoSetSize ++;
					ParetoSetInsertions[ParetoSetSize] = CandInsertions[NO_CandInsertions];
				}
  		    	  
				w = inst_tasks[LeftTasks[n]].inverse;

				NO_CandInsertions ++;
				CandInsertions[NO_CandInsertions].InsertedTask = w;
				CandInsertions[NO_CandInsertions].InsertRouteID = j;
				CandInsertions[NO_CandInsertions].InsertPos = k;
				CandInsertions[NO_CandInsertions].InsertCost = min_cost[inst_tasks[Routes1[j][k-1]].tail_node][inst_tasks[w].head_node]
				+min_cost[inst_tasks[w].tail_node][inst_tasks[Routes1[j][k]].head_node]
				-min_cost[inst_tasks[Routes1[j][k-1]].tail_node][inst_tasks[Routes1[j][k]].head_node];
		    
				CandInsertions[NO_CandInsertions].InsertVioLoad = IVLoad;
   			
				Out[0] = 0;
				Add = 1;
   			
				for (m = 1; m <= ParetoSetSize; m++)
				{
					if (CandInsertions[NO_CandInsertions].InsertCost > ParetoSetInsertions[m].InsertCost
						&& CandInsertions[NO_CandInsertions].InsertVioLoad > ParetoSetInsertions[m].InsertVioLoad)
					{
						Add = 0;
						break;
					}
					else if (CandInsertions[NO_CandInsertions].InsertCost < ParetoSetInsertions[m].InsertCost
						&& CandInsertions[NO_CandInsertions].InsertVioLoad < ParetoSetInsertions[m].InsertVioLoad)
					{
						Out[0] ++;
						Out[Out[0]] = m;
					}
				}
   			
				if (Add)
				{
					for (m = Out[0]; m > 0; m--)
					{
						for (l = Out[m]; l < ParetoSetSize; l++)
						{
							ParetoSetInsertions[l] = ParetoSetInsertions[l+1];
						}
						ParetoSetSize --;
					}
     			
					ParetoSetSize ++;
					ParetoSetInsertions[ParetoSetSize] = CandInsertions[NO_CandInsertions];
				}
			}
	  }
	  
	  //if (NO_CurrRoutes < vehicle_num)
	  //{
		  NO_CandInsertions ++;
		  CandInsertions[NO_CandInsertions].InsertedTask = LeftTasks[n];
		  CandInsertions[NO_CandInsertions].InsertRouteID = 0;
		  CandInsertions[NO_CandInsertions].InsertPos = 2;
		  CandInsertions[NO_CandInsertions].InsertCost = min_cost[DEPOT][inst_tasks[LeftTasks[n]].head_node]+
			  min_cost[inst_tasks[LeftTasks[n]].tail_node][DEPOT];
		  CandInsertions[NO_CandInsertions].InsertVioLoad = 0;
	  //}
  	  
	  Out[0] = 0;
	  Add = 1;

	  for (m = 1; m <= ParetoSetSize; m++)
	  {
		  if (CandInsertions[NO_CandInsertions].InsertCost > ParetoSetInsertions[m].InsertCost
			  && CandInsertions[NO_CandInsertions].InsertVioLoad > ParetoSetInsertions[m].InsertVioLoad)
		  {
			  Add = 0;
			  break;
		  }
		  else if (CandInsertions[NO_CandInsertions].InsertCost < ParetoSetInsertions[m].InsertCost
			  && CandInsertions[NO_CandInsertions].InsertVioLoad < ParetoSetInsertions[m].InsertVioLoad)
		  {
			  Out[0] ++;
			  Out[Out[0]] = m;
		  }
	  }
   		
	  if (Add)
	  {
		  for (m = Out[0]; m > 0; m--)
		  {
			  for (l = Out[m]; l < ParetoSetSize; l++)
			  {
				  ParetoSetInsertions[l] = ParetoSetInsertions[l+1];
			  }
			  ParetoSetSize --;
		  }
    			
		  ParetoSetSize ++;
		  ParetoSetInsertions[ParetoSetSize] = CandInsertions[NO_CandInsertions];
	  }
		
	  k = rand_choose(ParetoSetSize);
	  BestInsertion = ParetoSetInsertions[k];
		
	  if (BestInsertion.InsertRouteID == 0)
	  {
		  Routes1[0][0] ++;
		  Routes1[Routes1[0][0]][0] = 3;
		  Routes1[Routes1[0][0]][1] = 0;
		  Routes1[Routes1[0][0]][2] = BestInsertion.InsertedTask;
		  Routes1[Routes1[0][0]][3] = 0;
			
		  XCLds[0] ++;
		  XCLds[XCLds[0]] = inst_tasks[BestInsertion.InsertedTask].demand;
	  }
	  else
	  {
		  add_element(Routes1[BestInsertion.InsertRouteID], BestInsertion.InsertedTask, BestInsertion.InsertPos);
		  XCLds[BestInsertion.InsertRouteID] += inst_tasks[BestInsertion.InsertedTask].demand;
	  }
		
		/*for (i = 1; i <= LeftTasks[0]; i++)
		{
			if (LeftTasks[i] == BestInsertion.InsertedTask || LeftTasks[i] == inst_tasks[BestInsertion.InsertedTask].Inv)
				break;
		}
		delete_element(LeftTasks, i);*/
	}
	
	//if (Prt)
  //	printf("add ok\n");
	
	xed_child->sequence[0] = 1;
	for (i = 1; i <= Routes1[0][0]; i++)
	{
		if (Routes1[i][0] == 2)
			continue;
		
		xed_child->sequence[0] --;
		link_array(xed_child->sequence, Routes1[i]);
	}
	
	xed_child->total_cost = get_task_seq_total_cost(xed_child->sequence, inst_tasks);
	memcpy(xed_child->route_seg_load, XCLds, sizeof(XCLds));
	
	for (i = xed_child->route_seg_load[0]; i > 0; i--)
	{
		if (xed_child->route_seg_load[i] == 0)
			delete_element(xed_child->route_seg_load, i);
	}
	
	xed_child->total_vio_load = get_total_vio_load(xed_child->route_seg_load);
	get_route_seg_length(xed_child->route_seg_length, xed_child->sequence, inst_tasks);
	xed_child->max_length = max(xed_child->route_seg_length);
	
	int RouteLoad;
 	find_ele_positions(Positions, xed_child->sequence, 0);
 	for (i = 1; i < Positions[0]; i++)
 	{
 		RouteLoad = 0;
 		for (j = Positions[i]; j < Positions[i+1]; j++)
 		{
 			RouteLoad += inst_tasks[xed_child->sequence[j]].demand;
 		}
 		
 		if (RouteLoad != xed_child->route_seg_load[i])
 		{
 			printf("XChild Seq\n");
     	for (k = 1; k <= xed_child->sequence[0]; k++)
 	    {
     		printf("%d  ", xed_child->sequence[k]);
 	    }
     	printf("\n");
     	
     	printf("i = %d, RouteLoad = %d\n", i, RouteLoad);
     	
     	printf("XChild Loads\n");
     	for (k = 1; k <= xed_child->route_seg_load[0]; k++)
 	    {
     		printf("%d  ", xed_child->route_seg_load[k]);
 	    }
     	printf("\n");
     	exit(0);
 		}
 	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void lns_mut(individual *c, individual *p, const task *inst_tasks, double coef)
{
	indi_copy(c, p);
	
	//double coef;

	int imp = 1;
	while (imp)
	{
		imp = 0;

		if (c->total_vio_load == 0 && c->total_cost < min_f1)
			min_f1 = c->total_cost;
		if (c->total_vio_load == 0 && c->total_cost > max_f1)
			max_f1 = c->total_cost;
		if (c->total_vio_load == 0 && c->max_length < min_f2)
			min_f2 = c->max_length;
		if (c->total_vio_load == 0 && c->max_length > max_f2)
			max_f2 = c->max_length;

		norm_f1 = 1.0*(c->total_cost-min_f1)/(max_f1-min_f1);
		norm_f2 = 1.0*(c->max_length-min_f2)/(max_f2-min_f2);
		//norm_f1 = c->total_cost;
		//norm_f2 = c->max_length;
		c->fitness = coef*norm_f1+(1-coef)*norm_f2;

		individual tmp_indi;
		indi_copy(&tmp_indi, c);

		//printf("coef = %.20lf\n", coef);
		//printf("before:   tvl = %d, totalcost = %d, maxlength = %d, fit = %.20lf\n", c->total_vio_load, c->total_cost, c->max_length, c->fitness);
		//printf("tmp_indi: tvl = %d, totalcost = %d, maxlength = %d, fit = %.20lf\n", tmp_indi.total_vio_load, tmp_indi.total_cost, tmp_indi.max_length, tmp_indi.fitness);
		//printf("min_f1 = %d, min_f2 = %d, max_f1 = %d, max_f2 = %d\n", min_f1, min_f2, max_f1, max_f2);
		//print_one_dim_array(tmp_indi.sequence);

		lns(c, coef, 1, inst_tasks);

		//printf("after:    tvl = %d, totalcost = %d, maxlength = %d, fit = %.20lf\n", c->total_vio_load, c->total_cost, c->max_length, c->fitness);
		//printf("min_f1 = %d, min_f2 = %d, max_f1 = %d, max_f2 = %d\n", min_f1, min_f2, max_f1, max_f2);
		//print_one_dim_array(c->sequence);

		//norm_f1 = 1.0*(tmp_indi.total_cost-min_f1)/(max_f1-min_f1);
		//norm_f2 = 1.0*(tmp_indi.max_length-min_f2)/(max_f2-min_f2);
		//tmp_indi.fitness = coef*norm_f1+(1-coef)*norm_f2;
		//norm_f1 = 1.0*(c->total_cost-min_f1)/(max_f1-min_f1);
		//norm_f2 = 1.0*(c->max_length-min_f2)/(max_f2-min_f2);
		//c->fitness = coef*norm_f1+(1-coef)*norm_f2;

		if (indi_fit_cmp(c, &tmp_indi) == -1)
			imp = 1;
	}
	//printf("1over\n");

	int qual = 0;
	//if (c->fitness < 0.05)
		qual = 1;

	//printf("cfit = %lf, qual = %d\n", c->fitness, qual);

	imp = 1;
	while (imp && qual)
	{
		imp = 0;

		if (c->total_vio_load == 0 && c->total_cost < min_f1)
			min_f1 = c->total_cost;
		if (c->total_vio_load == 0 && c->total_cost > max_f1)
			max_f1 = c->total_cost;
		if (c->total_vio_load == 0 && c->max_length < min_f2)
			min_f2 = c->max_length;
		if (c->total_vio_load == 0 && c->max_length > max_f2)
			max_f2 = c->max_length;

		norm_f1 = 1.0*(c->total_cost-min_f1)/(max_f1-min_f1);
		norm_f2 = 1.0*(c->max_length-min_f2)/(max_f2-min_f2);
		//norm_f1 = c->total_cost;
		//norm_f2 = c->max_length;
		c->fitness = coef*norm_f1+(1-coef)*norm_f2;

		individual tmp_indi;
		indi_copy(&tmp_indi, c);

		//printf("coef = %lf\n", coef);
		//printf("before:   tvl = %d, totalcost = %d, maxlength = %d, fit = %.20lf\n", c->total_vio_load, c->total_cost, c->max_length, c->fitness);
		//printf("tmp_indi: tvl = %d, totalcost = %d, maxlength = %d, fit = %.20lf\n", tmp_indi.total_vio_load, tmp_indi.total_cost, tmp_indi.max_length, tmp_indi.fitness);
		//printf("min_f1 = %d, min_f2 = %d, max_f1 = %d, max_f2 = %d\n", min_f1, min_f2, max_f1, max_f2);

		lns(c, coef, 2, inst_tasks);

		//printf("after:    tvl = %d, totalcost = %d, maxlength = %d, fit = %.20lf\n", c->total_vio_load, c->total_cost, c->max_length, c->fitness);
		//printf("min_f1 = %d, min_f2 = %d, max_f1 = %d, max_f2 = %d\n", min_f1, min_f2, max_f1, max_f2);
		//print_one_dim_array(c->sequence);

		//norm_f1 = 1.0*(tmp_indi.total_cost-min_f1)/(max_f1-min_f1);
		//norm_f2 = 1.0*(tmp_indi.max_length-min_f2)/(max_f2-min_f2);
		//tmp_indi.fitness = coef*norm_f1+(1-coef)*norm_f2;
		//norm_f1 = 1.0*(c->total_cost-min_f1)/(max_f1-min_f1);
		//norm_f2 = 1.0*(c->max_length-min_f2)/(max_f2-min_f2);
		//c->fitness = coef*norm_f1+(1-coef)*norm_f2;

		if (indi_fit_cmp(c, &tmp_indi) == -1)
			imp = 1;
	}
	//printf("2over\n");

	imp = 1;
	while (imp && qual)
	{
		imp = 0;

		if (c->total_vio_load == 0 && c->total_cost < min_f1)
			min_f1 = c->total_cost;
		if (c->total_vio_load == 0 && c->total_cost > max_f1)
			max_f1 = c->total_cost;
		if (c->total_vio_load == 0 && c->max_length < min_f2)
			min_f2 = c->max_length;
		if (c->total_vio_load == 0 && c->max_length > max_f2)
			max_f2 = c->max_length;

		norm_f1 = 1.0*(c->total_cost-min_f1)/(max_f1-min_f1);
		norm_f2 = 1.0*(c->max_length-min_f2)/(max_f2-min_f2);
		//norm_f1 = c->total_cost;
		//norm_f2 = c->max_length;
		c->fitness = coef*norm_f1+(1-coef)*norm_f2;

		individual tmp_indi;
		indi_copy(&tmp_indi, c);

		//printf("before: totalcost = %d, maxlength = %d, min_f1 = %d, min_f2 = %d, fit = %lf\n", c->total_cost, c->max_length, min_f1, min_f2, c->fitness);
		//printf("max_f1 = %d, min_f1 = %d, max_f2 = %d, min_f2 = %d\n", max_f1, min_f1, max_f2, min_f2);

		lns(c, coef, 1, inst_tasks);

		//printf("after: totalcost = %d, maxlength = %d, min_f1 = %d, min_f2 = %d, fit = %lf\n", c->total_cost, c->max_length, min_f1, min_f2, c->fitness);
		//printf("max_f1 = %d, min_f1 = %d, max_f2 = %d, min_f2 = %d\n", max_f1, min_f1, max_f2, min_f2);
		//printf("%d\n", indi_fit_cmp(c, &tmp_indi));

		//norm_f1 = 1.0*(tmp_indi.total_cost-min_f1)/(max_f1-min_f1);
		//norm_f2 = 1.0*(tmp_indi.max_length-min_f2)/(max_f2-min_f2);
		//tmp_indi.fitness = coef*norm_f1+(1-coef)*norm_f2;
		//norm_f1 = 1.0*(c->total_cost-min_f1)/(max_f1-min_f1);
		//norm_f2 = 1.0*(c->max_length-min_f2)/(max_f2-min_f2);
		//c->fitness = coef*norm_f1+(1-coef)*norm_f2;

		if (indi_fit_cmp(c, &tmp_indi) == -1)
			imp = 1;
	}
	//printf("3over\n");

	/*if (c->total_vio_load == 0)
	{
		int positions[MAX_SEG_TAG_LENGTH];

		find_ele_positions(positions, c->sequence, 0);
		for (int i = 1; i < positions[0]; i++)
		{
			int route[MAX_TASK_ROUTE_LENGTH];
			int fh_route1[MAX_TASK_ROUTE_LENGTH], fh_route2[MAX_TASK_ROUTE_LENGTH], fh_route[MAX_TASK_ROUTE_LENGTH];
			copy_sub_array(route, c->sequence, positions[i], positions[i+1]);
			fred_heuristic(fh_route1, route, inst_tasks, 0);
			fred_heuristic(fh_route2, route, inst_tasks, 1);
			int cost1 = get_task_seq_total_cost(fh_route1, inst_tasks);
			int cost2 = get_task_seq_total_cost(fh_route2, inst_tasks);
			int cost3;
			if (cost1 < cost2)
			{
				memcpy(fh_route, fh_route1, sizeof(fh_route1));
				cost3 = cost1;
			}
			else
			{
				memcpy(fh_route, fh_route2, sizeof(fh_route2));
				cost3 = cost2;
			}
			int cost4 = get_task_seq_total_cost(route, inst_tasks);
			if (cost3 < cost4)
			{
				for (int j = positions[i]; j < positions[i+1]; j++)
				{
					c->sequence[j] = fh_route[j-positions[i]+1];
				}
				c->total_cost += cost3-cost4;
				c->route_seg_length[i] += cost3-cost4;
			}
		}
		c->max_length = max(c->route_seg_length);
	}*/
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void lns(individual *indi, double coef, int nsize, const task *inst_tasks)
{
	norm_f1 = 1.0*(indi->total_cost-min_f1)/(max_f1-min_f1);
	norm_f2 = 1.0*(indi->max_length-min_f2)/(max_f2-min_f2);
	//norm_f1 = indi->total_cost;
	//norm_f2 = indi->max_length;
	indi->fitness = coef*norm_f1+(1-coef)*norm_f2;

	if (nsize == 1) // traditional move operators, i.e., single insertion, double insertion, swap, etc.
	{
		move si_move, di_move, swap_move, next_move;
		next_move.add_vio_load = INF;
		next_move.fitness = INF;

		single_insertion(&si_move, indi, coef, inst_tasks);
		if (move_fit_cmp(&si_move, &next_move) == -1)
		{
			next_move = si_move;
		}
		if (mv_indi_fit_cmp(&next_move, indi) != -1)
		{
			double_insertion(&di_move, indi, coef, inst_tasks);
			if (move_fit_cmp(&di_move, &next_move) == -1)
			{
				next_move = di_move;
			}
			if (mv_indi_fit_cmp(&next_move, indi) != -1)
			{
				swap(&swap_move, indi, coef, inst_tasks);
				if (move_fit_cmp(&swap_move, &next_move) == -1)
				{
					next_move = swap_move;
				}
			}
		}

		//printf("in coef = %.20lf\n", coef);
		//printf("type = %d, task1 = %d, task2 = %d, orig_seg = %d, targ_seg = %d, orig_pos = %d, targ_pos = %d, add_vioload = %d, fitness = %.20lf\n",
		//	next_move.type, next_move.task1, next_move.task2, next_move.orig_seg, next_move.targ_seg, next_move.orig_pos, next_move.targ_pos,
		//	next_move.add_vio_load, next_move.fitness);

		//printf("move addvl = %d, totalcost = %d, maxlength = %d, fit = %.20lf\n", next_move.add_vio_load, next_move.total_cost, next_move.max_length, next_move.fitness);
		//printf("min_f1 = %d, min_f2 = %d, max_f1 = %d, max_f2 = %d\n", min_f1, min_f2, max_f1, max_f2);
		norm_f1 = 1.0*(next_move.total_cost-min_f1)/(max_f1-min_f1);
		norm_f2 = 1.0*(next_move.max_length-min_f2)/(max_f2-min_f2);
		//norm_f1 = next_move.total_cost;
		//norm_f2 = next_move.max_length;
		next_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
		//printf("norm_f1 = %lf, norm_f2 = %lf, fit = %.20lf\n", norm_f1, norm_f2, next_move.fitness);
		//printf("min_f1 = %d, min_f2 = %d, max_f1 = %d, max_f2 = %d\n", min_f1, min_f2, max_f1, max_f2);

		int orig_ptr, targ_ptr, seg_ptr1, seg_ptr2;
		orig_ptr = 0;
		targ_ptr = 0;
		seg_ptr1 = 0;
		seg_ptr2 = 0;
		for (int i = 1; i < indi->sequence[0]; i++)
		{
			if (indi->sequence[i] == 0)
			{
				if (seg_ptr1 < next_move.orig_seg)
					seg_ptr1 ++;
				if (seg_ptr2 < next_move.targ_seg)
					seg_ptr2 ++;
				if (seg_ptr1 == next_move.orig_seg && orig_ptr == 0)
					orig_ptr = i+next_move.orig_pos-1;
				if (seg_ptr2 == next_move.targ_seg && targ_ptr == 0)
					targ_ptr = i+next_move.targ_pos-1;
			}
			if (orig_ptr != 0 && targ_ptr != 0)
				break;
		}

		//printf("before\n");
		//print_one_dim_array(indi->sequence);
		//print_one_dim_array(indi->route_seg_length);
		//printf("totalcost = %d, maxlength = %d\n", indi->total_cost, indi->max_length);

		switch (next_move.type)
		{
		case SI:
			{
				delete_element(indi->sequence, orig_ptr);
				if (targ_ptr > orig_ptr)
					targ_ptr --;
				indi->route_seg_length[next_move.orig_seg] = next_move.orig_length;
				indi->route_seg_load[next_move.orig_seg] -= inst_tasks[next_move.task1].demand;
				if (next_move.targ_seg > indi->route_seg_length[0])
				{
					indi->sequence[0] ++;
					indi->sequence[indi->sequence[0]] = next_move.task1;
					indi->sequence[0] ++;
					indi->sequence[indi->sequence[0]] = 0;
					indi->route_seg_length[0] ++;
					indi->route_seg_length[indi->route_seg_length[0]] = next_move.targ_length;
					indi->route_seg_load[0] ++;
					indi->route_seg_load[indi->route_seg_load[0]] = inst_tasks[next_move.task1].demand;
				}
				else
				{
					add_element(indi->sequence, next_move.task1, targ_ptr);
					indi->route_seg_length[next_move.targ_seg] = next_move.targ_length;
					indi->route_seg_load[next_move.targ_seg] += inst_tasks[next_move.task1].demand;
				}
			}
			break;
		case DI:
			{
				delete_element(indi->sequence, orig_ptr+1);
				delete_element(indi->sequence, orig_ptr);
				if (targ_ptr > orig_ptr)
					targ_ptr -= 2;
				indi->route_seg_length[next_move.orig_seg] = next_move.orig_length;
				indi->route_seg_load[next_move.orig_seg] -= inst_tasks[next_move.task1].demand+inst_tasks[next_move.task2].demand;
				if (next_move.targ_seg > indi->route_seg_length[0])
				{
					indi->sequence[0] ++;
					indi->sequence[indi->sequence[0]] = next_move.task1;
					indi->sequence[0] ++;
					indi->sequence[indi->sequence[0]] = next_move.task2;
					indi->sequence[0] ++;
					indi->sequence[indi->sequence[0]] = 0;
					indi->route_seg_length[0] ++;
					indi->route_seg_length[indi->route_seg_length[0]] = next_move.targ_length;
					indi->route_seg_load[0] ++;
					indi->route_seg_load[indi->route_seg_load[0]] = inst_tasks[next_move.task1].demand+inst_tasks[next_move.task2].demand;
				}
				else
				{
					add_element(indi->sequence, next_move.task2, targ_ptr);
					add_element(indi->sequence, next_move.task1, targ_ptr);
					indi->route_seg_length[next_move.targ_seg] = next_move.targ_length;
					indi->route_seg_load[next_move.targ_seg] += inst_tasks[next_move.task1].demand+inst_tasks[next_move.task2].demand;
				}
			}
			break;
		case SWAP:
			{
				indi->sequence[targ_ptr] = next_move.task1;
				indi->sequence[orig_ptr] = next_move.task2;
				indi->route_seg_length[next_move.orig_seg] = next_move.orig_length;
				indi->route_seg_length[next_move.targ_seg] = next_move.targ_length;
				indi->route_seg_load[next_move.orig_seg] -= inst_tasks[next_move.task1].demand-inst_tasks[next_move.task2].demand;
				indi->route_seg_load[next_move.targ_seg] += inst_tasks[next_move.task1].demand-inst_tasks[next_move.task2].demand;
			}
			break;
		}
		indi->max_length = next_move.max_length;
		indi->total_cost = next_move.total_cost;
		indi->total_vio_load += next_move.add_vio_load;
		indi->fitness = next_move.fitness;

		if (indi->route_seg_load[next_move.orig_seg] == 0)
		{
			if (next_move.type == DI && next_move.orig_seg > next_move.targ_seg)
			{
				delete_element(indi->sequence, orig_ptr+1);
			}
			else
			{
				delete_element(indi->sequence, orig_ptr);
			}
			delete_element(indi->route_seg_length, next_move.orig_seg);
			delete_element(indi->route_seg_load, next_move.orig_seg);
		}

		//printf("after\n");
		//print_one_dim_array(indi->sequence);
		//print_one_dim_array(indi->route_seg_length);
		//printf("totalcost = %d, maxlength = %d\n", indi->total_cost, indi->max_length);
	}
	else
	{
		individual tmp_indi, next_indi;
		int task_routes[MAX_SEG_TAG_LENGTH][MAX_TASK_SEG_LENGTH];

		task_routes[0][0] = 1;
		task_routes[1][0] = 1;
		task_routes[1][1] = 0;
		for (int i = 2; i <= indi->sequence[0]; i++)
		{
			task_routes[task_routes[0][0]][0] ++;
			task_routes[task_routes[0][0]][task_routes[task_routes[0][0]][0]] = indi->sequence[i];

			if (indi->sequence[i] == 0 && i < indi->sequence[0])
			{
				task_routes[0][0] ++;
				task_routes[task_routes[0][0]][0] = 1;
				task_routes[task_routes[0][0]][1] = 0;
			}
		}

		int tmp_route_seg_length[MAX_SEG_TAG_LENGTH];
		memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));

		int tmp_routes[MAX_TASK_SEG_LENGTH], tmp;
		for (int i = 1; i < task_routes[0][0]; i++)
		{
			for (int j = i+1; j <= task_routes[0][0]; j++)
			{
				if (tmp_route_seg_length[j] > tmp_route_seg_length[i])
				{
					memcpy(tmp_routes, task_routes[i], sizeof(task_routes[i]));
					memcpy(task_routes[i], task_routes[j], sizeof(task_routes[j]));
					memcpy(task_routes[j], tmp_routes, sizeof(tmp_routes));
					tmp = tmp_route_seg_length[i];
					tmp_route_seg_length[i] = tmp_route_seg_length[j];
					tmp_route_seg_length[j] = tmp;
				}
			}
		}

		if (task_routes[0][0] < nsize)
			return;
		
		long long int ub_trial = 2;
		for (int i = 1; i < task_routes[0][0]; i++)
		{
			ub_trial *= 2;
		}
		//ub_trial = 2^task_routes[0][0];

		next_indi.total_vio_load = INF;

		int lns_routes[MAX_SEG_TAG_LENGTH];
		int maxcount = 100;
		int count = 0;
		for (int i = 0; i < ub_trial; i++)
		{
			lns_routes[0] = 0;
			for (int j = 1; j <= task_routes[0][0]; j++)
			{
				if (bit_one(i, j))
				{
					lns_routes[0] ++;
					lns_routes[lns_routes[0]] = j;
				}
			}

			if (lns_routes[0] != nsize)
				continue;

			count ++;

			int serve_mark[MAX_TASK_TAG_LENGTH];
			memset(serve_mark, 0, sizeof(serve_mark));
			serve_mark[0] = task_num;
			for (int j = 1; j <= lns_routes[0]; j++)
			{
				for (int k = 2; k < task_routes[lns_routes[j]][0]; k++)
				{
					serve_mark[task_routes[lns_routes[j]][k]] = 1;
					serve_mark[inst_tasks[task_routes[lns_routes[j]][k]].inverse] = 1;
				}
			}

			path_scanning(&tmp_indi, inst_tasks, serve_mark);

			for (int j = 1; j <= task_routes[0][0]; j++)
			{
				if (bit_one(i, j))
					continue;

				tmp_indi.sequence[0] --;
				link_array(tmp_indi.sequence, task_routes[j]);
				//tmp_indi.route_seg_load[0] ++;
				//tmp_indi.route_seg_load[tmp_indi.route_seg_load[0]] = indi->route_seg_load[j];
				//if (indi->route_seg_load[j] > capacity)
				//	tmp_indi.total_vio_load += indi->route_seg_load[j];
			}

			tmp_indi.total_cost = get_task_seq_total_cost(tmp_indi.sequence, inst_tasks);
			get_route_seg_length(tmp_indi.route_seg_length, tmp_indi.sequence, inst_tasks);
			tmp_indi.max_length = max(tmp_indi.route_seg_length);
			get_route_seg_load(tmp_indi.route_seg_load, tmp_indi.sequence, inst_tasks);
			tmp_indi.total_vio_load = get_total_vio_load(tmp_indi.route_seg_load);
			norm_f1 = 1.0*(tmp_indi.total_cost-min_f1)/(max_f1-min_f1);
			norm_f2 = 1.0*(tmp_indi.max_length-min_f2)/(max_f2-min_f2);
			//norm_f1 = tmp_indi.total_cost;
			//norm_f2 = tmp_indi.max_length;
			tmp_indi.fitness = coef*norm_f1+(1-coef)*norm_f2;

			if (indi_fit_cmp(&tmp_indi, &next_indi) == -1)
				indi_copy(&next_indi, &tmp_indi);

			norm_f1 = 1.0*(indi->total_cost-min_f1)/(max_f1-min_f1);
			norm_f2 = 1.0*(indi->max_length-min_f2)/(max_f2-min_f2);
			//norm_f1 = indi->total_cost;
			//norm_f2 = indi->max_length;
			indi->fitness = coef*norm_f1+(1-coef)*norm_f2;

			if (indi_fit_cmp(&next_indi, indi) == -1)
				break;

			//printf("count = %d\n", count);
			if (count >= maxcount)
				break;
		}

		if (indi_fit_cmp(&next_indi, indi) == -1)
			indi_copy(indi, &next_indi);

		norm_f1 = 1.0*(indi->total_cost-min_f1)/(max_f1-min_f1);
		norm_f2 = 1.0*(indi->max_length-min_f2)/(max_f2-min_f2);
		//norm_f1 = indi->total_cost;
		//norm_f2 = indi->max_length;
		indi->fitness = coef*norm_f1+(1-coef)*norm_f2;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void single_insertion(move *best_move, individual *indi, double coef, const task *inst_tasks)
{
	best_move->type = SI;
	best_move->add_vio_load = INF;

	int tmp_route_seg_length[MAX_SEG_TAG_LENGTH];

	int task_routes[MAX_SEG_TAG_LENGTH][MAX_TASK_SEG_LENGTH];
	task_routes[0][0] = 1;
	task_routes[1][0] = 1;
	task_routes[1][1] = 0;
	for (int i = 2; i <= indi->sequence[0]; i++)
	{
		task_routes[task_routes[0][0]][0] ++;
		task_routes[task_routes[0][0]][task_routes[task_routes[0][0]][0]] = indi->sequence[i];

		if (indi->sequence[i] == 0 && i < indi->sequence[0])
		{
			task_routes[0][0] ++;
			task_routes[task_routes[0][0]][0] = 1;
			task_routes[task_routes[0][0]][1] = 0;
		}
	}

	move tmp_move;
	for (int s1 = 1; s1 <= task_routes[0][0]; s1++)
	{
		tmp_move.orig_seg = s1;
		for (int i = 2; i < task_routes[s1][0]; i++)
		{
			tmp_move.orig_pos = i;
			for (int s2 = 1; s2 <= task_routes[0][0]+1; s2++) /* s2 > task_routes[0][0] --> create a new route */
			{
				if (s2 == s1)
					continue;

				tmp_move.targ_seg = s2;
				if (s2 > task_routes[0][0]/* && task_routes[0][0] < vehicle_num*/)
				{
					tmp_move.add_vio_load = 0;
					if (indi->route_seg_load[s1] > capacity)
						tmp_move.add_vio_load -= indi->route_seg_load[s1]-capacity;
					if (indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand > capacity)
						tmp_move.add_vio_load += indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand-capacity;

					tmp_move.task1 = task_routes[s1][i];
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] +=
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					tmp_route_seg_length[0] ++;
					tmp_route_seg_length[tmp_route_seg_length[0]] =
						min_cost[DEPOT][inst_tasks[tmp_move.task1].head_node]+inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][DEPOT];

					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[tmp_route_seg_length[0]];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[tmp_route_seg_length[0]]-indi->route_seg_length[s1];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}
					continue;
				}

				for (int j = 2; j <= task_routes[s2][0]; j++)
				{
					tmp_move.targ_pos = j;
					tmp_move.add_vio_load = 0;
					if (indi->route_seg_load[s1] > capacity)
						tmp_move.add_vio_load -= indi->route_seg_load[s1]-capacity;
					if (indi->route_seg_load[s2] > capacity)
						tmp_move.add_vio_load -= indi->route_seg_load[s2]-capacity;
					if (indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand > capacity)
						tmp_move.add_vio_load += indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand-capacity;
					if (indi->route_seg_load[s2]+inst_tasks[task_routes[s1][i]].demand > capacity)
						tmp_move.add_vio_load += indi->route_seg_load[s2]+inst_tasks[task_routes[s1][i]].demand-capacity;

					tmp_move.task1 = task_routes[s1][i];
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] += 
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					tmp_route_seg_length[s2] +=
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];
						
					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[s2];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
					//tmp_move.fitness = coef*tmp_move.total_cost+(1-coef)*tmp_move.max_length;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] += 
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					tmp_route_seg_length[s2] +=
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];
						
					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[s2];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
					//tmp_move.fitness = coef*tmp_move.total_cost+(1-coef)*tmp_move.max_length;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}
				}
			}
		}
	}

	//printf("task = %d, orig_route = %d, targ_route = %d, orig_load = %d, tar_load = %d, add_cost = %d, add_vioload = %d, add_fit = %lf\n", best_move->task1,
	//	best_move->orig_route, best_move->targ_route, best_move->orig_load, best_move->targ_load, best_move->add_cost, best_move->add_vio_load,
	//	best_move->add_fit);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void double_insertion(move *best_move, individual *indi, double coef, const task *inst_tasks)
{
	best_move->type = DI;
	best_move->add_vio_load = INF;

	int tmp_route_seg_length[MAX_SEG_TAG_LENGTH];

	int task_routes[MAX_SEG_TAG_LENGTH][MAX_TASK_SEG_LENGTH];
	task_routes[0][0] = 1;
	task_routes[1][0] = 1;
	task_routes[1][1] = 0;
	for (int i = 2; i <= indi->sequence[0]; i++)
	{
		task_routes[task_routes[0][0]][0] ++;
		task_routes[task_routes[0][0]][task_routes[task_routes[0][0]][0]] = indi->sequence[i];

		if (indi->sequence[i] == 0 && i < indi->sequence[0])
		{
			task_routes[0][0] ++;
			task_routes[task_routes[0][0]][0] = 1;
			task_routes[task_routes[0][0]][1] = 0;
		}
	}

	move tmp_move;
	for (int s1 = 1; s1 <= task_routes[0][0]; s1++)
	{
		if (task_routes[s1][0] < 4)
			continue;

		tmp_move.orig_seg = s1;
		for (int i = 2; i < task_routes[s1][0]-1; i++)
		{
			tmp_move.orig_pos = i;
			for (int s2 = 1; s2 <= task_routes[0][0]+1; s2++) /* s2 > task_routes[0][0] --> create a new route */
			{
				if (s2 == s1)
					continue;

				tmp_move.targ_seg = s2;
				if (s2 > task_routes[0][0]/* && task_routes[0][0] < vehicle_num*/)
				{
					if (task_routes[s1][0] > 4)
						continue;

					tmp_move.add_vio_load = 0;
					if (indi->route_seg_load[s1] > capacity)
						tmp_move.add_vio_load -= indi->route_seg_load[s1]-capacity;
					if (indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s1][i+1]].demand > capacity)
						tmp_move.add_vio_load += indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s1][i+1]].demand-capacity;

					tmp_move.task1 = task_routes[s1][i];
					tmp_move.task2 = task_routes[s1][i+1];
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] +=
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						inst_tasks[task_routes[s1][i+1]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node];
					tmp_route_seg_length[0] ++;
					tmp_route_seg_length[tmp_route_seg_length[0]] =
						min_cost[DEPOT][inst_tasks[tmp_move.task1].head_node]+inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+inst_tasks[tmp_move.task2].serv_cost+
						min_cost[inst_tasks[tmp_move.task2].tail_node][DEPOT];

					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[tmp_route_seg_length[0]];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[tmp_route_seg_length[0]]-indi->route_seg_length[s1];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
					//tmp_move.fitness = coef*tmp_move.total_cost+(1-coef)*tmp_move.max_length;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					tmp_move.task1 = task_routes[s1][i];
					tmp_move.task2 = inst_tasks[task_routes[s1][i+1]].inverse;
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] +=
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						inst_tasks[task_routes[s1][i+1]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node];
					tmp_route_seg_length[0] ++;
					tmp_route_seg_length[tmp_route_seg_length[0]] =
						min_cost[DEPOT][inst_tasks[tmp_move.task1].head_node]+inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+inst_tasks[tmp_move.task2].serv_cost+
						min_cost[inst_tasks[tmp_move.task2].tail_node][DEPOT];

					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[tmp_route_seg_length[0]];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[tmp_route_seg_length[0]]-indi->route_seg_length[s1];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
					//tmp_move.fitness = coef*tmp_move.total_cost+(1-coef)*tmp_move.max_length;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}
					continue;
				}

				for (int j = 2; j <= task_routes[s2][0]; j++) 
				{
					tmp_move.targ_pos = j;
					tmp_move.add_vio_load = 0;
					if (indi->route_seg_load[s1] > capacity)
						tmp_move.add_vio_load -= indi->route_seg_load[s1]-capacity;
					if (indi->route_seg_load[s2] > capacity)
						tmp_move.add_vio_load -= indi->route_seg_load[s2]-capacity;
					if (indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s1][i+1]].demand > capacity)
						tmp_move.add_vio_load += indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s1][i+1]].demand-capacity;
					if (indi->route_seg_load[s2]+inst_tasks[task_routes[s1][i]].demand+inst_tasks[task_routes[s1][i+1]].demand > capacity)
						tmp_move.add_vio_load += indi->route_seg_load[s2]+inst_tasks[task_routes[s1][i]].demand+inst_tasks[task_routes[s1][i+1]].demand-capacity;

					tmp_move.task1 = task_routes[s1][i];
					tmp_move.task2 = task_routes[s1][i+1];
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] += 
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						inst_tasks[task_routes[s1][i+1]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node];
					tmp_route_seg_length[s2] +=
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
						inst_tasks[tmp_move.task2].serv_cost+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];
						
					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[s2];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
					//tmp_move.fitness = coef*tmp_move.total_cost+(1-coef)*tmp_move.max_length;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					////

					tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
					tmp_move.task2 = task_routes[s1][i+1];
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] += 
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						inst_tasks[task_routes[s1][i+1]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node];
					tmp_route_seg_length[s2] +=
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
						inst_tasks[tmp_move.task2].serv_cost+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];
						
					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[s2];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
					//tmp_move.fitness = coef*tmp_move.total_cost+(1-coef)*tmp_move.max_length;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					////

					tmp_move.task1 = task_routes[s1][i];
					tmp_move.task2 = inst_tasks[task_routes[s1][i+1]].inverse;
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] += 
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						inst_tasks[task_routes[s1][i+1]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node];
					tmp_route_seg_length[s2] +=
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
						inst_tasks[tmp_move.task2].serv_cost+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];
						
					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[s2];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
					//tmp_move.fitness = coef*tmp_move.total_cost+(1-coef)*tmp_move.max_length;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					////

					tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
					tmp_move.task2 = inst_tasks[task_routes[s1][i+1]].inverse;
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] += 
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						inst_tasks[task_routes[s1][i+1]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i+1]].tail_node][inst_tasks[task_routes[s1][i+2]].head_node];
					tmp_route_seg_length[s2] +=
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[tmp_move.task2].head_node]+
						inst_tasks[tmp_move.task2].serv_cost+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node];
						
					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[s2];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
					//tmp_move.fitness = coef*tmp_move.total_cost+(1-coef)*tmp_move.max_length;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}
				}
			}
		}
	}

	//printf("task1 = %d, task2 = %d, orig_route = %d, targ_route = %d, orig_load = %d, tar_load = %d, add_cost = %d, add_vioload = %d, add_fit = %lf\n",
	//	best_move->task1, best_move->task2, best_move->orig_route, best_move->targ_route, best_move->orig_load, best_move->targ_load, best_move->add_cost,
	//	best_move->add_vio_load, best_move->add_fit);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void swap(move *best_move, individual *indi, double coef, const task *inst_tasks)
{
	best_move->type = SWAP;
	best_move->add_vio_load = INF;

	int tmp_route_seg_length[MAX_SEG_TAG_LENGTH];

	int task_routes[MAX_SEG_TAG_LENGTH][MAX_TASK_SEG_LENGTH];
	task_routes[0][0] = 1;
	task_routes[1][0] = 1;
	task_routes[1][1] = 0;
	for (int i = 2; i <= indi->sequence[0]; i++)
	{
		task_routes[task_routes[0][0]][0] ++;
		task_routes[task_routes[0][0]][task_routes[task_routes[0][0]][0]] = indi->sequence[i];

		if (indi->sequence[i] == 0 && i < indi->sequence[0])
		{
			task_routes[0][0] ++;
			task_routes[task_routes[0][0]][0] = 1;
			task_routes[task_routes[0][0]][1] = 0;
		}
	}

	move tmp_move;
	for (int s1 = 1; s1 < task_routes[0][0]; s1++)
	{
		tmp_move.orig_seg = s1;
		for (int i = 2; i < task_routes[s1][0]-1; i++)
		{
			tmp_move.orig_pos = i;
			for (int s2 = s1+1; s2 <= task_routes[0][0]; s2++)
			{
				tmp_move.targ_seg = s2;
				for (int j = 2; j < task_routes[s2][0]; j++)
				{
					tmp_move.targ_pos = j;

					tmp_move.add_vio_load = 0;
					if (indi->route_seg_load[s1] > capacity)
						tmp_move.add_vio_load -= indi->route_seg_load[s1]-capacity;
					if (indi->route_seg_load[s2] > capacity)
						tmp_move.add_vio_load -= indi->route_seg_load[s2]-capacity;
					if (indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand+inst_tasks[task_routes[s2][j]].demand > capacity)
						tmp_move.add_vio_load += indi->route_seg_load[s1]-inst_tasks[task_routes[s1][i]].demand+inst_tasks[task_routes[s2][j]].demand-capacity;
					if (indi->route_seg_load[s2]+inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s2][j]].demand > capacity)
						tmp_move.add_vio_load += indi->route_seg_load[s2]+inst_tasks[task_routes[s1][i]].demand-inst_tasks[task_routes[s2][j]].demand-capacity;

					tmp_move.task1 = task_routes[s1][i];
					tmp_move.task2 = task_routes[s2][j];
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] += 
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
						inst_tasks[tmp_move.task2].serv_cost+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					tmp_route_seg_length[s2] +=
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						inst_tasks[task_routes[s2][j]].serv_cost-
						min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];
						
					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[s2];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
					//tmp_move.fitness = coef*tmp_move.total_cost+(1-coef)*tmp_move.max_length;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					////

					tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
					tmp_move.task2 = task_routes[s2][j];
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] += 
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
						inst_tasks[tmp_move.task2].serv_cost+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					tmp_route_seg_length[s2] +=
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						inst_tasks[task_routes[s2][j]].serv_cost-
						min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];
						
					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[s2];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
					//tmp_move.fitness = coef*tmp_move.total_cost+(1-coef)*tmp_move.max_length;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					////

					tmp_move.task1 = task_routes[s1][i];
					tmp_move.task2 = inst_tasks[task_routes[s2][j]].inverse;
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] += 
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
						inst_tasks[tmp_move.task2].serv_cost+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					tmp_route_seg_length[s2] +=
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						inst_tasks[task_routes[s2][j]].serv_cost-
						min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];
						
					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[s2];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
					//tmp_move.fitness = coef*tmp_move.total_cost+(1-coef)*tmp_move.max_length;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}

					////

					tmp_move.task1 = inst_tasks[task_routes[s1][i]].inverse;
					tmp_move.task2 = inst_tasks[task_routes[s2][j]].inverse;
					memcpy(tmp_route_seg_length, indi->route_seg_length, sizeof(indi->route_seg_length));
					tmp_route_seg_length[s1] += 
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[tmp_move.task2].head_node]+
						inst_tasks[tmp_move.task2].serv_cost+
						min_cost[inst_tasks[tmp_move.task2].tail_node][inst_tasks[task_routes[s1][i+1]].head_node]-
						min_cost[inst_tasks[task_routes[s1][i-1]].tail_node][inst_tasks[task_routes[s1][i]].head_node]-
						inst_tasks[task_routes[s1][i]].serv_cost-
						min_cost[inst_tasks[task_routes[s1][i]].tail_node][inst_tasks[task_routes[s1][i+1]].head_node];
					tmp_route_seg_length[s2] +=
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[tmp_move.task1].head_node]+
						inst_tasks[tmp_move.task1].serv_cost+
						min_cost[inst_tasks[tmp_move.task1].tail_node][inst_tasks[task_routes[s2][j+1]].head_node]-
						min_cost[inst_tasks[task_routes[s2][j-1]].tail_node][inst_tasks[task_routes[s2][j]].head_node]-
						inst_tasks[task_routes[s2][j]].serv_cost-
						min_cost[inst_tasks[task_routes[s2][j]].tail_node][inst_tasks[task_routes[s2][j+1]].head_node];
						
					tmp_move.max_length = max(tmp_route_seg_length);

					tmp_move.orig_length = tmp_route_seg_length[s1];
					tmp_move.targ_length = tmp_route_seg_length[s2];
					tmp_move.total_cost = indi->total_cost+
						tmp_route_seg_length[s1]+tmp_route_seg_length[s2]-indi->route_seg_length[s1]-indi->route_seg_length[s2];

					norm_f1 = 1.0*(tmp_move.total_cost-min_f1)/(max_f1-min_f1);
					norm_f2 = 1.0*(tmp_move.max_length-min_f2)/(max_f2-min_f2);
					//norm_f1 = tmp_move.total_cost;
					//norm_f2 = tmp_move.max_length;
					tmp_move.fitness = coef*norm_f1+(1-coef)*norm_f2;
					//tmp_move.fitness = coef*tmp_move.total_cost+(1-coef)*tmp_move.max_length;

					if (move_fit_cmp(&tmp_move, best_move) == -1)
					{
						best_move->task1 = tmp_move.task1;
						best_move->task2 = tmp_move.task2;
						best_move->orig_seg = tmp_move.orig_seg;
						best_move->targ_seg = tmp_move.targ_seg;
						best_move->orig_pos = tmp_move.orig_pos;
						best_move->targ_pos = tmp_move.targ_pos;
						best_move->orig_length = tmp_move.orig_length;
						best_move->targ_length = tmp_move.targ_length;
						best_move->total_cost = tmp_move.total_cost;
						best_move->add_vio_load = tmp_move.add_vio_load;
						best_move->max_length = tmp_move.max_length;
						best_move->fitness = tmp_move.fitness;
					}
				}
			}
		}
	}

	//printf("task1 = %d, task2 = %d, orig_route = %d, targ_route = %d, orig_load = %d, tar_load = %d, add_cost = %d, add_vioload = %d, add_fit = %lf\n",
	//	best_move->task1, best_move->task2, best_move->orig_route, best_move->targ_route, best_move->orig_load, best_move->targ_load, best_move->add_cost,
	//	best_move->add_vio_load, best_move->add_fit);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int dominated_num[MAX_TOTALSIZE];
int dominating_indis[MAX_TOTALSIZE][MAX_TOTALSIZE];
int dom[MAX_TOTALSIZE][MAX_TOTALSIZE];
int S[MAX_TOTALSIZE];
int R[MAX_TOTALSIZE];
double dis[MAX_TOTALSIZE][MAX_TOTALSIZE];
double fitness[MAX_TOTALSIZE];
double D[MAX_TOTALSIZE];
individual tmp_indi[MAX_TOTALSIZE];
int idx[MAX_TOTALSIZE];
bool fla[MAX_TOTALSIZE];
double d[MAX_TOTALSIZE];
void domination_sort(individual *pop, int tmp_popsize)
{
	//for (int i = 0; i < tmp_popsize; i++)
	//{
	//	pop[i].rank = tmp_popsize+1;
	//}
	
	memset(dominated_num, 0, sizeof(dominated_num));
	for (int i = 0; i < MAX_TOTALSIZE; i++)
	{
		dominating_indis[i][0] = 0;
	}

	for (int i = 0; i < tmp_popsize-1; i++)
	{
		for (int j = i+1; j < tmp_popsize; j++)
		{
			int cmp = indi_cmp(&pop[i], &pop[j]);
			if (cmp == -1)
			{
				dominated_num[j] ++;
				dominating_indis[i][0] ++;
				dominating_indis[i][dominating_indis[i][0]] = j;
			}
			else if (cmp == 1)
			{
				dominated_num[i] ++;
				dominating_indis[j][0] ++;
				dominating_indis[j][dominating_indis[j][0]] = i;
			}
		}
	}

	int trial = 0;
	for (int i = 0; i < tmp_popsize; i++)
	{
		if (dominated_num[i] == 0)
		{
			pop[i].rank = 1;
			trial ++;
		}
	}
	//printf("trial: %d \n", trial);
	int curr_rank = 1;
	while (trial < tmp_popsize)
	{
		for (int i = 0; i < tmp_popsize; i++)
		{
			if (pop[i].rank == curr_rank)
			{
				for (int j = 1; j <= dominating_indis[i][0]; j++)
				{
					dominated_num[dominating_indis[i][j]] --;
					if (dominated_num[dominating_indis[i][j]] == 0)
					{
						pop[dominating_indis[i][j]].rank = curr_rank+1;
						trial ++;
					}
				}
			}
		}
		curr_rank ++;
	}

	individual tmp_indi;
	for (int i = 0; i < tmp_popsize-1; i++)
	{
		for (int j = i+1; j < tmp_popsize; j++)
		{
			if (pop[j].rank < pop[i].rank)
			{
				indi_copy(&tmp_indi, &pop[i]);
				indi_copy(&pop[i], &pop[j]);
				indi_copy(&pop[j], &tmp_indi);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int dominate(individual *po,int i,int j)
{
	//printf("i: max_length is %d, total_cost is %d\n", po[i].max_length, po[i].total_cost);
	//printf("j: max_length is %d, total_cost is %d\n", po[j].max_length, po[j].total_cost);
	if(po[i].max_length < po[j].max_length && po[i].total_cost < po[j].total_cost )
		return 1;
	if(po[i].max_length > po[j].max_length && po[i].total_cost > po[j].total_cost)
		return -1;
	return 0;
}

int dominate2(individual *po,int i,int j)
{
	//printf("i: max_length is %d, total_cost is %d\n", po[i].max_length, po[i].total_cost);
	//printf("j: max_length is %d, total_cost is %d\n", po[j].max_length, po[j].total_cost);
	if(po[i].max_length < po[j].max_length && po[i].total_cost < po[j].total_cost && po[i].total_vio_load <= po[j].total_vio_load)
		return 1;
	if(po[i].max_length > po[j].max_length && po[i].total_cost > po[j].total_cost && po[i].total_vio_load > po[j].total_vio_load)
		return -1;
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void swap(double* a, double* b) 
{
    double temp = *a;
    *a = *b;
    *b = temp;
}
// a是数组首地址，l和r分别是两端的索引,要求l<=r
double Lomuto_partition(double a[], int l, int r) 
{
    double pivot = a[l];
    int j = l;

    for(int i = l + 1; i <= r; ++i)
    {
        if( a[i] < pivot )
        {
            swap( &a[++j], &a[i]);
        }
    }

    swap(&a[l], &a[j]);
    return j;   // j是下标
}


// 返回第k小的元素值。
// 注意：k不是下标，表示第k个
double __quick_select(double a[], int l, int r, int k) 
{
    // s是分裂点，中轴的下标
    int s = Lomuto_partition(a, l, r);

    if(s == l+k-1)
        return a[s];
    else if( s > (l+k-1) )// 处理左边的部分
        return __quick_select(a, l, s-1, k);
    else // 处理右边的部分
        return __quick_select(a, s+1, r, k-(s-l+1)); 
}

// 快速选择主函数
double quick_select_min(double a[], int len, int k) 
{
    return __quick_select(a, 0, len-1, k);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int rand_one_number(int T,int popsize,int idx)
{
	int start = idx - floor(T / 2);
	int finish = idx + floor(T / 2);

	if (start < 0)
	{
		start = 0;
		finish = start + T - 1;
	}
	if (finish > popsize - 1)
	{
		finish = popsize - 1;
		start = finish - T + 1;
	}
   // printf("%d %d\n", start, finish);
	start = start + rand() % T;
    
    return start;
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SPEA2(individual *pop, int tmp_popsize,int * F,int * MMM)
{	// F:1 : 直接 Truncation
	if(F == 0)
	{
		for (int i = 0; i < tmp_popsize; i++)
		{
			for (int j = 0; j < tmp_popsize; j++)
			{
				dom[i][j] = 0;
			}
		}

		for (int i = 0; i < tmp_popsize-1; i++)
		{
			for (int j = i+1; j < tmp_popsize; j++)
			{
				if(dominate(pop, i, j) == 1)
					dom[i][j] = 1;
				else if(dominate(pop, i, j) == -1)
					dom[j][i] = 1;
			}
		}
		// printf("dom: \n");
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	for (int j = 0; j < tmp_popsize; j++)
		// 	{
		// 		printf("%d ",dom[i][j]);
		// 	}
		// 	printf("\n");
		// }


		for (int i = 0; i < tmp_popsize;i++)
		{
			S[i] = 0;
			for (int j = 0; j < tmp_popsize;j++)
			{
				S[i] = S[i] + dom[i][j];
			}
		}

		// printf("S: \n");
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	printf("%d ",S[i]);
			
		// }
		// printf("\n");

		for (int i = 0; i < tmp_popsize;i++)
		{
			R[i] = 0;
			for (int j = 0; j < tmp_popsize;j++)
			{
				if(dom[j][i] == 1)
					R[i] = R[i] + S[j];
			}
		}

		// printf("R: \n");
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	printf("%d ",R[i]);
			
		// }
		// printf("\n");
		

		int maxfm = 0, minfm=INF,maxft = 0, minft=INF;
		for (int i = 0; i < tmp_popsize;i++)
		{	
			if(pop[i].total_cost < minft)
				minft = pop[i].total_cost;

			if(pop[i].total_cost > maxft)
				maxft = pop[i].total_cost;
			
			if(pop[i].max_length < minfm)
				minfm = pop[i].max_length;

			if(pop[i].max_length > maxfm)
				maxfm = pop[i].max_length;
		}
		// printf("scale: \n");
		// printf("total: %d, %d \n",minft,maxft);
		// printf("max: %d, %d \n",minfm,maxfm);
		// printf("f1 f2: \n");
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	printf("%d %d \n",pop[i].max_length, pop[i].total_cost);
			
		// }
		// printf("\n");


		for (int i = 0; i < tmp_popsize; i++)
		{
			for (int j = 0; j < tmp_popsize; j++)
			{
				if(i == j)
					dis[i][j] = INF;
				else
				{
					double to ,ma;
					// a1 = (pop[i].max_length - pop[j].max_length) * (pop[i].max_length - pop[j].max_length) / ((maxfm - minfm) * (maxfm - minfm));
					// a2 = (pop[i].total_cost - pop[j].total_cost) * (pop[i].total_cost - pop[j].total_cost) / ((maxft - minft) * (maxft - minft));
					
					// a1 =  float(abs(pop[i].max_length - pop[j].max_length)) ;
					// a2 = float(abs(pop[i].total_cost - pop[j].total_cost));
					//dis[i][j] = sqrt(a1*a1 + a2*a2);
					
					// a1 =  float(abs(pop[i].max_length - pop[j].max_length)) / (maxfm - minfm + 1) ;
					// a2 = float(abs(pop[i].total_cost - pop[j].total_cost)) / (maxft - minft + 1);
					// dis[i][j] = a1 + a2;

					// printf("abs(%d - %d) / (%d) = %f \n", pop[i].max_length, pop[j].max_length,maxfm-minfm,a1 );
					// printf("abs(%d - %d) / (%d - %d) = %f \n", pop[i].total_cost, pop[j].total_cost,maxft,minft,a2 );
					

					to = std::max(pop[i].total_cost, pop[j].total_cost);
					ma = std::max(pop[i].max_length, pop[j].max_length);
					to = float(abs(pop[i].total_cost - to)/(maxft-minft));
					ma = float(abs(pop[i].max_length - ma)/(maxfm-minfm));
					dis[i][j] = float( sqrt(to*to + ma*ma) );


					// char avalue[100] = "avalue.dat";
					// FILE *fp;
					// fp = fopen(avalue, "w");
					// for (int i = 0; i < popsize;i++)
					// {	
					// 	fprintf(fp, "%f %f \n", a1, a2);
					// }
					// fprintf(fp, "\n");
					// fclose(fp);

				}
			}
			//printf(" \n");
		}

		// printf("dist: \n");
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	for (int j = 0; j < tmp_popsize; j++)
		// 	{
		// 		printf("%3f ",dis[i][j]);
		// 	}
		// 	printf("\n");
		// }

		//int k = int(floor(sqrt(tmp_popsize)));
		//int k = 1;
		for (int i = 0; i < tmp_popsize;i++)
		{
			//D[i] = float(1/(2+quick_select_min(dis[i], tmp_popsize, k)));
			D[i] = float(1/(2+quick_select_min(dis[i], tmp_popsize, 1)));
		}

		for (int i = 0; i < tmp_popsize;i++)
		{
			fitness[i] = R[i] + D[i];
			fla[i] = 0;
		}

		// printf("fitness: \n");
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	printf("%d + %f = %f ",R[i], D[i],fitness[i]);
			
		// }
		// printf("\n");

		int count = 0;
		for (int i = 0; i < tmp_popsize; i++)
		{
			if(fitness[i] < 1)
			{
				count ++;
				fla[i] = 1;
			}
		}
		// printf("count: %d\n", count);
		// printf("count: %d\n", tmp_popsize/2);
		if(count < tmp_popsize/2 + 1 )
		//if(count )
		{ // 直接截断

			for (int i = 0; i < tmp_popsize;i++)
			{
				idx[i] = i;
			}

			// for (int i = 0; i < tmp_popsize;i++)
			// {
			// 	printf("%f\n", fitness[i]);
			// }

			for (int i = 0; i < tmp_popsize-1; i++)
			{
				for (int j = i+1; j < tmp_popsize; j++)
				{
					if(fitness[i] > fitness[j])
					{
						int chang;
						chang = idx[i];
						idx[i] = idx[j];
						idx[j] = chang;

						double chang2;
						chang2 = fitness[i];
						fitness[i] = fitness[j];
						fitness[j] = chang2;
					}
				}
			}

			// for (int i = 0; i < tmp_popsize;i++)
			// {
			// 	printf("%f  %d\n", fitness[i], idx[i]);
			// }

			for (int i = 0; i < tmp_popsize; i++)
			{
				indi_copy(&tmp_indi[i], &pop[idx[i]]);
			}

			for (int i = 0; i < tmp_popsize;i++)
			{
				indi_copy( &pop[i],&tmp_indi[i]);	
			}
		}
		else 
		{ // Truncation
			// for (int i = 0; i < tmp_popsize;i++)
			// {
			// 	printf("%d  %d\n", pop[i].max_length, pop[i].total_cost);
			// }
			count = count - tmp_popsize/2;
			//printf("number of Truncation: %d \n", count);
			for (int it = 0; it < count; it++)
			{
				int min_indx, comp = INF;
				for (int i = 0; i < tmp_popsize; i++)
				{   // 计算 fitness<1 中的矩阵的最小值 储存在d[i]
					d[i] = INF;
					for (int j = 0; j < tmp_popsize; j++)
					{
						if(fla[i] == 1 && dis[i][j] < d[i])
							d[i] = dis[i][j];
					}

					if(fla[i] == 1 && d[i] < comp )
					{
						comp = d[i];
						min_indx = i;
					}

				}
				fla[min_indx] = 0;
			}

			int su;
			for (int i = 0; i < tmp_popsize; i++)
			{
				su = su + fla[i];
			}

			int jj1 = 0,jj2 = tmp_popsize/2;
			for (int i = 0; i < tmp_popsize; i++)
			{ 
				if(fla[i] == 1 )
				{
					indi_copy(&tmp_indi[jj1], &pop[i]);
					jj1 = jj1 + 1;
				}
				else 
				{
					indi_copy(&tmp_indi[jj2], &pop[i]);
					jj2 = jj2 + 1;
				}
			}
			// printf("jj1: %d ,jj2: %d \n", jj1, jj2);
			// int su2 = 0;
			// for (int i = 0; i < tmp_popsize;i++)
			// 	su2 += fla[i];
			// printf("su2: %d \n", su2);
			for (int i = 0; i < tmp_popsize;i++)
			{
				indi_copy( &pop[i],&tmp_indi[i]);	
			}
		}
	}
	else
	{
		int count = 0;
		for (int i = 0; i < tmp_popsize;i++)
		{
			//printf("%d %d \n", pop[i].max_length, pop[i].total_cost);
			fla[i] = 0;
			if(pop[i].rank == 1)
			{
				fla[i] = 1;
				count++;
			}
			//printf("%d %d %d \n", pop[i].max_length, pop[i].total_cost,fla[i]);
				
		}

		int maxfm = 0, minfm=INF,maxft = 0, minft=INF,maxfv = 0, minfv=INF;
		int total_inx = 0;
		for (int i = 0; i < tmp_popsize;i++)
		{	
			if(fla[i])
			{
				if(pop[i].total_cost < minft)
				{
					minft = pop[i].total_cost;
					total_inx = i;
				}

				if(pop[i].total_cost > maxft)
					maxft = pop[i].total_cost;
				
				if(pop[i].max_length < minfm)
					minfm = pop[i].max_length;

				if(pop[i].max_length > maxfm)
					maxfm = pop[i].max_length;

				if(pop[i].total_vio_load < minfv)
					minfv = pop[i].total_vio_load;

				if(pop[i].total_vio_load > maxfv)
					maxfv = pop[i].total_vio_load;
			}
		}


		
		for (int i = 0; i < tmp_popsize; i++)
		{
			for (int j = 0; j < tmp_popsize; j++)
			{
				if(i == j)
					dis[i][j] = INF;
				else if(fla[i] && fla[j])
				{
					double to,ma;
					// a1 = (pop[i].max_length - pop[j].max_length) * (pop[i].max_length - pop[j].max_length) / ((maxfm - minfm) * (maxfm - minfm));
					// a2 = (pop[i].total_cost - pop[j].total_cost) * (pop[i].total_cost - pop[j].total_cost) / ((maxft - minft) * (maxft - minft));
					
					// a1 =  float(abs(pop[i].max_length - pop[j].max_length)) ;
					// a2 = float(abs(pop[i].total_cost - pop[j].total_cost));
					//dis[i][j] = sqrt(a1*a1 + a2*a2);

					// a1 =  float(abs(pop[i].max_length - pop[j].max_length)) / (maxfm - minfm + 1) ;
					// a2 = float(abs(pop[i].total_cost - pop[j].total_cost)) / (maxft - minft + 1);
					// dis[i][j] = a1 + a2;

					//printf("abs(%d - %d) / (%d) = %f \n", pop[i].max_length, pop[j].max_length,maxfm-minfm,a1 );
					//printf("abs(%d - %d) / (%d - %d) = %f \n", pop[i].total_cost, pop[j].total_cost,maxft,minft,a2 );
					
					// char avalue[100] = "avalue.dat";
					// FILE *fp;
					// fp = fopen(avalue, "w");
					// for (int i = 0; i < popsize;i++)
					// {	
					// 	fprintf(fp, "%f %f \n", a1, a2);
					// }
					// fprintf(fp, "\n");
					// fclose(fp);
					to = std::max(pop[i].total_cost, pop[j].total_cost);
					ma = std::max(pop[i].max_length, pop[j].max_length);
					to = float(abs(pop[i].total_cost - to)/(maxft-minft));
					ma = float(abs(pop[i].max_length - ma)/(maxfm-minfm));
					dis[i][j] = float(sqrt(to*to+ma*ma));
				}
			}
			//printf(" \n");
		}
		count = count - tmp_popsize/2;
		//printf("times of Truncation: %d \n", count);
		// for (int it = 0; it < count; it++)
		// {
		// 	int min_indx, comp = INF;
		// 	int min_indx2, comp2 = INF; // 第二小的值
		// 	for (int i = 0; i < tmp_popsize; i++)
		// 	{   // 计算 fitness<1 中的矩阵的最小值 储存在d[i]
		// 		d[i] = INF;
		// 		if(fla[i] > 0)
		// 		{
		// 			for (int j = 0; j < tmp_popsize; j++)
		// 			{
		// 				if(fla[i] == 1 && dis[i][j] < d[i])
		// 					d[i] = dis[i][j];
		// 			}

		// 			if( d[i] < comp )
		// 			{
		// 				// 最小变为第二小，更新最小
		// 				comp = d[i];
		// 				comp2 = comp;
		// 				min_indx2 = min_indx;
		// 				min_indx = i;

		// 			}
		// 			else if(d[i] >= comp || d[i] > comp2)
		// 			{
		// 				// 只更新第二个
		// 				min_indx2 = i;
		// 				comp2 = d[i];
		// 			}
		// 		}

		// 	}
		// 	if(min_indx != total_inx)
		// 		fla[min_indx] = 0;
		// 	else
		// 		fla[min_indx2] = 0;
		// }

		// int su;
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	su = su + fla[i];
		// }
		int iind[300],co[300];
		for (int i = 0; i < 300;i++)
			iind[i] = 0;
		int num = 1,po = 0;
		
		for (int it = 0; it < count+1; it++)
		{	
			co[0] = 0;
			po = 0;
			num = 1;
			int ti;
			for (int i = 0;i < tmp_popsize;i++)
			{
				if(fla[i] == 1)
				{
					for (int j = i + 1;j < tmp_popsize;j++)
					{
						if(fla[j] == 1)
						{
							if (pop[i].total_cost == pop[j].total_cost && pop[i].max_length == pop[j].max_length)
							{
								num++;
							}
							else
							{
								co[po] = i; // 标记第i个簇的在pop中每一个簇的第一个下标
								iind[po] = num; //标记第i个簇的个数
								po++;
								num = 1;
							}
							break;
						}
					}
					ti = i;
				}
				
			}
			
				co[po] = ti;
				iind[po] = num;
				po++;
				//printf("**************\n");
				int mma = 0,mma_i;
				for (int iii = 0; iii < po;iii++)
				{
					if(mma < iind[iii])
					{
						mma_i = co[iii];
						mma = iind[iii];
						*MMM = mma;
					}
				}
			if(it < count)
			{
				fla[mma_i] = 0;
			}
		}
		
		*F = po;
		//printf("po: %d\n", po);

		// printf("\n");
		// for (int i = 0; i < tmp_popsize;i++)
		// {
		// 	//printf("%d %d \n", pop[i].max_length, pop[i].total_cost);
		// 	printf("%d : %d %d %d \n",i, pop[i].max_length, pop[i].total_cost,fla[i]);
				
		// }
		// for (int i = 0; i < po ;i++)
		// {
		// 	//printf("%d %d \n", pop[i].max_length, pop[i].total_cost);
		// 	printf("%d %d \n", co[i], iind[i]);
				
		// }
		// printf("**************\n");

		int jj1 = 0,jj2 = tmp_popsize/2;
		for (int i = 0; i < tmp_popsize; i++)
		{ 
			if(fla[i] == 1 )
			{
				indi_copy(&tmp_indi[jj1], &pop[i]);
				jj1 = jj1 + 1;
			}
			else 
			{
				indi_copy(&tmp_indi[jj2], &pop[i]);
				jj2 = jj2 + 1;
			}
		}
		// printf("jj1: %d ,jj2: %d \n", jj1, jj2);
		// int su2 = 0;
		// for (int i = 0; i < tmp_popsize;i++)
		// 	su2 += fla[i];
		// printf("su2: %d \n", su2);
		for (int i = 0; i < tmp_popsize;i++)
		{
			indi_copy( &pop[i],&tmp_indi[i]);	
		}


	}

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void SPEA3(individual *pop, int tmp_popsize,int F)
{	// F:1 : 直接 Truncation
	if(F == 0)
	{
		for (int i = 0; i < tmp_popsize; i++)
		{
			for (int j = 0; j < tmp_popsize; j++)
			{
				dom[i][j] = 0;
			}
		}


		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	printf("%d %d %d \n",pop[i].max_length,pop[i].total_cost,pop[i].total_vio_load);
			
		// }
		// printf("\n");


		for (int i = 0; i < tmp_popsize-1; i++)
		{
			for (int j = i+1; j < tmp_popsize; j++)
			{
				if(dominate2(pop, i, j) == 1)
					dom[i][j] = 1;
				else if(dominate2(pop, i, j) == -1)
					dom[j][i] = 1;
			}
		}
		// printf("dom: \n");
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	for (int j = 0; j < tmp_popsize; j++)
		// 	{
		// 		printf("%d ",dom[i][j]);
		// 	}
		// 	printf("\n");
		// }


		for (int i = 0; i < tmp_popsize;i++)
		{
			S[i] = 0;
			for (int j = 0; j < tmp_popsize;j++)
			{
				S[i] = S[i] + dom[i][j];
			}
		}

		// printf("S: \n");
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	printf("%d ",S[i]);
			
		// }
		// printf("\n");

		for (int i = 0; i < tmp_popsize;i++)
		{
			R[i] = 0;
			for (int j = 0; j < tmp_popsize;j++)
			{
				if(dom[j][i] == 1)
					R[i] = R[i] + S[j];
			}
		}

		// printf("R: \n");
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	printf("%d ",R[i]);
			
		// }
		// printf("\n");
		

		int maxfm = 0, minfm=INF,maxft = 0, minft=INF,maxfv = 0, minfv=INF;
		//int total_fir = 0, total_secd = 0, maxft_secd = 0;
		for (int i = 0; i < tmp_popsize;i++)
		{	
			if(pop[i].total_cost < minft)
				minft = pop[i].total_cost;

			if(pop[i].total_cost > maxft)
				maxft = pop[i].total_cost;
			
			if(pop[i].max_length < minfm)
				minfm = pop[i].max_length;

			if(pop[i].max_length > maxfm)
				maxfm = pop[i].max_length;

			if(pop[i].total_vio_load < minfv)
				minfv = pop[i].total_vio_load;

			if(pop[i].total_vio_load > maxfv)
				maxfv = pop[i].total_vio_load;
		}
		// printf("scale: \n");
		// printf("total: %d, %d \n",minft,maxft);
		// printf("max: %d, %d \n",minfm,maxfm);
		// printf("f1 f2: \n");
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	printf("%d %d \n",pop[i].max_length, pop[i].total_cost);
			
		// }
		// printf("\n");


		for (int i = 0; i < tmp_popsize; i++)
		{
			for (int j = 0; j < tmp_popsize; j++)
			{
				if(i == j)
					dis[i][j] = INF;
				else
				{
					double vi,to ,ma;
					// a1 = (pop[i].max_length - pop[j].max_length) * (pop[i].max_length - pop[j].max_length) / ((maxfm - minfm) * (maxfm - minfm));
					// a2 = (pop[i].total_cost - pop[j].total_cost) * (pop[i].total_cost - pop[j].total_cost) / ((maxft - minft) * (maxft - minft));
					
					// a1 =  float(abs(pop[i].max_length - pop[j].max_length)) ;
					// a2 = float(abs(pop[i].total_cost - pop[j].total_cost));
					//dis[i][j] = sqrt(a1*a1 + a2*a2);
					
					// a1 =  float(abs(pop[i].max_length - pop[j].max_length)) / (maxfm - minfm + 1) ;
					// a2 = float(abs(pop[i].total_cost - pop[j].total_cost)) / (maxft - minft + 1);
					// dis[i][j] = a1 + a2;

					// printf("abs(%d - %d) / (%d) = %f \n", pop[i].max_length, pop[j].max_length,maxfm-minfm,a1 );
					// printf("abs(%d - %d) / (%d - %d) = %f \n", pop[i].total_cost, pop[j].total_cost,maxft,minft,a2 );
					

					to = std::max(pop[i].total_cost, pop[j].total_cost);
					ma = std::max(pop[i].max_length, pop[j].max_length);
					vi = std::max(pop[i].total_vio_load, pop[j].total_vio_load);
					to = float(abs(pop[i].total_cost - to)/(maxft-minft+1));
					ma = float(abs(pop[i].max_length - ma)/(maxfm-minfm+1));
					vi = float(abs(pop[i].total_vio_load- vi)/(maxfv-minfv+1));
					dis[i][j] = float( sqrt(to*to + ma*ma + vi*vi) );



					// char avalue[100] = "avalue.dat";
					// FILE *fp;
					// fp = fopen(avalue, "w");
					// for (int i = 0; i < popsize;i++)
					// {	
					// 	fprintf(fp, "%f %f \n", a1, a2);
					// }
					// fprintf(fp, "\n");
					// fclose(fp);

				}
			}
			//printf(" \n");
		}

		// printf("dist: \n");
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	for (int j = 0; j < tmp_popsize; j++)
		// 	{
		// 		printf("%3f ",dis[i][j]);
		// 	}
		// 	printf("\n");
		// }

		//int k = int(floor(sqrt(tmp_popsize)));
		//int k = 1;
		for (int i = 0; i < tmp_popsize;i++)
		{
			//D[i] = float(1/(2+quick_select_min(dis[i], tmp_popsize, k)));
			D[i] = float(1/(2+quick_select_min(dis[i], tmp_popsize, 1)));
		}

		for (int i = 0; i < tmp_popsize;i++)
		{
			fitness[i] = R[i] + D[i];
			fla[i] = 0;
		}

		// printf("fitness: \n");
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	printf("%d + %f = %f ",R[i], D[i],fitness[i]);
			
		// }
		// printf("\n");

		int count = 0;
		for (int i = 0; i < tmp_popsize; i++)
		{
			if(fitness[i] < 1)
			{
				count ++;
				fla[i] = 1;
			}
		}
		// printf("count: %d\n", count);
		// printf("count: %d\n", tmp_popsize/2);
		if(count < tmp_popsize/2 + 1 )
		//if(count )
		{ // 直接截断

			for (int i = 0; i < tmp_popsize;i++)
			{
				idx[i] = i;
			}

			// for (int i = 0; i < tmp_popsize;i++)
			// {
			// 	printf("%f\n", fitness[i]);
			// }

			for (int i = 0; i < tmp_popsize-1; i++)
			{
				for (int j = i+1; j < tmp_popsize; j++)
				{
					if(fitness[i] > fitness[j])
					{
						int chang;
						chang = idx[i];
						idx[i] = idx[j];
						idx[j] = chang;

						double chang2;
						chang2 = fitness[i];
						fitness[i] = fitness[j];
						fitness[j] = chang2;
					}
				}
			}

			// for (int i = 0; i < tmp_popsize;i++)
			// {
			// 	printf("%f  %d\n", fitness[i], idx[i]);
			// }

			for (int i = 0; i < tmp_popsize; i++)
			{
				indi_copy(&tmp_indi[i], &pop[idx[i]]);
			}

			for (int i = 0; i < tmp_popsize;i++)
			{
				indi_copy( &pop[i],&tmp_indi[i]);	
			}
		}
		else 
		{ // Truncation
			// for (int i = 0; i < tmp_popsize;i++)
			// {
			// 	printf("%d  %d\n", pop[i].max_length, pop[i].total_cost);
			// }
			count = count - tmp_popsize/2;
			//printf("number of Truncation: %d \n", count);
			for (int it = 0; it < count; it++)
			{
				int min_indx, comp = INF;
				for (int i = 0; i < tmp_popsize; i++)
				{   // 计算 fitness<1 中的矩阵的最小值 储存在d[i]
					d[i] = INF;
					for (int j = 0; j < tmp_popsize; j++)
					{
						if(fla[i] == 1 && dis[i][j] < d[i])
							d[i] = dis[i][j];
					}

					if(fla[i] == 1 && d[i] < comp )
					{
						comp = d[i];
						min_indx = i;
					}

				}
				fla[min_indx] = 0;
			}

			int su;
			for (int i = 0; i < tmp_popsize; i++)
			{
				su = su + fla[i];
			}

			int jj1 = 0,jj2 = tmp_popsize/2;
			for (int i = 0; i < tmp_popsize; i++)
			{ 
				if(fla[i] == 1 )
				{
					indi_copy(&tmp_indi[jj1], &pop[i]);
					jj1 = jj1 + 1;
				}
				else 
				{
					indi_copy(&tmp_indi[jj2], &pop[i]);
					jj2 = jj2 + 1;
				}
			}
			// printf("jj1: %d ,jj2: %d \n", jj1, jj2);
			// int su2 = 0;
			// for (int i = 0; i < tmp_popsize;i++)
			// 	su2 += fla[i];
			// printf("su2: %d \n", su2);
			for (int i = 0; i < tmp_popsize;i++)
			{
				indi_copy( &pop[i],&tmp_indi[i]);	
			}
		}
	}
	else
	{

		int count = 0;
		for (int i = 0; i < tmp_popsize;i++)
		{
			fla[i] = 0;
			if(pop[i].rank == 1)
			{
				fla[i] = 1;
				count++;
			}
				
		}

		int maxfm = 0, minfm=INF,maxft = 0, minft=INF,maxfv = 0, minfv=INF;
		int total_inx = 0;
		for (int i = 0; i < tmp_popsize;i++)
		{	
			if(fla[i])
			{
				if(pop[i].total_cost < minft)
				{
					minft = pop[i].total_cost;
					total_inx = i;
				}

				if(pop[i].total_cost > maxft)
					maxft = pop[i].total_cost;
				
				if(pop[i].max_length < minfm)
					minfm = pop[i].max_length;

				if(pop[i].max_length > maxfm)
					maxfm = pop[i].max_length;

				if(pop[i].total_vio_load < minfv)
					minfv = pop[i].total_vio_load;

				if(pop[i].total_vio_load > maxfv)
					maxfv = pop[i].total_vio_load;
			}
		}

		
		for (int i = 0; i < tmp_popsize; i++)
		{
			for (int j = 0; j < tmp_popsize; j++)
			{
				if(i == j)
					dis[i][j] = INF;
				else if(fla[i] && fla[j])
				{
					double to,ma,vi;
					// a1 = (pop[i].max_length - pop[j].max_length) * (pop[i].max_length - pop[j].max_length) / ((maxfm - minfm) * (maxfm - minfm));
					// a2 = (pop[i].total_cost - pop[j].total_cost) * (pop[i].total_cost - pop[j].total_cost) / ((maxft - minft) * (maxft - minft));
					
					// a1 =  float(abs(pop[i].max_length - pop[j].max_length)) ;
					// a2 = float(abs(pop[i].total_cost - pop[j].total_cost));
					//dis[i][j] = sqrt(a1*a1 + a2*a2);

					// a1 =  float(abs(pop[i].max_length - pop[j].max_length)) / (maxfm - minfm + 1) ;
					// a2 = float(abs(pop[i].total_cost - pop[j].total_cost)) / (maxft - minft + 1);
					// dis[i][j] = a1 + a2;

					//printf("abs(%d - %d) / (%d) = %f \n", pop[i].max_length, pop[j].max_length,maxfm-minfm,a1 );
					//printf("abs(%d - %d) / (%d - %d) = %f \n", pop[i].total_cost, pop[j].total_cost,maxft,minft,a2 );
					
					// char avalue[100] = "avalue.dat";
					// FILE *fp;
					// fp = fopen(avalue, "w");
					// for (int i = 0; i < popsize;i++)
					// {	
					// 	fprintf(fp, "%f %f \n", a1, a2);
					// }
					// fprintf(fp, "\n");
					// fclose(fp);

					to = std::max(pop[i].total_cost, pop[j].total_cost);
					ma = std::max(pop[i].max_length, pop[j].max_length);
					vi = std::max(pop[i].total_vio_load, pop[j].total_vio_load);
					to = float(abs(pop[i].total_cost - to)/(maxft-minft));
					ma = float(abs(pop[i].max_length - ma)/(maxfm-minfm));
					vi = float(abs(pop[i].total_vio_load- vi)/(maxfv-minfv));
					dis[i][j] = float( sqrt(to*to + ma*ma + vi*vi) );
				}
			}
			//printf(" \n");
		}
		count = count - tmp_popsize/2;
		//printf("times of Truncation: %d \n", count);
		for (int it = 0; it < count; it++)
		{
			int min_indx, comp = INF;   // 最小的值
			int min_indx2, comp2 = INF; // 第二小的值
			for (int i = 0; i < tmp_popsize; i++)
			{   // 计算 fitness<1 中的矩阵的最小值 储存在d[i]
				d[i] = INF;
				if(fla[i] == 1)
				{
					for (int j = 0; j < tmp_popsize; j++)
					{
						if(fla[i] == 1 && dis[i][j] < d[i])
							d[i] = dis[i][j];
					}

					if( d[i] < comp )
					{
						// 最小变为第二小，更新最小
						comp = d[i];
						comp2 = comp;
						min_indx2 = min_indx;
						min_indx = i;

					}
					else if(d[i] >= comp || d[i] > comp2)
					{
						// 只更新第二个
						min_indx2 = i;
						comp2 = d[i];
					}
				}

			}
			if(min_indx != total_inx)
				fla[min_indx] = 0;
			else
				fla[min_indx2] = 0;
		}

		// int su;
		// for (int i = 0; i < tmp_popsize; i++)
		// {
		// 	su = su + fla[i];
		// }

		int jj1 = 0,jj2 = tmp_popsize/2;
		for (int i = 0; i < tmp_popsize; i++)
		{ 
			if(fla[i] == 1 )
			{
				indi_copy(&tmp_indi[jj1], &pop[i]);
				jj1 = jj1 + 1;
			}
			else 
			{
				indi_copy(&tmp_indi[jj2], &pop[i]);
				jj2 = jj2 + 1;
			}
		}
		// printf("jj1: %d ,jj2: %d \n", jj1, jj2);
		// int su2 = 0;
		// for (int i = 0; i < tmp_popsize;i++)
		// 	su2 += fla[i];
		// printf("su2: %d \n", su2);
		for (int i = 0; i < tmp_popsize;i++)
		{
			indi_copy( &pop[i],&tmp_indi[i]);	
		}


	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void front_crowdness_sort(individual *pop, int rank, int tmp_popsize)
{
	int rank_indis[MAX_TOTALSIZE];
	rank_indis[0] = 0;
	int start, finish;
	for (start = 0; start < tmp_popsize; start++)
	{
		if (pop[start].rank == rank)
			break;
	}
	for (finish = start; finish < tmp_popsize; finish++)
	{
		if (finish == tmp_popsize-1 || pop[finish+1].rank > rank)
			break;
	}
	for (int i = start; i <= finish; i++)
	{
		rank_indis[0] ++;
		rank_indis[rank_indis[0]] = i;
	}

	for (int i = 1; i <= rank_indis[0]; i++)
	{
		pop[rank_indis[i]].crowdness = 0;
	}

	if (rank_indis[0] == 1)
	{
		pop[rank_indis[1]].crowdness = INF;
		return;
	}

	int maxf, minf;

	/* sort by total_cost */

	for (int i = 1; i < rank_indis[0]; i++)
	{
		for (int j = i+1; j <= rank_indis[0]; j++)
		{
			if (pop[rank_indis[j]].total_cost < pop[rank_indis[i]].total_cost)
			{
				int tmp = rank_indis[i];
				rank_indis[i] = rank_indis[j];
				rank_indis[j] = tmp;
			}
		}
	}

	minf = pop[rank_indis[1]].total_cost;
	maxf = pop[rank_indis[rank_indis[0]]].total_cost;
	pop[rank_indis[1]].crowdness = INF;
	pop[rank_indis[rank_indis[0]]].crowdness = INF;

	for (int i = 2; i < rank_indis[0]; i++)
	{
		pop[rank_indis[i]].crowdness += 1.0*(pop[rank_indis[i+1]].total_cost-pop[rank_indis[i-1]].total_cost)/(maxf-minf);
	}

	/* sort by makespan */

	for (int i = 1; i < rank_indis[0]; i++)
	{
		for (int j = i+1; j <= rank_indis[0]; j++)
		{
			if (pop[rank_indis[j]].max_length < pop[rank_indis[i]].max_length)
			{
				int tmp = rank_indis[i];
				rank_indis[i] = rank_indis[j];
				rank_indis[j] = tmp;
			}
		}
	}

	minf = pop[rank_indis[1]].max_length;
	maxf = pop[rank_indis[rank_indis[0]]].max_length;
	pop[rank_indis[1]].crowdness = INF;
	pop[rank_indis[rank_indis[0]]].crowdness = INF;

	for (int i = 2; i < rank_indis[0]; i++)
	{
		if (pop[rank_indis[i]].crowdness < INF)
			pop[rank_indis[i]].crowdness += 1.0*(pop[rank_indis[i+1]].max_length-pop[rank_indis[i-1]].max_length)/(maxf-minf);
	}

	individual tmp_indi;

	for (int i = start; i < finish; i++)
	{
		for (int j = i+1; j <= finish; j++)
		{
			if (pop[j].crowdness > pop[i].crowdness)
			{
				indi_copy(&tmp_indi, &pop[i]);
				indi_copy(&pop[i], &pop[j]);
				indi_copy(&pop[j], &tmp_indi);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int move_fit_cmp(move *mv1, move *mv2)
// return -1 if mv1 is better, and 1 if mv2 is better
{
	if (mv1->add_vio_load < mv2->add_vio_load)
		return -1;

	if (mv1->add_vio_load > mv2->add_vio_load)
		return 1;

	if ((mv1->total_cost == mv2->total_cost && mv1->max_length < mv2->max_length) ||
		(mv1->max_length == mv2->max_length && mv1->total_cost < mv2->total_cost))
		return -1;

	if ((mv1->total_cost == mv2->total_cost && mv1->max_length > mv2->max_length) ||
		(mv1->max_length == mv2->max_length && mv1->total_cost > mv2->total_cost))
		return 1;

	if (mv1->fitness < mv2->fitness)
		return -1;

	if (mv1->fitness > mv2->fitness)
		return 1;

	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int indi_fit_cmp(individual *indi1, individual *indi2)
{
	if (indi1->total_vio_load < indi2->total_vio_load)
		return -1;

	if (indi1->total_vio_load > indi2->total_vio_load)
		return 1;

	if ((indi1->total_cost == indi2->total_cost && indi1->max_length < indi2->max_length) ||
		(indi1->max_length == indi2->max_length && indi1->total_cost < indi2->total_cost))
		return -1;

	if ((indi1->total_cost == indi2->total_cost && indi1->max_length > indi2->max_length) ||
		(indi1->max_length == indi2->max_length && indi1->total_cost > indi2->total_cost))
		return 1;

	//printf("fit1 = %.20lf, fit2 = %.20lf\n", indi1->fitness, indi2->fitness);
	if (indi1->fitness < indi2->fitness)
		return -1;

	if (indi1->fitness > indi2->fitness)
		return 1;

	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int mv_indi_fit_cmp(move *mv, individual *indi)
{
	if (mv->add_vio_load < 0)
		return -1;

	if (mv->add_vio_load > 0)
		return 1;

	if ((mv->total_cost == indi->total_cost && mv->max_length < indi->max_length) ||
		(mv->max_length == indi->max_length && mv->total_cost < indi->total_cost))
		return -1;

	if ((mv->total_cost == indi->total_cost && mv->max_length > indi->max_length) ||
		(mv->max_length == indi->max_length && mv->total_cost > indi->total_cost))
		return 1;

	if (mv->fitness < indi->fitness)
		return -1;

	if (mv->fitness > indi->fitness)
		return 1;

	return 0;
}