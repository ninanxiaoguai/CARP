
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mod_dijkstra()
{
	int i, j, k, m, minimum;
	
	for (i = 1; i <= vertex_num; i++)
	{
		for (j = 1; j <= vertex_num; j++)
		{
			if (j == i)
				continue;
			
			shortest_path[i][j][0] = 1;
			shortest_path[i][j][1] = i;
			min_cost[i][j] = INF;
		}
	}

	int mark[MAX_NODE_TAG_LENGTH], dist[MAX_NODE_TAG_LENGTH], dist1[MAX_NODE_TAG_LENGTH], nearest_neighbor[MAX_NODE_TAG_LENGTH];

	for (i = 1; i <= vertex_num; i++)
	{
		mark[i] = 1;
		
		for (j = 1; j <= vertex_num; j++)
		{
			if (j == i)
				continue;
			
			mark[j] = 0;
			dist[j] = trav_cost[i][j];
			dist1[j] = dist[j];
		}
		
		for (k = 1; k < vertex_num; k++)
		{
			minimum = INF;
			nearest_neighbor[0] = 0;
			
			for (j = 1; j <= vertex_num; j++)
			{
				if (mark[j])
					continue;
					
				if (dist1[j] == INF)
					continue;
					
				if (dist1[j] < minimum)
					minimum = dist1[j];    // 找到已i为起点的边的最小值
			}
			
			if (minimum == INF)
				continue;
			
			for (j = 1; j <= vertex_num; j++)
			{
				if (mark[j])
					continue;
					
				if (dist1[j] == minimum)
				{
					nearest_neighbor[0] ++;
					nearest_neighbor[nearest_neighbor[0]] = j;
				}
			}
			
			int v = nearest_neighbor[1];
			dist1[v] = INF;
			mark[v] = 1;
			
			if (shortest_path[i][v][0] == 0 || (shortest_path[i][v][0] > 0 && shortest_path[i][v][shortest_path[i][v][0]] != v))
			{
				shortest_path[i][v][0] ++;
				shortest_path[i][v][shortest_path[i][v][0]] = v;
			}
				
			for (j = 1; j <= vertex_num; j++)
			{
				if (mark[j])
					continue;
					
				if (minimum+trav_cost[v][j] < dist[j])
				{
					dist[j] = minimum+trav_cost[v][j];
					dist1[j] = minimum+trav_cost[v][j];
					for (m = 0; m <= shortest_path[i][v][0]; m++)
					{
						shortest_path[i][j][m] = shortest_path[i][v][m];
					}
				}
			}
			
			for (j = 1; j <= vertex_num; j++)
			{
				if (j == i)
					continue;
				
				min_cost[i][j] = dist[j];
			}
		}		
	}
	
	for (i = 1; i <= vertex_num; i++)
	{
		for (j = 1; j <= vertex_num; j++)
		{
			if (shortest_path[i][j][0] == 1)
				shortest_path[i][j][0] = 0;
		}
	}
	
	for (i = 1; i <= vertex_num; i++)
	{
		shortest_path[i][i][0] = 1;
		shortest_path[i][i][1] = i;
		min_cost[i][i] = 0;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void path_scanning(individual *ps_indi, const task *inst_tasks, int *serve_mark)
{
	int i, k;
	int serve_task_num = 0;
	for (i = req_edge_num+1; i <= task_num; i++)
	{
		if (serve_mark[i])
			serve_task_num ++;
	}
	
	int load, trial, mindist;

	int unserved_task[MAX_TASK_TAG_LENGTH], candi_task[MAX_TASK_TAG_LENGTH], nearest_task[MAX_TASK_TAG_LENGTH];
	int nearest_isol_task[MAX_TASK_TAG_LENGTH], nearest_inci_task[MAX_TASK_TAG_LENGTH], sel_task[MAX_TASK_TAG_LENGTH];
	int current_task, next_task;

	int positions[MAX_TASK_SEQ_LENGTH];
	
	individual tmp_indi1, tmp_indi2, tmp_indi3, tmp_indi4, tmp_indi5;

	ps_indi->total_cost = INF;

	int dep_dist[MAX_TASK_TAG_LENGTH];
	double yield[MAX_TASK_TAG_LENGTH];
	int max_dep_dist, min_dep_dist;
	double max_yield, min_yield;
	
	for (i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		dep_dist[i] = min_cost[inst_tasks[i].tail_node][DEPOT]; // depot 到每个任务的的trav_cost
		yield[i] = 1.0*inst_tasks[i].demand/inst_tasks[i].serv_cost;
	}
	
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    tmp_indi1.sequence[0] = 1;
	tmp_indi1.sequence[1] = 0;
	tmp_indi1.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi1.sequence[tmp_indi1.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= capacity-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			tmp_indi1.sequence[0] ++;
	  		tmp_indi1.sequence[tmp_indi1.sequence[0]] = 0;
		  	tmp_indi1.route_seg_load[0] ++;
			tmp_indi1.route_seg_load[tmp_indi1.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		max_dep_dist = -1;
		sel_task[0] = 0;
		for (i = 1; i <= nearest_isol_task[0]; i++)
		{
			if (dep_dist[nearest_isol_task[i]] > max_dep_dist)
			{
				max_dep_dist = dep_dist[nearest_isol_task[i]];
				sel_task[0] = 1;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
			else if (dep_dist[nearest_isol_task[i]] == max_dep_dist)
			{
				sel_task[0] ++;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi1.sequence[0] ++;
		tmp_indi1.sequence[tmp_indi1.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi1.sequence[0] ++;
	tmp_indi1.sequence[tmp_indi1.sequence[0]] = 0;
	tmp_indi1.route_seg_load[0] ++;
	tmp_indi1.route_seg_load[tmp_indi1.route_seg_load[0]] = load;
		
	tmp_indi1.total_cost = get_task_seq_total_cost(tmp_indi1.sequence, inst_tasks);
	tmp_indi1.total_vio_load = get_total_vio_load(tmp_indi1.route_seg_load);
	
	if (tmp_indi1.total_cost < ps_indi->total_cost)
	{
		memcpy(ps_indi->sequence, tmp_indi1.sequence, (tmp_indi1.sequence[0]+1)*sizeof(int));
		memcpy(ps_indi->route_seg_load, tmp_indi1.route_seg_load, (tmp_indi1.route_seg_load[0]+1)*sizeof(int));
		ps_indi->total_cost = tmp_indi1.total_cost;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	tmp_indi2.sequence[0] = 1;
	tmp_indi2.sequence[1] = 0;
	tmp_indi2.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi2.sequence[tmp_indi2.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= capacity-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
			tmp_indi2.sequence[0] ++;
	  		tmp_indi2.sequence[tmp_indi2.sequence[0]] = 0;
		  	tmp_indi2.route_seg_load[0] ++;
			tmp_indi2.route_seg_load[tmp_indi2.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		min_dep_dist = INF;
		sel_task[0] = 0;
		for (i = 1; i <= nearest_isol_task[0]; i++)
		{
			if (dep_dist[nearest_isol_task[i]] < min_dep_dist)
			{
				min_dep_dist = dep_dist[nearest_isol_task[i]];
				sel_task[0] = 1;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
			else if (dep_dist[nearest_isol_task[i]] == min_dep_dist)
			{
				sel_task[0] ++;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi2.sequence[0] ++;
		tmp_indi2.sequence[tmp_indi2.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi2.sequence[0] ++;
	tmp_indi2.sequence[tmp_indi2.sequence[0]] = 0;
	tmp_indi2.route_seg_load[0] ++;
	tmp_indi2.route_seg_load[tmp_indi2.route_seg_load[0]] = load;
		
	tmp_indi2.total_cost = get_task_seq_total_cost(tmp_indi2.sequence, inst_tasks);
	tmp_indi2.total_vio_load = get_total_vio_load(tmp_indi2.route_seg_load);
	
	if (tmp_indi2.total_cost < ps_indi->total_cost)
	{
		memcpy(ps_indi->sequence, tmp_indi2.sequence, (tmp_indi2.sequence[0]+1)*sizeof(int));
		memcpy(ps_indi->route_seg_load, tmp_indi2.route_seg_load, (tmp_indi2.route_seg_load[0]+1)*sizeof(int));
		ps_indi->total_cost = tmp_indi2.total_cost;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	tmp_indi3.sequence[0] = 1;
	tmp_indi3.sequence[1] = 0;
	tmp_indi3.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi3.sequence[tmp_indi3.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= capacity-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			tmp_indi3.sequence[0] ++;
	  		tmp_indi3.sequence[tmp_indi3.sequence[0]] = 0;
		  	tmp_indi3.route_seg_load[0] ++;
			tmp_indi3.route_seg_load[tmp_indi3.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		max_yield = -1;
		sel_task[0] = 0;
		for (i = 1; i <= nearest_isol_task[0]; i++)
		{
			if (yield[nearest_isol_task[i]] > max_yield)
			{
				max_yield = yield[nearest_isol_task[i]];
				sel_task[0] = 1;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
			else if (yield[nearest_isol_task[i]] == max_yield)
			{
				sel_task[0] ++;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi3.sequence[0] ++;
		tmp_indi3.sequence[tmp_indi3.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi3.sequence[0] ++;
	tmp_indi3.sequence[tmp_indi3.sequence[0]] = 0;
	tmp_indi3.route_seg_load[0] ++;
	tmp_indi3.route_seg_load[tmp_indi3.route_seg_load[0]] = load;
		
	tmp_indi3.total_cost = get_task_seq_total_cost(tmp_indi3.sequence, inst_tasks);
	tmp_indi3.total_vio_load = get_total_vio_load(tmp_indi3.route_seg_load);
	
	if (tmp_indi3.total_cost < ps_indi->total_cost)
	{
		memcpy(ps_indi->sequence, tmp_indi3.sequence, (tmp_indi3.sequence[0]+1)*sizeof(int));
		memcpy(ps_indi->route_seg_load, tmp_indi3.route_seg_load, (tmp_indi3.route_seg_load[0]+1)*sizeof(int));
		ps_indi->total_cost = tmp_indi3.total_cost;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	tmp_indi4.sequence[0] = 1;
	tmp_indi4.sequence[1] = 0;
	tmp_indi4.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi4.sequence[tmp_indi4.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= capacity-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			tmp_indi4.sequence[0] ++;
	  		tmp_indi4.sequence[tmp_indi4.sequence[0]] = 0;
		  	tmp_indi4.route_seg_load[0] ++;
			tmp_indi4.route_seg_load[tmp_indi4.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		min_yield = INF;
		sel_task[0] = 0;
		for (i = 1; i <= nearest_isol_task[0]; i++)
		{
			if (yield[nearest_isol_task[i]] < min_yield)
			{
				min_yield = yield[nearest_isol_task[i]];
				sel_task[0] = 1;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
			else if (yield[nearest_isol_task[i]] == min_yield)
			{
				sel_task[0] ++;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi4.sequence[0] ++;
		tmp_indi4.sequence[tmp_indi4.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi4.sequence[0] ++;
	tmp_indi4.sequence[tmp_indi4.sequence[0]] = 0;
	tmp_indi4.route_seg_load[0] ++;
	tmp_indi4.route_seg_load[tmp_indi4.route_seg_load[0]] = load;
		
	tmp_indi4.total_cost = get_task_seq_total_cost(tmp_indi4.sequence, inst_tasks);
	tmp_indi4.total_vio_load = get_total_vio_load(tmp_indi4.route_seg_load);
	
	if (tmp_indi4.total_cost < ps_indi->total_cost)
	{
		memcpy(ps_indi->sequence, tmp_indi4.sequence, (tmp_indi4.sequence[0]+1)*sizeof(int));
		memcpy(ps_indi->route_seg_load, tmp_indi4.route_seg_load, (tmp_indi4.route_seg_load[0]+1)*sizeof(int));
		ps_indi->total_cost = tmp_indi4.total_cost;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    tmp_indi5.sequence[0] = 1;
	tmp_indi5.sequence[1] = 0;
	tmp_indi5.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi5.sequence[tmp_indi5.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= capacity-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			tmp_indi5.sequence[0] ++;
	  		tmp_indi5.sequence[tmp_indi5.sequence[0]] = 0;
		  	tmp_indi5.route_seg_load[0] ++;
			tmp_indi5.route_seg_load[tmp_indi5.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		if (load < capacity/2)
		{
			max_dep_dist = -1;
			sel_task[0] = 0;
			for (i = 1; i <= nearest_isol_task[0]; i++)
			{
				if (dep_dist[nearest_isol_task[i]] > max_dep_dist)
				{
					max_dep_dist = dep_dist[nearest_isol_task[i]];
					sel_task[0] = 1;
					sel_task[sel_task[0]] = nearest_isol_task[i];
				}
				else if (dep_dist[nearest_isol_task[i]] == max_dep_dist)
				{
					sel_task[0] ++;
					sel_task[sel_task[0]] = nearest_isol_task[i];
				}
			}
		}
		else
		{
			min_dep_dist = INF;
			sel_task[0] = 0;
			for (i = 1; i <= nearest_isol_task[0]; i++)
			{
				if (dep_dist[nearest_isol_task[i]] < min_dep_dist)
				{
					min_dep_dist = dep_dist[nearest_isol_task[i]];
					sel_task[0] = 1;
					sel_task[sel_task[0]] = nearest_isol_task[i];
				}
				else if (dep_dist[nearest_isol_task[i]] == min_dep_dist)
				{
					sel_task[0] ++;
					sel_task[sel_task[0]] = nearest_isol_task[i];
				}
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi5.sequence[0] ++;
		tmp_indi5.sequence[tmp_indi5.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi5.sequence[0] ++;
	tmp_indi5.sequence[tmp_indi5.sequence[0]] = 0;
	tmp_indi5.route_seg_load[0] ++;
	tmp_indi5.route_seg_load[tmp_indi5.route_seg_load[0]] = load;
	
	tmp_indi5.total_cost = get_task_seq_total_cost(tmp_indi5.sequence, inst_tasks);
	tmp_indi5.total_vio_load = get_total_vio_load(tmp_indi5.route_seg_load);
	
	if (tmp_indi5.total_cost < ps_indi->total_cost)
	{
		memcpy(ps_indi->sequence, tmp_indi5.sequence, (tmp_indi5.sequence[0]+1)*sizeof(int));
		memcpy(ps_indi->route_seg_load, tmp_indi5.route_seg_load, (tmp_indi5.route_seg_load[0]+1)*sizeof(int));
		ps_indi->total_cost = tmp_indi5.total_cost;
	}

	ps_indi->total_vio_load = 0;
	get_route_seg_length(ps_indi->route_seg_length, ps_indi->sequence, inst_tasks);
	ps_indi->max_length = max(ps_indi->route_seg_length);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void rand_scanning(individual *rs_indi, const task *inst_tasks, int *serve_mark)
{
	int i, k;
	int serve_task_num = 0;
	for (i = req_edge_num+1; i <= task_num; i++)
	{   // 统计 serve_mark 中 1 的个数
		// 即为 edge 的个数(51)
		if (serve_mark[i])
			serve_task_num ++;
	}

	int load, trial, mindist;

	int unserved_task[MAX_TASK_TAG_LENGTH], candi_task[MAX_TASK_TAG_LENGTH], nearest_task[MAX_TASK_TAG_LENGTH];
	int current_task, next_task;

	int positions[MAX_TASK_SEQ_LENGTH];

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	rs_indi->sequence[0] = 1;
	rs_indi->sequence[1] = 0;
	rs_indi->route_seg_load[0] = 0;
	// task_num = 2*serve_task_num = 
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{ // 记录 serve_mark 中为 1 的的个数,储存在[0]中
	  // 并给 serve_mark 中[1~个数]分别给一个序号，与下标相同
		if (!serve_mark[i])
			continue;

		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}

	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = rs_indi->sequence[rs_indi->sequence[0]];

		candi_task[0] = 0;

		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= capacity-load)
			{// 如果有边的cost超过了capacity，那么这个边就不能考虑了
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
			//candi_task[0]: 目前所考虑的task的个数
			// candi_task[1~#]：所考虑的边的编号
		}

		if (candi_task[0] == 0)
		{  // 重新开启一个新路线
			rs_indi->sequence[0] ++;
			rs_indi->sequence[rs_indi->sequence[0]] = 0;
			//储存当前这一次的路径满足多少的demand
			rs_indi->route_seg_load[0] ++;
			rs_indi->route_seg_load[rs_indi->route_seg_load[0]] = load;
			load = 0;
			// 注意trial并没有增加
			continue;
		}

		mindist = INF;
		nearest_task[0] = 0;

		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{ // nearest_task[0] 中储存以 current_task的tail_node 为起点，距哪个edges的端点所消耗的 trav最小的有几个
			  // nearest_task[1~nearest_task[0]]，这几个的edges的标号(edges 因为是双向，所以为 2*#edges)
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}

		//k = rand_choose(nearest_task[0]);
		//next_task = nearest_task[k];

		k = rand_choose(candi_task[0]); // 在最近的这几个边里随机选一个， k为 记录边ID的数组candi_task的下标
		next_task = candi_task[k];

		trial ++;
		rs_indi->sequence[0] ++;
		rs_indi->sequence[rs_indi->sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;   // 加上next_task 后 load的值

		//查询在 unserved_task 中有 next_task 这个标号的个数(记录在[0])
		//并记录在[1~个数]里，记录unserved_task数组中与next_task相同的下标
		find_ele_positions(positions, unserved_task, next_task); 
		delete_element(unserved_task, positions[1]);

		if (inst_tasks[next_task].inverse > 0)
		{ // 如果 next_task 改边为edges即为双向边，那么就要再删一次
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
		//printf("No. leaving is %d\n",unserved_task[0]);
	}
	
	rs_indi->sequence[0] ++;
	rs_indi->sequence[rs_indi->sequence[0]] = 0;
	rs_indi->route_seg_load[0] ++;
	rs_indi->route_seg_load[rs_indi->route_seg_load[0]] = load;

	// printf("No. leaving is %d\n",rs_indi->route_seg_load[0]);
	// printf("\n");
	// 至此该个体的路线储存在sequence，用 0 分开，
	rs_indi->total_cost = get_task_seq_total_cost(rs_indi->sequence, inst_tasks);
	rs_indi->total_vio_load = get_total_vio_load(rs_indi->route_seg_load);
	// 储存每一个路线的cost 总和即为 rs_indi->total_cost
	get_route_seg_length(rs_indi->route_seg_length, rs_indi->sequence, inst_tasks);
	rs_indi->max_length = max(rs_indi->route_seg_length);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct route
{
	int perm[MAX_TASK_ROUTE_LENGTH];
	int cost;
	int load;
} route;

route task_routes[MAX_TASK_SEQ_LENGTH], tmp_route;

void augment_merge(individual *am_indi, const task *inst_tasks, int *serve_mark)
{
	int serve_task_num = 0;
	for (int i = req_edge_num+1; i <= task_num; i++)
	{
		if (serve_mark[i])
			serve_task_num ++;
	}

	/* generating phase */

	for (int i = 1; i <= req_edge_num+req_arc_num; i++)
	{
		task_routes[i].perm[0] = 3;
		task_routes[i].perm[1] = 0;
		task_routes[i].perm[2] = i+req_edge_num;
		task_routes[i].perm[3] = 0;
		task_routes[i].cost = get_task_seq_total_cost(task_routes[i].perm, inst_tasks);
		task_routes[i].load = inst_tasks[i+req_edge_num].demand;
	}

	for (int i = 1; i < req_edge_num+req_arc_num; i++)
	{
		for (int j = i+1; j <= req_edge_num+req_arc_num; j++)
		{
			if (task_routes[j].cost > task_routes[i].cost)
			{
				tmp_route = task_routes[i];
				task_routes[i] = task_routes[j];
				task_routes[j] = tmp_route;
			}
		}
	}

	/* augment phase */

	for (int i = 1; i < req_edge_num+req_arc_num; i++)
	{
		if (task_routes[i].perm[0] == 0)
			continue;

		for (int j = i+1; j <= req_edge_num+req_arc_num; j++)
		{
			if (task_routes[j].perm[0] == 0)
				continue;
			if (task_routes[j].load+task_routes[i].load > capacity)
				continue;

			for (int k = 1; k < task_routes[i].perm[0]; k++)
			{
				if (min_cost[inst_tasks[task_routes[i].perm[k]].tail_node][inst_tasks[task_routes[i].perm[k+1]].head_node] ==
					min_cost[inst_tasks[task_routes[i].perm[k]].tail_node][inst_tasks[task_routes[j].perm[2]].head_node]+
					inst_tasks[task_routes[j].perm[2]].dead_cost+
					min_cost[inst_tasks[task_routes[j].perm[2]].tail_node][inst_tasks[task_routes[i].perm[k+1]].head_node])
				{
					add_element(task_routes[i].perm, task_routes[j].perm[2], k+1);
					task_routes[i].load += task_routes[j].load;
					task_routes[j].perm[0] = 0;
					break;
				}
				else if (min_cost[inst_tasks[task_routes[i].perm[k]].tail_node][inst_tasks[task_routes[i].perm[k+1]].head_node] ==
					min_cost[inst_tasks[task_routes[i].perm[k]].tail_node][inst_tasks[task_routes[j].perm[2]].tail_node]+
					inst_tasks[task_routes[j].perm[2]].dead_cost+
					min_cost[inst_tasks[task_routes[j].perm[2]].head_node][inst_tasks[task_routes[i].perm[k+1]].head_node])
				{
					add_element(task_routes[i].perm, inst_tasks[task_routes[j].perm[2]].inverse, k+1);
					task_routes[i].load += task_routes[j].load;
					task_routes[j].perm[0] = 0;
					break;
				}
			}
		}
	}

	/* merge phase */

	typedef struct merge
	{
		int route1;
		int route2;
		int reverse;
		int saving;
	} merge;

	merge tmp_merge, best_merge;

	while(1)
	{
		best_merge.saving = 0;
		for (int i = 1; i <= req_edge_num+req_arc_num; i++)
		{
			if (task_routes[i].perm[0] == 0)
				continue;

			tmp_merge.route1 = i;
			for (int j = 1; j <= req_edge_num+req_arc_num; j++)
			{
				if (j == i)
					continue;
				if (task_routes[j].perm[0] == 0)
					continue;
				if (task_routes[j].load+task_routes[i].load > capacity)
					continue;

				tmp_merge.route2 = j;

				int flag = 1;
				for (int k = 2; k < task_routes[j].perm[0]; k++)
				{
					if (inst_tasks[task_routes[j].perm[k]].inverse == 0)
					{
						flag = 0;
						break;
					}
				}

				tmp_merge.reverse = 0;
				tmp_merge.saving =
					min_cost[inst_tasks[task_routes[i].perm[task_routes[i].perm[0]-1]].tail_node][DEPOT]+
					min_cost[DEPOT][inst_tasks[task_routes[j].perm[2]].head_node]-
					min_cost[inst_tasks[task_routes[i].perm[task_routes[i].perm[0]-1]].tail_node][inst_tasks[task_routes[j].perm[2]].head_node];

				if (tmp_merge.saving > best_merge.saving)
					best_merge = tmp_merge;

				if (flag)
				{
					tmp_merge.reverse = 1;
					tmp_merge.saving =
						min_cost[inst_tasks[task_routes[i].perm[task_routes[i].perm[0]-1]].tail_node][DEPOT]+
						min_cost[DEPOT][inst_tasks[task_routes[j].perm[task_routes[j].perm[0]-1]].tail_node]-
						min_cost[inst_tasks[task_routes[i].perm[task_routes[i].perm[0]-1]].tail_node][inst_tasks[task_routes[j].perm[task_routes[j].perm[0]-1]].tail_node];

					if (tmp_merge.saving > best_merge.saving)
						best_merge = tmp_merge;
				}
			}
		}

		if (best_merge.saving == 0)
			break;

		if (best_merge.reverse)
		{
			task_routes[best_merge.route1].perm[0] --;
			for (int i = task_routes[best_merge.route2].perm[0]-1; i > 0; i--)
			{
				task_routes[best_merge.route1].perm[0] ++;
				task_routes[best_merge.route1].perm[task_routes[best_merge.route1].perm[0]] = inst_tasks[task_routes[best_merge.route2].perm[i]].inverse;
			}
		}
		else
		{
			task_routes[best_merge.route1].perm[0] --;
			for (int i = 2; i <= task_routes[best_merge.route2].perm[0]; i++)
			{
				task_routes[best_merge.route1].perm[0] ++;
				task_routes[best_merge.route1].perm[task_routes[best_merge.route1].perm[0]] = task_routes[best_merge.route2].perm[i];
			}
		}

		task_routes[best_merge.route1].cost += task_routes[best_merge.route2].cost-best_merge.saving;
		task_routes[best_merge.route1].load += task_routes[best_merge.route2].load;
		task_routes[best_merge.route2].perm[0] = 0;
	}

	am_indi->sequence[0] = 1;
	am_indi->route_seg_length[0] = 0;
	am_indi->route_seg_load[0] = 0;
	am_indi->total_cost = 0;
	am_indi->max_length = 0;
	for (int i = 1; i <= req_edge_num+req_arc_num; i++)
	{
		if (task_routes[i].perm[0] == 0)
			continue;

		am_indi->sequence[0] --;
		link_array(am_indi->sequence, task_routes[i].perm);
		am_indi->route_seg_length[0] ++;
		am_indi->route_seg_length[am_indi->route_seg_length[0]] = task_routes[i].cost;
		am_indi->route_seg_load[0] ++;
		am_indi->route_seg_load[am_indi->route_seg_load[0]] = task_routes[i].load;
		am_indi->total_cost += task_routes[i].cost;
		if (task_routes[i].cost > am_indi->max_length)
			am_indi->max_length = task_routes[i].cost;
	}

	//print_one_dim_array(am_indi->sequence);
	//printf("%d, %d, %d\n", am_indi->total_cost, am_indi->max_length, am_indi->total_vio_load);
	//printf("%d\n", get_task_seq_total_cost(am_indi->sequence, inst_tasks));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ulusoy_heuristic(individual *uh_indi, const task *inst_tasks, int *serve_mark)
{
	int i, k;
	int serve_task_num = 0;
	for (i = req_edge_num+1; i <= task_num; i++)
	{
		if (serve_mark[i])
			serve_task_num ++;
	}

	int max_load = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		max_load += inst_tasks[i].demand;
	}
	
	int load, trial, mindist;

	int unserved_task[MAX_TASK_TAG_LENGTH], candi_task[MAX_TASK_TAG_LENGTH], nearest_task[MAX_TASK_TAG_LENGTH];
	int nearest_isol_task[MAX_TASK_TAG_LENGTH], nearest_inci_task[MAX_TASK_TAG_LENGTH], sel_task[MAX_TASK_TAG_LENGTH];
	int current_task, next_task;

	int positions[MAX_TASK_SEQ_LENGTH];
	
	individual tmp_indi1, tmp_indi2, tmp_indi3, tmp_indi4, tmp_indi5;
	chromosome tmp_chro1, tmp_chro2, tmp_chro3, tmp_chro4, tmp_chro5, uh_chro;

	uh_chro.total_cost = INF;

	int dep_dist[MAX_TASK_TAG_LENGTH];
	double yield[MAX_TASK_TAG_LENGTH];
	int max_dep_dist, min_dep_dist;
	double max_yield, min_yield;
	
	for (i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		dep_dist[i] = min_cost[inst_tasks[i].tail_node][DEPOT];
		yield[i] = 1.0*inst_tasks[i].demand/inst_tasks[i].serv_cost;
	}
	
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    tmp_indi1.sequence[0] = 1;
	tmp_indi1.sequence[1] = 0;
	tmp_indi1.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi1.sequence[tmp_indi1.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= max_load-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			tmp_indi1.sequence[0] ++;
	  		tmp_indi1.sequence[tmp_indi1.sequence[0]] = 0;
		  	tmp_indi1.route_seg_load[0] ++;
			tmp_indi1.route_seg_load[tmp_indi1.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		max_dep_dist = -1;
		sel_task[0] = 0;
		for (i = 1; i <= nearest_isol_task[0]; i++)
		{
			if (dep_dist[nearest_isol_task[i]] > max_dep_dist)
			{
				max_dep_dist = dep_dist[nearest_isol_task[i]];
				sel_task[0] = 1;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
			else if (dep_dist[nearest_isol_task[i]] == max_dep_dist)
			{
				sel_task[0] ++;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi1.sequence[0] ++;
		tmp_indi1.sequence[tmp_indi1.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi1.sequence[0] ++;
	tmp_indi1.sequence[tmp_indi1.sequence[0]] = 0;
	tmp_indi1.route_seg_load[0] ++;
	tmp_indi1.route_seg_load[tmp_indi1.route_seg_load[0]] = load;

	memcpy(tmp_chro1.sequence, tmp_indi1.sequence, sizeof(tmp_indi1.sequence));
	remove_task_seq_delimiters(tmp_chro1.sequence);
	split(&tmp_indi1, &tmp_chro1, inst_tasks);
	
	if (tmp_chro1.total_cost < uh_chro.total_cost)
	{
		memcpy(uh_chro.sequence, tmp_chro1.sequence, (tmp_chro1.sequence[0]+1)*sizeof(int));
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	tmp_indi2.sequence[0] = 1;
	tmp_indi2.sequence[1] = 0;
	tmp_indi2.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi2.sequence[tmp_indi2.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= max_load-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
			tmp_indi2.sequence[0] ++;
	  		tmp_indi2.sequence[tmp_indi2.sequence[0]] = 0;
		  	tmp_indi2.route_seg_load[0] ++;
			tmp_indi2.route_seg_load[tmp_indi2.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		min_dep_dist = INF;
		sel_task[0] = 0;
		for (i = 1; i <= nearest_isol_task[0]; i++)
		{
			if (dep_dist[nearest_isol_task[i]] < min_dep_dist)
			{
				min_dep_dist = dep_dist[nearest_isol_task[i]];
				sel_task[0] = 1;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
			else if (dep_dist[nearest_isol_task[i]] == min_dep_dist)
			{
				sel_task[0] ++;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi2.sequence[0] ++;
		tmp_indi2.sequence[tmp_indi2.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi2.sequence[0] ++;
	tmp_indi2.sequence[tmp_indi2.sequence[0]] = 0;
	tmp_indi2.route_seg_load[0] ++;
	tmp_indi2.route_seg_load[tmp_indi2.route_seg_load[0]] = load;
		
	memcpy(tmp_chro2.sequence, tmp_indi2.sequence, sizeof(tmp_indi2.sequence));
	remove_task_seq_delimiters(tmp_chro2.sequence);
	split(&tmp_indi2, &tmp_chro2, inst_tasks);
	
	if (tmp_chro2.total_cost < uh_chro.total_cost)
	{
		memcpy(uh_chro.sequence, tmp_chro2.sequence, (tmp_chro2.sequence[0]+1)*sizeof(int));
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	tmp_indi3.sequence[0] = 1;
	tmp_indi3.sequence[1] = 0;
	tmp_indi3.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi3.sequence[tmp_indi3.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= max_load-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			tmp_indi3.sequence[0] ++;
	  		tmp_indi3.sequence[tmp_indi3.sequence[0]] = 0;
		  	tmp_indi3.route_seg_load[0] ++;
			tmp_indi3.route_seg_load[tmp_indi3.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		max_yield = -1;
		sel_task[0] = 0;
		for (i = 1; i <= nearest_isol_task[0]; i++)
		{
			if (yield[nearest_isol_task[i]] > max_yield)
			{
				max_yield = yield[nearest_isol_task[i]];
				sel_task[0] = 1;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
			else if (yield[nearest_isol_task[i]] == max_yield)
			{
				sel_task[0] ++;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi3.sequence[0] ++;
		tmp_indi3.sequence[tmp_indi3.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi3.sequence[0] ++;
	tmp_indi3.sequence[tmp_indi3.sequence[0]] = 0;
	tmp_indi3.route_seg_load[0] ++;
	tmp_indi3.route_seg_load[tmp_indi3.route_seg_load[0]] = load;
		
	memcpy(tmp_chro3.sequence, tmp_indi3.sequence, sizeof(tmp_indi3.sequence));
	remove_task_seq_delimiters(tmp_chro3.sequence);
	split(&tmp_indi3, &tmp_chro3, inst_tasks);
	
	if (tmp_chro3.total_cost < uh_chro.total_cost)
	{
		memcpy(uh_chro.sequence, tmp_chro3.sequence, (tmp_chro3.sequence[0]+1)*sizeof(int));
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	tmp_indi4.sequence[0] = 1;
	tmp_indi4.sequence[1] = 0;
	tmp_indi4.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi4.sequence[tmp_indi4.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= max_load-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			tmp_indi4.sequence[0] ++;
	  		tmp_indi4.sequence[tmp_indi4.sequence[0]] = 0;
		  	tmp_indi4.route_seg_load[0] ++;
			tmp_indi4.route_seg_load[tmp_indi4.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		min_yield = INF;
		sel_task[0] = 0;
		for (i = 1; i <= nearest_isol_task[0]; i++)
		{
			if (yield[nearest_isol_task[i]] < min_yield)
			{
				min_yield = yield[nearest_isol_task[i]];
				sel_task[0] = 1;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
			else if (yield[nearest_isol_task[i]] == min_yield)
			{
				sel_task[0] ++;
				sel_task[sel_task[0]] = nearest_isol_task[i];
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi4.sequence[0] ++;
		tmp_indi4.sequence[tmp_indi4.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi4.sequence[0] ++;
	tmp_indi4.sequence[tmp_indi4.sequence[0]] = 0;
	tmp_indi4.route_seg_load[0] ++;
	tmp_indi4.route_seg_load[tmp_indi4.route_seg_load[0]] = load;
		
	memcpy(tmp_chro4.sequence, tmp_indi4.sequence, sizeof(tmp_indi4.sequence));
	remove_task_seq_delimiters(tmp_chro4.sequence);
	split(&tmp_indi4, &tmp_chro4, inst_tasks);
	
	if (tmp_chro4.total_cost < uh_chro.total_cost)
	{
		memcpy(uh_chro.sequence, tmp_chro4.sequence, (tmp_chro4.sequence[0]+1)*sizeof(int));
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    tmp_indi5.sequence[0] = 1;
	tmp_indi5.sequence[1] = 0;
	tmp_indi5.route_seg_load[0] = 0;
	
	unserved_task[0] = 0;
	for(i = 1; i <= task_num; i++)
	{
		if (!serve_mark[i])
			continue;
		
		unserved_task[0] ++;
		unserved_task[unserved_task[0]] = i;
	}
	
	load = 0;
	trial = 0;
	while (trial < serve_task_num)
	{
		current_task = tmp_indi5.sequence[tmp_indi5.sequence[0]];
		
		candi_task[0] = 0;
		
		for (i = 1; i <= unserved_task[0]; i++)
		{
			if (inst_tasks[unserved_task[i]].demand <= max_load-load)
			{
				candi_task[0] ++;
				candi_task[candi_task[0]] = unserved_task[i];
			}
		}
		
		if (candi_task[0] == 0)
		{
  			tmp_indi5.sequence[0] ++;
	  		tmp_indi5.sequence[tmp_indi5.sequence[0]] = 0;
		  	tmp_indi5.route_seg_load[0] ++;
			tmp_indi5.route_seg_load[tmp_indi5.route_seg_load[0]] = load;
			load = 0;
			continue;
		}
		
		mindist = INF;
		nearest_task[0] = 0;
		
		for (i = 1; i <= candi_task[0]; i++)
		{
			if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] < mindist)
			{
				mindist = min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node];
				nearest_task[0] = 1;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
			else if (min_cost[inst_tasks[current_task].tail_node][inst_tasks[candi_task[i]].head_node] == mindist)
			{
				nearest_task[0] ++;
				nearest_task[nearest_task[0]] = candi_task[i];
			}
		}
		
		nearest_isol_task[0] = 0;
		nearest_inci_task[0] = 0;
		
		for (i = 1; i <= nearest_task[0]; i++)
		{
			if (inst_tasks[nearest_task[i]].tail_node == 1)
			{
				nearest_inci_task[0] ++;
				nearest_inci_task[nearest_inci_task[0]] = nearest_task[i];
			}
			else
			{
				nearest_isol_task[0] ++;
				nearest_isol_task[nearest_isol_task[0]] = nearest_task[i];
			}
		}
		
		if (nearest_isol_task[0] == 0)
		{
			memcpy(nearest_isol_task, nearest_inci_task, (nearest_inci_task[0]+1)*sizeof(int));
		}
		
		if (load < max_load/2)
		{
			max_dep_dist = -1;
			sel_task[0] = 0;
			for (i = 1; i <= nearest_isol_task[0]; i++)
			{
				if (dep_dist[nearest_isol_task[i]] > max_dep_dist)
				{
					max_dep_dist = dep_dist[nearest_isol_task[i]];
					sel_task[0] = 1;
					sel_task[sel_task[0]] = nearest_isol_task[i];
				}
				else if (dep_dist[nearest_isol_task[i]] == max_dep_dist)
				{
					sel_task[0] ++;
					sel_task[sel_task[0]] = nearest_isol_task[i];
				}
			}
		}
		else
		{
			min_dep_dist = INF;
			sel_task[0] = 0;
			for (i = 1; i <= nearest_isol_task[0]; i++)
			{
				if (dep_dist[nearest_isol_task[i]] < min_dep_dist)
				{
					min_dep_dist = dep_dist[nearest_isol_task[i]];
					sel_task[0] = 1;
					sel_task[sel_task[0]] = nearest_isol_task[i];
				}
				else if (dep_dist[nearest_isol_task[i]] == min_dep_dist)
				{
					sel_task[0] ++;
					sel_task[sel_task[0]] = nearest_isol_task[i];
				}
			}
		}
		
		//k = rand_choose(sel_task[0]);
		k = 1;
		next_task = sel_task[k];
		
		trial ++;
		tmp_indi5.sequence[0] ++;
		tmp_indi5.sequence[tmp_indi5.sequence[0]] = next_task;
		load += inst_tasks[next_task].demand;
		find_ele_positions(positions, unserved_task, next_task);
		delete_element(unserved_task, positions[1]);
		
		if (inst_tasks[next_task].inverse > 0)
		{
			find_ele_positions(positions, unserved_task, inst_tasks[next_task].inverse);
			delete_element(unserved_task, positions[1]);
		}
	}
	
	tmp_indi5.sequence[0] ++;
	tmp_indi5.sequence[tmp_indi5.sequence[0]] = 0;
	tmp_indi5.route_seg_load[0] ++;
	tmp_indi5.route_seg_load[tmp_indi5.route_seg_load[0]] = load;
	
	memcpy(tmp_chro5.sequence, tmp_indi5.sequence, sizeof(tmp_indi5.sequence));
	remove_task_seq_delimiters(tmp_chro5.sequence);
	split(&tmp_indi5, &tmp_chro5, inst_tasks);
	
	if (tmp_chro5.total_cost < uh_chro.total_cost)
	{
		memcpy(uh_chro.sequence, tmp_chro5.sequence, (tmp_chro5.sequence[0]+1)*sizeof(int));
	}

	split(uh_indi, &uh_chro, inst_tasks);
	uh_indi->total_vio_load = 0;
	get_route_seg_length(uh_indi->route_seg_length, uh_indi->sequence, inst_tasks);
	uh_indi->max_length = max(uh_indi->route_seg_length);
}