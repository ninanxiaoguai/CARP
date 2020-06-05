
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int* rand_perm(int num)
{
	int *a = (int *)malloc((num+1)*sizeof(int));
	int *left_ele = (int *)malloc((num+1)*sizeof(int));
	left_ele[0] = num;
	for (int i = 1; i <= num; i++)
	{
		left_ele[i] = i;
	}

	a[0] = num;
	for (int i = 1; i <= num; i++)
	{
		int k = rand_choose(left_ele[0]);
		a[i] = left_ele[k];
		delete_element(left_ele, k);
	}

	free(left_ele);

	return a;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int rand_choose(int num)
{
	int k = rand()%num;

	k++;

	return k;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_one_dim_array(int *a)
{
	for (int i = 1; i <= a[0]; i++)
	{
		printf("%d ", a[i]);
	}
	printf("\n");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_two_dim_matrix(int **a, int row, int col)
{
	for (int i = 1; i <= row; i++)
	{
		for (int j = 1; j <= col; j++)
		{
			printf("%d ", a[i][j]);
		}
		printf("\n");
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void delete_element(int *a, int k)
{
	if (k < 1 || k > a[0])
	{
		printf("the deleting position is wrong!\n");
		exit(0);
	}

	for (int i = k; i < a[0]; i++)
	{
		a[i] = a[i+1];
	}
	a[0] --;
}

int delete_element(individual *indis, int indis_size, int k)
{
	if (k < 0 || k > indis_size)
	{
		printf("the deleting position is wrong!\n");
		exit(0);
	}

	for (int i = k; i < indis_size; i++)
	{
		indi_copy(&indis[i], &indis[i+1]);
	}
	int tmp = indis_size;
	tmp --;
	return tmp;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void add_element(int *a, int e, int k)
{
	if (k < 1 || k > a[0]+1)
	{
		printf("the inserting position is wrong!\n");
		exit(0);
	}

	a[0] ++;
	for (int i = a[0]; i > k; i--)
	{
		a[i] = a[i-1];
	}
	a[k] = e;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void reverse_direction(int *a, int k1, int k2)
{
	if (k1 > k2 || k1 < 1 || k2 > a[0])
	{
		printf("the reversing positions are wrong!\n");
		exit(0);
	}
	
	int half = (k2-k1+1)/2;

	for (int i = k1; i < k1+half; i++)
	{
		int tmp = a[i];
		a[i] = a[k1+k2-i];
		a[k1+k2-i] = tmp;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void find_ele_positions(int *positions, int *a, int e)
{
	positions[0] = 0;
	for (int i = 1; i <= a[0]; i++)
	{
		if (a[i] == e)
		{
			positions[0] ++;
			positions[positions[0]] = i;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void copy_sub_array(int *b, int *a, int k1, int k2)
{
	b[0] = k2-k1+1;
	for (int i = k1; i <= k2; i++)
	{
		b[i-k1+1] = a[i];
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void insert_array(int *b, int *a, int k)
{
	if (k < 1 || k > b[0]+1)
	{
		printf("the inserting position is wrong!\n");
		exit(0);
	}

	for (int i = b[0]; i >= k; i--)
	{
		b[i+a[0]] = b[i];
	}

	for (int i = 1; i <= a[0]; i++)
	{
		b[k+i-1] = a[i];
	}

	b[0] += a[0];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void replace_sub_array(int *b, int *a, int k1, int k2)
{
	if (k1 < k2)
	{
		printf("the replace positions are wrong!\n");
		exit(0);
	}

	for (int i = k2+1; i <= b[0]; i++)
	{
		a[a[0]+i-k2] = b[i];
	}

	b[0] += a[0]-(k2-k1+1);
	for (int i = 1; i <= a[0]+b[0]-k2; i++)
	{
		b[i+k1-1] = a[i];
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void link_array(int *b, int *a)
{
	for (int i = 1; i <= a[0]; i++)
	{
		b[b[0]+i] = a[i];
	}
	b[0] += a[0];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int equal(int *a, int *b)
{
	if (a[0] != b[0])
		return 0;
	
	for (int i = 1; i < a[0]; i++)
	{
		if (a[i] != b[i])
			return 0;
	}

	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int min(int *a)
{
	int min_val = INF;
	for (int i = 1; i <= a[0]; i++)
	{
		if (a[i] < min_val)
			min_val = a[i];
	}
	
	return min_val;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int max(int *a)
{
	int max_val = -INF;
	for (int i = 1; i <= a[0]; i++)
	{
		if (a[i] > max_val)
			max_val = a[i];
	}
	
	return max_val;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int set_bit(int num, int k)
{
	int set_num, unit;
	unit = 1;
	set_num = num | (unit<<(k-1));
	return set_num;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int unset_bit(int num, int k)
{
	if (!bit_one(num, k))
	{
		printf("num = %d, k = %d\n", num, k);
		return 0;
	}
	int set_num, zero;
	zero = 1;
	set_num = num & ~(zero<<(k-1));
	return set_num;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int bit_one(int num, int k)
{
	int check, unit;
	unit = 1;
	check = num & (unit<<(k-1));
	return (check != 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int find_arc(int head_node, int tail_node, const arc *inst_arcs)
{
	for (int i = 1; i < total_arc_num; i++)
	{
		if (inst_arcs[i].head_node == head_node && inst_arcs[i].tail_node == tail_node)
			return i;
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int find_task(int head_node, int tail_node, const task *inst_tasks)
{
	for (int i = 1; i <= task_num; i++)
	{
		if (inst_tasks[i].head_node == head_node && inst_tasks[i].tail_node == tail_node)
			return i;
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int indi_cmp(individual *indi1, individual *indi2)
// The two compared individuals are with the same amount of constraint violation, i.e., indi1->total_vio_load = indi2->total_vio_load.
{
	// compare "total_vio_load" first
	if (indi1->total_vio_load < indi2->total_vio_load)
		return -1;

	if (indi1->total_vio_load > indi2->total_vio_load)
		return 1;

	int better = 0, equal = 0, worse = 0;
	if (indi1->max_length < indi2->max_length)
	{
		better ++;
	}
	else if (indi1->max_length == indi2->max_length)
	{
		equal ++;
	}
	else
	{
		worse ++;
	}
	if (indi1->total_cost < indi2->total_cost)
	{
		better ++;
	}
	else if (indi1->total_cost == indi2->total_cost)
	{
		equal ++;
	}
	else
	{
		worse ++;
	}

	if (worse == 0 && better > 0) // indi1 dominates indi2
	{
		return -1;
	}
	else if (better == 0 && worse > 0) // indi1 is dominated by indi2
	{
		return 1;
	}
	else if (better == 0 && worse == 0) // indi1 equals indi2
	{
		return -2;
	}
	else // indi1 and indi2 are nondominated by each other
	{
		return 0;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 判断 indi 是否与nondominated_indis内的前nondominated_indis_size个构成非支配关系，进行删减添加，返回当前的非支配解集
int  modify_nondominated_indis(individual nondominated_indis[MAX_NONDOMINATED_NUM], int nondominated_indis_size, individual *indi)
{
	int tmp_size = nondominated_indis_size;
	int dominated = 0;
	for (int i = nondominated_indis_size-1; i >= 0; i--)
	{
		if (indi_cmp(&nondominated_indis[i], indi) == -1 || indi_cmp(&nondominated_indis[i], indi) == -2) 
		/* indi is dominated */
		{
			dominated = 1;
			break;
		}
		else if (indi_cmp(&nondominated_indis[i], indi) == 1)
		{ // nondominated_indis[i] 被 indi 支配，因此在列表中删掉它
			tmp_size = delete_element(nondominated_indis, tmp_size, i);
		}
	}
	if (!dominated)
	{// 添加 indi
		//printf("tmp_size = %d\n", tmp_size);
		//printf("%d, %d, %d\n", indi->max_length, indi->total_cost, indi->total_vio_load);
		indi_copy(&nondominated_indis[tmp_size], indi);
		tmp_size ++;
	}

	max_ndtotalcost = 0;
	min_ndtotalcost = INF;
	max_ndmaxlength = 0;
	min_ndmaxlength = INF;
	for (int i = 0; i < tmp_size; i++)
	{
		if (nondominated_indis[i].max_length > max_ndmaxlength)
			max_ndmaxlength = nondominated_indis[i].max_length;
		if (nondominated_indis[i].max_length < min_ndmaxlength)
			min_ndmaxlength = nondominated_indis[i].max_length;
		if (nondominated_indis[i].total_cost > max_ndtotalcost)
			max_ndtotalcost = nondominated_indis[i].total_cost;
		if (nondominated_indis[i].total_cost < min_ndtotalcost)
			min_ndtotalcost = nondominated_indis[i].total_cost;
	}

	return tmp_size;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void remove_task_seq_delimiters(int *task_seq)
{
	for (int i = task_seq[0]; i > 0; i--)
	{
		if (task_seq[i] == 0)
			delete_element(task_seq, i);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void split(individual *indi, chromosome *chro, const task *inst_tasks)
{
	int V[MAX_TASK_SEQ_LENGTH], P[MAX_TASK_SEQ_LENGTH];
	V[0] = 0;
	P[0] = 0;

	for (int i = 1; i <= chro->sequence[0]; i++)
	{
		V[i] = INF;
	}

	for (int i = 1; i <= chro->sequence[0]; i++)
	{
		int load = 0, cost = 0;
		int j = i;

		while (j <= chro->sequence[0] && load <= capacity)
		{
			load += inst_tasks[chro->sequence[j]].demand;

			if (j == i)
			{
				cost = min_cost[DEPOT][inst_tasks[chro->sequence[j]].head_node]+inst_tasks[chro->sequence[j]].serv_cost
					+min_cost[inst_tasks[chro->sequence[j]].tail_node][DEPOT];
			}
			else
			{
				cost += min_cost[inst_tasks[chro->sequence[j-1]].tail_node][inst_tasks[chro->sequence[j]].head_node]
				+inst_tasks[chro->sequence[j]].serv_cost+min_cost[inst_tasks[chro->sequence[j]].tail_node][DEPOT]
				-min_cost[inst_tasks[chro->sequence[j-1]].tail_node][DEPOT];
			}

			if (load <= capacity)
			{
				int V_new = V[i-1]+cost;

				if (V_new < V[j])
				{
					V[j] = V_new;
					P[j] = i-1;
				}

				j ++;
			}
		}
	}

	indi->sequence[0] = 1;
	indi->sequence[1] = 0;
	int j = chro->sequence[0];
	int ptr = P[j];

	while (ptr > 0)
	{
		for (int k = ptr+1; k <= j; k++)
		{
			indi->sequence[0] ++;
			indi->sequence[indi->sequence[0]] = chro->sequence[k];
		}

		indi->sequence[0] ++;
		indi->sequence[indi->sequence[0]] = 0;

		j = ptr;
		ptr = P[j];
	}

	for (int k = 1; k <= j; k++)
	{
		indi->sequence[0] ++;
		indi->sequence[indi->sequence[0]] = chro->sequence[k];
	}

	indi->sequence[0] ++;
	indi->sequence[indi->sequence[0]] = 0;

	get_route_seg_length(indi->route_seg_length, indi->sequence, inst_tasks);
	get_route_seg_load(indi->route_seg_load, indi->sequence, inst_tasks);
	indi->total_cost = V[chro->sequence[0]];
	chro->total_cost = indi->total_cost;
	indi->max_length = 0;
	for (int i = 1; i <= indi->route_seg_length[0]; i++)
	{
		if (indi->route_seg_length[i] > indi->max_length)
			indi->max_length = indi->route_seg_length[i];
	}
	chro->max_length = indi->max_length;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int fleet_limited_split(int *split_task_seq, int *one_task_seq, const task *inst_tasks)
{
	int aux_cost[MAX_TASK_SEQ_LENGTH][MAX_TASK_ROUTE_LENGTH], opt_cost[MAX_TASK_SEQ_LENGTH][MAX_ROUTE_TAG_LENGTH];
	int pred_task[MAX_TASK_SEQ_LENGTH][MAX_ROUTE_TAG_LENGTH], tmp_task_seq[MAX_TASK_SEQ_LENGTH];

	for (int i = 0; i <= one_task_seq[0]; i++)
	{
		for (int j = 0; j <= one_task_seq[0]; j++)
		{
			aux_cost[i][j] = INF;
		}
	}

	for (int i = 0; i <= one_task_seq[0]; i++)
	{
		for (int j = 0; j <= one_task_seq[0]; j++)
		{
			printf("%d ", aux_cost[i][j]);
		}
		printf("\n");
	}

	for (int i = 0; i <= one_task_seq[0]; i++)
	{
		for (int j = 0; j <= vehicle_num; j++)
		{
			opt_cost[i][j] = INF;
		}
	}
	opt_cost[0][0] = 0;

	for (int i = 0; i <= one_task_seq[0]; i++)
	{
		int cost;
		int load = 0;
		int j = i;

		while (j <= one_task_seq[0])
		{
			load += inst_tasks[one_task_seq[j]].demand;

			if (load >capacity)
				break;

			if (j == i)
			{
				cost = min_cost[DEPOT][inst_tasks[one_task_seq[j]].head_node]+inst_tasks[one_task_seq[j]].serv_cost
					+min_cost[inst_tasks[one_task_seq[j]].tail_node][DEPOT];
			}
			else
			{
				cost += min_cost[inst_tasks[one_task_seq[j-1]].tail_node][inst_tasks[one_task_seq[j]].head_node]
				+inst_tasks[one_task_seq[j]].serv_cost+min_cost[inst_tasks[one_task_seq[j]].tail_node][DEPOT]
				-min_cost[inst_tasks[one_task_seq[j-1]].tail_node][DEPOT];
			}

			aux_cost[i][j] = cost;

			j ++;
		}
	}

	for (int i = 1; i <= vehicle_num; i++)
	{
		for (int j = 1; j <= one_task_seq[0]; j++)
		{
			for (int k = 0; k < j; k++)
			{
				if (opt_cost[k][i-1] < INF && opt_cost[k][i-1]+aux_cost[k+1][j] < opt_cost[j][i])
				{
					opt_cost[j][i] = opt_cost[k][i-1]+aux_cost[k+1][j];
					pred_task[j][i] = k;
				}
			}
		}
	}

	int all_opt_num;
	int all_opt_cost = INF;
	for (int i = 1; i <= vehicle_num; i++)
	{
		if (opt_cost[one_task_seq[0]][i] < all_opt_cost)
		{
			all_opt_cost = opt_cost[one_task_seq[0]][i];
			all_opt_num = i;
		}
	}

	if (all_opt_cost == INF)
	{
		split_task_seq[0] = 0;
		return all_opt_cost;
	}

	tmp_task_seq[0] = 1;
	tmp_task_seq[tmp_task_seq[0]] = one_task_seq[0];
	while (1)
	{
		tmp_task_seq[0] ++;
		tmp_task_seq[tmp_task_seq[0]] = pred_task[tmp_task_seq[tmp_task_seq[0]-1]][all_opt_num];
		all_opt_num --;

		if (tmp_task_seq[tmp_task_seq[0]] == 0)
			break;
	}

	split_task_seq[0] = 1;
	split_task_seq[1] = 0;

	for (int i = tmp_task_seq[0]; i > 1; i--)
	{
		for (int j = tmp_task_seq[i]+1; j <= tmp_task_seq[i-1]; j++)
		{
			split_task_seq[0] ++;
			split_task_seq[split_task_seq[0]] = one_task_seq[j];
		}

		split_task_seq[0] ++;
		split_task_seq[split_task_seq[0]] = 0;
	}

	return all_opt_cost;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
int legal(int *task_seq, const task *inst_tasks)
{
	for (int i = 1; i < task_seq[0]; i++)
	{
		if (task_seq[i] == 0 && task_seq[i+1] == 0)
			return 0;
	}

	int tmp_task_seq[MAX_TASK_SEQ_LENGTH];
	memcpy(tmp_task_seq, task_seq, sizeof(task_seq));
	remove_task_seq_delimiters(tmp_task_seq);

	if (tmp_task_seq[0] < req_edge_num+req_arc_num)
		return 0;

	for (int i = 1; i < tmp_task_seq[0]; i++)
	{
		for (int j = i+1; j <= tmp_task_seq[0]; j++)
		{
			if (tmp_task_seq[j] == tmp_task_seq[i] || tmp_task_seq[j] == inst_tasks[tmp_task_seq[i]].inverse)
				return 0;
		}
	}

	//free(tmp_task_seq);

	return 1;
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int get_task_seq_total_cost(int *task_seq, const task *inst_tasks)
{
	int total_cost =  0;
	for (int i = 1; i < task_seq[0]; i++)
	{
		total_cost += min_cost[inst_tasks[task_seq[i]].tail_node][inst_tasks[task_seq[i+1]].head_node]+inst_tasks[task_seq[i]].serv_cost;
	}

	return total_cost;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void get_route_seg_load(int *route_seg_load, int *task_seq, const task *inst_tasks)
{
	route_seg_load[0] = 1;
	route_seg_load[1] = 0;
	for (int i = 2; i < task_seq[0]; i++)
	{
		if (task_seq[i] == 0)
		{
			route_seg_load[0] ++;
			route_seg_load[route_seg_load[0]] = 0;
			continue;
		}

		route_seg_load[route_seg_load[0]] += inst_tasks[task_seq[i]].demand;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int get_total_vio_load(int *route_seg_load)
{
	int total_vio_load = 0;
	for (int i = 1; i <= route_seg_load[0]; i++)
	{
		if (route_seg_load[i] > capacity)
			total_vio_load += route_seg_load[i]-capacity;
	}
	
	return total_vio_load;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void get_route_seg_length(int *route_seg_length, int *task_seq, const task *inst_tasks)
{
	route_seg_length[0] = 1;
	route_seg_length[1] = 0;
	for (int i = 2; i <= task_seq[0]; i++)
	{
		route_seg_length[route_seg_length[0]] += min_cost[inst_tasks[task_seq[i-1]].tail_node][inst_tasks[task_seq[i]].head_node]+inst_tasks[task_seq[i]].serv_cost;
		if (task_seq[i] == 0)
		{
			route_seg_length[0] ++;
			route_seg_length[route_seg_length[0]] = 0;
		}
	}
	route_seg_length[0] --;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void get_assignment(int *assignment, int *task_seq, const task *inst_tasks)
{
	assignment[0] = req_edge_num+req_arc_num;
	int curr_route = 1;

	for (int i = 2; i <= task_seq[0]; i++)
	{
		if (task_seq[i] == 0)
		{
			curr_route ++;
			continue;
		}

		int tmp_task = task_seq[i];
    	if (tmp_task > req_edge_num)
	    	tmp_task -= req_edge_num;
    
	    assignment[tmp_task] = curr_route;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void indi_copy(individual *target, individual *source)
{
	memcpy(target->sequence, source->sequence, (source->sequence[0]+1)*sizeof(int));
	memcpy(target->route_seg_load, source->route_seg_load, (source->route_seg_load[0]+1)*sizeof(int));
	memcpy(target->route_seg_length, source->route_seg_length, (source->route_seg_length[0]+1)*sizeof(int));
	target->total_cost = source->total_cost;
	target->max_length = source->max_length;
	target->total_vio_load = source->total_vio_load;
	target->fitness = source->fitness;
	target->rank = source->rank;
	target->crowdness = source->crowdness;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*void TransTaskSeqToOrderedList(int (*OrderedList)[250], int *TaskSeq, Task *ARPTask)
{
	int i, j, k;
	
	OrderedList[0][0] = 1;
	
	for (i = 2; i <= TaskSeq[0]; i++)
	{
		if (TaskSeq[i] == 0)
		{
 			for (j = 1; j <= ShortestPath[ARPTask[TaskSeq[i-1]].Head][ARPTask[TaskSeq[i]].Tail][0]; j++)
  		{
	  		OrderedList[OrderedList[0][0]][0] ++;
		  	OrderedList[OrderedList[0][0]][OrderedList[OrderedList[0][0]][0]] =
		  	ShortestPath[ARPTask[TaskSeq[i-1]].Head][ARPTask[TaskSeq[i]].Tail][j];
		  }
			
			OrderedList[0][0] ++;
			OrderedList[OrderedList[0][0]][0] = 0;
		}
		else
		{
			for (j = 1; j <= ShortestPath[ARPTask[TaskSeq[i-1]].Head][ARPTask[TaskSeq[i]].Tail][0]; j++)
	 		{
	  		OrderedList[OrderedList[0][0]][0] ++;
		  	OrderedList[OrderedList[0][0]][OrderedList[OrderedList[0][0]][0]] =
		  	ShortestPath[ARPTask[TaskSeq[i-1]].Head][ARPTask[TaskSeq[i]].Tail][j];
		  }
		}
	}
	
	OrderedList[0][0] --;
}

void TransTaskRouteToOrderedRoute(int *OrderedRoute, int *TaskRoute, Task *ARPTask)
{
	int i, j, k;
	
	OrderedRoute[0] = 0;
	
	for (i = 2; i <= TaskRoute[0]; i++)
	{
		for (j = 1; j <= ShortestPath[ARPTask[TaskRoute[i-1]].Head][ARPTask[TaskRoute[i]].Tail][0]; j++)
		{
	 		OrderedRoute[0] ++;
	  	OrderedRoute[OrderedRoute[0]] = ShortestPath[ARPTask[TaskRoute[i-1]].Head][ARPTask[TaskRoute[i]].Tail][j];
	  }
	}
}*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*int CalcOrderedRouteTotalCost(int *OrderedRoute)
{
	int i;
	int TotalCost = 0;
	
	for (i = 1; i < OrderedRoute[0]; i++)
	{
		TotalCost += Cost[OrderedRoute[i]][OrderedRoute[i+1]];
	}
	
	return TotalCost;
}*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*void GetROccurrence(int *ROccurrence, int *OrderedRoute, int p1, int p2, Task *ARPTask)
{
	int i;
	int TmpTask;
	
	for (i = p1; i < p2; i++)
  {
    TmpTask = FindTask(OrderedRoute[i], OrderedRoute[i+1], ARPTask);
    
    if (TmpTask == 0)
    	continue;
    
    if (TmpTask > NO_ReqEdge)
    	TmpTask -= NO_ReqEdge;
    
    ROccurrence[TmpTask] ++;
  }
}*/
