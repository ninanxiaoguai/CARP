
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"

int ad_matrix[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int tad_matrix[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int match_matrix[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int min_span_tree[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int rep_node1[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int rep_node2[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
graph piece_graph;
graph deg_graph;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void fred_heuristic(int *fh_route, int *orig_route, const task *inst_tasks, int eta)
{
	//printf("orig_route\n");
	//print_one_dim_array(orig_route);
	memset(ad_matrix, 0, sizeof(ad_matrix));
	memset(tad_matrix, 0, sizeof(tad_matrix));

	int orig_deg[MAX_NODE_TAG_LENGTH];
	memset(orig_deg, 0, sizeof(orig_deg));
	for (int i = 1; i < vertex_num; i++)
	{
		for (int j = i+1; j <= vertex_num; j++)
		{
			if (trav_cost[i][j] < INF)
			{
				orig_deg[i] ++;
				orig_deg[j] ++;
			}
		}
	}

	typedef struct node_char
	{
		int inci_stat;
		int mark_stat;
		int scan_stat;
		int neighbors[MAX_NODE_TAG_LENGTH-1];
	} node_char;

	node_char g_node_char[MAX_NODE_TAG_LENGTH];

	for (int i = 1; i <= vertex_num; i++)
	{
		g_node_char[i].inci_stat = 0;
		g_node_char[i].neighbors[0] = 0;
	}

	g_node_char[1].inci_stat = 1;
	g_node_char[1].mark_stat = 0;
	g_node_char[1].scan_stat = 0;

	for (int i = 2; i < orig_route[0]; i++)
	{
		g_node_char[inst_tasks[orig_route[i]].head_node].inci_stat = 1;
		g_node_char[inst_tasks[orig_route[i]].head_node].mark_stat = 0;
		g_node_char[inst_tasks[orig_route[i]].head_node].scan_stat = 0;
		g_node_char[inst_tasks[orig_route[i]].head_node].neighbors[0] ++;
		g_node_char[inst_tasks[orig_route[i]].head_node].neighbors[g_node_char[inst_tasks[orig_route[i]].head_node].neighbors[0]] =
			inst_tasks[orig_route[i]].tail_node;

		g_node_char[inst_tasks[orig_route[i]].tail_node].inci_stat = 1;
		g_node_char[inst_tasks[orig_route[i]].tail_node].mark_stat = 0;
		g_node_char[inst_tasks[orig_route[i]].tail_node].scan_stat = 0;
		g_node_char[inst_tasks[orig_route[i]].tail_node].neighbors[0] ++;
		g_node_char[inst_tasks[orig_route[i]].tail_node].neighbors[g_node_char[inst_tasks[orig_route[i]].tail_node].neighbors[0]] =
			inst_tasks[orig_route[i]].head_node;

		ad_matrix[inst_tasks[orig_route[i]].head_node][inst_tasks[orig_route[i]].tail_node] = 1;
		ad_matrix[inst_tasks[orig_route[i]].tail_node][inst_tasks[orig_route[i]].head_node] = 1;
		tad_matrix[inst_tasks[orig_route[i]].head_node][inst_tasks[orig_route[i]].tail_node] = 1;
		tad_matrix[inst_tasks[orig_route[i]].tail_node][inst_tasks[orig_route[i]].head_node] = 1;
	}

	//for (int i = 2; i < orig_route[0]; i++)
	//{
	//	printf("%d, %d\n", inst_tasks[orig_route[i]].head_node, inst_tasks[orig_route[i]].tail_node);
	//}
	//print_two_dim_matrix(ad_matrix, vertex_num, vertex_num);

	/* finding connected pieces */
	int curr_mark = 0;
	while (1)
	{
		int sel_node;
		int flag = 0;
		for (int i = 1; i <= vertex_num; i++)
		{
			if (g_node_char[i].inci_stat && g_node_char[i].mark_stat != 0 && g_node_char[i].scan_stat == 0)
			{
				sel_node = i;
				flag = 1;
				break;
			}
		}

		if (!flag)
		{
			flag = 0;

			for (int i = 1; i <= vertex_num; i++)
			{
				if (g_node_char[i].inci_stat && g_node_char[i].mark_stat == 0)
				{
					sel_node = i;
					curr_mark ++;
					g_node_char[sel_node].mark_stat = curr_mark;
					flag = 1;
					break;
				}
			}

			if (!flag)
				break;
		}

		for (int i = 1; i <= g_node_char[sel_node].neighbors[0]; i++)
		{
			g_node_char[g_node_char[sel_node].neighbors[i]].mark_stat = curr_mark;
		}
		g_node_char[sel_node].scan_stat = 1;
	}

	//for (int i = 1; i <= vertex_num; i++)
	//{
	//	printf("%d, %d, %d\n", i, g_node_char[i].inci_stat, g_node_char[i].mark_stat);
	//}

	if (curr_mark > 1)
	{		
		piece_graph.node_num = curr_mark;
		for (int i = 1; i <= curr_mark; i++)
		{
			for (int j = 1; j <= curr_mark; j++)
			{
				piece_graph.weight_matrix[i][j] = INF;
			}
		}

		for (int i = 1; i < vertex_num; i++)
		{
			if (!g_node_char[i].inci_stat)
				continue;

			for (int j = i+1; j <= vertex_num; j++)
			{
				if (!g_node_char[j].inci_stat)
					continue;

				if (g_node_char[i].mark_stat == g_node_char[j].mark_stat)
					continue;

				int dist = min_cost[i][j]+eta*(orig_deg[i]+orig_deg[j]-4);
				if (dist < piece_graph.weight_matrix[g_node_char[i].mark_stat][g_node_char[j].mark_stat])
				{
					piece_graph.weight_matrix[g_node_char[i].mark_stat][g_node_char[j].mark_stat] = dist;
					rep_node1[g_node_char[i].mark_stat][g_node_char[j].mark_stat] = i;
					rep_node2[g_node_char[i].mark_stat][g_node_char[j].mark_stat] = j;
					piece_graph.weight_matrix[g_node_char[j].mark_stat][g_node_char[i].mark_stat] = dist;
					rep_node1[g_node_char[j].mark_stat][g_node_char[i].mark_stat] = j;
					rep_node2[g_node_char[j].mark_stat][g_node_char[i].mark_stat] = i;
				}
			}
		}

		get_min_span_tree();

		for (int i = 1; i < curr_mark; i++)
		{
			for (int j = i+1; j <= curr_mark; j++)
			{
				if (min_span_tree[i][j])
				{
					ad_matrix[rep_node1[i][j]][rep_node2[i][j]] ++;
					ad_matrix[rep_node2[i][j]][rep_node1[i][j]] ++;
				}
			}
		}
	}

	/* finding odd nodes */
	int deg[MAX_NODE_TAG_LENGTH], odd_nodes[MAX_NODE_TAG_LENGTH];
	memset(deg, 0, sizeof(deg));
	odd_nodes[0] = 0;

	for (int i = 1; i < vertex_num; i++)
	{
		for (int j = i+1; j <= vertex_num; j++)
		{
			deg[i] += ad_matrix[i][j];
			deg[j] += ad_matrix[i][j];
		}
	}

	for (int i = 1; i <= vertex_num; i++)
	{
		if (deg[i]%2)
		{
			odd_nodes[0] ++;
			odd_nodes[odd_nodes[0]] = i;
		}
	}

	/* even the graph */

	deg_graph.node_num = odd_nodes[0];
	for (int i = 1; i <= odd_nodes[0]; i++)
	{
		for (int j = 1; j <= odd_nodes[0]; j++)
		{
			deg_graph.weight_matrix[i][j] = min_cost[odd_nodes[i]][odd_nodes[j]];
		}
	}

	even_graph();

	for (int i = 1; i < odd_nodes[0]; i++)
	{
		for (int j = i+1; j <= odd_nodes[0]; j++)
		{
			if (match_matrix[i][j])
			{
				ad_matrix[odd_nodes[i]][odd_nodes[j]] ++;
				ad_matrix[odd_nodes[j]][odd_nodes[i]] ++;
			}
		}
	}

	/* get Euler route */

	int Euler_route[MAX_NODE_ROUTE_LENGTH];
	get_Euler_route(Euler_route, vertex_num);

	//printf("euler_route\n");
	//print_one_dim_array(Euler_route);

	/* transform to task seq */

	fh_route[0] = 1;
	fh_route[1] = 0;

	for (int i = 1; i < Euler_route[0]; i++)
	{
		if (tad_matrix[Euler_route[i]][Euler_route[i+1]])
		{
			tad_matrix[Euler_route[i]][Euler_route[i+1]] --;
			tad_matrix[Euler_route[i+1]][Euler_route[i]] --;

			fh_route[0] ++;
			fh_route[fh_route[0]] = find_task(Euler_route[i], Euler_route[i+1], inst_tasks);
		}
	}

	fh_route[0] ++;
	fh_route[fh_route[0]] = 0;

	//printf("fh_route\n");
	//print_one_dim_array(fh_route);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void get_min_span_tree()
{
	memset(min_span_tree, 0, sizeof(min_span_tree));

	int node_visited[MAX_NODE_TAG_LENGTH], edge_color[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
	memset(node_visited, 0, sizeof(node_visited));
	memset(edge_color, 0, sizeof(edge_color));

	/* start from node 1 */
	node_visited[1] = 1;

	for (int i = 1; i <= piece_graph.node_num; i++)
	{
		if (piece_graph.weight_matrix[1][i] < INF)
		{
			edge_color[1][i] = 1;
			edge_color[i][1] = 1;
		}
	}

	int min_head, min_tail;
	for (int n = 0; n < piece_graph.node_num; n++)
	{
		int min_weight = INF;
		for (int i = 1; i <= piece_graph.node_num; i++)
		{
			for (int j = 1; j <= piece_graph.node_num; j++)
			{
				if (edge_color[i][j] == 1)
				{
					if (piece_graph.weight_matrix[i][j] < min_weight)
					{
						min_weight = piece_graph.weight_matrix[i][j];
						min_head = i;
						min_tail = j;
					}
				}
			}
		}

		if (!node_visited[min_head])
		{
			node_visited[min_head] = 1;
			for (int i = 1; i <= piece_graph.node_num; i++)
			{
				if (piece_graph.weight_matrix[i][min_head] < INF && !node_visited[i])
				{
					edge_color[i][min_head] = 1;
					edge_color[min_head][i] = 1;
				}
			}
		}
		else
		{
			node_visited[min_tail] = 1;
			for (int i = 1; i <= piece_graph.node_num; i++)
			{
				if (piece_graph.weight_matrix[i][min_tail] < INF && !node_visited[i])
				{
					edge_color[i][min_tail] = 1;
					edge_color[min_tail][i] = 1;
				}
			}
		}

		min_span_tree[min_head][min_tail] = 1;
		min_span_tree[min_tail][min_head] = 1;
		edge_color[min_head][min_tail] = 2;
		edge_color[min_tail][min_head] = 2;

		for (int i = 1; i <= piece_graph.node_num; i++)
		{
			if (node_visited[i])
				continue;

			min_weight = INF;
			for (int j = 1; j <= piece_graph.node_num; j++)
			{
				if (piece_graph.weight_matrix[i][j] < INF && node_visited[j])
				{
					edge_color[i][j] = 0;
					edge_color[j][i] = 0;

					if (piece_graph.weight_matrix[i][j] < min_weight)
					{
						min_weight = piece_graph.weight_matrix[i][j];
						min_head = i;
						min_tail = j;
					}
				}
			}

			edge_color[min_head][min_tail] = 1;
			edge_color[min_tail][min_head] = 1;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void even_graph()
{
	memset(match_matrix, 0, sizeof(match_matrix));

	int left_nodes[MAX_NODE_TAG_LENGTH];
	left_nodes[0] = deg_graph.node_num;
	for (int i = 1; i <= deg_graph.node_num; i++)
	{
		left_nodes[i] = i;
	}

	matching best_matching;

	if (deg_graph.node_num <= 6)
	{
		matching curr_matching;
		curr_matching.matching_size = 0;
		curr_matching.total_weight = 0;

		enum_matching(&best_matching, &curr_matching, left_nodes);

		for (int i = 0; i < best_matching.matching_size; i++)
    	{
	    	match_matrix[best_matching.links[i].head_node][best_matching.links[i].tail_node] ++;
		    match_matrix[best_matching.links[i].tail_node][best_matching.links[i].head_node] ++;
	    }
	}
	else
	{
		best_matching.matching_size = 0;
		best_matching.total_weight = 0;

		for (int n = 0; n < deg_graph.node_num/2; n++)
		{
			edge_char tmp_link, best_link;
			best_link.cost = INF;
			int best_i, best_j;

			for (int i = 1; i < left_nodes[0]; i++)
			{
				for (int j = i+1; j <= left_nodes[0]; j++)
				{
					tmp_link.head_node = left_nodes[i];
					tmp_link.tail_node = left_nodes[j];
					tmp_link.cost = deg_graph.weight_matrix[left_nodes[i]][left_nodes[j]];

					if (tmp_link.cost < best_link.cost)
					{
						best_link = tmp_link;
						best_i = i;
						best_j = j;
					}
				}
			}

			best_matching.links[best_matching.matching_size] = best_link;
			best_matching.matching_size ++;
			best_matching.total_weight += best_link.cost;

			delete_element(left_nodes, best_j);
			delete_element(left_nodes, best_i);
		}

		for (int i = 0; i < best_matching.matching_size; i++)
    	{
	    	match_matrix[best_matching.links[i].head_node][best_matching.links[i].tail_node] ++;
		    match_matrix[best_matching.links[i].tail_node][best_matching.links[i].head_node] ++;
	    }
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void enum_matching(matching *best_matching, matching *curr_matching, int *left_nodes)
{
	if (left_nodes[0] == 0)
	{
		best_matching->matching_size = curr_matching->matching_size;
		best_matching->total_weight = curr_matching->total_weight;
		memcpy(best_matching->links, curr_matching->links, curr_matching->matching_size*sizeof(edge_char));
		//for (int i = 0; i < curr_matching->matching_size; i++)
		//{
		//	best_matching->links[i] = curr_matching->links[i];
		//}
	}
	else
	{
		int best_total_weight = INF;

		for (int i = 1; i < left_nodes[0]; i++)
		{
			for (int j = i+1; j <= left_nodes[0]; j++)
			{
				int node1, node2;
				node1 = left_nodes[i];
				node2 = left_nodes[j];
				curr_matching->links[curr_matching->matching_size].head_node = node1;
				curr_matching->links[curr_matching->matching_size].tail_node = node2;
				curr_matching->matching_size ++;
				curr_matching->total_weight += deg_graph.weight_matrix[node1][node2];
				delete_element(left_nodes, j);
				delete_element(left_nodes, i);

				matching tmp_best_matching;
				enum_matching(&tmp_best_matching, curr_matching, left_nodes);

				if (tmp_best_matching.total_weight < best_total_weight)
				{
					best_matching->matching_size = tmp_best_matching.matching_size;
            		best_matching->total_weight = tmp_best_matching.total_weight;
					memcpy(best_matching->links, tmp_best_matching.links, tmp_best_matching.matching_size*sizeof(edge_char));
            		//for (int i = 0; i < tmp_best_matching->matching_size; i++)
            		//{
            		//	best_matching->links[i] = tmp_best_matching->links[i];
            		//}
				}

				curr_matching->matching_size --;
				curr_matching->total_weight -= deg_graph.weight_matrix[node1][node2];
				add_element(left_nodes, node1, i);
				add_element(left_nodes, node2, j);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void get_Euler_route(int *Euler_route, int node_num)
{
	int dep_pos;

	Euler_route[0] = 0;

	augment_Euler_route(Euler_route, 1, node_num);

	for (dep_pos = 1; dep_pos <= Euler_route[0]; dep_pos++)
	{
		if (Euler_route[dep_pos] == DEPOT)
			break;
	}

	for (int i = 1; i <= dep_pos; i++)
	{
		Euler_route[Euler_route[0]+i] = Euler_route[i];
	}

	for (int i = 1; i <= Euler_route[0]; i++)
	{
		Euler_route[i] = Euler_route[dep_pos+i-1];
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void augment_Euler_route(int *Euler_route, int curr_node, int node_num)
{
	int neighbors[MAX_NODE_TAG_LENGTH];
	neighbors[0] = 0;

	for (int i = 1; i <= node_num; i++)
	{
		if (ad_matrix[curr_node][i])
		{
			neighbors[0] ++;
			neighbors[neighbors[0]] = i;
		}
	}

	if (neighbors[0] == 0)
	{
		Euler_route[0] ++;
		Euler_route[Euler_route[0]] = curr_node;
	}
	else
	{
		while (neighbors[0] > 0)
		{
			int k = neighbors[0];

			ad_matrix[curr_node][neighbors[k]] --;
			ad_matrix[neighbors[k]][curr_node] --;

			augment_Euler_route(Euler_route, neighbors[k], node_num);

			neighbors[0] = 0;
			for (int i = 1; i <= node_num; i++)
			{
				if (ad_matrix[curr_node][i])
				{
					neighbors[0] ++;
					neighbors[neighbors[0]] = i;
				}
			}
		}

		Euler_route[0] ++;
		Euler_route[Euler_route[0]] = curr_node;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void post_opt(individual *indi, const task *inst_tasks)
{
	int positions[MAX_TASK_SEQ_LENGTH];
	find_ele_positions(positions, indi->sequence, 0);
	for (int p = 1; p < positions[0]; p++)
	{
		if (positions[p+1]-positions[p] < 5)
			continue;

		int route[MAX_TASK_ROUTE_LENGTH];
		int fh_route1[MAX_TASK_ROUTE_LENGTH], fh_route2[MAX_TASK_ROUTE_LENGTH], fh_route[MAX_TASK_ROUTE_LENGTH];
		copy_sub_array(route, indi->sequence, positions[p], positions[p+1]);
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
			for (int j = positions[p]; j < positions[p+1]; j++)
			{
				indi->sequence[j] = fh_route[j-positions[p]+1];
			}
			indi->total_cost += cost3-cost4;
		}
	}

	get_route_seg_length(indi->route_seg_length, indi->sequence, inst_tasks);
	indi->max_length = max(indi->route_seg_length);
}
