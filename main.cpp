#include "functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>


#include "arrayoperations.cpp"
#include "heuristic.cpp"
#include "initialization.cpp"
#include "searchoperators.cpp"

using namespace std;


int vertex_num;
int req_edge_num;
int req_arc_num;
int nonreq_edge_num;
int nonreq_arc_num;
int task_num;
int total_arc_num;
int vehicle_num;
int capacity;
int lower_bound;

int DEPOT;

int trav_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int serve_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int shortest_path[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
int min_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];

individual nondominated_indis[MAX_NONDOMINATED_NUM];
individual pop[MAX_TOTALSIZE];
individual inter_pop[MAX_TOTALSIZE];
individual sub_pop[MAX_TOTALSIZE];
individual total_pop[MAX_TOTALSIZE];
individual best_one;
int popsize = 60;
int nondominated_num;
int max_ndtotalcost;
int min_ndtotalcost;
int max_ndmaxlength;
int min_ndmaxlength;

int exc_min_f1;
int exc_min_f2;
int exc_max_f1;
int exc_max_f2; // f1 = total_cost, f2 = max_length
int min_f1;
int max_f1;
int min_f2;
int max_f2;
double norm_f1;
double norm_f2;
int change_m;
int MMM;
char dummy_string[50];
char curr_file[50];
char filename_string[50];
int get_edge[1005][1005];

// 0 ~ 22  gdb
// 23 ~ 46 eglese
// 47 ~ 80 bccm
const char *input_files[] = {
	"DATA\\gdb\\gdb1.dat",
	"DATA\\gdb\\gdb2.dat",
	"DATA\\gdb\\gdb3.dat",
	"DATA\\gdb\\gdb4.dat",
	"DATA\\gdb\\gdb5.dat",
	"DATA\\gdb\\gdb6.dat",
	"DATA\\gdb\\gdb7.dat",
	"DATA\\gdb\\gdb8.dat",
	"DATA\\gdb\\gdb9.dat",
	"DATA\\gdb\\gdb10.dat",
	"DATA\\gdb\\gdb11.dat",
	"DATA\\gdb\\gdb12.dat",
	"DATA\\gdb\\gdb13.dat",
	"DATA\\gdb\\gdb14.dat",
	"DATA\\gdb\\gdb15.dat",
	"DATA\\gdb\\gdb16.dat",
	"DATA\\gdb\\gdb17.dat",
	"DATA\\gdb\\gdb18.dat",
	"DATA\\gdb\\gdb19.dat",
	"DATA\\gdb\\gdb20.dat",
	"DATA\\gdb\\gdb21.dat",
	"DATA\\gdb\\gdb22.dat",
	"DATA\\gdb\\gdb23.dat",

	"DATA\\eglese\\egl-e1-A.dat",
	"DATA\\eglese\\egl-e1-B.dat",
	"DATA\\eglese\\egl-e1-C.dat",
	"DATA\\eglese\\egl-e2-A.dat",
	"DATA\\eglese\\egl-e2-B.dat",
	"DATA\\eglese\\egl-e2-C.dat",
	"DATA\\eglese\\egl-e3-A.dat",
	"DATA\\eglese\\egl-e3-B.dat",
	"DATA\\eglese\\egl-e3-C.dat",
	"DATA\\eglese\\egl-e4-A.dat",
	"DATA\\eglese\\egl-e4-B.dat",
	"DATA\\eglese\\egl-e4-C.dat",
	"DATA\\eglese\\egl-s1-A.dat",
	"DATA\\eglese\\egl-s1-B.dat",
	"DATA\\eglese\\egl-s1-C.dat",
	"DATA\\eglese\\egl-s2-A.dat",
	"DATA\\eglese\\egl-s2-B.dat",
	"DATA\\eglese\\egl-s2-C.dat",
	"DATA\\eglese\\egl-s3-A.dat",
	"DATA\\eglese\\egl-s3-B.dat",
	"DATA\\eglese\\egl-s3-C.dat",
	"DATA\\eglese\\egl-s4-A.dat",
	"DATA\\eglese\\egl-s4-B.dat",
	"DATA\\eglese\\egl-s4-C.dat",

	"DATA\\bccm\\1A.dat",
	"DATA\\bccm\\1B.dat",
	"DATA\\bccm\\1C.dat",
	"DATA\\bccm\\2A.dat",
	"DATA\\bccm\\2B.dat",
	"DATA\\bccm\\2C.dat",
	"DATA\\bccm\\3A.dat",
	"DATA\\bccm\\3B.dat",
	"DATA\\bccm\\3C.dat",
	"DATA\\bccm\\4A.dat",
	"DATA\\bccm\\4B.dat",
	"DATA\\bccm\\4C.dat",
	"DATA\\bccm\\4D.dat",
	"DATA\\bccm\\5A.dat",
	"DATA\\bccm\\5B.dat",
	"DATA\\bccm\\5C.dat",
	"DATA\\bccm\\5D.dat",
	"DATA\\bccm\\6A.dat",
	"DATA\\bccm\\6B.dat",
	"DATA\\bccm\\6C.dat",
	"DATA\\bccm\\7A.dat",
	"DATA\\bccm\\7B.dat",
	"DATA\\bccm\\7C.dat",
	"DATA\\bccm\\8A.dat",
	"DATA\\bccm\\8B.dat",
	"DATA\\bccm\\8C.dat",
	"DATA\\bccm\\9A.dat",
	"DATA\\bccm\\9B.dat",
	"DATA\\bccm\\9C.dat",
	"DATA\\bccm\\9D.dat",
	"DATA\\bccm\\10A.dat",
	"DATA\\bccm\\10B.dat",
	"DATA\\bccm\\10C.dat",
	"DATA\\bccm\\10D.dat",
};


task inst_tasks[MAX_TASK_TAG_LENGTH];
arc inst_arcs[MAX_ARCS_TAG_LENGTH];

string solution_to_string(individual a) {
    vector <string> prefix; // depot, route_id, served_demand, cost of route, #traveled edges
    vector <string> ans;
    int cost_so_far = 0;
    int served_demand = 0;
    int pre_vertex = DEPOT;
    int route_id = 1;
    int num_of_edges = 1;
    bool flag = false;
    string data;
    int tot_cost = 0;
    for (int i = 2; i <= a.sequence[0]; i++) {
        if (a.sequence[i] == 0) {
            string temp = "(D,0,";
            temp += to_string(DEPOT) + ",";
            temp += to_string(DEPOT) + ",";
            temp += to_string(served_demand) + ",";
            cost_so_far += min_cost[DEPOT][pre_vertex];
            tot_cost += cost_so_far;
            temp += to_string(cost_so_far) + ")";
            ans.push_back(temp);
            pre_vertex = DEPOT;

            string temp_prefix;
            temp_prefix += to_string(DEPOT) + " ";
            temp_prefix += to_string(route_id++) + " ";
            temp_prefix += to_string(served_demand) + " ";
            temp_prefix += to_string(cost_so_far) + " ";
            temp_prefix += to_string(num_of_edges) + " ";
            prefix.push_back(temp_prefix);
            served_demand = 0;
            cost_so_far = 0;
            num_of_edges = 0;

            data += temp_prefix;
            for(auto s:ans) {
                data += s;
            }
            //cout<<data;
            data += "\n";
            ans.clear();
        } else {
            auto S = inst_tasks[a.sequence[i]];
            num_of_edges++;
            while (shortest_path[pre_vertex][S.head_node][0] != 1) {
                num_of_edges++;
                int now_vertex = shortest_path[pre_vertex][S.head_node][2];
                int edge_id = get_edge[pre_vertex][now_vertex];
                string temp = "(P," + to_string(edge_id) + ",";
                temp += to_string(pre_vertex) + ",";
                temp += to_string(now_vertex) + ",";
                temp += to_string(served_demand) + ",";
                cost_so_far += min_cost[now_vertex][pre_vertex];
                temp += to_string(cost_so_far) + ")";
//                ans += temp;
                ans.push_back(temp);
                pre_vertex = now_vertex;
            }
            num_of_edges++;
            int now_vertex = shortest_path[pre_vertex][S.head_node][1];
            int edge_id = get_edge[pre_vertex][now_vertex];
            string temp = "(P," + to_string(edge_id) + ",";
            temp += to_string(pre_vertex) + ",";
            temp += to_string(now_vertex) + ",";
            temp += to_string(served_demand) + ",";
            cost_so_far += min_cost[now_vertex][pre_vertex];
            temp += to_string(cost_so_far) + ")";
            ans.push_back(temp);
            pre_vertex = now_vertex;

            edge_id = get_edge[S.head_node][S.tail_node];
            temp = "(S," + to_string(edge_id) + ",";
            temp += to_string(S.head_node) + ",";
            temp += to_string(S.tail_node) + ",";
            served_demand += S.demand;
            cost_so_far += S.serv_cost;
            temp += to_string(served_demand) + ",";
            temp += to_string(cost_so_far) + ")";
            ans.push_back(temp);
            pre_vertex = S.tail_node;
        }
    }
    //printf("tot cost\n",tot_cost);
    return data;
}

// 0 ~ 22  gdb
// 23 ~ 46 eglese
// 47 ~ 80 bccm
int main(void)

{
	int Run = 1;
	int function_num = 1;
	for (int num_f = 0; num_f < function_num; num_f++)
	{
		change_m = popsize; // 聚类中最大的个数

		FILE *fp;
		fp = fopen(input_files[num_f], "r");
		//printf("%s", input_files[num_f]);
		req_arc_num = 0;
		nonreq_arc_num = 0;

		while (1)
		{
			fscanf(fp, "%s", dummy_string);
			//printf("dummy_string = %s\n", dummy_string);

			if (strcmp(dummy_string, "VERTICES") == 0)
			{
				fscanf(fp, "%s", dummy_string);
				fscanf(fp, "%d", &vertex_num);
			}
			else if (strcmp(dummy_string, "ARISTAS_REQ") == 0)
			{
				fscanf(fp, "%s", dummy_string);
				fscanf(fp, "%d", &req_edge_num);
			}
			else if (strcmp(dummy_string, "ARISTAS_NOREQ") == 0)
			{
				fscanf(fp, "%s", dummy_string);
				fscanf(fp, "%d", &nonreq_edge_num);
			}
			else if (strcmp(dummy_string, "VEHICULOS") == 0)
			{
				fscanf(fp, "%s", dummy_string);
				fscanf(fp, "%d", &vehicle_num);
			}
			else if (strcmp(dummy_string, "CAPACIDAD") == 0)
			{
				fscanf(fp, "%s", dummy_string);
				fscanf(fp, "%d", &capacity);
			}
			else if (strcmp(dummy_string, "LISTA_ARISTAS_REQ") == 0)
			{
				task_num = 2 * req_edge_num + req_arc_num;
				total_arc_num = task_num + 2 * nonreq_edge_num + nonreq_arc_num;

				fscanf(fp, "%s", dummy_string);

				for (int i = 1; i <= req_edge_num; i++)
				{
					fscanf(fp, "%s", dummy_string);
					fscanf(fp, "%d,", &inst_tasks[i].head_node);
					fscanf(fp, "%d)", &inst_tasks[i].tail_node);
					fscanf(fp, "%s", dummy_string);
					fscanf(fp, "%d", &inst_tasks[i].serv_cost);
					fscanf(fp, "%s", dummy_string);
					fscanf(fp, "%d", &inst_tasks[i].demand);
					inst_tasks[i].dead_cost = inst_tasks[i].serv_cost;
					inst_tasks[i].inverse = i + req_edge_num;

					inst_tasks[i + req_edge_num].head_node = inst_tasks[i].tail_node;
					inst_tasks[i + req_edge_num].tail_node = inst_tasks[i].head_node;
					inst_tasks[i + req_edge_num].dead_cost = inst_tasks[i].dead_cost;
					inst_tasks[i + req_edge_num].serv_cost = inst_tasks[i].serv_cost;
					inst_tasks[i + req_edge_num].demand = inst_tasks[i].demand;
					inst_tasks[i + req_edge_num].inverse = i;

					inst_arcs[i].head_node = inst_tasks[i].head_node;
					inst_arcs[i].tail_node = inst_tasks[i].tail_node;
					inst_arcs[i].trav_cost = inst_tasks[i].dead_cost;
					inst_arcs[i + req_edge_num].head_node = inst_arcs[i].tail_node;
					inst_arcs[i + req_edge_num].tail_node = inst_arcs[i].head_node;
					inst_arcs[i + req_edge_num].trav_cost = inst_arcs[i].trav_cost;
				
					get_edge[inst_tasks[i].head_node][inst_tasks[i].tail_node] = i;
                	get_edge[inst_tasks[i].tail_node][inst_tasks[i].head_node] = i;
				}
			}
			else if (strcmp(dummy_string, "LISTA_ARISTAS_NOREQ") == 0)
			{
				fscanf(fp, "%s", dummy_string);

				for (int i = task_num + 1; i <= task_num + nonreq_edge_num; i++)
				{
					fscanf(fp, "%s", dummy_string);
					fscanf(fp, "%d,", &inst_arcs[i].head_node);
					fscanf(fp, "%d)", &inst_arcs[i].tail_node);
					fscanf(fp, "%s", dummy_string);
					fscanf(fp, "%d", &inst_arcs[i].trav_cost);

					inst_arcs[i + nonreq_edge_num].head_node = inst_arcs[i].tail_node;
					inst_arcs[i + nonreq_edge_num].tail_node = inst_arcs[i].head_node;
					inst_arcs[i + nonreq_edge_num].trav_cost = inst_arcs[i].trav_cost;
				}
			}
			else if (strcmp(dummy_string, "DEPOSITO") == 0)
			{
				fscanf(fp, "%s", dummy_string);
				fscanf(fp, "%d", &DEPOT);
				break;
			}
			else if (strcmp(dummy_string, "NOMBRE") == 0)
			{
				fscanf(fp, "%s", dummy_string);
				fscanf(fp, "%s", filename_string);
				//printf("Problem:%s\n", filename_string);
			}
		}
		fclose(fp);

		inst_tasks[DUMMY_CYCLE].tail_node = DEPOT;
		inst_tasks[DUMMY_CYCLE].head_node = DEPOT;
		inst_tasks[DUMMY_CYCLE].dead_cost = 0;
		inst_tasks[DUMMY_CYCLE].serv_cost = 0;
		inst_tasks[DUMMY_CYCLE].demand = 0;
		inst_tasks[DUMMY_CYCLE].inverse = DUMMY_CYCLE;
		inst_arcs[DUMMY_CYCLE].tail_node = 1;
		inst_arcs[DUMMY_CYCLE].head_node = 1;
		inst_arcs[DUMMY_CYCLE].trav_cost = 0;

		for (int i = 1; i <= vertex_num; i++)
		{
			for (int j = 1; j <= vertex_num; j++)
			{
				trav_cost[i][j] = INF;
				serve_cost[i][j] = 0;
			}
		}

		// total_arc_num：两倍的edges + 两倍的arcs
		for (int i = 1; i <= total_arc_num; i++)
		{ //total_arc_num = 196
			trav_cost[inst_arcs[i].tail_node][inst_arcs[i].head_node] = inst_arcs[i].trav_cost;
		}

		// task_num：两倍的edges
		for (int i = 1; i <= task_num; i++)
		{ // task_num = 102
			serve_cost[inst_tasks[i].tail_node][inst_tasks[i].head_node] = inst_tasks[i].serv_cost;
		}

		mod_dijkstra();

		/* get output_file */
		//printf("id:%d  %s\n ", num_f, filename_string);
		for (int run = 0; run < Run; run++)
		{
			char wei[20] = ".dat";
			char str[25];
			char path[20] = "res\\";
			char path2[20] = "opt\\";
			char pro_id[20];
			itoa(run, str, 10);
			itoa(num_f, pro_id, 10);
		    //char output_file[40] = "NSGA_ALL.dat";
			char output_file[40] = "SPEA_ALL.dat";
			// nsga
			char nsga[100] = "_nsga_nondominated";
			//char output_nsga_all[100] = "nsga_all.dat";
			char *output_nsga = (char *)malloc(strlen(nsga) + strlen(path) + strlen(pro_id) + strlen(str) + strlen(wei));
			strcpy(output_nsga, path);
			strcat(output_nsga, pro_id);
			strcat(output_nsga, nsga);
			strcat(output_nsga, str);
			strcat(output_nsga, wei);
			//printf("---%s---\n", output_nsga);

			char nsga_opt1[100] = "_nsga_opt";
			char *opt_nsga = (char *)malloc(strlen(nsga_opt1) + strlen(path2) + strlen(pro_id) + strlen(str) + strlen(wei));
			strcpy(opt_nsga, path2);
			strcat(opt_nsga, pro_id);
			strcat(opt_nsga, nsga_opt1);
			strcat(opt_nsga, str);
			strcat(opt_nsga, wei);
			//printf("---%s---\n", opt_nsga);

			// spea
			char spea[100] = "_spea_nondominated";
			//char output_spea_all[100] = "spea_all.dat";
			char *output_spea = (char *)malloc(strlen(spea) + strlen(path) + strlen(pro_id) + strlen(str) + strlen(wei));
			strcpy(output_spea, path);
			strcat(output_spea, pro_id);
			strcat(output_spea, spea);
			strcat(output_spea, str);
			strcat(output_spea, wei);
			//printf("---%s---\n", output_spea);

			char spea_opt1[100] = "_spea_opt";
			char *opt_spea = (char *)malloc(strlen(spea_opt1) + strlen(path2) + strlen(pro_id) + strlen(str) + strlen(wei));
			strcpy(opt_spea, path2);
			strcat(opt_spea, pro_id);
			strcat(opt_spea, spea_opt1);
			strcat(opt_spea, str);
			strcat(opt_spea, wei);
			//printf("---%s---\n", opt_spea);

			// mix
			char mix[100] = "_mix_nondominated";
			//char output_mix_all[100] = "mix_all.dat";
			char *output_mix = (char *)malloc(strlen(mix) + strlen(path) + strlen(pro_id) + strlen(mix) + strlen(wei));
			strcpy(output_mix, path);
			strcat(output_mix, pro_id);
			strcat(output_mix, mix);
			strcat(output_mix, str);
			strcat(output_mix, wei);
			//printf("---%s---\n", output_mix);

			char mix_opt1[100] = "_mix_opt";
			char *opt_mix = (char *)malloc(strlen(mix_opt1) + strlen(path2) + strlen(pro_id) + strlen(str) + strlen(wei));
			strcpy(opt_mix, path2);
			strcat(opt_mix, pro_id);
			strcat(opt_mix, mix_opt1);
			strcat(opt_mix, str);
			strcat(opt_mix, wei);
			//printf("---%s---\n", opt_mix);

			/* begin to run */
			//int flag = 1; // 1:nsga     2:spea2     3:spea2 + nsga
			printf("ID: %d.  Problem:%s\nRun: %d \n",num_f, filename_string,run);
			printf("Generation: ");
			for (int flag = 1; flag < 4; flag++)
			{
				
				
				MMM = NHSIZE;
				if (flag == 1 || flag == 2)
					continue;

				/* random seed */
				int tm;
				tm = time(NULL);
				srand(tm);

				if (flag == 1)
				{
					fp = fopen(opt_nsga, "w");
					fprintf(fp, "\n");
					fclose(fp);
					fp = fopen(output_file, "w");
					fprintf(fp, "\n");
					fclose(fp);
					
				}
				else if (flag == 2)
				{
					fp = fopen(opt_spea, "w");
					fprintf(fp, "\n");
					fclose(fp);
				}
				else if (flag == 3)
				{
					fp = fopen(opt_mix, "w");
					fprintf(fp, "\n");
					fclose(fp);
				}

				clock_t start = clock();

				//printf("tm = %d\n", tm);

				// fp = fopen(output_file, "w");
				// fprintf(fp, "The random seed is %d.\n", tm);
				// fprintf(fp, "There are %d vertices, %d tasks, and the capacities of vehicles is %d.\n\n", vertex_num, req_edge_num + req_arc_num, capacity);
				// fclose(fp);

				/* got LBs of makespan */
				double lamda[MAX_TOTALSIZE];
				for (int i = 0; i < popsize; i++)
				{
					lamda[i] = 1.0 * i / (popsize - 1);
					//lamda[i] = 1.0*(i+1)/(popsize+1);
				}
				int neighborhoods[MAX_TOTALSIZE][MAX_TOTALSIZE];
				for (int i = 0; i < popsize; i++)
				{
					neighborhoods[i][0] = 0;
				}

				// neighborhoods：储存lamda的的邻居编号
				for (int i = 0; i < popsize; i++)
				{
					int start = i - NHSIZE / 2;
					int finish = i + NHSIZE / 2;

					if (start < 0)
					{
						start = 0;
						finish = start + NHSIZE - 1;
					}
					if (finish > popsize - 1)
					{
						finish = popsize - 1;
						start = finish - NHSIZE + 1;
					}

					for (int j = start; j <= finish; j++)
					{
						neighborhoods[i][0]++;
						neighborhoods[i][neighborhoods[i][0]] = j;
					}
				}

				/* initialization */

				individual elite_indi[3], tmp_indi;
				int serve_mark[MAX_TASK_TAG_LENGTH];
				nondominated_num = 0;

				int used;

				// path_scanning
				memset(serve_mark, 0, sizeof(serve_mark));
				serve_mark[0] = task_num;
				for (int i = 1; i <= task_num; i++)
				{
					serve_mark[i] = 1;
				}
				path_scanning(&elite_indi[0], inst_tasks, serve_mark);
				elite_indi[0].total_vio_load = 0;

				// augment_merge
				memset(serve_mark, 0, sizeof(serve_mark));
				serve_mark[0] = task_num;
				for (int i = 1; i <= task_num; i++)
				{
					serve_mark[i] = 1;
				}
				augment_merge(&elite_indi[1], inst_tasks, serve_mark);
				elite_indi[1].total_vio_load = 0;

				// ulusoy_heuristic
				memset(serve_mark, 0, sizeof(serve_mark));
				serve_mark[0] = task_num;
				for (int i = 1; i <= task_num; i++)
				{
					serve_mark[i] = 1;
				}
				ulusoy_heuristic(&elite_indi[2], inst_tasks, serve_mark);
				elite_indi[2].total_vio_load = 0;

				// Algorithm 1: The assignment of representatives
				for (int i = 0; i < 2; i++)
				{
					for (int j = i + 1; j < 3; j++)
					{
						if (elite_indi[j].max_length < elite_indi[i].max_length ||
							(elite_indi[j].max_length == elite_indi[i].max_length && elite_indi[j].total_cost > elite_indi[i].total_cost))
						{
							indi_copy(&tmp_indi, &elite_indi[i]);
							indi_copy(&elite_indi[i], &elite_indi[j]);
							indi_copy(&elite_indi[j], &tmp_indi);
						}
					}
				}

				int tmp_popsize = 0;
				individual init_indi;
				while (tmp_popsize < 3)
				{
					indi_copy(&init_indi, &elite_indi[tmp_popsize]);
					tmp_popsize++;
					indi_copy(&pop[tmp_popsize - 1], &init_indi);
					nondominated_num = modify_nondominated_indis(nondominated_indis, nondominated_num, &init_indi);
				}
				while (tmp_popsize < popsize)
				{
					int trial = 0;
					while (trial < M_trial)
					{
						memset(serve_mark, 0, sizeof(serve_mark));
						serve_mark[0] = task_num;
						for (int i = 1; i <= task_num; i++)
						{
							serve_mark[i] = 1;
						}

						rand_scanning(&init_indi, inst_tasks, serve_mark);

						used = 0;
						for (int i = 0; i < tmp_popsize; i++)
						{
							if (init_indi.max_length == pop[i].max_length && init_indi.total_cost == pop[i].total_cost)
							{
								used = 1; // 由于pop中的前三个都是有之前的三个函数初始化，如果有重复了需要重新生成
								break;
							}
						}

						if (!used)
							break;
					}

					if (trial == M_trial && used)
					{
						popsize = tmp_popsize;
						break;
					}

					tmp_popsize++;
					indi_copy(&pop[tmp_popsize - 1], &init_indi);
					nondominated_num = modify_nondominated_indis(nondominated_indis, nondominated_num, &init_indi);
				} // 至此已经完成初始化

				// for (int i = 0; i < popsize;i++)
				// {
				// 	printf(" the number of ver is %d\n",pop[i].sequence[0]);
				// }
				//printf("nondominated\n");
				// for (int i = 0; i < nondominated_num; i++)
				// {
				// 	//printf("(%d, %d, %d)\n", nondominated_indis[i].max_length, nondominated_indis[i].total_cost, nondominated_indis[i].total_vio_load);
				// }

				// 更新 以下四个变量
				min_f1 = INF;
				max_f1 = 0;
				min_f2 = INF;
				max_f2 = 0;

				for (int i = 0; i < popsize; i++)
				{
					//if (pop[i].total_vio_load > 0)
					//	continue;

					if (pop[i].total_vio_load == 0 && pop[i].total_cost < min_f1)
						min_f1 = pop[i].total_cost;
					if (pop[i].total_vio_load == 0 && pop[i].total_cost > max_f1)
						max_f1 = pop[i].total_cost;
					if (pop[i].total_vio_load == 0 && pop[i].max_length < min_f2)
						min_f2 = pop[i].max_length;
					if (pop[i].total_vio_load == 0 && pop[i].max_length > max_f2)
						max_f2 = pop[i].max_length;
				}
				// end 更新

				/* searching phase */

				int ite = 0, wite = 0;
				individual parent1, parent2, xed_child, mted_child, child;

				//int offsize = popsize; // the number of the generated offsprings
				//int totalsize = popsize+offsize; // the number of the sum of the current population and the offsprings

				while (ite < M_ite)
				{
					if (MMM > NHSIZE)
					{
						if (MMM % 2 == 0)
						{
							if (MMM + 1 < popsize)
							{
								MMM++;
							}
							else
							{
								MMM--;
							}
						}
						if (MMM + 5 < popsize)
						{
							MMM += 5;
						}
					}
					else
					{
						MMM = NHSIZE;
					}
					
					//printf("min_f1 = %d, min_f2 = %d, max_f1 = %d, max_f2 = %d\n", min_f1, min_f2, max_f1, max_f2);
					ite++;
					wite++;
					
					// sorting of pop
					// Algorithm 1: The assignment of representatives
					for (int i = 0; i < popsize - 1; i++)
					{
						for (int j = i + 1; j < popsize; j++)
						{
							if (pop[j].max_length < pop[i].max_length ||
								(pop[j].max_length == pop[i].max_length && pop[j].total_cost > pop[i].total_cost))
							{
								indi_copy(&tmp_indi, &pop[i]);
								indi_copy(&pop[i], &pop[j]);
								indi_copy(&pop[j], &tmp_indi);
							}
						}
					}
					// fp = fopen(output_file, "a");
					// for (int ii = 0; ii < popsize - 1; ii++)
					// {
					// 	fprintf(fp, "%d %d %d\n", pop[ii].total_cost, pop[ii].max_length,pop[ii].total_vio_load);
					// }
					// //fprintf(fp, "\n");
					// fclose(fp);
					// for (int i = 0; i < popsize; i++)
					// 	printf("max_length: %d, total_cost: %d\n", pop[i].max_length, pop[i].total_cost);

					/* generate offsprings */
					for (int i = 0; i < popsize; i++) // for each individual, i.e., each subproblem
					{
						child.total_cost = 0; // child does not exist

						int k1, k2, pid1, pid2;
						int candi[MAX_TOTALSIZE + 1];
						memcpy(candi, neighborhoods[i], sizeof(neighborhoods[i]));

						///////////////////////////
						// k1 = rand_choose(candi[0]); // 随机在neighbor里随机两个个体编号 k1 k2
						// while (1)
						// {
						// 	k2 = rand_choose(candi[0]);
						// 	if (k2 != k1)
						// 		break;
						// }

						// pid1 = candi[k1];
						// pid2 = candi[k2];
						/////////////////////////////////////
						k1 = rand_one_number(MMM, popsize, i);

						while (1)
						{
							k2 = rand_one_number(MMM, popsize, i);
							if (k2 != k1)
								break;
						}

						pid1 = k1;
						pid2 = k2;
						///////////////////////////////////
						//printf("k1 = %d, k2 = %d, popsize = %d\n", k1, k2, popsize);
						indi_copy(&parent1, &pop[pid1]);
						indi_copy(&parent2, &pop[pid2]);

						//printf("ite = %d, i = %d\n", ite, i);

						SBX(&xed_child, &parent1, &parent2, inst_tasks); // crossover

						//printf("xed totalcost = %d, maxlength = %d, tvl = %d\n", xed_child.total_cost, xed_child.max_length, xed_child.total_vio_load);

						used = 0;
						for (int j = 1; j <= neighborhoods[i][0]; j++)
						{
							if (xed_child.max_length == pop[neighborhoods[i][j]].max_length && xed_child.total_cost == pop[neighborhoods[i][j]].total_cost &&
								xed_child.total_vio_load == pop[neighborhoods[i][j]].total_vio_load)
							{ // 当前的父代中已经有了这个后代的解，重复了
								used = 1;
								break;
							}
						}

						if (!used)
						{ // 如果不重复，赋给 child
							indi_copy(&child, &xed_child);
						}

						//printf("ite = %d, i = %d\n", ite, i);

						double random = 1.0 * rand() / RAND_MAX;
						//printf("iter: %d mutation, rand(): %f,  Pro: %f \n", ite, random, M_PROB * popsize / change_m);
						if (random < M_PROB * popsize / change_m)
						//if (ite%10 == 0)
						{
							//printf("iter: %d mutation, rand(): %f,  Pro: %f \n", ite, random, M_PROB * popsize / change_m);
							// 变异算子 mted_child为结果，输入为xed_child
							//printf("before totalcost = %d, maxlength = %d, tvl = %d\n", xed_child.total_cost, xed_child.max_length, xed_child.total_vio_load);
							lns_mut(&mted_child, &xed_child, inst_tasks, lamda[i]);
							//printf("after totalcost = %d, maxlength = %d, tvl = %d\n", mted_child.total_cost, mted_child.max_length, mted_child.total_vio_load);

							used = 0;
							for (int j = 1; j <= neighborhoods[i][0]; j++)
							{ // 当前的父代中已经有了这个后代的解，重复了
								if (mted_child.max_length == pop[neighborhoods[i][j]].max_length && mted_child.total_cost == pop[neighborhoods[i][j]].total_cost &&
									mted_child.total_vio_load == pop[neighborhoods[i][j]].total_vio_load)
								{
									used = 1;
									break;
								}
							}

							if (!used)
							{ // 如果不重复，赋给 child
								indi_copy(&child, &mted_child);
							}
						}

						//indi_copy(&child, &mted_child);
						// 更新min_f1 max_f1 min_f2 max_f2
						if (child.total_vio_load == 0 && child.total_cost > 0)
						{ // 判断 indi 是否与nondominated_indis内的前nondominated_indis_size个构成非支配关系，进行删减添加，返回当前的非支配解集
							nondominated_num = modify_nondominated_indis(nondominated_indis, nondominated_num, &child);

							if (child.total_cost < min_f1)
								min_f1 = child.total_cost;
							if (child.total_cost > max_f1)
								max_f1 = child.total_cost;
							if (child.max_length < min_f2)
								min_f2 = child.max_length;
							if (child.max_length > max_f2)
								max_f2 = child.max_length;
						}

						if (child.total_cost == 0)
						{ //说明未产生成功该个体，因此保留该父代
							indi_copy(&inter_pop[i], &pop[i]);
							//printf("in Iter : %d, pop[%d] is not replaced\n", ite, i);
						}
						else
						{
							indi_copy(&inter_pop[i], &child);
						}
					} //至此 inter_pop 为产生的后代的集合

					// sorting of pop
					// 对产生的后代集合排序即，
					//Algorithm 1: The assignment of representatives
					for (int i = 0; i < popsize - 1; i++)
					{
						for (int j = i + 1; j < popsize; j++)
						{
							if (inter_pop[j].max_length < inter_pop[i].max_length ||
								(inter_pop[j].max_length == inter_pop[i].max_length && inter_pop[j].total_cost > inter_pop[i].total_cost))
							{
								indi_copy(&tmp_indi, &inter_pop[i]);
								indi_copy(&inter_pop[i], &inter_pop[j]);
								indi_copy(&inter_pop[j], &tmp_indi);
							}
						}
					}

					//把 父代 + 子代 放到 total_pop 中
					for (int i = 0; i < popsize; i++)
					{
						indi_copy(&total_pop[i], &pop[i]);
						indi_copy(&total_pop[i + popsize], &inter_pop[i]);
					}
					// 对total_pop根据非支配排序 rank(越小越好)，从小到大排序

					if (flag == 1)
					{
						domination_sort(total_pop, 2 * popsize);
						// for (int i = 0; i < popsize * 2; i++)
						// 	printf("No. %d  rank is %d \n", i,total_pop[i].rank);

						// 如果在截断的过程中，边界的个体rank相同，
						// 将边界的rank找到，并根据crowding进行选择
						if (total_pop[popsize - 1].rank == total_pop[popsize].rank)
						{
							front_crowdness_sort(total_pop, total_pop[popsize - 1].rank, 2 * popsize);
						}
					}
					else if (flag == 2)
					{
						domination_sort(total_pop, 2 * popsize);
						int rank1 = 0;
						for (int i = 0; i < tmp_popsize; i++)
							if (total_pop[i].rank == 1)
								rank1 += 1;
						//printf("rank1_sum = %d\n ", rank1);
						//SPEA2(total_pop, 2 * popsize, 0);
						SPEA3(total_pop, 2 * popsize, 0);
						//domination_sort(total_pop, 2 * popsize);
					}
					else if (flag == 3)
					{
						domination_sort(total_pop, 2 * popsize);
						int rank1 = 0;
						for (int i = 0; i < 2 * popsize; i++)
							if (total_pop[i].rank == 1)
								rank1 += 1;

						//printf("rank1_sum = %d, rand: %d \n", rank1, ra);
						if (rank1 > popsize)
						{
							if (rank1 > popsize)
							{
								//printf("run_partial_spea: %d\n",run_spea);
								//SPEA3(total_pop, 2 * popsize, 1);
								change_m = 1;
								SPEA2(total_pop, 2 * popsize, &change_m, &MMM);
								//printf("change_m: %d\n", change_m);
							}
							else
							{
								SPEA3(total_pop, 2 * popsize, 0);
							}
							//domination_sort(total_pop, 2 * popsize);
						}
						else
						{
							//printf("nsga\n");
							//domination_sort(total_pop, 2 * popsize);
							if (total_pop[popsize - 1].rank == total_pop[popsize].rank)
							{
								front_crowdness_sort(total_pop, total_pop[popsize - 1].rank, 2 * popsize);
							}
						}
					}

					// 更新当前种群
					for (int i = 0; i < popsize; i++)
					{
						indi_copy(&pop[i], &total_pop[i]);
					}
					// 在储存的非支配解集中 根据total_cost进行从小到大排序
					for (int i = 0; i < nondominated_num - 1; i++)
					{
						for (int j = i + 1; j < nondominated_num; j++)
						{
							if (nondominated_indis[j].total_cost < nondominated_indis[i].total_cost)
							{
								indi_copy(&tmp_indi, &nondominated_indis[i]);
								indi_copy(&nondominated_indis[i], &nondominated_indis[j]);
								indi_copy(&nondominated_indis[j], &tmp_indi);
							}
						}
					}
					int costt = INF;
					for (int i = 0; i < nondominated_num; i++)
					{
						if (costt > nondominated_indis[i].total_cost)
						{
							costt = nondominated_indis[i].total_cost;
							indi_copy(&best_one, &nondominated_indis[i]);
						}
					}
					
					
					if(ite%60 == 0)
						printf(" %d ",ite);

					

				}

				int costt = INF;
				for (int i = 0; i < nondominated_num; i++)
				{
					if (costt > nondominated_indis[i].total_cost)
					{
						costt = nondominated_indis[i].total_cost;
						indi_copy(&best_one, &nondominated_indis[i]);
					}
				}

				clock_t finish = clock();
				if (flag == 1)
				{ // nsga2
					fp = fopen(output_nsga, "w");
					for (int i = 0; i < nondominated_num; i++)
						fprintf(fp, "%d %d\n", nondominated_indis[i].max_length, nondominated_indis[i].total_cost);

					fclose(fp);
					printf("NSGA2 is done! time_cost: %f \n", (double)(finish - start) / CLOCKS_PER_SEC);
					printf("total cost:  %d \n\n", best_one.total_cost);
				}
				else if (flag == 2)
				{ //spea2
					fp = fopen(output_spea, "w");
					for (int i = 0; i < nondominated_num; i++)
						fprintf(fp, "%d %d\n", nondominated_indis[i].max_length, nondominated_indis[i].total_cost);

					fclose(fp);
					printf("SPEA2 is done! time_cost: %f \n", (double)(finish - start) / CLOCKS_PER_SEC);
					printf("total cost:  %d \n\n", best_one.total_cost);
				}
				else if (flag == 3)
				{ //mix
					fp = fopen(output_mix, "w");
					for (int i = 0; i < nondominated_num; i++)
						fprintf(fp, "%d %d\n", nondominated_indis[i].max_length, nondominated_indis[i].total_cost);

					fclose(fp);
					printf("\nmix is done! time_cost: %f \n", (double)(finish - start) / CLOCKS_PER_SEC);
					printf("Total cost: %d. Max_cost: %d \n\n", best_one.total_cost, best_one.max_length);
				}

				string best_one_route;
				best_one_route =  solution_to_string(best_one);
				//cout << best_one_route << endl;
				if (flag == 1)
					{
						fp = fopen(opt_nsga, "a");
						fprintf(fp,"gen: %d.  max_No.cluster: %d total: %d, max: %d\n ", ite,MMM, best_one.total_cost, best_one.max_length);
						fclose(fp);
					}
					else if (flag == 2)
					{
						fp = fopen(opt_spea, "a");
						fprintf(fp,"gen: %d.  max_No.cluster: %d total: %d, max: %d\n ", ite,MMM, best_one.total_cost, best_one.max_length);
						fclose(fp);
					}
					else if (flag == 3)
					{
						fp = fopen(opt_mix, "a");
						fprintf(fp,"total: %d, max: %d\n", best_one.total_cost, best_one.max_length);
						fprintf(fp, "%s", best_one_route.c_str());
						fclose(fp);
					}

			}

			
		}
	}
}