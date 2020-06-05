
#ifndef _TEST_H
#define _TEST_H

#define INF 2100000000
#define DUMMY_CYCLE 0
#define MAX_TASK_NUM 1000
#define MAX_NODE_NUM 300
#define MAX_ARCS_TAG_LENGTH 1001
#define MAX_TASK_TAG_LENGTH 1001
#define MAX_NODE_TAG_LENGTH 300
#define MAX_ROUTE_TAG_LENGTH 50
#define MAX_SEG_TAG_LENGTH 50
#define MAX_TASK_SEG_LENGTH 550
#define MAX_NODE_SEG_LENGTH 1000
#define MAX_TASK_ROUTE_LENGTH 550
#define MAX_NODE_ROUTE_LENGTH 1000
#define MAX_TASK_SEQ_LENGTH 550
#define MAX_FACS_TAG_LENGTH 5
#define SI 1
#define DI 2
#define SWAP 3

#define MAX_POPSIZE 200
#define MAX_TOTALSIZE 200
#define MAX_NONDOMINATED_NUM 1000

#define M_trial 10
#define POP_SIZE 150
#define M_PROB 0.1
#define M_ite 900
#define M_wite 600

#define NHSIZE 9 // size of neighborhood of the weighted vector

#define MAX_NSIZE 10 // upper bound of nsize

extern int vertex_num;
extern int req_edge_num;
extern int req_arc_num;
extern int nonreq_edge_num;
extern int nonreq_arc_num;
extern int task_num;
extern int total_arc_num;
extern int vehicle_num;
extern int capacity;
extern int lower_bound;

extern int DEPOT;

extern int trav_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
extern int serve_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
extern int shortest_path[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
extern int min_cost[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];

/* randomness */

int* rand_perm(int num);
// generate a random permutation from 1 to num

int rand_choose(int num);
// choose a number randomly between 1 and num (num must be not too large)

/* operations for arrays */

void print_one_dim_array(int *a);

void print_two_dim_matrix(int **a, int row, int col);

void delete_element(int *a, int k);
// delete the kth element from array a

void add_element(int *a, int e, int k);
// add element e to the kth position of array a

void reverse_direction(int *a, int k1, int k2);
// reverse the direction of array a from position k1 to k2

void find_ele_positions(int *positions, int *a, int e);
// find the positions that element e occurs in array a

void copy_sub_array(int *b, int *a, int k1, int k2);
// copy sub-array a[k1:k2] to array b

void insert_array(int *b, int *a, int k);
// insert array a into position k of array b

void replace_sub_array(int *b, int *a, int k1, int k2);
// replace the sub-array b[k1:k2] with array a

void link_array(int *b, int *a);
// link array a behind array b

int equal(int *a, int *b);
// if array a equals array b, return 1, otherwise return 0

int min(int *a);
// get minimal value within array a

int max(int *a);
// get maximal value within array a

int set_bit(int num, int k);
// set kth bit of num to 1

int unset_bit(int num, int k);
// set kth bit of num to 0

int bit_one(int num, int k);
// return 1 if the kth bit of num is 1, otherwise return 0

/* operations for task sequences (a special array) */

typedef struct task
{
	int head_node;
	int tail_node;
	int dead_cost;
	int serv_cost;
	int demand;
	int inverse;
} task;

typedef struct arc
{
	int tail_node;
	int head_node;
	int trav_cost;
} arc;

typedef struct individual
{
	int sequence[MAX_TASK_SEQ_LENGTH];
	int route_seg_load[MAX_SEG_TAG_LENGTH];
	int route_seg_length[MAX_SEG_TAG_LENGTH];
	int total_cost;
	int max_length;
	int total_vio_load;
	double fitness;
	int rank;
	double crowdness;
} individual;

typedef struct chromosome
{
	int sequence[MAX_TASK_SEQ_LENGTH];
	int total_cost;
	int max_length;
	int rank;
	double crowdness;
} chromosome;

extern individual nondominated_indis[MAX_NONDOMINATED_NUM];
extern individual pop[MAX_POPSIZE];
extern individual inter_pop[MAX_TOTALSIZE];
extern int popsize;
extern int nondominated_num;
extern int max_ndtotalcost;
extern int min_ndtotalcost;
extern int max_ndmaxlength;
extern int min_ndmaxlength;

extern int exc_min_f1;
extern int exc_min_f2;
extern int exc_max_f1;
extern int exc_max_f2; // f1 = total_cost, f2 = max_length
extern int min_f1;
extern int max_f1;
extern int min_f2;
extern int max_f2;
extern double norm_f1;
extern double norm_f2;
extern int fvs[5];

int delete_element(individual *indis, int indis_size, int k);
// delete the kth individual from the individual set indis with its size as indis_size

int find_arc(int head_node, int tail_node, const arc *inst_arcs);
// find arc index according to head_node and tail_node

int find_task(int head_node, int tail_node, const task *inst_tasks);
// find task index according to head_node and tail_node

int indi_cmp(individual *indi1, individual *indi2);
// compare indi1 and indi2 with pareto domination comparison

int modify_nondominated_indis(individual nondominated_indis[MAX_NONDOMINATED_NUM], int nondominated_indis_size, individual *indi);
// modify the nodominated individual set with indi, and re-calculate the range of the objectives with the nondominated solutions

void remove_task_seq_delimiters(int *task_seq);
// remove delimiters '0' from task_seq

void split(individual *indi, chromosome *chro, const task *inst_tasks);
// split operator

int fleet_limited_split(int *split_task_seq, int *one_task_seq, const task *inst_tasks);
// split operator subject to vehicle_num

int legal(int *task_seq, const task *inst_tasks);
// if task_seq is legal, return 1, otherwise return 0

int get_task_seq_total_cost(int *task_seq, const task *inst_tasks);
// get total cost of task_seq

int get_total_vio_load(int *route_loads);
// get total_vio_load according to route_loads

void get_route_seg_length(int *route_seg_length, int *task_seq, const task *inst_tasks);
// get length of each route according to task_seq

void get_route_seg_load(int *route_seg_load, int *task_seq, const task *inst_tasks);
// get load of each route according to task_seq

void get_assignment(int *assignment, int *task_seq, const task *inst_tasks);
// get assignment of each task according to task_seq

void indi_copy(individual *target, individual *source);
// copy source to target

/* operations for task sequence <--> ordered list */

//void trans_task_seq_to_ordered_list(int (*ordered_list)[250], int *task_seq, const task *inst_tasks);

//void trans_task_route_to_node_route(int *node_route, int *task_route, const task *inst_tasks);

//int CalcOrderedRouteTotalCost(int *OrderedRoute);

//void GetROccurrence(int *ROccurrence, int *OrderedRoute, int p1, int p2, Task *ARPTask);

/* initialization and preprocessing functions */

void mod_dijkstra();

void path_scanning(individual *ps_indi, const task *inst_tasks, int *serve_mark);

void rand_scanning(individual *rs_indi, const task *inst_tasks, int *serve_mark);

void augment_merge(individual *am_indi, const task *inst_tasks, int *serve_mark);

void ulusoy_heuristic(individual *uh_indi, const task *inst_tasks, int *serve_mark);

/* Frederickson's heuristic */

typedef struct graph
{
	int node_num; // from 1 to node_num
	int weight_matrix[MAX_NODE_TAG_LENGTH][MAX_NODE_TAG_LENGTH];
} graph;

typedef struct edge_char
{
	int head_node;
	int tail_node;
	int cost;
} edge_char;

typedef struct matching
{
	int matching_size;
	edge_char links[MAX_NODE_NUM/2+1];
	int total_weight;
} matching;

void fred_heuristic(int *fh_seg, int *orig_seg, const task *inst_tasks, int eta);

void get_min_span_tree();

void even_graph();

void enum_matching(matching *best_matching, matching *curr_matching, int *left_nodes);

void get_Euler_route(int *Euler_route, int node_num);

void augment_Euler_route(int *Euler_route, int curr_node, int node_num);

void post_opt(individual &indi, const task *inst_tasks);
// applying fred_heuristic one each route of indi

/* search operators */

typedef struct move
{
	int type;
	int task1;
	int task2;
	int orig_seg;
	int targ_seg;
	int orig_pos;
	int targ_pos;
	int orig_length;
	int targ_length;
	int total_cost;
	int add_vio_load;
	int max_length;
	double fitness;
} move;

void rand_selection(individual *p1, individual *p2, individual *pop);

void tour_selection(individual *p1, individual *p2, individual *pop);

void SBX(individual *xed_child, individual *p1, individual *p2, const task *inst_tasks);

void lns_mut(individual *c, individual *p, const task *inst_tasks, double coef);

void lns(individual *indi, double coef, int nsize, const task *inst_tasks);

void single_insertion(move *best_move, individual *curr_indi, double coef, const task *inst_tasks);

void double_insertion(move *best_move, individual *curr_indi, double coef, const task *inst_tasks);

void swap(move *best_move, individual *curr_indi, double coef, const task *inst_tasks);

void global_repair_operator(individual *indi, const task *inst_tasks);

void domination_sort(individual *pop, int tmp_popsize);

void front_crowdness_sort(individual *pop, int rank, int tmp_popsize);

void SPEA2(individual *pop, int tmp_popsize,int * f, int * mMMMm);

void SPEA3(individual *pop, int tmp_popsize,int f);

int rand_one_number(int T, int popsize, int idx);

int move_fit_cmp(move *mv1, move *mv2);

int indi_fit_cmp(individual *indi1, individual *indi2);

int mv_indi_fit_cmp(move *mv, individual *indi);

#endif