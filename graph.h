#ifndef GRAPH_H_HHWU
#define GRAPH_H_HHWU

#include <sys/statfs.h>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <algorithm>
#include <set>
#include <queue>
#include <iostream>
#include "Timer.h"
#include "io.h"

using namespace std;
typedef int vertex_id;

const int BLK_SZ = 1024*1024;
const int FILE_NUM = 1024;
const int max_level=1023;
const int SZ_PTR=sizeof(void*);
const int LABEL_INNER_NODE_NUMBER= 1024*1024;
const int LABEL_INNER_EDGE_NUMBER= LABEL_INNER_NODE_NUMBER *6;

struct edgepair {
	vertex_id node1;
	vertex_id node2;
	int weight;
};

const int SZ_VERTEX = sizeof(int);
const int VERTEX_PER_BLK = BLK_SZ/SZ_VERTEX;
const int SZ_OFFSET = sizeof(long);
const int OFFSET_PER_BLK = BLK_SZ/SZ_OFFSET;
const int SZ_EDGE = sizeof(edgepair);
const int EDGE_PER_BLK = BLK_SZ/SZ_EDGE;

struct noderange{
	vertex_id nodeid;
	vertex_id in_degree;
	vertex_id out_degree;
	int begin;
	int r_begin;
};

bool compare(const noderange &a, const noderange &b)
{
	long a_in = a.in_degree;
	long a_out = a.out_degree;
	long b_in = b.in_degree;
	long b_out = b.out_degree;
	
	if((a_in * a_out) == (b_in * b_out))
		return a.nodeid < b.nodeid;
	return (a_in * a_out) < (b_in * b_out);
}

bool compare_edge(const edgepair& a, const edgepair& b)
{
	if(a.node1==b.node1){
		if(a.node2 == b.node2)
			return a.weight < b.weight;
		else
			return a.node2 < b.node2;

	}
	return a.node1 < b.node1;
}

bool compare_r_edge(const edgepair &a, const edgepair&b)
{
	if(a.node2 == b.node2){
		if(a.node1 == b.node1)
			return a.weight < b.weight;
		else
			return a.node1 < b.node1;
	}
	return a.node2 < b.node2;
}


struct minvertex {
	vertex_id nodeid;
	vertex_id in_degree;
	vertex_id out_degree;
	int fileindex;
	friend bool operator < (const minvertex &a, const minvertex &b)
	{
		long a_in = a.in_degree;
		long a_out = a.out_degree;
		long b_in = b.in_degree;
		long b_out = b.out_degree;	
		if((a_in * a_out) == (b_in * b_out)){
			return a.nodeid > b.nodeid;
		}
		return (a_in * a_out > b_in * b_out);  //The smaller one has higher priority
	}
};

struct minedge {
	vertex_id node1;
	vertex_id node2;
	int weight;
	int fileindex;
	friend bool operator < (const minedge &a, const minedge &b)
	{
		if(a.node1 == b.node1)
			return a.node2 > b.node2;
		return a.node1 > b.node1;
	}
};

struct r_minedge{
	vertex_id node1;
	vertex_id node2;
	int weight;
	int fileindex;
	friend bool operator < (const r_minedge &a, const r_minedge &b)
	{
		if(a.node2 == b.node2)
			return a.node1 > b.node1;
		return a.node2 > b.node2;
	}
};

struct weightedge {	
	vertex_id nodeid;
	int weight;
};

bool compare_weightedge(const weightedge& a, const weightedge& b)
{
	return a.nodeid < b.nodeid;
}

struct minvid {
	vertex_id nodeid;
	int levelindex;
};

bool compare_vid(const minvid& a, const minvid& b)
{
	return a.nodeid < b.nodeid;
}

class Graph
{
public:
    Graph(); 
	void set_vertex_num(int v_num); 
	void set_mem(int mem_size);
	void initial(); 
	void write_graph(int &ptr_node_index, vector <noderange> & indexNode, vector <weightedge> & e, vector <weightedge> & r_e, WriteBuffer & write_io); // write the nodes and edges by write_io
	void write_edge(int & ptr_edge_index, vector <edgepair> & edge, WriteBuffer & write_io); // write the edges into file by write_io
	void sortFile(const char* filePath1, const char* filePath2); //sort the graph according to the vertex degree
	void computeIS(const char* filePath1, const char* filePath2, const char* filePath3); //computing independent set, getting the new graph
	
	void tri_loop(int & ptr_node_index, vector <noderange> & indexNode, vector <weightedge> & e, vector <weightedge> & r_e, vector <int> & out_curr, vector <int> & in_curr, vector <int> & new_out_degree, vector <int> & new_in_degree, vector <int> & pos_node_inner, vector <bool> & trin, vector <bool> & r_trin, vector <noderange> & indexNode_inner, vector <weightedge> & e_inner, vector <weightedge> & r_e_inner, int & ptr_node_index_inner, int & ptr_e_inner, int & ptr_r_e_inner); //when the memory for the inner loop is full, applying triangle inequality 
	void inner_loop(int & ptr_node_index, int & ptr_e, int & ptr_r_e, vector <noderange> & indexNode, vector <weightedge> & e, vector <weightedge> & r_e, vector <int> & out_curr, vector <int> & in_curr, vector <int> & new_out_degree, vector <int> & new_in_degree, ReadBuffer & read_inner, WriteBuffer & write_io, long & node_number_total, long & edge_number_total, long & r_edge_number_total); // when the memory for outer loop is full, call the function to load inner loop 
	void triangle_inequality(const char* filePath1, const char* filePath2); // applying triangle inequality to make the shrink the graph size 
	
	void labelInit(int &u, int &u_out_deg, int &u_in_deg, vector <weightedge> & neighbor, WriteBuffer & write_label); // the initial labels when doing independent set
	bool label_loop(int & ptr_node_index, vector <noderange> & indexNode, vector <weightedge> & e, vector <weightedge> & r_e, vector < vector <weightedge> > & e2, vector < vector <weightedge> > & r_e2, vector <int> & curr, vector <int> & r_curr, vector <int> & pos_node_inner, int &ptr_node_index_inner, int &ptr_e_inner, int &ptr_r_e_inner, vector <noderange> & indexNode_inner, vector <weightedge> & e_inner, vector <weightedge> & r_e_inner, long & label_num); // when the memory for the inner loop is full, updating the labels
	int labelUpdate_inner_loop(int & ptr_node_index, int & ptr_e, int & ptr_r_e, vector <noderange> & indexNode, vector <weightedge> & e, vector <weightedge> & r_e, vector < vector <weightedge> > & e2, vector < vector <weightedge> > & r_e2, vector <int> & curr, vector <int> & r_curr, FILE *sortfile, int lev, WriteBuffer &write_level,  WriteBuffer &write_label, long & curr_offset, vector <long> & offset); // when the memory for outer loop is full, call the function to update the labels
	void labelUpdate(int lev, WriteBuffer &write_level, WriteBuffer &write_label, long & curr_offset, vector <long> & offset); // update the labels of level lev
	void topdownLabel(const char* filePath1, const char* filePath2); //doing top down labeling process
	void labelinttobyte(const char* filePath1, const char* filePath2);
	void print_statistics(const char* filePath); //write out the statistics of the remaining graph

	void run(const char* filePath, bool flag_c, bool flag_k, int round_number); 

public:
    long V, E, mem;
    long max_node_m, max_edge_m;
    long INNER_NODE_NUMBER, INNER_EDGE_NUMBER;
    long  max_node_label_m, max_edge_label_m;
    FILE ** level;
    int round;
    long v_num_old, e_num_old;
    long v_num_new, e_num_new;
    long label_number_total;
    double ratio;
    long limit_label_m;
    bool stopis;
    int max_file_num;
};


Graph::Graph()
{
    mem = 4096;
}


void Graph::set_vertex_num(int v_num)
{
    V = v_num;
}


void Graph::set_mem(int mem_size)
{
    mem = mem_size;
}


void Graph::initial()
{
	stopis= false;
    level = (FILE**)malloc((max_level+1)*SZ_PTR);
    
    max_node_m = mem*1024*2;  
	max_edge_m = max_node_m*6;
	
	INNER_NODE_NUMBER = mem*512;
	INNER_EDGE_NUMBER = INNER_NODE_NUMBER*6;
	
	max_node_label_m= mem*256;
	max_edge_label_m= max_node_label_m*12;
	
	limit_label_m = max_edge_label_m*32;
	
	
	struct statfs diskInfo;  
    statfs(".", &diskInfo);  
    unsigned long long blocksize = diskInfo.f_bsize;   
    unsigned long long totalsize = blocksize * diskInfo.f_blocks;
    unsigned long long availableDisk = diskInfo.f_bavail * blocksize;  //available disk size
    unsigned long long num = availableDisk>>30; 
    
    max_file_num=num;
    if(max_file_num >= 256){
    	max_file_num = 256;
    }
    else if(max_file_num < 88){
    	max_file_num = 88;
    }
    
}

// the start point of running
void Graph::run(const char* filePath, bool flag_c, bool flag_k, int round_number)
{
	Timer tt;
	tt.start();

	initial();
		
	char SORT_INPUT_NAME[500], REMAIN_GRAPH_NAME[500], GRAPH_TRIANGLE_NAME[500], LABEL_FILE_NAME[500], FINAL_LABEL_FILE_NAME[500], OFFSET_FILE_NAME[500], GRAPH_GK_info[500];
	strcpy(SORT_INPUT_NAME,filePath);
	strcat(SORT_INPUT_NAME,".sort");
	strcpy(REMAIN_GRAPH_NAME, filePath);
	strcat(REMAIN_GRAPH_NAME,".remain");
	strcpy(GRAPH_TRIANGLE_NAME, filePath);
	strcat(GRAPH_TRIANGLE_NAME,".gk");
	strcpy(LABEL_FILE_NAME, filePath);
	strcat(LABEL_FILE_NAME,".label_temp");
	strcpy(FINAL_LABEL_FILE_NAME, filePath);
	strcat(FINAL_LABEL_FILE_NAME,".label");
	strcpy(OFFSET_FILE_NAME, filePath);
	strcat(OFFSET_FILE_NAME,".offset");
	strcpy(GRAPH_GK_info, filePath);
	strcat(GRAPH_GK_info,".info");
	
	ratio=0.0;
	round = 0;
	cout<<"round: "<<round+1<<endl;
	sortFile(filePath, SORT_INPUT_NAME);
	computeIS(filePath, SORT_INPUT_NAME, REMAIN_GRAPH_NAME);	
	triangle_inequality(REMAIN_GRAPH_NAME, GRAPH_TRIANGLE_NAME);
	
	if(flag_c){ // complete the round
		for (round=1; ; round++){
			cout<<"round: "<<round+1<<endl;
			sortFile(GRAPH_TRIANGLE_NAME, SORT_INPUT_NAME);
			computeIS(GRAPH_TRIANGLE_NAME, SORT_INPUT_NAME, REMAIN_GRAPH_NAME);	
			if(stopis)
				break;
			triangle_inequality(REMAIN_GRAPH_NAME, GRAPH_TRIANGLE_NAME);
		}
	}
	else if(flag_k){ // set the number of levels
		for (round=1; round<round_number; round++){
			cout<<"round: "<<round+1<<endl;
			sortFile(GRAPH_TRIANGLE_NAME, SORT_INPUT_NAME);
			computeIS(GRAPH_TRIANGLE_NAME, SORT_INPUT_NAME, REMAIN_GRAPH_NAME);	
			if(stopis)
				break;
			triangle_inequality(REMAIN_GRAPH_NAME, GRAPH_TRIANGLE_NAME);
		}
	}
	else { //automatically
		for (round=1; ratio<=0.95; round++){
			cout<<"round: "<<round+1<<endl;
			sortFile(GRAPH_TRIANGLE_NAME, SORT_INPUT_NAME);
			computeIS(GRAPH_TRIANGLE_NAME, SORT_INPUT_NAME, REMAIN_GRAPH_NAME);	
			if(stopis)
				break;
			triangle_inequality(REMAIN_GRAPH_NAME, GRAPH_TRIANGLE_NAME);
		}
	}
		
	
	topdownLabel(LABEL_FILE_NAME, OFFSET_FILE_NAME);
	labelinttobyte(LABEL_FILE_NAME, FINAL_LABEL_FILE_NAME);
	print_statistics(GRAPH_GK_info);
	
	remove(SORT_INPUT_NAME);
	remove(REMAIN_GRAPH_NAME);
	remove(LABEL_FILE_NAME);
	
	tt.stop();
	cout<<"The total time: " << tt.GetRuntime()  << " seconds" << endl;
} 

// write the nodes in indexNode and edges in e by write_io
void Graph::write_graph(int &ptr_node_index, vector <noderange> & indexNode, vector <weightedge> & e, vector <weightedge> & r_e, WriteBuffer & write_io)
{
	for(int i = 0; i < ptr_node_index; i++)
	{
		write_io.write(&indexNode[i].nodeid, 1);
		write_io.write(&indexNode[i].out_degree, 1);
		write_io.write(&indexNode[i].in_degree, 1);

		for(int j=0; j<indexNode[i].out_degree; j++){
			write_io.write(&e[j+indexNode[i].begin].nodeid, 1);
			write_io.write(&e[j+indexNode[i].begin].weight, 1);
		}
		for(int j=0; j<indexNode[i].in_degree; j++){
			write_io.write(&r_e[j+indexNode[i].r_begin].nodeid, 1);
			write_io.write(&r_e[j+indexNode[i].r_begin].weight, 1);
		}
	}
	write_io.flush();
	ptr_node_index = 0;
}

// write the edges in edge into file by write_io
void Graph::write_edge(int & ptr_edge_index, vector <edgepair> & edge, WriteBuffer & write_io)
{	
	edgepair previous;
	previous.node1 = previous.node2 = -1;

	for(int i=0; i<ptr_edge_index; i++){
		if(edge[i].node1 != previous.node1 || edge[i].node2 != previous.node2){
			write_io.write(&edge[i].node1, 1);
			write_io.write(&edge[i].node2, 1);
			write_io.write(&edge[i].weight, 1);
			
			previous.node1=edge[i].node1;
			previous.node2=edge[i].node2;
		}
	}
	write_io.flush();
	ptr_edge_index = 0;
}

//sort the graph according to the vertex degree
void Graph::sortFile(const char* filePath1, const char* filePath2)
{
	Timer tt;
	tt.start();

	FILE* file = fopen(filePath1,"rb");
	FILE* sortfile = fopen(filePath2,"wb");
	
	int ptr_file_number = 0;
	FILE **seed = (FILE **)malloc(FILE_NUM*SZ_PTR);
	
	vector <noderange> indexNode; 
	vector <weightedge> e;
	vector <weightedge> r_e;
	
	indexNode.resize(max_node_m+1);
	e.resize(max_edge_m+1);
	r_e.resize(max_edge_m+1);
	

	int ptr_node_index = 0, ptr_e = 0, ptr_r_e = 0;
	long node_number_total =0, edge_number_total = 0, r_edge_number_total = 0;

	ReadBuffer read_io(file);	
	int u, u_out_deg, u_in_deg;
	int v, weight;
	while(!read_io.isend)
	{
		read_io.read(&u);
		read_io.read(&u_out_deg);
		read_io.read(&u_in_deg);
		if(ptr_node_index >= max_node_m || (ptr_e + u_out_deg) >= max_edge_m || (ptr_r_e + u_in_deg) >= max_edge_m) 
		{
			std::sort(indexNode.begin(), indexNode.begin()+ptr_node_index, compare);//sort according to the node degree 
			seed[ptr_file_number] = tmpfile();
			if(seed[ptr_file_number]==NULL){
				cout << "ERROR: cannot create temp file" <<endl;
				exit(1);
			}
			
			WriteBuffer write_io(seed[ptr_file_number]);
			write_graph(ptr_node_index, indexNode, e, r_e, write_io);
			write_io.flush();
			ptr_e = 0;
			ptr_r_e = 0;
			ptr_file_number++;
		}
		indexNode[ptr_node_index].nodeid = u;
		indexNode[ptr_node_index].out_degree = u_out_deg;
		indexNode[ptr_node_index].in_degree = u_in_deg;
		indexNode[ptr_node_index].begin = ptr_e;
		indexNode[ptr_node_index].r_begin = ptr_r_e;
		ptr_node_index++;
		
		node_number_total++;
		edge_number_total += u_out_deg;
		r_edge_number_total += u_in_deg;

		//vertex store in indexNode, edge store in e
		for(int i = 0; i < (u_out_deg+u_in_deg); i++){
			read_io.read(&v);
			read_io.read(&weight);
			
			if(i < u_out_deg)
			{
				e[ptr_e].nodeid = v;
		        e[ptr_e].weight = weight; 
		        ptr_e++;
			}
			else
			{
				r_e[ptr_r_e].nodeid = v;
				r_e[ptr_r_e].weight = weight;
				ptr_r_e++;
			}
		}
	}
	
	if(ptr_node_index > 0)
	{ //for the last part of the file
		std::sort(indexNode.begin(), indexNode.begin()+ptr_node_index, compare);
		seed[ptr_file_number] = tmpfile();
		if(seed[ptr_file_number]==NULL){
			cout << "ERROR: cannot create temp file" <<endl;
			exit(1);
		}
		WriteBuffer write_io(seed[ptr_file_number]);
		write_graph(ptr_node_index, indexNode, e, r_e, write_io);
		write_io.flush();
		ptr_e = 0;
		ptr_r_e = 0;
		ptr_file_number++;
	}

	fclose(file);
	
//	cout<<"edge_number_total:" <<edge_number_total <<endl;
//	cout<<"r_edge_number_total: " <<r_edge_number_total <<endl;
//	cout<<"node_number_total:" << node_number_total << endl;
//	cout<<"ptr_file_number:" <<ptr_file_number <<endl;
	
	v_num_old = node_number_total;
	e_num_old = edge_number_total;
	
	
	//merge Sort
	priority_queue<minvertex> a;
	minvertex new_vertex;
	ReadBuffer seed_io[ptr_file_number+1];
	for(int i = 0; i < ptr_file_number; i++){
		seed_io[i].open(seed[i]);
		seed_io[i].read(&new_vertex.nodeid);
		seed_io[i].read(&new_vertex.out_degree);
		seed_io[i].read(&new_vertex.in_degree);
		new_vertex.fileindex = i;
		a.push(new_vertex);
	}
	
	WriteBuffer write_io_1(sortfile);
	node_number_total = 0;
	edge_number_total = 0;
	r_edge_number_total = 0;
	
	minvertex minone;
	while(!a.empty())
	{
		minone = a.top();
		a.pop();
		
		write_io_1.write(&minone.nodeid, 1);
		write_io_1.write(&minone.out_degree, 1);
		write_io_1.write(&minone.in_degree, 1);
		
		node_number_total++;
		edge_number_total += minone.out_degree;
		r_edge_number_total += minone.in_degree;

		for(int i = 0; i < minone.out_degree+minone.in_degree; i++)
		{
			seed_io[minone.fileindex].read(&v);
			seed_io[minone.fileindex].read(&weight);
			
			write_io_1.write(&v, 1);
			write_io_1.write(&weight, 1);
		}
		
		if(!seed_io[minone.fileindex].isend)
		{
			seed_io[minone.fileindex].read(&new_vertex.nodeid);
			seed_io[minone.fileindex].read(&new_vertex.out_degree);
			seed_io[minone.fileindex].read(&new_vertex.in_degree);
			
			new_vertex.fileindex = minone.fileindex;
			a.push(new_vertex);
 		}
	}
	write_io_1.flush();
	
	for(int i = 0; i < ptr_file_number; i++){
		fclose(seed[i]);
	}
	
	fclose(sortfile);
	tt.stop();
	cout<<"The time used for sorting:" << tt.GetRuntime() << " seconds" << endl;
}

//computing independent set, getting the new graph
void Graph::computeIS(const char* filePath1, const char* filePath2, const char* filePath3)
{
	Timer tt;
	tt.start();
	
//	cout<<"Start computing independent set:" <<endl;
	FILE* file = fopen(filePath1,"rb");
	FILE* sortfile = fopen(filePath2,"rb");
	FILE* remainfile = fopen(filePath3,"wb");
	
	vector <bool> mark, deleteable; 
	vector <edgepair> edge;
	vector <edgepair> r_edge;
	vector <weightedge> neighbor;
	
	mark.resize(V+1, false);
	deleteable.resize(V+1, true);
	neighbor.resize(2*V+1);
	
	ReadBuffer read_io(sortfile);
	int ptr_edge_index=0, ptr_r_edge_index = 0, mark_node_count=0;
	
	FILE ** seed = (FILE**)malloc(FILE_NUM*SZ_PTR);
	FILE ** r_seed = (FILE**)malloc(FILE_NUM*SZ_PTR);
	int ptr_file_number = 0;
	
	int u, u_out_deg, u_in_deg;
	int v, weight;
	edgepair new_edgepair;
	while(!read_io.isend) //according to the sortfile, get the independent set
	{
		read_io.read(&u);
		read_io.read(&u_out_deg);
		read_io.read(&u_in_deg);
		
		for(int i = 0; i < u_out_deg+u_in_deg; i++){
			read_io.read(&neighbor[i].nodeid);
			read_io.read(&neighbor[i].weight);
		}

		if(deleteable[u] && !mark[u] )
		{
			if(u_in_deg*u_out_deg + ptr_edge_index >= max_edge_m)
			{
				sort(edge.begin(), edge.begin()+ptr_edge_index, compare_edge);
				seed[ptr_file_number] = tmpfile();
				WriteBuffer write_io(seed[ptr_file_number]);
				write_edge(ptr_edge_index, edge, write_io);
				edge.clear();
				
				sort(r_edge.begin(), r_edge.begin()+ptr_r_edge_index, compare_r_edge);
				r_seed[ptr_file_number] = tmpfile();
				WriteBuffer write_io_r(r_seed[ptr_file_number]);
				write_edge(ptr_r_edge_index, r_edge, write_io_r);
				r_edge.clear();
				
				ptr_file_number ++;
				
				if(ptr_file_number>=max_file_num){
					fprintf(stderr,"ERROR: Too Many Temp File, out of disk size\n");
					exit(1);
				}
			}
			
			mark[u]=true;
			mark_node_count++;//delete node

			for(int i=0; i<u_out_deg+u_in_deg; i++)
				deleteable[neighbor[i].nodeid]=false;

			for(int j=0; j<u_out_deg; j++){
				for(int k=u_out_deg; k<u_out_deg+u_in_deg; k++){
					if(neighbor[k].nodeid!=neighbor[j].nodeid){
						new_edgepair.node1=neighbor[k].nodeid; // in node
						new_edgepair.node2=neighbor[j].nodeid; // out node
						new_edgepair.weight=neighbor[j].weight+neighbor[k].weight;
						edge.push_back(new_edgepair);
						r_edge.push_back(new_edgepair);
						
						ptr_edge_index++;
						ptr_r_edge_index++;
					}
				}
			}
		}
	}
	
	if(ptr_edge_index>0){
		std::sort(edge.begin(), edge.begin()+ptr_edge_index, compare_edge);
		seed[ptr_file_number] = tmpfile();
		WriteBuffer write_io(seed[ptr_file_number]);
		write_edge(ptr_edge_index, edge, write_io);
		edge.clear();
		
		sort(r_edge.begin(), r_edge.begin()+ptr_r_edge_index, compare_r_edge);
		r_seed[ptr_file_number] = tmpfile();
		WriteBuffer write_io_r(r_seed[ptr_file_number]);
		write_edge(ptr_r_edge_index, r_edge, write_io_r);
		r_edge.clear();
		ptr_file_number ++;
		
		if(ptr_file_number>=max_file_num){
			fprintf(stderr,"ERROR: Too Many Temp File, out of disk size\n");
			exit(1);
		}
	}
	vector<bool>().swap(deleteable);
	fclose(sortfile);
	
//	cout<<"ptr_file_number:" <<ptr_file_number <<endl;
//	cout<<"Getting augmenting edges:" <<endl;
	// merge sort of the augmenting edges
	priority_queue<minedge> a;
	minedge new_edge;
	ReadBuffer seed_io[ptr_file_number];	
	FILE * tempfile = tmpfile();
	WriteBuffer write_io(tempfile);
	
	for(int i=0;i<ptr_file_number;i++)
	{
		seed_io[i].open(seed[i]);
		seed_io[i].read(&new_edge.node1);
		seed_io[i].read(&new_edge.node2);
		seed_io[i].read(&new_edge.weight);
		new_edge.fileindex = i;
		a.push(new_edge);
	}
	
	int min_file_index=0;
	minedge minone_edge;
	while (!a.empty()) {	
		minone_edge = a.top();
		a.pop();
		
		write_io.write(&minone_edge.node1, 1);
		write_io.write(&minone_edge.node2, 1);
		write_io.write(&minone_edge.weight, 1);
		
		min_file_index = minone_edge.fileindex;

		if(!seed_io[minone_edge.fileindex].isend)
		{
			seed_io[minone_edge.fileindex].read(&new_edge.node1);
			seed_io[minone_edge.fileindex].read(&new_edge.node2);
			seed_io[minone_edge.fileindex].read(&new_edge.weight);
			
			new_edge.fileindex = minone_edge.fileindex;
			a.push(new_edge);
 		}
	}
	write_io.flush();
	vector<edgepair>().swap(edge);
	
	for(int i = 0; i < ptr_file_number; i++){
		fclose(seed[i]);
	}
	
	// merge sort reverse
	priority_queue<r_minedge> r_a;
	r_minedge r_new_edge;
	ReadBuffer r_seed_io[ptr_file_number];
	
	FILE * r_tempfile = tmpfile();
	WriteBuffer r_write_io(r_tempfile);	
	for(int i=0;i<ptr_file_number;i++)
	{
		r_seed_io[i].open(r_seed[i]);
		r_seed_io[i].read(&r_new_edge.node1);
		r_seed_io[i].read(&r_new_edge.node2);
		r_seed_io[i].read(&r_new_edge.weight);
		r_new_edge.fileindex = i;
		r_a.push(r_new_edge);
	}
	
	int r_min_file_index=0;
	r_minedge r_minone_edge;
	while (!r_a.empty()) {	
		r_minone_edge = r_a.top();
		r_a.pop();
		
		r_write_io.write(&r_minone_edge.node1, 1);
		r_write_io.write(&r_minone_edge.node2, 1);
		r_write_io.write(&r_minone_edge.weight, 1);
		
		r_min_file_index = r_minone_edge.fileindex;

		if(!r_seed_io[r_minone_edge.fileindex].isend)
		{
			r_seed_io[r_minone_edge.fileindex].read(&r_new_edge.node1);
			r_seed_io[r_minone_edge.fileindex].read(&r_new_edge.node2);
			r_seed_io[r_minone_edge.fileindex].read(&r_new_edge.weight);
			
			r_new_edge.fileindex = r_minone_edge.fileindex;
			r_a.push(r_new_edge);
 		}
	}
	r_write_io.flush();
	vector<edgepair>().swap(r_edge);
	
	for(int i = 0; i < ptr_file_number; i++){
		fclose(r_seed[i]);
	}
	
//	cout<<"Generating the new graph:" <<endl;
	//scan original graph to get new graph		
	rewind(tempfile);
	rewind(r_tempfile);
	WriteBuffer write_io_1(remainfile);
	
	vector <int> out_neighborlist;
	vector <int> in_neighborlist;	
	out_neighborlist.resize(V+1, -1);
	in_neighborlist.resize(V+1, -1);
	
	vector <int> pos;
	pos.resize(V+1, -1);
	vector <int> r_pos;
	r_pos.resize(V+1, -1);
	
	int ptr_node_index=0, ptr_e=0, ptr_r_e = 0;
	edgepair * edgebuff = (edgepair *)malloc(BLK_SZ);
	edgepair * r_edgebuff = (edgepair *)malloc(BLK_SZ);
	int edge_read = fread(edgebuff,SZ_EDGE,EDGE_PER_BLK,tempfile);
	int r_edge_read = fread(r_edgebuff,SZ_EDGE,EDGE_PER_BLK,r_tempfile);

	int ptr_edge_buff=0, ptr_r_edge_buff = 0, remarknumber=0;
	long node_number_total=0, edge_number_total=0, r_edge_number_total  = 0;

	ReadBuffer read_io_1(file);

	level[round] = tmpfile(); //tmpfile()to store different level
	WriteBuffer write_label(level[round]);

	vector <noderange> indexNode;
	indexNode.resize(max_node_m+1);
	
	vector <weightedge> e;
	e.resize(max_edge_m+1);
	
	vector <weightedge> r_e;
	r_e.resize(max_edge_m+1);

	while(!read_io_1.isend)
	{
		read_io_1.read(&u);
		read_io_1.read(&u_out_deg);
		read_io_1.read(&u_in_deg);

		if(!mark[u]) //u is not in independent set
		{
			if(ptr_node_index >= max_node_m || (ptr_e+16*u_out_deg)>= max_edge_m || (ptr_r_e+16*u_in_deg)>=max_edge_m ) //need to flush out the edges
			{
				node_number_total+=ptr_node_index;
				for(int i=0; i<ptr_node_index; i++){
					sort(e.begin()+indexNode[i].begin, e.begin()+indexNode[i].begin+indexNode[i].out_degree, compare_weightedge);
					edge_number_total += indexNode[i].out_degree;
				}
				for(int i=0; i<ptr_node_index; i++){
					sort(r_e.begin()+indexNode[i].r_begin, r_e.begin()+indexNode[i].r_begin+indexNode[i].in_degree, compare_weightedge);
					r_edge_number_total += indexNode[i].in_degree;
				}
				
				write_graph(ptr_node_index, indexNode, e, r_e, write_io_1);
				write_io_1.flush();
				ptr_e = 0;
				ptr_r_e = 0;
			}

			indexNode[ptr_node_index].nodeid = u;
			indexNode[ptr_node_index].out_degree = u_out_deg;
			indexNode[ptr_node_index].in_degree = u_in_deg;
			indexNode[ptr_node_index].begin = ptr_e;
			indexNode[ptr_node_index].r_begin = ptr_r_e;
			ptr_node_index++;

			for(int i = 0; i < u_out_deg; i++)
			{
				read_io_1.read(&v);
				read_io_1.read(&weight);
				if(!mark[v])
				{
					in_neighborlist[v] = u; //u -> v
			        pos[v] = ptr_e;
                    e[ptr_e].nodeid = v;
                    e[ptr_e].weight = weight;
                    ptr_e++;
				}
				else 
					indexNode[ptr_node_index-1].out_degree--;
			}
			
			for(int i = 0; i < u_in_deg; i++)
			{
				read_io_1.read(&v);
				read_io_1.read(&weight);
				if(!mark[v])
				{
					out_neighborlist[v] = u; //pay attention: v -> u
					r_pos[v] = ptr_r_e;
					r_e[ptr_r_e].nodeid = v;
					r_e[ptr_r_e].weight = weight;
					ptr_r_e++;
				}
				else
					indexNode[ptr_node_index-1].in_degree--;
			}

		}
		else
		{
			for(int i=0;i<u_in_deg+u_out_deg;i++){
				read_io_1.read(&neighbor[i].nodeid);
				read_io_1.read(&neighbor[i].weight);
			}
			remarknumber++;
			labelInit(u, u_out_deg, u_in_deg, neighbor, write_label); //the inital labels in level L[round]
		}

		if(!mark[u]){
			while (!feof(tempfile) || ptr_edge_buff!=edge_read) {
				if(edgebuff[ptr_edge_buff].node1>u) 
					break;
				if(edgebuff[ptr_edge_buff].node1 == u){
					if(in_neighborlist[edgebuff[ptr_edge_buff].node2] != u)
					{
						in_neighborlist[edgebuff[ptr_edge_buff].node2] = u;
						indexNode[ptr_node_index-1].out_degree++;
						e[ptr_e].nodeid=edgebuff[ptr_edge_buff].node2;
						e[ptr_e].weight=edgebuff[ptr_edge_buff].weight;
						pos[edgebuff[ptr_edge_buff].node2]=ptr_e;
						ptr_e++;
						
						if(ptr_e>=max_edge_m){  //need to flush out the edges
							ptr_node_index--;
							int last_ptr_node = ptr_node_index; 
							
							node_number_total+=ptr_node_index;
							for(int i=0; i<ptr_node_index; i++){
								sort(e.begin()+indexNode[i].begin, e.begin()+indexNode[i].begin+indexNode[i].out_degree, compare_weightedge);
								edge_number_total += indexNode[i].out_degree;
							}
							for(int i=0; i<ptr_node_index; i++){
								sort(r_e.begin()+indexNode[i].r_begin, r_e.begin()+indexNode[i].r_begin+indexNode[i].in_degree, compare_weightedge);
								r_edge_number_total += indexNode[i].in_degree;
							}
				
							write_graph(ptr_node_index, indexNode, e, r_e, write_io_1);
							write_io_1.flush();
							ptr_e = 0;
							ptr_r_e = 0;
							
							if(indexNode[last_ptr_node].out_degree >= max_edge_m-1 || indexNode[last_ptr_node].in_degree >= max_edge_m-1){
								cout<<"ERROR: Too High Degree Vertex"<<endl;
								exit(1);
							}
							
							for(int i=indexNode[last_ptr_node].begin; i<indexNode[last_ptr_node].begin+indexNode[last_ptr_node].out_degree; i++){
								e[ptr_e].nodeid=e[i].nodeid;
								e[ptr_e].weight=e[i].weight;
								pos[e[i].nodeid]=ptr_e;
								ptr_e++;
							}
							
							for(int i=indexNode[last_ptr_node].r_begin; i<indexNode[last_ptr_node].r_begin+indexNode[last_ptr_node].in_degree; i++){
								r_e[ptr_r_e].nodeid=r_e[i].nodeid;
								r_e[ptr_r_e].weight=r_e[i].weight;
								r_pos[r_e[i].nodeid]=ptr_r_e;
								ptr_r_e++;
							}
							
							ptr_node_index=1;
						}
					}
					else if(edgebuff[ptr_edge_buff].weight < e[pos[edgebuff[ptr_edge_buff].node2]].weight)
						e[pos[edgebuff[ptr_edge_buff].node2]].weight=edgebuff[ptr_edge_buff].weight;
				}
				ptr_edge_buff++;
				if(ptr_edge_buff==edge_read){
					if(!feof(tempfile)){
						edge_read = fread(edgebuff,SZ_EDGE,EDGE_PER_BLK,tempfile);
						ptr_edge_buff=0;
					}
				}
			}
			
			while (!feof(r_tempfile) || ptr_r_edge_buff!=r_edge_read) {
				if(r_edgebuff[ptr_r_edge_buff].node2>u) 
					break;
				if(r_edgebuff[ptr_r_edge_buff].node2 == u){
					if(out_neighborlist[r_edgebuff[ptr_r_edge_buff].node1] != u)
					{
						out_neighborlist[r_edgebuff[ptr_r_edge_buff].node1] = u;
						indexNode[ptr_node_index-1].in_degree++;
						r_e[ptr_r_e].nodeid=r_edgebuff[ptr_r_edge_buff].node1; //
						r_e[ptr_r_e].weight=r_edgebuff[ptr_r_edge_buff].weight;
						r_pos[r_edgebuff[ptr_r_edge_buff].node1]=ptr_r_e;
						ptr_r_e++;
						
						if(ptr_r_e>=max_edge_m){  //need to flush out the edges
							ptr_node_index--;
							int last_ptr_node = ptr_node_index; 
							
							node_number_total+=ptr_node_index;
							for(int i=0; i<ptr_node_index; i++){
								sort(e.begin()+indexNode[i].begin, e.begin()+indexNode[i].begin+indexNode[i].out_degree, compare_weightedge);
								edge_number_total += indexNode[i].out_degree;
							}
							for(int i=0; i<ptr_node_index; i++){
								sort(r_e.begin()+indexNode[i].r_begin, r_e.begin()+indexNode[i].r_begin+indexNode[i].in_degree, compare_weightedge);
								r_edge_number_total += indexNode[i].in_degree;
							}
				
							write_graph(ptr_node_index, indexNode, e, r_e, write_io_1);
							write_io_1.flush();
							ptr_e = 0;
							ptr_r_e = 0;
							
							if(indexNode[last_ptr_node].out_degree >= max_edge_m-1 || indexNode[last_ptr_node].in_degree >= max_edge_m-1){
								cout<<"ERROR: Too High Degree Vertex"<<endl;
								exit(1);
							}
							
							for(int i=indexNode[last_ptr_node].begin; i<indexNode[last_ptr_node].begin+indexNode[last_ptr_node].out_degree; i++){
								e[ptr_e].nodeid=e[i].nodeid;
								e[ptr_e].weight=e[i].weight;
								pos[e[i].nodeid]=ptr_e;
								ptr_e++;
							}
							
							for(int i=indexNode[last_ptr_node].r_begin; i<indexNode[last_ptr_node].r_begin+indexNode[last_ptr_node].in_degree; i++){
								r_e[ptr_r_e].nodeid=r_e[i].nodeid;
								r_e[ptr_r_e].weight=r_e[i].weight;
								r_pos[r_e[i].nodeid]=ptr_r_e;
								ptr_r_e++;
							}
							
							ptr_node_index=1;
						}
					}
					else if(r_edgebuff[ptr_r_edge_buff].weight < r_e[r_pos[r_edgebuff[ptr_r_edge_buff].node1]].weight)
						r_e[r_pos[r_edgebuff[ptr_r_edge_buff].node1]].weight=r_edgebuff[ptr_r_edge_buff].weight;
				}
				ptr_r_edge_buff++;
				if(ptr_r_edge_buff==r_edge_read){
					if(!feof(r_tempfile)){
						r_edge_read = fread(r_edgebuff,SZ_EDGE,EDGE_PER_BLK,r_tempfile);
						ptr_r_edge_buff=0;
					}
				}
			}
			
		}
		
	}
	
	//for the last part
	if(ptr_node_index > 0){
		node_number_total+=ptr_node_index;
		for(int i=0; i<ptr_node_index; i++){
			sort(e.begin()+indexNode[i].begin, e.begin()+indexNode[i].begin+indexNode[i].out_degree, compare_weightedge);
			edge_number_total += indexNode[i].out_degree;
		}
		for(int i=0; i<ptr_node_index; i++){
			sort(r_e.begin()+indexNode[i].r_begin, r_e.begin()+indexNode[i].r_begin+indexNode[i].in_degree, compare_weightedge);
			r_edge_number_total += indexNode[i].in_degree;
		}
		
		write_graph(ptr_node_index, indexNode, e, r_e, write_io_1);
		write_io_1.flush();
		ptr_e = 0;
		ptr_r_e = 0;
	}
	
	fclose(tempfile);
	fclose(r_tempfile);
	fclose(file);
	fclose(remainfile);
	
	if(node_number_total==0)
		stopis=true;
	else 
		stopis=false;
	
	cout<<"The number of deleted nodes:" <<remarknumber <<endl;
	cout<<"Edge number in total:" << edge_number_total <<endl;
//	cout << "Reverse edge number in total: " << r_edge_number_total << endl;
	cout<<"Node number in total:" << node_number_total << endl;
	
	tt.stop();
	cout<<"The time used for getting independent set and adding augmented edges:" << tt.GetRuntime() << " seconds" << endl;
}

//when the memory for the inner loop is full, applying triangle inequality 
void Graph::tri_loop(int & ptr_node_index, vector <noderange> & indexNode, vector <weightedge> & e, vector <weightedge> & r_e, vector <int> & out_curr, vector <int> & in_curr, vector <int> & new_out_degree, vector <int> & new_in_degree, vector <int> & pos_node_inner, vector <bool> & trin, vector <bool> & r_trin, vector <noderange> & indexNode_inner, vector <weightedge> & e_inner, vector <weightedge> & r_e_inner, int & ptr_node_index_inner, int & ptr_e_inner, int & ptr_r_e_inner)
{
	int lastid=indexNode_inner[ptr_node_index_inner-1].nodeid;
	for(int i = 0; i<ptr_node_index; i++){ //u
		while(out_curr[i]<indexNode[i].begin+indexNode[i].out_degree) {
			if(e[out_curr[i]].nodeid > lastid) 
				break; //handle lastid
			if(pos_node_inner[e[out_curr[i]].nodeid] != -1) {
				int k=indexNode_inner[pos_node_inner[e[out_curr[i]].nodeid]].begin;
				int end=k +indexNode_inner[pos_node_inner[e[out_curr[i]].nodeid]].out_degree;
				for(int j=indexNode[i].begin; j<indexNode[i].begin+indexNode[i].out_degree; j++){
					if(trin[j]){ // the edge not be removed by triangle inequality 
						while(k<end){
							if(e_inner[k].nodeid>e[j].nodeid)
								break;	
							if(e_inner[k].nodeid==e[j].nodeid){  // triangle inequality
								if(e[j].weight>=e_inner[k].weight+e[out_curr[i]].weight){
									trin[j]=false;	
									new_out_degree[i]--;	
								}
							}
							k++;	
						}

						if(k==end)
							break;
					}
				}
			}
			out_curr[i]++;
		}
	}
	
	for(int i = 0; i<ptr_node_index; i++){ //u
		while(in_curr[i]<indexNode[i].r_begin+indexNode[i].in_degree) {
			if(r_e[in_curr[i]].nodeid > lastid) 
				break; //handle lastid
			if(pos_node_inner[r_e[in_curr[i]].nodeid] != -1) {
				int k=indexNode_inner[pos_node_inner[r_e[in_curr[i]].nodeid]].r_begin;
				int end=k +indexNode_inner[pos_node_inner[r_e[in_curr[i]].nodeid]].in_degree;
				for(int j=indexNode[i].r_begin; j<indexNode[i].r_begin+indexNode[i].in_degree; j++){
					if(r_trin[j]){ // the edge not be removed by triangle inequality 
						while(k<end){
							if(r_e_inner[k].nodeid>r_e[j].nodeid)
								break;	
							if(r_e_inner[k].nodeid==r_e[j].nodeid){  // triangle inequality
								if(r_e[j].weight>=r_e_inner[k].weight+r_e[in_curr[i]].weight){
									r_trin[j]=false;	
									new_in_degree[i]--;	
								}
							}
							k++;	
						}

						if(k==end)
							break;
					}
				}
			}
			in_curr[i]++;
		}
	}
		
	
	
	ptr_node_index_inner = 0;
	ptr_e_inner = 0;
	ptr_r_e_inner = 0;
}


// when the memory for outer loop is full, call the function to load inner loop 
void Graph::inner_loop(int & ptr_node_index, int & ptr_e, int & ptr_r_e, vector <noderange> & indexNode, vector <weightedge> & e, vector <weightedge> & r_e, vector <int> & out_curr, vector <int> & in_curr, vector <int> & new_out_degree, vector <int> & new_in_degree, ReadBuffer & read_inner, WriteBuffer & write_io, long & node_number_total, long & edge_number_total, long & r_edge_number_total)
{
	vector <int> pos_node_inner;
	pos_node_inner.resize(V+1, -1);
	vector <bool> trin;
	trin.resize(max_edge_m+1, true);
	vector <bool> r_trin;
	r_trin.resize(max_edge_m+1, true);

	vector <noderange> indexNode_inner;
	indexNode_inner.resize(INNER_NODE_NUMBER+1);	
	vector <weightedge> e_inner;
	e_inner.resize(INNER_EDGE_NUMBER+1);
	vector <weightedge> r_e_inner;
	r_e_inner.resize(INNER_EDGE_NUMBER+1);

	int ptr_node_index_inner = 0, ptr_e_inner=0, ptr_r_e_inner = 0;
	int u_inner, u_out_deg_inner, u_in_deg_inner;
	
	while(!read_inner.isend){
		read_inner.read(&u_inner);
		read_inner.read(&u_out_deg_inner);
		read_inner.read(&u_in_deg_inner);
		if(ptr_node_index_inner >= INNER_NODE_NUMBER || (ptr_e_inner+u_out_deg_inner)>= INNER_EDGE_NUMBER || (ptr_r_e_inner+u_in_deg_inner)>= INNER_EDGE_NUMBER){ //the memory for the inner loop is full
			tri_loop(ptr_node_index, indexNode, e, r_e, out_curr, in_curr, new_out_degree, new_in_degree, pos_node_inner, trin, r_trin, indexNode_inner, e_inner, r_e_inner, ptr_node_index_inner, ptr_e_inner, ptr_r_e_inner);
		}

		indexNode_inner[ptr_node_index_inner].nodeid = u_inner;
		indexNode_inner[ptr_node_index_inner].out_degree = u_out_deg_inner;
		indexNode_inner[ptr_node_index_inner].in_degree = u_in_deg_inner;
		indexNode_inner[ptr_node_index_inner].begin = ptr_e_inner;
		indexNode_inner[ptr_node_index_inner].r_begin = ptr_r_e_inner;
		pos_node_inner[u_inner] = ptr_node_index_inner;
		ptr_node_index_inner++;

		for(int i=0;i<u_out_deg_inner;i++){
			read_inner.read(&e_inner[ptr_e_inner].nodeid);
			read_inner.read(&e_inner[ptr_e_inner].weight);
			ptr_e_inner++;
		}
		for(int i=0;i<u_in_deg_inner;i++){
			read_inner.read(&r_e_inner[ptr_r_e_inner].nodeid);
			read_inner.read(&r_e_inner[ptr_r_e_inner].weight);
			ptr_r_e_inner++;
		}
	}
	if(ptr_node_index_inner > 0){ //handle the last part
		tri_loop(ptr_node_index, indexNode, e, r_e, out_curr, in_curr, new_out_degree, new_in_degree, pos_node_inner, trin, r_trin, indexNode_inner, e_inner, r_e_inner, ptr_node_index_inner, ptr_e_inner, ptr_r_e_inner);
	}

	//wirte out the remaing edges after applying triangle inequality
	for(int i=0; i<ptr_node_index; i++){
		node_number_total++;
		int ptr = 0;
		write_io.write(&indexNode[i].nodeid, 1);
		write_io.write(&new_out_degree[i], 1);
		write_io.write(&new_in_degree[i], 1);
		edge_number_total+= new_out_degree[i];
		r_edge_number_total+=new_in_degree[i];
		for(int j=indexNode[i].begin; j<indexNode[i].begin+indexNode[i].out_degree; j++){
			if(trin[j]){ 
				write_io.write(&e[j].nodeid, 1);
				write_io.write(&e[j].weight, 1);
			}
		}
		for(int j=indexNode[i].r_begin; j<indexNode[i].r_begin+indexNode[i].in_degree; j++){
			if(r_trin[j]){ 
				write_io.write(&r_e[j].nodeid, 1);
				write_io.write(&r_e[j].weight, 1);
			}
		}
	}
	write_io.flush();

	ptr_node_index = 0;
	ptr_e = 0;
	ptr_r_e = 0;
}

// applying triangle inequality to make the shrink the graph size 
void Graph::triangle_inequality(const char* filePath1, const char* filePath2)
{
	Timer tt;
	tt.start();
	
//	cout<<"Start the triangle inequality:" <<endl;
	FILE* remain = fopen(filePath1,"rb");
	FILE* remain2 = fopen(filePath1,"rb");
	FILE* newgraph = fopen(filePath2,"wb");

	ReadBuffer read_io(remain);
	WriteBuffer write_io(newgraph);
	
	vector <int> new_out_degree;
	new_out_degree.resize(max_node_m+1);
	vector <int> new_in_degree;
	new_in_degree.resize(max_node_m+1);
	vector <int> out_curr;
	out_curr.resize(max_node_m+1);
	vector <int> in_curr;
	in_curr.resize(max_node_m+1);
	
	vector <noderange> indexNode;
	indexNode.resize(max_node_m+1);	
	vector <weightedge> e;
	e.resize(max_edge_m+1);
	vector <weightedge> r_e;
	r_e.resize(max_edge_m+1);
	
	int ptr_node_index=0, ptr_e=0, ptr_r_e = 0;
	long edge_number_total = 0, r_edge_number_total = 0, node_number_total = 0;
	
	int u, u_out_deg, u_in_deg;
	int v, weight;
	while(!read_io.isend)
	{
		read_io.read(&u);
		read_io.read(&u_out_deg);
		read_io.read(&u_in_deg);
		if(ptr_node_index >= max_node_m || (ptr_e+u_out_deg)>= max_edge_m || (ptr_r_e+u_in_deg) >= max_edge_m ){ //the memory for outer loop is full
			ReadBuffer read_io_2(remain2);
			inner_loop(ptr_node_index, ptr_e, ptr_r_e, indexNode, e, r_e, out_curr, in_curr, new_out_degree, new_in_degree, read_io_2, write_io, node_number_total, edge_number_total, r_edge_number_total);
		}

		indexNode[ptr_node_index].nodeid = u;
		indexNode[ptr_node_index].out_degree = u_out_deg;
		indexNode[ptr_node_index].in_degree = u_in_deg;
		indexNode[ptr_node_index].begin = ptr_e;
		indexNode[ptr_node_index].r_begin = ptr_r_e;
		new_out_degree[ptr_node_index] = u_out_deg;
		new_in_degree[ptr_node_index] = u_in_deg;
		out_curr[ptr_node_index] = ptr_e;
		in_curr[ptr_node_index] = ptr_r_e;
		ptr_node_index++;

		for(int i=0;i<u_out_deg;i++){
			read_io.read(&e[ptr_e].nodeid);
			read_io.read(&e[ptr_e].weight);
			ptr_e++;
		}
		
		for(int i=0;i<u_in_deg;i++){
			read_io.read(&r_e[ptr_r_e].nodeid);
			read_io.read(&r_e[ptr_r_e].weight);
			ptr_r_e++;
		}
	}
	if(ptr_node_index>0){ // handle the last part
		ReadBuffer read_io_2(remain2);
		inner_loop(ptr_node_index, ptr_e, ptr_r_e, indexNode, e, r_e, out_curr, in_curr, new_out_degree, new_in_degree, read_io_2, write_io, node_number_total, edge_number_total, r_edge_number_total);
	}
	
	fclose(newgraph);
	fclose(remain);
	fclose(remain2);
	
	cout<<"After the round, The remaining graph:" <<endl;
	cout<<"Edge number in total:" << edge_number_total <<endl;
//	cout << "Reverse edge number in total: " << r_edge_number_total << endl;
	cout<<"Node number in total:" << node_number_total << endl;
	
	v_num_new = node_number_total;
	e_num_new = edge_number_total;
	ratio = (e_num_new+v_num_new)*1.0/(e_num_old+v_num_old); //get the ratio comparing with previous graph
	
//	cout<<"Ratio:" << ratio <<endl;
	tt.stop();
	cout<<"The time used for triangle inequality:" << tt.GetRuntime() << " seconds" << endl;
}

//the inital labels, it is called when computing independent set 
void Graph::labelInit(int &u, int &u_out_deg, int &u_in_deg, vector <weightedge> & neighbor, WriteBuffer & write_label)
{
	int u_out_deg_new = u_out_deg+1;
	int u_in_deg_new = u_in_deg+1;
	write_label.write(&u, 1);
	write_label.write(&u_out_deg_new, 1);
	write_label.write(&u_in_deg_new, 1);

	bool finish = false;
	int dis = 0;
	int x = u;
	for(int i = 0; i < u_out_deg; i++)
	{
		if(!finish)
		{
			if(neighbor[i].nodeid > u)
			{
				finish = true;
				write_label.write(&x, 1);
				write_label.write(&dis, 1);
			}
		}
		write_label.write(&neighbor[i].nodeid, 1);
		write_label.write(&neighbor[i].weight, 1);
	}
	if(!finish)
	{
		finish = true;
		write_label.write(&x, 1);
		write_label.write(&dis, 1);
	}
	
	finish = false;
	for(int i = u_out_deg; i<u_in_deg+u_out_deg; i++)
	{
		if(!finish)
		{
			if(neighbor[i].nodeid > u)
			{
				finish = true;
				write_label.write(&x, 1);
				write_label.write(&dis, 1);
			}
		}
		write_label.write(&neighbor[i].nodeid, 1);
		write_label.write(&neighbor[i].weight, 1);
	}
	if(!finish)
	{
		finish = true;
		write_label.write(&x, 1);
		write_label.write(&dis, 1);
	}
	
	write_label.flush();
}

// when the memory for the inner loop is full, updating the labels
bool Graph::label_loop(int & ptr_node_index, vector <noderange> & indexNode, vector <weightedge> & e, vector <weightedge> & r_e, vector < vector <weightedge> > & e2, vector < vector <weightedge> > & r_e2, vector <int> & curr, vector <int> & r_curr, vector <int> & pos_node_inner, int &ptr_node_index_inner, int &ptr_e_inner, int &ptr_r_e_inner, vector <noderange> & indexNode_inner, vector <weightedge> & e_inner, vector <weightedge> & r_e_inner, long & label_num)
{
	int lastid=indexNode_inner[ptr_node_index_inner-1].nodeid;
	for(int i=0; i<ptr_node_index; i++) {      //out_label
		while(curr[i]<indexNode[i].begin+indexNode[i].out_degree) {
			if(e[curr[i]].nodeid > lastid) 
				break;
			if(pos_node_inner[e[curr[i]].nodeid] != -1) { //e[curr[i]].nodeid is found
				int k=indexNode_inner[pos_node_inner[e[curr[i]].nodeid]].begin;
				int end=k+indexNode_inner[pos_node_inner[e[curr[i]].nodeid]].out_degree;
				int size=(int)e2[i].size();
				for(int j=0; j<size; j++){
					while(k<end){
						if(e_inner[k].nodeid>e2[i][j].nodeid) 
							break;
						if(e_inner[k].nodeid==e2[i][j].nodeid){
							if(e2[i][j].weight>=e_inner[k].weight+e[curr[i]].weight)
								e2[i][j].weight=e_inner[k].weight+e[curr[i]].weight;
						}
						if(e_inner[k].nodeid<e2[i][j].nodeid){  // add new label
							weightedge temp;
							temp.nodeid=e_inner[k].nodeid;
							temp.weight=e_inner[k].weight+e[curr[i]].weight;
							e2[i].push_back(temp);
							
							label_num++;
							if(label_num >= limit_label_m) // memory is full
								return true;
						}
						k++;
					}
					if(k==end) break;
				}
				while(k<end){ // add new label
					weightedge temp;
					temp.nodeid=e_inner[k].nodeid;
					temp.weight=e_inner[k].weight+e[curr[i]].weight;
					e2[i].push_back(temp);
					k++;
					
					label_num++;
					if(label_num >= limit_label_m) // memory is full
						return true;
				}
				sort(e2[i].begin(), e2[i].end(), compare_weightedge);
			}
			curr[i]++;
		}
	}


	for(int i=0; i<ptr_node_index; i++) {      //in_label
		while(r_curr[i]<indexNode[i].r_begin+indexNode[i].in_degree) {
			if(r_e[r_curr[i]].nodeid > lastid) 
				break;
			if(pos_node_inner[r_e[r_curr[i]].nodeid] != -1) { //e[curr[i]].nodeid is found
				int k=indexNode_inner[pos_node_inner[r_e[r_curr[i]].nodeid]].r_begin;
				int end=k+indexNode_inner[pos_node_inner[r_e[r_curr[i]].nodeid]].in_degree;
				int size=(int)r_e2[i].size();
				for(int j=0; j<size; j++){
					while(k<end){
						if(r_e_inner[k].nodeid>r_e2[i][j].nodeid) 
							break;
						if(r_e_inner[k].nodeid==r_e2[i][j].nodeid){
							if(r_e2[i][j].weight>=r_e_inner[k].weight+r_e[r_curr[i]].weight)
								r_e2[i][j].weight=r_e_inner[k].weight+r_e[r_curr[i]].weight;
						}
						if(r_e_inner[k].nodeid<r_e2[i][j].nodeid){  // add new label
							weightedge temp;
							temp.nodeid=r_e_inner[k].nodeid;
							temp.weight=r_e_inner[k].weight+r_e[r_curr[i]].weight;
							r_e2[i].push_back(temp);
							
							label_num++;
							if(label_num >= limit_label_m) // memory is full
								return true;
						}
						k++;
					}
					if(k==end) break;
				}
				while(k<end){ // add new label
					weightedge temp;
					temp.nodeid=r_e_inner[k].nodeid;
					temp.weight=r_e_inner[k].weight+r_e[r_curr[i]].weight;
					r_e2[i].push_back(temp);
					k++;
					
					label_num++;
					if(label_num >= limit_label_m) // memory is full
						return true;
				}
				sort(r_e2[i].begin(), r_e2[i].end(), compare_weightedge);
			}
			r_curr[i]++;
		}
	}

	ptr_node_index_inner=0;
	ptr_e_inner=0;
	ptr_r_e_inner=0;
	
	return false;
}


int Graph::labelUpdate_inner_loop(int & ptr_node_index, int & ptr_e, int & ptr_r_e, vector <noderange> & indexNode, vector <weightedge> & e, vector <weightedge> & r_e, vector < vector <weightedge> > & e2, vector < vector <weightedge> > & r_e2, vector <int> & curr, vector <int> & r_curr, FILE *sortfile, int lev, WriteBuffer &write_level,  WriteBuffer &write_label, long & curr_offset, vector <long> & offset)
{	
	vector <int> pos_node_inner;
	pos_node_inner.resize(V+1, -1);
	
	vector <noderange> indexNode_inner;
	indexNode_inner.resize(LABEL_INNER_NODE_NUMBER+1);	
	vector <weightedge> e_inner;
	e_inner.resize(LABEL_INNER_EDGE_NUMBER+1);
	vector <weightedge> r_e_inner;
	r_e_inner.resize(LABEL_INNER_EDGE_NUMBER+1);
	
	ReadBuffer read_inner(sortfile);
	int ptr_node_index_inner=0, ptr_e_inner=0, ptr_r_e_inner = 0;
	int u_inner, u_out_deg_inner, u_in_deg_inner;
	
	long label_num=ptr_e + ptr_r_e;
	bool label_exceed_mem = false;
	while(!read_inner.isend){
		read_inner.read(&u_inner);
		read_inner.read(&u_out_deg_inner);
		read_inner.read(&u_in_deg_inner);
		if(ptr_node_index_inner >= LABEL_INNER_NODE_NUMBER || (ptr_e_inner+u_out_deg_inner)>= LABEL_INNER_EDGE_NUMBER || (ptr_r_e_inner+u_in_deg_inner)>= LABEL_INNER_EDGE_NUMBER) { //the memory for the inner loop is full
			label_exceed_mem=label_loop(ptr_node_index, indexNode, e, r_e, e2, r_e2, curr, r_curr, pos_node_inner, ptr_node_index_inner, ptr_e_inner, ptr_r_e_inner, indexNode_inner, e_inner, r_e_inner, label_num);
			if(label_exceed_mem) 
				return 0;	
		}

		indexNode_inner[ptr_node_index_inner].nodeid = u_inner;
		indexNode_inner[ptr_node_index_inner].out_degree = u_out_deg_inner;
		indexNode_inner[ptr_node_index_inner].in_degree = u_in_deg_inner;
		indexNode_inner[ptr_node_index_inner].begin = ptr_e_inner;
		indexNode_inner[ptr_node_index_inner].r_begin = ptr_r_e_inner;
		pos_node_inner[u_inner] = ptr_node_index_inner;
		ptr_node_index_inner++;		

		for(int i=0;i<u_out_deg_inner;i++){
			read_inner.read(&e_inner[ptr_e_inner].nodeid);
			read_inner.read(&e_inner[ptr_e_inner].weight);
			ptr_e_inner++;
		}
		for(int i=0;i<u_in_deg_inner;i++){
			read_inner.read(&r_e_inner[ptr_r_e_inner].nodeid);
			read_inner.read(&r_e_inner[ptr_r_e_inner].weight);
			ptr_r_e_inner++;
		}
	}
	if(ptr_node_index_inner>0){ //handle the last part
		label_exceed_mem=label_loop(ptr_node_index, indexNode, e, r_e, e2, r_e2, curr, r_curr, pos_node_inner, ptr_node_index_inner, ptr_e_inner, ptr_r_e_inner, indexNode_inner, e_inner, r_e_inner, label_num);
		if(label_exceed_mem)
			return 0;
	}

	//write the labels into labelfile and a new file for level[lev]
	int size=0, r_size=0;
	for(int i=0; i<ptr_node_index; i++){  
		size = (int) e2[i].size();
		r_size = (int) r_e2[i].size();
		write_label.write(&indexNode[i].nodeid, 1);
		write_label.write(&size, 1);
		write_label.write(&r_size, 1);
		write_level.write(&indexNode[i].nodeid, 1);
		write_level.write(&size, 1);
		write_level.write(&r_size, 1);
		
		offset[indexNode[i].nodeid]=curr_offset;
		curr_offset += 12+5*(size+r_size);  // u:4, u_out_deg:4, u_in_deg:4, v: 4, weight: 1

		for(int j=0; j<size; j++){
			write_label.write(&e2[i][j].nodeid, 1);
			write_label.write(&e2[i][j].weight, 1);
			write_level.write(&e2[i][j].nodeid, 1);
			write_level.write(&e2[i][j].weight, 1);
		}
		for(int j=0; j<r_size; j++){
			write_label.write(&r_e2[i][j].nodeid, 1);
			write_label.write(&r_e2[i][j].weight, 1);
			write_level.write(&r_e2[i][j].nodeid, 1);
			write_level.write(&r_e2[i][j].weight, 1);
		}
	}
	write_label.flush();
	write_level.flush();
	for(int i=0; i<ptr_node_index; i++){
		vector<weightedge>().swap(e2[i]);
		vector<weightedge>().swap(r_e2[i]);
	}
	
	int num = ptr_node_index;
	ptr_node_index = 0;
	ptr_e = 0;
	ptr_r_e = 0;
	
	return num;
}

// update the labels for level lev
void Graph::labelUpdate(int lev, WriteBuffer &write_level, WriteBuffer &write_label, long & curr_offset, vector <long> & offset)
{
	
	FILE * sortfile; 
	vector <minvid> vidarray;
	vidarray.resize(round+1);
	
	int u, u_out_deg, u_in_deg;
	int v, weight;
	
	if(lev+1!=round-1) //merge and sort label L[lev+1...round-1] node into vidarrary
	{
		sortfile = tmpfile();
		WriteBuffer write_io(sortfile);
		ReadBuffer read_seed[round+1];

		for(int i=lev+1; i<round; i++){
			read_seed[i].open(level[i]);
			read_seed[i].read(&vidarray[i].nodeid);
			vidarray[i].levelindex=i;
		}
		sort(vidarray.begin()+lev+1, vidarray.begin()+round, compare_vid);

		int level_size = round-lev-1, min_level_index=0;
		while(level_size!=0){
			min_level_index=vidarray[lev+1].levelindex;
			u = vidarray[lev+1].nodeid;
			read_seed[min_level_index].read(&u_out_deg);
			read_seed[min_level_index].read(&u_in_deg);
			
			write_io.write(&u, 1);
			write_io.write(&u_out_deg, 1);
			write_io.write(&u_in_deg, 1);
			
			for(int i=0; i<u_out_deg+u_in_deg; i++){
				read_seed[min_level_index].read(&v);
				read_seed[min_level_index].read(&weight);
				write_io.write(&v, 1);
				write_io.write(&weight, 1);
			}

			if(!read_seed[min_level_index].isend){ //sort
				read_seed[min_level_index].read(&vidarray[lev+1].nodeid);
				vidarray[lev+1].levelindex=min_level_index;			
				sort(vidarray.begin()+lev+1, vidarray.begin()+lev+1+level_size, compare_vid);
			}
			else {
				level_size--;
				vidarray.erase(vidarray.begin()+lev+1);
			}
		}
		write_io.flush();
	}
	else
		sortfile = level[lev+1];

	vector <noderange> indexNode;
	indexNode.resize(max_node_label_m+1);	
	vector <weightedge> e;
	e.resize(max_edge_label_m+1);
	vector <weightedge> r_e;
	r_e.resize(max_edge_label_m+1);
	vector < vector <weightedge> > e2;
	e2.resize(max_node_label_m+1);
	vector < vector <weightedge> > r_e2;
	r_e2.resize(max_node_label_m+1);
	int ptr_node_index=0, ptr_e=0, ptr_r_e = 0;
	
	vector <int> curr;
	curr.resize(max_node_label_m+1);
	vector <int> r_curr;
	r_curr.resize(max_node_label_m+1);

	ReadBuffer read_lev(level[lev]);
	
	vector <noderange> indexNode_pack;
	indexNode_pack.resize(max_node_label_m+1);	
	vector <weightedge> e_pack;
	e_pack.resize(max_edge_label_m+1); 
	vector <weightedge> r_e_pack;
	r_e_pack.resize(max_edge_label_m+1); 
	int ptr_node_index_pack=0;
	
	int ptr_total=0, ptr_start=0;
	int processed_num=0;
	while(!read_lev.isend){
		read_lev.read(&u);
		read_lev.read(&u_out_deg);
		read_lev.read(&u_in_deg);
		if(ptr_node_index >= max_node_label_m || (ptr_e+u_out_deg)>= max_edge_label_m || (ptr_e+u_in_deg)>= max_edge_label_m) { // when the memory for outer loop is full
			indexNode_pack = indexNode;
			e_pack = e;
			r_e_pack = r_e;
			ptr_node_index_pack = ptr_node_index;
			ptr_total=0;
			processed_num=0;
			
			while(ptr_total < ptr_node_index_pack){
				processed_num=labelUpdate_inner_loop(ptr_node_index, ptr_e, ptr_r_e, indexNode, e, r_e, e2, r_e2, curr, r_curr, sortfile, lev, write_level, write_label, curr_offset, offset);
				
				if(processed_num==0){  //out of memory, retry
		//			cout<<"processed_num:" <<processed_num<<endl;
				
					for(int i=0; i<max_node_label_m; i++){
						vector<weightedge>().swap(e2[i]);
						vector<weightedge>().swap(r_e2[i]);
					}
					max_node_label_m = max_node_label_m/2;
					
					ptr_node_index=0;
					ptr_e=0;
					ptr_r_e = 0;
					
					for(int i=ptr_start; i<ptr_node_index_pack; i++){
						if(ptr_node_index >= max_node_label_m || (ptr_e+indexNode_pack[i].out_degree)>= max_edge_label_m || (ptr_r_e+indexNode_pack[i].in_degree)>= max_edge_label_m){
							if(ptr_node_index<=1){
								fprintf(stderr,"ERROR: Too High Degree");
								exit(1);
							} 
							break;	
						}
						indexNode[ptr_node_index].nodeid = indexNode_pack[i].nodeid; 
						indexNode[ptr_node_index].out_degree = indexNode_pack[i].out_degree; 
						indexNode[ptr_node_index].in_degree = indexNode_pack[i].in_degree; 
						indexNode[ptr_node_index].begin = ptr_e;
						indexNode[ptr_node_index].r_begin = ptr_r_e;
						curr[ptr_node_index] = ptr_e;
						r_curr[ptr_node_index] = ptr_r_e;
						ptr_node_index ++;
				
						for(int j=indexNode_pack[i].begin; j<indexNode_pack[i].begin+indexNode_pack[i].out_degree; j++){
							e[ptr_e]=e_pack[j];
							e2[ptr_node_index-1].push_back(e[ptr_e]);
							ptr_e++;
						}
						
						for(int j=indexNode_pack[i].begin; j<indexNode_pack[i].r_begin+indexNode_pack[i].in_degree; j++){
							r_e[ptr_r_e]=r_e_pack[j];
							r_e2[ptr_node_index-1].push_back(r_e[ptr_r_e]);
							ptr_r_e++;
						}
					}
					
					
				}
				else {
					ptr_total+=processed_num;
					ptr_start+=processed_num;
					
					ptr_node_index=0;
					ptr_e=0;
					ptr_r_e = 0;
					
					for(int i=ptr_start; i<ptr_node_index_pack; i++){
						if(ptr_node_index >= max_node_label_m || (ptr_e+indexNode_pack[i].out_degree)>= max_edge_label_m || (ptr_r_e+indexNode_pack[i].in_degree)>= max_edge_label_m){
							if(ptr_node_index<=1){
								fprintf(stderr,"ERROR: Too High Degree");
								exit(1);
							} 
							break;	
						}
						indexNode[ptr_node_index].nodeid = indexNode_pack[i].nodeid; 
						indexNode[ptr_node_index].out_degree = indexNode_pack[i].out_degree; 
						indexNode[ptr_node_index].in_degree = indexNode_pack[i].in_degree; 
						indexNode[ptr_node_index].begin = ptr_e;
						indexNode[ptr_node_index].r_begin = ptr_r_e;
						curr[ptr_node_index] = ptr_e;
						r_curr[ptr_node_index] = ptr_r_e;
						ptr_node_index ++;
				
						for(int j=indexNode_pack[i].begin; j<indexNode_pack[i].begin+indexNode_pack[i].out_degree; j++){
							e[ptr_e]=e_pack[j];
							e2[ptr_node_index-1].push_back(e[ptr_e]);
							ptr_e++;
						}
						
						for(int j=indexNode_pack[i].begin; j<indexNode_pack[i].r_begin+indexNode_pack[i].in_degree; j++){
							r_e[ptr_r_e]=r_e_pack[j];
							r_e2[ptr_node_index-1].push_back(r_e[ptr_r_e]);
							ptr_r_e++;
						}
					}
					
				}
			}	
		}
	
		indexNode[ptr_node_index].nodeid = u;
		indexNode[ptr_node_index].out_degree = u_out_deg;
		indexNode[ptr_node_index].in_degree = u_in_deg;
		indexNode[ptr_node_index].begin = ptr_e;
		indexNode[ptr_node_index].r_begin = ptr_r_e;
		curr[ptr_node_index] = ptr_e;
		r_curr[ptr_node_index] = ptr_r_e;
		ptr_node_index++;
		for(int i=0;i<u_out_deg;i++){
			read_lev.read(&e[ptr_e].nodeid);
			read_lev.read(&e[ptr_e].weight);
			e2[ptr_node_index-1].push_back(e[ptr_e]);
			ptr_e++;
		}
		for(int i=0;i<u_in_deg;i++){
			read_lev.read(&r_e[ptr_r_e].nodeid);
			read_lev.read(&r_e[ptr_r_e].weight);
			r_e2[ptr_node_index-1].push_back(r_e[ptr_r_e]);
			ptr_r_e++;
		}
		
	}
	if(ptr_node_index>0){ //handle the last part
		indexNode_pack = indexNode;
		e_pack = e;
		r_e_pack = r_e;
		ptr_node_index_pack = ptr_node_index;
		ptr_total=0;
		processed_num=0;
		
		while(ptr_total < ptr_node_index_pack){
			processed_num=labelUpdate_inner_loop(ptr_node_index, ptr_e, ptr_r_e, indexNode, e, r_e, e2, r_e2, curr, r_curr, sortfile, lev, write_level, write_label, curr_offset, offset);
			
			if(processed_num==0){  //out of memory, retry
		//		cout<<"processed_num:" <<processed_num<<endl;
			
				for(int i=0; i<max_node_label_m; i++){
					vector<weightedge>().swap(e2[i]);
					vector<weightedge>().swap(r_e2[i]);
				}
				max_node_label_m = max_node_label_m/2;
				
				ptr_node_index=0;
				ptr_e=0;
				ptr_r_e = 0;
				
				for(int i=ptr_start; i<ptr_node_index_pack; i++){
					if(ptr_node_index >= max_node_label_m || (ptr_e+indexNode_pack[i].out_degree)>= max_edge_label_m || (ptr_r_e+indexNode_pack[i].in_degree)>= max_edge_label_m){
						if(ptr_node_index<=1){
							fprintf(stderr,"ERROR: Too High Degree");
							exit(1);
						} 
						break;	
					}
					indexNode[ptr_node_index].nodeid = indexNode_pack[i].nodeid; 
					indexNode[ptr_node_index].out_degree = indexNode_pack[i].out_degree; 
					indexNode[ptr_node_index].in_degree = indexNode_pack[i].in_degree; 
					indexNode[ptr_node_index].begin = ptr_e;
					indexNode[ptr_node_index].r_begin = ptr_r_e;
					curr[ptr_node_index] = ptr_e;
					r_curr[ptr_node_index] = ptr_r_e;
					ptr_node_index ++;
			
					for(int j=indexNode_pack[i].begin; j<indexNode_pack[i].begin+indexNode_pack[i].out_degree; j++){
						e[ptr_e]=e_pack[j];
						e2[ptr_node_index-1].push_back(e[ptr_e]);
						ptr_e++;
					}
					
					for(int j=indexNode_pack[i].begin; j<indexNode_pack[i].r_begin+indexNode_pack[i].in_degree; j++){
						r_e[ptr_r_e]=r_e_pack[j];
						r_e2[ptr_node_index-1].push_back(r_e[ptr_r_e]);
						ptr_r_e++;
					}
				}
				
				
			}
			else {
				ptr_total+=processed_num;
				ptr_start+=processed_num;
				
				ptr_node_index=0;
				ptr_e=0;
				ptr_r_e = 0;
				
				for(int i=ptr_start; i<ptr_node_index_pack; i++){
					if(ptr_node_index >= max_node_label_m || (ptr_e+indexNode_pack[i].out_degree)>= max_edge_label_m || (ptr_r_e+indexNode_pack[i].in_degree)>= max_edge_label_m){
						if(ptr_node_index<=1){
							fprintf(stderr,"ERROR: Too High Degree");
							exit(1);
						} 
						break;	
					}
					indexNode[ptr_node_index].nodeid = indexNode_pack[i].nodeid; 
					indexNode[ptr_node_index].out_degree = indexNode_pack[i].out_degree; 
					indexNode[ptr_node_index].in_degree = indexNode_pack[i].in_degree; 
					indexNode[ptr_node_index].begin = ptr_e;
					indexNode[ptr_node_index].r_begin = ptr_r_e;
					curr[ptr_node_index] = ptr_e;
					r_curr[ptr_node_index] = ptr_r_e;
					ptr_node_index ++;
			
					for(int j=indexNode_pack[i].begin; j<indexNode_pack[i].begin+indexNode_pack[i].out_degree; j++){
						e[ptr_e]=e_pack[j];
						e2[ptr_node_index-1].push_back(e[ptr_e]);
						ptr_e++;
					}
					
					for(int j=indexNode_pack[i].begin; j<indexNode_pack[i].r_begin+indexNode_pack[i].in_degree; j++){
						r_e[ptr_r_e]=r_e_pack[j];
						r_e2[ptr_node_index-1].push_back(r_e[ptr_r_e]);
						ptr_r_e++;
					}
				}
				
			}
		}	
	}
	
	if(lev+1!=round-1)
		fclose(sortfile);
}

//doing top down labeling 
void Graph::topdownLabel(const char* filePath1, const char* filePath2) //labelfile, offsetfile
{
	Timer tt;
	tt.start();

	vector <long> offset;
	offset.resize(V+1, -1);
	
	long curr_offset=0;

	ReadBuffer read_io(level[round-1]); //level[h]->labelFile
	FILE * labelFile = fopen(filePath1,"wb");
	WriteBuffer write_label(labelFile);

	int u, u_out_deg, u_in_deg;
	int v, weight;
	while(!read_io.isend){
		read_io.read(&u);
		read_io.read(&u_out_deg);
		read_io.read(&u_in_deg);

		offset[u]=curr_offset;
		curr_offset += 12+5*(u_out_deg+u_in_deg); // u:4, u_out_deg:4, u_in_deg:4, v: 4, weight: 1
		
		write_label.write(&u, 1);
		write_label.write(&u_out_deg, 1);
		write_label.write(&u_in_deg, 1);

		for(int i=0; i<u_out_deg+u_in_deg; i++){
			read_io.read(&v);
			read_io.read(&weight);
			write_label.write(&v, 1);
			write_label.write(&weight, 1);
		}
	}
	write_label.flush();

	FILE * new_levelfile;
	if(round >= 2){	
		for(int i=round-2; i>=0; i--){ // from h-1 to 1
			new_levelfile = tmpfile();
			WriteBuffer write_level(new_levelfile);
			cout<<"update label i:" << i << endl;
			labelUpdate(i, write_level, write_label, curr_offset, offset); //update the labels for level i
			fclose(level[i]);
			level[i] = new_levelfile;
		}
		for(int i=round-1; i>=0; i--){
			fclose(level[i]);
		}
	}
	write_label.flush();
	fclose(labelFile);

	//record the label offset in offset file 
	FILE * offsetFile = fopen(filePath2,"wb");
	long * write_buff_offset = (long *) malloc(BLK_SZ);
	int ptr_write_buff_offset = 0;
	for(int i=0; i<V; i++){
		write_buff_offset[ptr_write_buff_offset++]=offset[i];
		if(ptr_write_buff_offset==OFFSET_PER_BLK){
			fwrite(write_buff_offset,SZ_OFFSET,ptr_write_buff_offset,offsetFile);
			ptr_write_buff_offset = 0;
		}
	}
	
	if(ptr_write_buff_offset>0){
		fwrite(write_buff_offset,SZ_OFFSET,ptr_write_buff_offset,offsetFile);
		ptr_write_buff_offset = 0;
	}		
	fclose(offsetFile);
	free(write_buff_offset);

	tt.stop();
	cout<<"The time used for labeling:" << tt.GetRuntime() << " seconds" << endl;
}

//change the label file from int to byte
void Graph::labelinttobyte(const char* filePath1, const char* filePath2) //labelfile, finallabelfile
{
	FILE * intfile = fopen(filePath1,"rb");
	FILE * bytefile = fopen(filePath2,"wb");
	
	label_number_total=0;
	
	char * write_buff = (char *)malloc(BLK_SZ); 
	char * read_buff_1 = (char *)malloc(BLK_SZ);
	vertex_id * buff = (vertex_id *) read_buff_1;
	
	int num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,intfile);
	int ptr_buff = 0;
	
	int num_write=0;
	while(!feof(intfile) || ptr_buff!=num_read){
		int u = buff[ptr_buff];
		ptr_buff++;
		if(ptr_buff==num_read){
			num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,intfile);
			ptr_buff = 0;
		}
		
		memcpy(write_buff+num_write, &u, sizeof(u));
		num_write += sizeof(u);
		if(num_write > BLK_SZ-sizeof(int)){
			fwrite(write_buff,1,num_write,bytefile);
			num_write = 0;
		}
		
		int u_out_deg = buff[ptr_buff];
		ptr_buff++;
		if(ptr_buff==num_read){
			num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,intfile);
			ptr_buff = 0;
		}
		
		memcpy(write_buff+num_write, &u_out_deg, sizeof(u_out_deg));
		num_write += sizeof(u_out_deg);
		if(num_write > BLK_SZ-sizeof(int)){
			fwrite(write_buff,1,num_write,bytefile);
			num_write = 0;
		}
		
		int u_in_deg = buff[ptr_buff];
		ptr_buff++;
		if(ptr_buff==num_read){
			num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,intfile);
			ptr_buff = 0;
		}
		
		memcpy(write_buff+num_write, &u_in_deg, sizeof(u_in_deg));
		num_write += sizeof(u_in_deg);
		if(num_write > BLK_SZ-sizeof(int)){
			fwrite(write_buff,1,num_write,bytefile);
			num_write = 0;
		}
		
		label_number_total += u_out_deg+u_in_deg;

		for(int j=0; j<u_out_deg+u_in_deg; j++){
			int v = buff[ptr_buff];
			ptr_buff++;
			if(ptr_buff==num_read){
				num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,intfile);
				ptr_buff = 0;
			}
			
			memcpy(write_buff+num_write, &v, sizeof(v));
			num_write += sizeof(v);
			if(num_write > BLK_SZ-sizeof(char)){
				fwrite(write_buff,1,num_write,bytefile);
				num_write = 0;
			}

			char distance = (char) buff[ptr_buff];
			ptr_buff++;
			if(ptr_buff==num_read){
				num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,intfile);
				ptr_buff = 0; 
			}
			
			memcpy(write_buff+num_write, &distance, sizeof(distance));
			num_write += sizeof(distance);
			if(num_write > BLK_SZ-sizeof(int)){
				fwrite(write_buff,1,num_write,bytefile);
				num_write = 0;
			}

		}		
	}

		
	if(num_write>0){
		fwrite(write_buff,1,num_write,bytefile);
		num_write = 0;
	}
	
	fclose(bytefile);
	fclose(intfile);
	free(write_buff);
	free(read_buff_1);
}

//print out the statistics for remaing graph
void Graph::print_statistics(const char* filePath)
{
	FILE * graph_gk=fopen(filePath, "w");
	fprintf(graph_gk, "%ld %ld %ld %ld\n", V, v_num_new, e_num_new, label_number_total);
	fclose(graph_gk);
}

#endif

