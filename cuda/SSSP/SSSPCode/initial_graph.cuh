#ifndef INITIAL_GRAPH_HPP
#define INITIAL_GRAPH_HPP

#include <vector>

#include "user_specified_structures.cuh"

class neighbor {
public:
	Edge edgeValue;
	unsigned int srcIndex;
};

class initial_vertex {
public:
	Vertex vertexValue;
	Vertex_static VertexValueStatic;
	std::vector<neighbor> nbrs;
	initial_vertex():
		nbrs(0)
	{}
	Vertex& get_vertex_ref() {
		return vertexValue;
	}
};

typedef struct
{
    unsigned int src;
    unsigned int dst;
    int weight;
} SdwEdge;

#endif	//	INITIAL_GRAPH_HPP
