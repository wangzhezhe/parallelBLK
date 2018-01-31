
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "genresult.cuh"
#include <vector>
#include <algorithm>

using namespace std;

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

struct DataItem
{
    int start;
    int end;
    int weight;
    int flag;
    struct DataItem *next;
};

typedef struct HeadNode
{
    int maxweight;
    int flag;
    vector<struct DataItem *> edgeV;
} HeadNode;

int hashCode(int start)
{
    return start;
}

void vmodifymatColumn(int *cIndex, int *newvindex, int len)
{
    int i = 0;
    int val;
    for (i = 0; i < len; i++)
    {

        val = cIndex[i];
        cIndex[i] = newvindex[val];
    }
    return;
}

void modifymatColumn(MatrixInfo *mat, struct DataItem *hashCindex[], int *oldvtonew)
{
    int i;
    struct DataItem *entry = NULL;
    int cIndex, elemIndex, newcIndex;
    for (i = 0; i < mat->N; i++)
    {
        entry = hashCindex[i];

        while (entry != NULL)
        {
            cIndex = entry->start;
            elemIndex = entry->end;
            newcIndex = oldvtonew[cIndex];
            mat->cIndex[elemIndex] = newcIndex;
            //printf("elei %d oldc %d newc %d\n",elemIndex,cIndex,newcIndex);
            entry = entry->next;
        }
    }
}

int vsearch(vector<struct HeadNode *> &graph, int start, int end)
{
    //get the hash
    int hashIndex = start;

    if (graph.size() == 0)
    {
        return 0;
    }

    //struct DataItem *entry = (struct DataItem *)malloc(sizeof(struct DataItem));
    //entry = hashArray[hashIndex];

    if (graph[hashIndex]->edgeV.size() == 0)
    {
        return 0;
    }

    //range vector find the same one

    int i;
    int len = graph[hashIndex]->edgeV.size();
    for (i = 0; i < len; i++)
    {
        if (graph[hashIndex]->edgeV[i]->start == start && graph[hashIndex]->edgeV[i]->end == end)
        {

            graph[hashIndex]->edgeV[i]->weight++;
            if (graph[hashIndex]->edgeV[i]->weight > graph[hashIndex]->maxweight)
            {
                graph[hashIndex]->maxweight = graph[hashIndex]->edgeV[i]->weight;
            }
            return 1;
        }
    }

    return 0;
}

struct DataItem *search(struct DataItem *hashArray[], int start, int end)
{
    //get the hash
    int hashIndex = hashCode(start);

    if (hashArray == NULL)
    {
        return NULL;
    }

    struct DataItem *entry = (struct DataItem *)malloc(sizeof(struct DataItem));
    entry = hashArray[hashIndex];

    if (entry == NULL)
    {
        return NULL;
    }

    while (entry != NULL)
    {

        if (entry->start == start && entry->end == end)
        {
            return entry;
        }
        else
        {
            entry = entry->next;
        }
    }

    return NULL;
}

void vinsert(vector<struct HeadNode *> &graph, int start, int end, struct DataItem *maxEdge)
{
    //if exist, weight ++
    //struct DataItem *item = NULL;

    //printf("start vsearch s %d e %d\n", start, end);

    int update = vsearch(graph, start, end);

    if (update == 1)
    {
        return;
    }

    struct DataItem *newitem = (struct DataItem *)malloc(sizeof(struct DataItem));

    newitem->start = start;
    newitem->end = end;
    //default init weight is 1
    newitem->weight = 1;
    //init flag of every edge is 0
    newitem->flag = 0;
    newitem->next = NULL;
    //get the hash
    int hashIndex = hashCode(start);

    graph[hashIndex]->edgeV.push_back(newitem);

    if (newitem->weight > graph[hashIndex]->maxweight)
    {
        graph[hashIndex]->maxweight = newitem->weight;
    }
}

void insert(struct DataItem *hashArray[], int start, int end, struct DataItem *maxEdge)
{

    //if exist, weight ++
    struct DataItem *item = NULL;

    item = search(hashArray, start, end);

    if (item != NULL)
    {
        item->weight = item->weight + 1;
        if (maxEdge != NULL && item->weight > maxEdge->weight)
        {

            maxEdge->start = item->start;
            maxEdge->end = item->end;
            maxEdge->weight = item->weight;
        }
        return;
    }

    struct DataItem *newitem = (struct DataItem *)malloc(sizeof(struct DataItem));

    newitem->start = start;
    newitem->end = end;
    //default init weight is 1
    newitem->weight = 1;
    //init flag of every edge is 0
    newitem->flag = 0;
    newitem->next = NULL;
    //get the hash
    int hashIndex = hashCode(start);

    struct DataItem *entry = hashArray[hashIndex];
    if (entry == NULL)
    {
        hashArray[hashIndex] = newitem;
    }
    else
    {
        //move in array until an empty or deleted cell
        while (entry->next != NULL)
        {

            entry = entry->next;
        }
        entry->next = newitem;
    }

    if (maxEdge != NULL && newitem->weight > maxEdge->weight)
    {

        maxEdge->start = newitem->start;
        maxEdge->end = newitem->end;
        maxEdge->weight = newitem->weight;
    }
}

bool cmparray(struct DataItem *v1, struct DataItem *v2) //注意：本函数的参数的类型一定要与vector中元素的类型一致
{
    return v1->weight > v2->weight; //升序排列
}

bool cmphead(struct HeadNode *v1, struct HeadNode *v2) //注意：本函数的参数的类型一定要与vector中元素的类型一致
{
    return v1->maxweight > v2->maxweight; //升序排列
}

void sortgraph(vector<struct HeadNode *> &graph)
{
    sort(graph.begin(), graph.end(), cmphead);

    int len = graph.size();
    int i = 0;
    for (i = 0; i < len; i++)
    {
        if (graph[i]->edgeV.empty() != true)
        {
            sort(graph[i]->edgeV.begin(), graph[i]->edgeV.end(), cmparray);
        }
    }
    return;
}

void vdisplay(vector<struct HeadNode *> &graph)
{
    int i = 0, j = 0;
    for (i = 0; i < graph.size(); i++)
    {
        if (graph[i]->edgeV.empty() != true)
        {
            printf("index %d maxweight %d\n", i, graph[i]->maxweight);
        }

        for (j = 0; j < graph[i]->edgeV.size(); j++)
        {
            printf("s %d e %d w %d\n",
                   graph[i]->edgeV[j]->start, graph[i]->edgeV[j]->end, graph[i]->edgeV[j]->weight);
        }
    }
}

void display(struct DataItem *hashArray[], int maxrow)
{
    int i = 0;

    for (i = 0; i < maxrow; i++)
    {
        struct DataItem *entry = hashArray[i];

        if (entry != NULL)
        {
            while (entry != NULL)
            {
                printf("(s %d,e %d,w %d)\n", entry->start, entry->end, entry->weight);
                entry = entry->next;
            }
        }
        else
        {
            printf("~~\n");
        }
    }

    printf("\n");
}

void matTraverse(
    MatrixInfo *oldmat,
    MatrixInfo *newmat,
    int *newrownum,
    int *rowchange)
{
    int *vflag = (int *)calloc((oldmat->N + 1), sizeof(int));
    int oldrow, newrow;
    int i;
    for (i = 0; i < oldmat->nz; i++)
    {
        oldrow = oldmat->rIndex[i];
        if (i == 0)
        {
            newrow = 0;
        }
        else
        {
            if (oldmat->rIndex[i] == oldmat->rIndex[i - 1])
            {
                newrow = newmat->rIndex[i - 1];
            }
            else
            {
                newrow = newmat->rIndex[i - 1] + 1;
            }
        }

        newmat->rIndex[i] = newrow;
        rowchange[newrow] = oldrow;
    }
    *newrownum = newrow;
    return;
}

void addVbasedWarp(int *vaccess,
                   vector<struct HeadNode *> &graph,
                   int spanSize,
                   int workperWarp,
                   int warpNum,
                   struct DataItem *maxEdge)
{
    int warpid, i, j, s;
    int spanNum = workperWarp / spanSize;
    int start, end;
    int sver, ever;
    for (warpid = 0; warpid < warpNum; warpid++)
    {
        for (s = 0; s < spanNum; s++)
        {
            start = warpid * workperWarp + spanSize * s;
            end = min(start + spanSize, (warpid + 1) * workperWarp);
            //printf("start %d end %d\n", start, end);
            for (i = start; i < end; i++)
            {
                for (j = i + 1; j < end; j++)
                {
                    sver = min(vaccess[i], vaccess[j]);
                    ever = max(vaccess[i], vaccess[j]);
                    //printf("insert s %d e %d\n", i, j, sver, ever);
                    if (sver != ever)
                    {
                        vinsert(graph, sver, ever, maxEdge);
                    }
                }
            }
        }
    }
    return;
}

struct DataItem *findMax(int index, int eleN, struct DataItem *graph[], struct DataItem *reversegraph[])
{
    int i = 0;
    struct DataItem *maxE = (struct DataItem *)malloc(sizeof(struct DataItem));
    maxE->weight = 0;

    for (i = 0; i < eleN; i++)
    {
        //traverse graph
        struct DataItem *entry = graph[i];
        while (entry != NULL)
        {
            //printf("search entry s %d e %d f %d\n", entry->start, entry->end, entry->flag);
            if (entry->weight > maxE->weight && entry->flag == 0)
            {
                maxE->start = entry->start;
                maxE->end = entry->end;
                maxE->weight = entry->weight;
            }
            entry = entry->next;
        }
        //traverse reverseg
        entry = reversegraph[i];
        while (entry != NULL)
        {
            if (entry->weight > maxE->weight && entry->flag == 0)
            {
                maxE->start = entry->start;
                maxE->end = entry->end;
                maxE->weight = entry->weight;
            }
            entry = entry->next;
        }
    }
    //printf("maxe s %d e %d w %d\n", maxE->start, maxE->end, maxE->weight);
    return maxE;
}

struct DataItem *vgetMax(vector<struct HeadNode *> &graph, int index, int *vflag)
{

    int i = 0, j = 0;
    int lena = graph.size();
    int lenb;

    for (i = index; i < lena; i++)
    {
        lenb = graph[i]->edgeV.size();
        if (graph[i]->flag == 1)
        {
            continue;
        }
        //printf("i %d lena %d lenb %d\n", i, lena, lenb);
        for (j = 0; j < lenb; j++)
        {
            //wht there are seg f?

            if (graph[i]->edgeV[j]->flag != 1)
            {
                if (vflag[j] == 1 &&  vflag[i] == 1)
                {
                    graph[i]->edgeV[j]->flag = 1;
                    continue;
                }
                graph[i]->edgeV[j]->flag = 1;
                return graph[i]->edgeV[j];
            }
        }
        graph[i]->flag = 1;
    }

    return NULL;
}

void vgetReorderedVindex(
    vector<struct HeadNode *> &graph,
    int *reorderv,
    int spanSize,
    int vecNum)
{
    int count = 0;
    int *vflag = (int *)calloc(vecNum, sizeof(int));
    //struct DataItem *newMax = graph[0]->edgeV[0];
    struct DataItem *newMax = vgetMax(graph, 0, vflag);
    int countinspan = 0;
    int s, e;
    struct DataItem *maxe1 = NULL;
    struct DataItem *maxe2 = NULL;
    struct DataItem *globalMax = NULL;
    while (count < vecNum)
    {
        //if (count % 10000 == 0)
        if (count % 10000 == 0 || count > 36000)
        {
            printf("curr count %d ns %d ne %d nw %d\n", count, newMax->start, newMax->end, newMax->weight);
        }
        if (newMax == NULL)
        {
            printf("failed, newMax is NULL\n");
            exit(1);
        }

        s = newMax->start;
        e = newMax->end;
        //printf("curr s %d e %d flag %d\n", s, e, newMax->flag);
        newMax->flag = 1;
        if (vflag[s] == 0)
        {
            //the reordered vector
            reorderv[count] = s;
            count++;
            countinspan++;
            vflag[s] = 1;
            //printf("add index %d into vector \n", s);
        }

        if (vflag[e] == 0)
        {
            reorderv[count] = e;
            count++;
            countinspan++;
            vflag[e] = 1;
            //printf("add index %d into vector\n", e);
        }
        if (countinspan >= spanSize)
        {
            //update new node, continue

            if (globalMax != NULL)
            {
                //printf("curr max start %d \n", globalMax->start);
                newMax = vgetMax(graph, globalMax->start, vflag);
            }
            if (newMax == NULL)
            {
                newMax = vgetMax(graph, 0, vflag);
            }

            //newMax = vgetMax(graph, 0);
            //printf("get span s %d e %d\n", newMax->start, newMax->end);
            countinspan = 0;
            continue;
        }
        maxe1 = vgetMax(graph, s, vflag);
        maxe2 = vgetMax(graph, e, vflag);

        if (maxe1 == NULL && maxe2 != NULL)
        {
            newMax = maxe2;
        }
        else if (maxe1 != NULL && maxe2 == NULL)
        {
            newMax = maxe1;
        }
        else if (maxe1 != NULL && maxe2 != NULL)
        {
            if (maxe1->weight > maxe2->weight)
            {
                newMax = maxe1;
            }
            else
            {
                newMax = maxe2;
            }
        }
        else if (maxe1 == NULL && maxe2 == NULL)
        {
            //search from fist elem
            if (globalMax != NULL)
            {
                //printf("curr max start %d \n", globalMax->start);
                newMax = vgetMax(graph, globalMax->start, vflag);
            }
            if (newMax == NULL)
            {
                newMax = vgetMax(graph, 0, vflag);
            }

            globalMax = newMax;
        }
    }
}

void getReorderedVindex(
    struct DataItem *graph[],
    struct DataItem *reversegraph[],
    struct DataItem *maxEdge,
    int *reorderv,
    int vecNum)
{

    struct DataItem *newMax = maxEdge;
    int count = 0;
    int *vflag = (int *)calloc(vecNum, sizeof(int) + 1);
    while (count < vecNum)
    {
        printf("curr count %d\n", count);
        int s = newMax->start;
        int e = newMax->end;
        //printf("curr s %d e %d\n", s, e);
        if (vflag[s] == 0)
        {
            //the reordered vector
            reorderv[count] = s;
            count++;
            vflag[s] = 1;
            //printf("add index %d into vector \n", s);
        }
        if (vflag[e] == 0)
        {
            reorderv[count] = e;
            count++;
            vflag[e] = 1;
            //printf("add index %d into vector\n", e);
        }

        struct DataItem *maxInGraph = search(graph, s, e);
        struct DataItem *maxInrevGraph = search(reversegraph, e, s);
        maxInGraph->flag = 1;
        maxInrevGraph->flag = 1;

        struct DataItem *maxE = NULL;
        //go through both graph and reversegraph to find the max one
        maxE = findMax(s, vecNum, graph, reversegraph);

        //printf("find max(%d %d %d)\n",
        //       maxE->start, maxE->end, maxE->weight);

        if (maxE->weight == 0)
        {
            printf("bug for get reorderedv");
        }
        else
        {
            newMax = maxE;
        }

        //printf("newmax %d %d %d\n",newMax->start,newMax->end,newMax->weight);
    }
}

/*
int main()
{
    vector<struct HeadNode *> HashArray;
    int i;
    for (i = 0; i < 11; i++)
    {

        vector<struct DataItem *> ev;
        struct DataItem *item = (struct DataItem *)malloc(sizeof(struct DataItem));
        item->start = i + 1;
        item->end = i + 1;
        item->weight = 1;
        item->flag = 0;
        item->next = NULL;
        ev.push_back(item);

        HeadNode *node = new (HeadNode);
        node->edgeV = ev;
        node->maxweight = i;

        HashArray.push_back(node);

        printf("size %d maxweight %d\n", HashArray.size(), HashArray[i]->maxweight);

        if (i > 1)
        {
            printf("detect 1 size %d maxweight %d\n", HashArray.size(), HashArray[1]->maxweight);
        }
    }

    int count = HashArray.size();

    for (i = 0; i < count; i++)
    {
        if (HashArray[i]->edgeV[0] != NULL)
        {
            printf("maxweight %d weight %d start %d end %d\n",
                   HashArray[i]->maxweight, HashArray[i]->edgeV[0]->weight, HashArray[i]->edgeV[0]->start, HashArray[i]->edgeV[0]->end);
        }
    }

  

    sort(HashArray.begin(), HashArray.end(), cmp);

    printf("count %d \n", count);
    for (i = 0; i < count; i++)
    {
        if (HashArray[i]->edgeV[0] != NULL && HashArray[i] != NULL)
        {
            printf("maxweight %d weight %d start %d end %d\n",
                   HashArray[i]->maxweight, HashArray[i]->edgeV[0]->weight, HashArray[i]->edgeV[0]->start, HashArray[i]->edgeV[0]->end);
        }
    }
    
}
*/
