#include <bits/stdc++.h>

using namespace std;
class Vertex;
class Edge;
class Face;
class DCEL;

/**
 * @brief Vertex class for denoting vertices on a 2-D plane
 */
class Vertex
{
public:
    /**
     * x coordinate of vertex
     */
    float x;
    /**
     * y coordinate of vertex
     */
    float y;
    /**
     * outgoing edge from the vertex
     */
    Edge *outgoing;
    /**
     * @brief Construct a new Vertex object
     *
     */
    Vertex()
    {
    }
    /**
     * @brief Construct a new Vertex object
     *
     * @param xp x coordinate of vertex
     * @param yp y coordinate of vertex
     */
    Vertex(float xp, float yp)
    {
        x = xp;
        y = yp;
    }
    /**
     * @brief sets the outgoing edge of the vertex
     *
     * @param out outgoing edge of the vertex
     */
    void setedge(Edge *out)
    {
        outgoing = out;
    }
};
/**
 * @brief Edge class for denoting edges in our 2D plane
 *
 */
class Edge
{
public:
    /**
     * @brief denotes the starting vertex of a directed edge
     *
     */
    Vertex *start;
    /**
     * @brief denotes ending vertex of directed edge
     *
     */
    Vertex *end;
    /**
     * @brief denotes twin of edge which is basically
     * an edge joining the same vertices in the opposite direction
     *
     */
    Edge *twin;
    /**
     * @brief denotes the next directed edge in the same polygon
     *
     */
    Edge *next;
    /**
     * @brief denotes the previous directed edge in the same polygon
     *
     */
    Edge *prev;
    /**
     * @brief denotes the face to which the directed edge belongs
     * in this case the face to the right of the edge
     *
     */
    Face *rightface;
    /**
     * @brief Construct a new Edge object
     *
     */
    Edge()
    {
    }
    /**
     * @brief Construct a new Edge object
     *
     * @param st start vertex
     * @param en end vertex
     */
    Edge(Vertex *st, Vertex *en)
    {
        start = st;
        end = en;
    }
    /**
     * @brief sets the twin edge
     *
     * @param tw twin edge
     */
    void settwin(Edge *tw)
    {
        twin = tw;
    }
    /**
     * @brief sets the next edge
     *
     * @param n next edge
     */
    void setnext(Edge *n)
    {
        next = n;
    }
    /**
     * @brief sets the previous edge
     *
     * @param p previous edge
     */
    void setprev(Edge *p)
    {
        prev = p;
    }
    /**
     * @brief sets the right face
     *
     * @param r right face
     */
    void setface(Face *r)
    {
        rightface = r;
    }
};

/**
 * @brief Face class used to define faces or polygons in 2D space
 *
 */
class Face
{
public:
    /**
     * @brief denotes one of the edges belonging to the polygon
     *
     */
    Edge *adjacentedge;
    /**
     * @brief Construct a new Face object
     *
     */
    Face()
    {
    }
    /**
     * @brief Construct a new Face object
     *
     * @param adj adjacent edge
     */
    Face(Edge *adj)
    {
        adjacentedge = adj;
    }
};

/**
 * @brief DCEL class or the doubly connected edge list class
 * it is used to denote all the shapes, edges and vertices present in the 2D plane
 *
 */
class DCEL
{
public:
    /**
     * @brief denotes a vector containing all the edges in the DCEL
     *
     */
    vector<Edge *> edges;
    /**
     * @brief denotes a vector of all the vertices in the DCEL
     *
     */
    vector<Vertex *> vertices;
    /**
     * @brief denotes a vector of all the faces in the DCEL
     *
     */
    vector<Face *> faces;
    bool inside(vector<Edge *> polygon, Vertex *point);
    bool notch(Edge *e1, Edge *e2);
    int countEdges(Face *f);
    void boundingbox(Vertex &minlimits, Vertex &maxlimits, vector<Edge *> vertices);
    vector<Edge *> split(Vertex *v1, Vertex *v2, Face *polygon);
    vector<Edge *> dissolveedge(Edge *dummyedge, Face *polygon);
    void merger();
    void mp1();

};

/**
 * @brief checks if given vertex/point is inside a polygon/face
 * for a given convex polygon, a point lies inside the polygon if
 * the cross product of every directed edge and
 * vector from start vertex of that edge to the point
 * all have the same sign.
 * @param polygon is of type Face
 * @param point is of type Vertex
 * @return true returned when point is inside polygon
 * @return false returned when point is outside polygon
 */
bool DCEL::inside(vector<Edge *> polygon, Vertex *point)
{
    // cross product works for convex polygons
    int n = polygon.size();
    for (int i = 0; i < n; i++)
    {
        float x1 = polygon[i]->start->x;
        float y1 = polygon[i]->start->y;
        float x2 = polygon[i]->end->x;
        float y2 = polygon[i]->end->y;
        float v1i = x2 - x1;
        float v1j = y2 - y1;
        float v2i = point->x - x1;
        float v2j = point->y - y1;
        float crossproduct = (v1i * v2j) - (v2i * v1j);
        if (crossproduct > 0.0)
            return false;
    }
    return true;
}

/**
 * @brief checks if given two edges they form a notch or reflex angle.
 * it works by calculating the cross product of the two directed edges.
 * if it is positive/pointing out of the plane that means it is a reflex angle
 * else it is not
 *
 * @param e1 edge 1
 * @param e2 edge 2
 * @return true if it is a notch
 * @return false if it is not a notch
 */
bool DCEL::notch(Edge *e1, Edge *e2)
{
    float e1i = e1->end->x - e1->start->x;
    float e1j = e1->end->y - e1->start->y;
    float e2i = e2->end->x - e2->start->x;
    float e2j = e2->end->y - e2->start->y;
    float crossproduct = (e1i * e2j) - (e2i * e1j);
    if (crossproduct > 0.0)
        return true;
    return false;
}

/**
 * @brief counts the number of edges given a polygon/face
 * iterates over all edges starting with the adjacent edge of the face until it reaches back the same
 * @param f face/polygon
 * @return int number of edges in polygon
 */
int DCEL::countEdges(Face *f)
{
    Edge *e = f->adjacentedge;
    int n = 1;
    Edge *cur = e->next;
    while (cur != e)
    {
        n++;
        cur = cur->next;
    }
    return n;
}

/**
 * @brief calculates the bounding box given a set of vertices
 * iterates through each of the vertices and calculates the maximum and minimum x and y coordinates
 * to create a box that tightly encapsulates all these vertices
 *
 * @param minlimits max x and max y coordinates are stored here and returned
 * @param maxlimits min x and min y coordinates are stored here and returned
 * @param vertices vector of vertices passed as input
 */
void DCEL::boundingbox(Vertex &minlimits, Vertex &maxlimits, vector<Edge *> vertices)
{
    float xmin = (float)INT_MAX;
    float ymin = (float)INT_MAX;
    float xmax = (float)INT_MIN;
    float ymax = (float)INT_MIN;
    for (auto &v : vertices)
    {
        xmin = min(xmin, v->start->x);
        ymin = min(ymin, v->start->y);
        xmax = max(xmax, v->start->x);
        ymax = max(ymax, v->start->y);
    }
    minlimits.x = xmin;
    minlimits.y = ymin;
    maxlimits.x = xmax;
    maxlimits.y = ymax;
}

/**

* @brief Splits a polygon by adding a new pair of edges (edge and its twin) between two vertices and returns the edges of the resulting split polygon.
* @param v1 The first vertex of the edge to be added.
* @param v2 The second vertex of the edge to be added.
* @param polygon The polygon to be split.
* @param dcel A pointer to the DCEL data structure.
* @return A vector of pointers to the edges of the split polygon.
*/
vector<Edge *> DCEL::split(Vertex *v1, Vertex *v2, Face *polygon)
{
    Edge *dummyedge;
    Face *dummyface;
    Edge *dummyadjacent;
    Edge *dummytwin;
    dummyedge = new Edge(v1, v2);
    dummyface = new Face(dummyedge);
    dummytwin = new Edge(v2, v1);

    Edge *cur = polygon->adjacentedge;
    Edge *l0;
    Edge *lback;

    vector<Edge *> L;

    do
    {
        if (cur->start == v2)
            break;
        cur = cur->next;
    } while (cur != polygon->adjacentedge);

    do
    {
        L.push_back(cur);
        if (cur->end == v1)
            break;
        cur = cur->next;
    } while (cur != polygon->adjacentedge);

    dummyedge->next = L[0];
    dummyedge->prev = L.back();
    dummytwin->next = L.back()->next;
    dummytwin->prev = L[0]->prev;

    L[0]->prev->next = dummytwin;
    L.back()->next->prev = dummytwin;
    L.back()->next = dummyedge;
    L[0]->prev = dummyedge;

    dummyedge->twin = dummytwin;
    dummytwin->twin = dummyedge;
    dummytwin->rightface = polygon;

    for (auto &e : L)
    {
        e->rightface = dummyface;
    }

    dummyedge->rightface = dummyface;
    polygon->adjacentedge = dummytwin->next;
    L.push_back(dummyedge);

    edges.push_back(dummyedge);
    edges.push_back(dummytwin);
    faces.push_back(dummyface);
    return L;
}

/**

*@brief Dissolves an edge in a DCEL.
* This function dissolves an edge in a DCEL data structure by creating a new polygon
* with the same adjacent edges as the old polygon on the left side of the dissolved edge.
* @param dummyedge The dummy edge to be dissolved.
* @param polygon The polygon containing the dummy edge to be dissolved.
* @param dcel The DCEL data structure.
* @return A vector of edges representing the new polygon created after dissolving the edge.
*/

vector<Edge *> DCEL::dissolveedge(Edge *dummyedge, Face *polygon)
{
    vector<Edge *> L;
    Edge *dummytwin = dummyedge->twin;
    Face *dummyface = dummyedge->rightface;
    Edge *cur = dummyedge->next;
    while (cur != dummyedge)
    {
        L.push_back(cur);
        cur = cur->next;
    }
    polygon->adjacentedge = L[0];
    for (auto &e : L)
    {
        e->rightface = polygon;
    }

    L[0]->prev = dummytwin->prev;
    L.back()->next = dummytwin->next;
    L.back()->next->prev = L.back();
    L[0]->prev->next = L[0];

    for (vector<Edge *>::iterator e = edges.begin(); e != edges.end(); e++)
    {
        if (*e == dummyedge || *e == dummytwin)
        {
            edges.erase(e);
            e--;
        }
    }
    for (vector<Face *>::iterator f = faces.begin(); f < faces.end(); f++)
    {
        if (*f == dummyface)
        {
            faces.erase(f);
            // f--;
        }
    }
    return L;
}

/**
 * @brief Merges polygons by deleting edges which do not lead to non-convex polygons.
 * It iterates through all the new edges created and checks if removing it causes a notch
 * @param dcel A pointer to the DCEL data structure.
 */
void DCEL::merger()
{
    int n = vertices.size();
    Edge *e;
    vector<Edge *> partitionedges;
    copy(edges.begin() + 2 * n, edges.end(), back_inserter(partitionedges));
    for (int eidx = 0; eidx < partitionedges.size(); eidx++)
    {
        Edge *e = partitionedges[eidx];
        if ((notch(e->next, e->twin->prev)) && (!notch(e->prev, e->twin->next)))
        {
            dissolveedge(e->twin, e->rightface);
            eidx--;
            for (vector<Edge *>::iterator eit = partitionedges.begin(); eit != partitionedges.end(); eit++)
            {
                if (*eit == e || *eit == e->twin)
                {
                    partitionedges.erase(eit);
                    eit--;
                }
            }
        }
    }
}

/**
 * @brief Rotates a vector of DCEL edges by a specified number of positions.
 * @param dcel_edges The vector of DCEL edges to be rotated.
 * @param index The index from which rotation should be performed.
 * @param num_rotate The number of positions by which the vector should be rotated.
 * @return A vector of DCEL edges rotated by the specified number of positions.
 */
vector<Edge *> rotate_array(vector<Edge *> dcel_edges, int index, int num_rotate)
{
    while (num_rotate--)
    {
        Edge *temp = dcel_edges[index];
        for (int i = index; i > 0; i--)
        {
            dcel_edges[i] = dcel_edges[i - 1];
        }

        dcel_edges[0] = temp;

        Edge *temp2 = dcel_edges.back();
        for (int i = dcel_edges.size() - 1; i > index; i--)
        {
            dcel_edges[i] = dcel_edges[i - 1];
        }

        dcel_edges[index + 1] = temp2;
    }

    return dcel_edges;
}

/**
 * @brief This mp1 algorithm performs the decomposition of a simple polygon into convex polygons as per the MP1 algorithm described in the research paper: "Fernández, J., Cánovas, L., & Pelegrı́n, B. (2000). Algorithms for the decomposition of a polygon into convex polygons. European Journal of Operational Research, 121(2), 330-342"
 * @param dcel A pointer to the DCEL data structure representing the polygonal subdivision.
 * @return void
 */
void DCEL::mp1()
{
    Face *polygon = faces[1];
    vector<Edge *> L;
    while (countEdges(polygon) > 3)
    {
        L.clear();
        L.push_back(polygon->adjacentedge);
        Edge *curedge = L[0]->next;
        Edge *eto0;
        while (curedge != L[0])
        {
            eto0 = new Edge(curedge->end, L[0]->start);
            if (!notch(curedge->prev, curedge))
            {
                if (!notch(curedge, eto0))
                {
                    if (!notch(eto0, L[0]))
                    {
                        L.push_back(curedge);
                        curedge = curedge->next;
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
            else
            {
                break;
            }
        }
        if (curedge == L[0])
        {

            return;
        }

        // backtracking
        while (true)
        {
            if (L.size() == 1)
            {
                polygon->adjacentedge = polygon->adjacentedge->next;
                break;
            }
            if (L.size() == 0)
            {
                polygon->adjacentedge = polygon->adjacentedge->next;
                break;
            }

            L = split(L.back()->end, L[0]->start, polygon);
            Vertex minl, maxl;
            boundingbox(minl, maxl, L);
            Edge *it = polygon->adjacentedge->next;
            bool flag;
            while (it != polygon->adjacentedge->prev)
            {
                flag = false;
                if (it->start->x <= maxl.x && it->start->x >= minl.x)
                {
                    if (it->start->y <= maxl.y && it->start->y >= minl.y)
                    {
                        if (notch(it->prev, it))
                        {
                            if (inside(L, it->start))
                            {

                                L = dissolveedge(L.back(), polygon);
                                L.pop_back();
                                flag = true;
                                break;
                            }
                        }
                    }
                }
                it = it->next;
            }
            // no points inside new polygon
            if (!flag)
            {
                break;
            }
        }
    }
}

/**
 * @brief this is the main function of our program
 * argument passed to the main function is the number of vertices (n).
 * and a file named generate_test2.txt with all the vertices in clockwise order listed - has to be present in the same folder
 */

int main(int argc, char *argv[])
{

    clock_t start, end;
    start = clock();
    int x = atoi(argv[1]);

    for (int rot = 0; rot < x; rot++)
    {

        FILE *fp = fopen("generate_test2.txt", "r");
        if (fp == NULL)
            return 0;

        char *line = NULL;
        size_t len = 0;

        int n;

        fscanf(fp, "%d", &n);

        vector<Vertex *> input_vertices(n);
        vector<Edge *> input_edges;
        for (int i = 0; i < n; i++)
        {
            float xp, yp;
            fscanf(fp, "%f %f", &xp, &yp);
            input_vertices[i] = new Vertex(xp, yp);
        }

        fclose(fp);

        for (int i = 0; i < n; i++)
        {
            input_edges.push_back(new Edge(input_vertices[i], input_vertices[(i + 1) % n]));
        }

        input_edges.push_back(new Edge(input_vertices[0], input_vertices[n - 1]));

        for (int i = n - 1; i > 0; i--)
        {
            input_edges.push_back(new Edge(input_vertices[i], input_vertices[(i - 1 + n) % n]));
        }
        for (int i = 0; i < n; i++)
        {
            input_edges[i]->setnext(input_edges[(i + 1) % n]);
            input_edges[(i + 1) % n]->setprev(input_edges[i]);
            input_edges[i]->settwin(input_edges[2 * n - 1 - i]);
            input_edges[2 * n - 1 - i]->settwin(input_edges[i]);
        }

        input_edges[2 * n - 1]->setnext(input_edges[n]);
        input_edges[n]->setprev(input_edges[2 * n - 1]);

        for (int i = n; i < 2 * n - 1; i++)
        {
            input_edges[i]->setnext(input_edges[i + 1]);
            input_edges[i + 1]->setprev(input_edges[i]);
        }

        vector<Face *> faces;

        vector<Edge *> auxilary_edges = rotate_array(input_edges, n - 1, rot);

        faces.push_back(new Face(auxilary_edges[0]));
        faces.push_back(new Face(auxilary_edges[n]));

        DCEL *dcel = new DCEL();
        dcel->edges = auxilary_edges;
        dcel->faces = faces;
        dcel->vertices = input_vertices;

        swap(dcel->faces[0], dcel->faces[1]);
        dcel->mp1();
        dcel->merger();

        ofstream myfile;
        myfile.open("plot_answer.txt");

        myfile << dcel->faces.size() << endl;
        for (auto &a : dcel->faces)
        {
            Edge *starte = a->adjacentedge;
            Edge *cure = starte->next;
            myfile << dcel->countEdges(a) << endl;
            myfile << starte->start->x << " " << starte->start->y << "\n";
            while (cure != starte)
            {
                myfile << cure->start->x << " " << cure->start->y << "\n";
                cure = cure->next;
            }
            myfile << cure->start->x << " " << cure->start->y << "\n";
        }

        myfile.close();

        string for_call = "python3 to_test.py ";
        for_call += to_string(n);
        for_call += " ";
        for_call += to_string(rot);

        system(for_call.c_str());
    }
    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    ofstream myfile1;
    myfile1.open("times.txt", ios_base::app);
    myfile1 << fixed << time_taken << setprecision(5) << " ";
    myfile1.close();
    return 0;
}
