#include "hexMesh.h"

using namespace std;

/**
  Boundary information
 */
struct Bdry {
    string name;
    IntVector index;
};

/**
  Mesh generator application
 */
int main(int argc,char* argv[]) {
    using namespace Mesh;
    using namespace Util;
    Vertices key_points;
    vector<Bdry> patches;
    MergeObject bMerge;
    string str;
    string default_name;
    string i_file_name(argv[1]);
    string o_file_name;
    bool Import = false;
    char c;

    /*command line arguments*/
    for(int i = 1;i < argc;i++) {
        if(!strcmp(argv[i],"-i")) {
            Import = true;
            i++;
            i_file_name = argv[i];
        } else if(!strcmp(argv[i],"-o")) {
            i++;
            o_file_name = argv[i];
        } else if(!strcmp(argv[i],"-h")) {
            std::cout << "Usage:\n"
                << "  ./mesh <inputfile> <Options>\n"
                << "Options:\n"
                << "  -o     --  Output file name including extension\n"
                << "  -i     --  Input file name including extension\n"
                << "  -h     --  Display this message\n\n";
            return 0;
        } 
    }

    /** Import or Generate grid */
    if(Import) {
        str = i_file_name;
        std::string ext = str.substr(str.size() - 4);
        if(ext == ".txt") {
            std::ifstream is(str);
            gMesh.readTextMesh(is);
        } else if(ext == ".bin") {
            Util::ifstream_bin is(str);
            gMesh.readTextMesh(is);
        } else if(ext == ".bin") {
            std::ifstream is(str);
            gMesh.readMshMesh(is);
        }
        gBCS = gCells.size();
    } else {
        /*input stream*/
        ifstream input(i_file_name);

        /*read key points*/
        if(Util::nextc(input))
            input >> key_points;

        while((c = Util::nextc(input)) != 0) {
            char symbol;
            if(isdigit(c)) {
                /*read indices to key points*/
                IntVector index;
                input >> index;

                Vertices corners(index.size(),Vector(0));
                forEach(corners,i)
                    corners[i] = key_points[index[i]];

                IntVector n;
                Int type;
                vector<Scalar> s(12,Scalar(1));
                vector<Int> t(12);

                input >> str;
                if(!compare(str,"linear")) {
                    input >> n;
                    type = LINEAR;
                    t.assign(12,type);
                } else { 
                    if(!compare(str,"geometric")) type = GEOMETRIC;
                    else if(!compare(str,"wall")) type = WALL;
                    else if(!compare(str,"mixed")) type = MIXED;
                    else return 1;

                    //read divisions
                    vector<Scalar> ts(s);
                    vector<Int> tt(t);

                    Int sz;
                    input >> n;
                    input >> sz >> symbol;

                    if(type == MIXED) {
                        for(Int i = 0;i < sz ;i++) {
                            input >> symbol;
                            switch(symbol) {
                                case 'l':
                                case 'L':
                                    type = LINEAR;
                                    break;
                                case 'g':
                                case 'G':
                                    type = GEOMETRIC;
                                    break;
                                case 'w':
                                case 'W':
                                    type = WALL;
                                    break;
                            }
                            tt[i] = type;
                            input >> ts[i];
                        }
                    } else {
                        for(Int i = 0;i < sz ;i++) 
                            input >> ts[i];
                        tt.assign(12,type);
                    }

                    input >> symbol;

                    //assign to each side
                    Int r = 12 / sz;
                    for(Int i = 0;i < sz;i++) {
                        for(Int j = 0;j < r;j++) {
                            if(i * r + j < 12) {
                                s[i * r + j] = ts[i];
                                t[i * r + j] = tt[i];
                            }
                        }
                    }
                }

                //curved edges
                static const int sides[12][2] = {
                    {0,1}, {3,2}, {7,6}, {4,5},
                    {0,3}, {1,2}, {5,6}, {4,7},
                    {0,4}, {1,5}, {2,6}, {3,7}
                };
                vector<Edge> edges(12);
                Sphere sphere;
                bool has_sphere = false;
                for(Int i = 0;i < 12;i++) {
                    edges[i].v[0] = corners[sides[i][0]];
                    edges[i].v[1] = corners[sides[i][1]];
                }

                if((c = Util::nextc(input)) && (c == 'e' || c == 's')) {
                    Int sz,side,key;
                    input >> str;
                    if(!compare(str,"edges")) {
                        input >> sz >> symbol;
                        for(Int i = 0;i < sz;i++) {
                            input >> str >> side >> key;    
                            Edge& e = edges[side];
                            e.v[2] = key_points[key];
                            if(!compare(str,"arc")) {
                                e.type = ARC;
                            } else if(!compare(str,"cosine")) {
                                e.type = COSINE;
                            } else if(!compare(str,"quad")) {
                                e.type = QUAD;
                            } else {
                                e.type = NONE;
                            }
                        }
                        input >> symbol;
                    } else if(!compare(str,"sphere")) {
                        Scalar radiusi, radiuso;
                        input >> key >> radiusi >> radiuso;
                        sphere.radiusi = radiusi;
                        sphere.radiuso = radiuso;
                        sphere.center = key_points[key];
                        has_sphere = true;
                    } else {
                        Bdry b;
                        b.name = str;
                        while((c = Util::nextc(input)) && isdigit(c)) {
                            input >> b.index;
                            patches.push_back(b);
                        }
                    }
                }

                //generate mesh
                MeshObject mo;
                hexMesh(&n[0],&s[0],&t[0],&corners[0],&edges[0],mo,
                        has_sphere ? &sphere : 0);
                merge(gMesh,bMerge,mo);
            } else {
                /*read boundaries*/
                Bdry b;
                input >> b.name;
                if(b.name == "default") {
                    input >> default_name;
                } else {
                    while((c = Util::nextc(input)) && isdigit(c)) {
                        input >> b.index;
                        patches.push_back(b);
                    }
                }
            }
        }
        /*merge boundary & internals*/
        merge(gMesh,bMerge);

        /*boundaries*/
        forEach(patches,i) {
            IntVector list;
            IntVector& b = patches[i].index;

            //pick faces on same plane as patch
            Vector N = (key_points[b[1]] - key_points[b[0]])
                     ^ (key_points[b[2]] - key_points[b[0]]);
            N /= mag(N);
            forEach(gMesh.mPatches,j) {
                Patch& p = gMesh.mPatches[j];
                Vector H = (p.C - key_points[b[0]]);
                Scalar d1 = mag(N ^ p.N);
                Scalar d2 = sqrt(mag(N & H));
                if(d1 <= 10e-4 && d2 <= 10e-4) {
                    for(Int k = p.from;k < p.to;k++)
                        list.push_back(k);
                }
            }

            if(!list.empty()) {
                IntVector& gB = gBoundaries[patches[i].name.c_str()];
                auto it = find(gB.begin(),gB.end(),list[0]);
                if(it == gB.end()) {
                    forEach(list,j)
                        gB.push_back(list[j]);
                }
            }
        }
        /*default specified*/
        if(!default_name.empty()) {
            IntVector faceInB;
            faceInB.assign(gFacets.size(),0);
            forEachIt(gBoundaries,it) {
                IntVector& gB = it->second; 
                forEach(gB,j)
                    faceInB[gB[j]] = 1;
            }

            IntVector& gB = gBoundaries[default_name.c_str()];
            forEachS(gFacets,i,gMesh.mNF) {
                if(!faceInB[i])
                    gB.push_back(i);
            }
        }
    }

    /*write the grid in requested format*/
    std::string ext = o_file_name.substr(o_file_name.size() - 4);
    if(ext == ".txt") {
        ofstream os(o_file_name);
        os << gMesh;
    } else if(ext == ".bin") {
        Util::ofstream_bin os(o_file_name);
        os << gMesh;
    } else if(ext == ".msh") {
        Mesh::gMesh.addBoundaryCells();
        Mesh::gMesh.calcGeometry(0);
        ofstream os(o_file_name);
        gMesh.writeMshMesh(os);
    } else {
        cerr << "Use proper extension for output file name." << endl;
        return 1;
    }
    return 0;
}
