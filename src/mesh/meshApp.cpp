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
    Vertices corners;
    vector<Bdry> patches;
    MergeObject bMerge;
    string str;
    string default_name;
    char* i_file_name = argv[1];
    char* e_file_name = 0;
    bool Import = false;
    bool Export = false;
    char c;

    /*command line arguments*/
    for(int i = 1;i < argc;i++) {
        if(!strcmp(argv[i],"-i")) {
            i++;
            Import = true;
            i_file_name = argv[i];
        } else if(!strcmp(argv[i],"-o")) {
            i++;
            Export = true;
            e_file_name = argv[i];
        } else if(!strcmp(argv[i],"-h")) {
            std::cout << "Usage:\n"
                << "  ./mesh <inputfile> <Options>\n"
                << "Options:\n"
                << "  -i     --  Import from Fluent .msh file\n"
                << "  -o     --  Export to Fluent .msh file format\n"
                << "  -h     --  Display this message\n\n";
            return 0;
        } 
    }

    /*export to msh file format*/
    if(Export) {
        ofstream output(e_file_name);
        if(Import) str = i_file_name;
        else str = "grid";
        Mesh::gMeshName = str;
        Mesh::gMesh.readMesh();
        Mesh::gMesh.addBoundaryCells();
        Mesh::gMesh.calcGeometry();

        output << hex;
        gMesh.writeMshMesh(output);
        output << dec;
        return 0;
    }

    /*input stream*/
    ifstream input(i_file_name);

    /*import*/
    if(Import) {
        input >> hex;
        gMesh.readMshMesh(input);
        input >> dec;

        gMesh.writeMesh(cout);
        return 0;
    }

    /*read key points*/
    if(Util::nextc(input))
        input >> corners;

    while((c = Util::nextc(input)) != 0) {
        char symbol;
        if(isdigit(c)) {
            /*read indices to corners*/
            IntVector index;
            input >> index;

            Vertices v(index.size(),Vector(0));
            forEach(v,i)
                v[i] = corners[index[i]];

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
            for(Int i = 0;i < 12;i++) {
                edges[i].v[0] = v[sides[i][0]];
                edges[i].v[1] = v[sides[i][1]];
            }

            if((c = Util::nextc(input)) && (c == 'e')) {
                Int sz,side,key;
                input >> str;
                if(!compare(str,"edges")) {
                    input >> sz >> symbol;
                    for(Int i = 0;i < sz;i++) {
                        input >> str >> side >> key;    
                        Edge& e = edges[side];
                        e.v[2] = corners[key];
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
            hexMesh(&n[0],&s[0],&t[0],&v[0],&edges[0],mo);
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
        Vector N = (corners[b[1]] - corners[b[0]]) ^ (corners[b[2]] - corners[b[0]]);
        N /= mag(N);

        forEach(gMesh.mPatches,j) {
            Patch& p = gMesh.mPatches[j];
            Vector H = (p.C - corners[b[0]]);
            Scalar d1 = mag(N ^ p.N);
            Scalar d2 = sqrt(mag(N & H));
            // if(Mesh::pointInPolygon(corners,b,p.C))
            if(d1 <= 10e-4 && d2 <= 10e-4) {
                for(Int k = p.from;k < p.to;k++)
                    list.push_back(k);
            }
        }

        if(!list.empty()) {
            IntVector& gB = gBoundaries[patches[i].name.c_str()];
            IntVector::iterator it = find(gB.begin(),gB.end(),list[0]);
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
        forEachIt(Boundaries,gBoundaries,it) {
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
    /*write it*/
    gMesh.writeMesh(cout);
    return 0;
}
