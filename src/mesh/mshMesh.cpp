#include "mesh.h"

using namespace std;

/**
  Read ANSYS ascii mesh (.msh format)
 */
void Mesh::MeshObject::readMshMesh(std::istream& is) {
    char symbol,c;
    int id,ND,zone,findex,lindex,
        type,bctype,ftype,etype,
        node_start = 0,facet_start = 0;
    map<int,string> bnames;

    is >> hex;

    /*read id*/
    while((c = Util::nextc(is)) != 0) {
        int braces = 1;
        is >> symbol >> id;
        switch(id) {
            case 0x0: 
                do{ is >> c; } while(c != ')');
                break;
            case 0x2:
                is >> ND >> symbol;
                break;
            case 0x10:
                is >> symbol >> zone;
                is >> findex >> lindex;
                is >> type >> ND;
                is >> symbol >> symbol;
                if(zone != 0) {
                    Vertex v;
                    for(int i = findex;i <= lindex;i++) {
                        is >> v;
                        mVertices.push_back(v);
                    }
                    is >> symbol >> symbol;
                } else {
                    node_start = findex;
                }
                break;
            case 0x12:
                is >> symbol >> zone;
                is >> findex >> lindex;
                is >> type;
                if((c = Util::nextc(is)) == ')');
                else is >> etype;
                is >> symbol >> symbol;

                while(symbol == '(') {
                    do{ is >> c; } while(c != ')');
                    is >> symbol;
                }
                if(zone == 0) {
                    mCells.resize(lindex);
                }
                break;
            case 0x13:
                is >> symbol >> zone;
                is >> findex >> lindex;
                is >> bctype;
                if((c = Util::nextc(is)) == ')') ftype = bctype;
                else is >> ftype;
                is >> symbol >> symbol;

                if(zone != 0) {
                    std::stringstream name;
                    name << "zone" << zone;
                    IntVector& gB = mBoundaries[name.str().c_str()];

                    Facet f;
                    int n,c0,c1,k;
                    for(int i = findex;i <= lindex;i++) {
                        f.clear();
                        is >> n;
                        for(int j = 0;j < n;j++) {
                            is >> k;
                            f.push_back(k - node_start);
                        }
                        mFacets.push_back(f);
                        gB.push_back(i - facet_start);

                        is >> c0 >> c1;
                        if(c0 == 0) {
                            mFOC.push_back(Constants::MAX_INT);
                        } else {
                            mFOC.push_back(c0 - 1);
                        }
                        if(c1 == 0) {
                            mFNC.push_back(Constants::MAX_INT);
                        } else {
                            mFNC.push_back(c1 - 1);
                        }
                    }
                    is >> symbol >> symbol;
                } else {
                    facet_start = findex;
                }
                break;
            case 0x39:
            case 0x45:
                is >> symbol >> dec >> zone >> hex;
                {
                    string str;
                    char buf[64];
                    is >> str;
                    int i = 0;
                    do{ is >> c; } while(c != ')' && (buf[i++] = c));
                    buf[i] = 0;
                    bnames[zone] = buf;
                }
            default:
                while((c = Util::nextc(is))) {
                    is >> c;
                    if(c == '(') braces++;
                    else if(c == ')') {
                        braces--;
                        if(!braces) break;
                    }
                }
                break;
        }
    }
    /*rename*/
    forEachIt(bnames,it) {
        std::stringstream name;
        name << "zone" << dec << it->first;
        auto it1 = mBoundaries.find(name.str().c_str());
        if(it1 != mBoundaries.end()) {
            mBoundaries[it->second] = it1->second; 
            mBoundaries.erase(it1);
        }
    }
    /*add cells*/
    Int co,cn;
    forEach(mFacets,i) {
        co = mFOC[i];
        cn = mFNC[i];
        if(co != Constants::MAX_INT)
            mCells[co].push_back(i);
        if(cn != Constants::MAX_INT)
            mCells[cn].push_back(i);
    }

    is >> dec;
}
/**
  Write ANSYS ascii mesh (.msh format)
 */
void Mesh::MeshObject::writeMshMesh(std::ostream& os) {
    os << hex;
    os << "(0 \"ASCII msh file\")" << endl << endl;
    os << "(0 \"Dimension:\")" << endl;
    os << "(2 3)" << endl << endl;

    //vertices
    os << "(0 \"Vertices:\")" << endl;
    os << "(10 (0 1 " << mVertices.size() << " 0 3))" << endl << endl;
    os << "(10 (1 1 " << mVertices.size() << " 1 3)" << endl;
    os << "(" << endl;
    os.precision(10);
    forEach(mVertices,i)
        os << scientific << mVertices[i] << endl;
    os << "))" << endl << endl;

    //facets
    os << "(0 \"Facets:\")" << endl;
    os << "(13 (0 1 " << mFacets.size() << " 0 0))" << endl << endl;

    Int zone = 1;
    Int start = 1;

    //internal
    Int nInternal = mFacets.size() - (mCells.size() - mBCS);
    os << "(0 \"Internal faces:\")" << endl;
    os << "(39 (" << dec << zone << hex << " interior " << "interior-1" 
        << ")())" << endl;
    os << "(13 (" << zone << " " << start << " " 
        << nInternal << " 2 0)" << endl;
    zone++;
    start += nInternal;
    os << "(" << endl;
    forEach(mFacets,f) {
        if(mFNC[f] >= mBCS)
            continue;
        Facet& mf = mFacets[f];
        os << mf.size() << " ";
        forEach(mf,j)
            os << mf[j] + 1 << " ";
        if(mReversed[f])
            os << mFOC[f] + 1 << " " << mFNC[f] + 1 << endl;
        else
            os << mFNC[f] + 1 << " " << mFOC[f] + 1 << endl;  
    }
    os << "))" << endl << endl;

    //boundary
    forEachIt(mBoundaries,it) {
        const IntVector& fvec = it->second; 
        os << "(0 \"" << it->first << "\")" << endl;

        string bname = "pressure-inlet";
        Int bid = 4;
        if(it->first.find("WALL") != std::string::npos) {
            bname = "wall";
            bid = 3;
        }
        os << "(39 (" << dec << zone << hex << " " << bname << " " << it->first 
            << ")())" << endl;
        os << "(13 (" << zone << " " << start << " " 
            << start + fvec.size() - 1 << " " << bid << " 0)" << endl;
        zone++;
        start += fvec.size();

        os << "(" << endl;
        forEach(fvec,i) {
            Int f = fvec[i];
            Facet& mf = mFacets[f];
            os << mf.size() << " ";
            forEach(mf,j)
                os << mf[j] + 1 << " ";
            if(mReversed[f])
                os << mFOC[f] + 1 << " 0" << endl;
            else
                os << "0 " << mFOC[f] + 1 << endl;
        }
        os << "))" << endl << endl;
    }

    //cells
    os << "(0 \"Cells:\")" << endl;
    os << "(12 (0 1 " << mBCS << " 0 0))" << endl;
    os << "(12 (1 1 " << mBCS << " 1 0)(" << endl;
    for(Int i = 0;i < mBCS;i++)
        os << "4 ";
    os << endl << ")())" << endl;

    os << dec;
}
