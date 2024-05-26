#ifndef HMSH_HMSH_H
#define HMSH_HMSH_H
#include "typedefs.h"
#include "topology.h"
#include <queue>

namespace pddg {
struct Hmsh;
struct Vert;
struct Edge;
struct Face;
struct Half;
struct AdjVH;
struct AdjFH;
struct AdjLH;
template<typename N> struct AdjIter;

struct Elem { int id; Hmsh* m; };

struct Face : Elem {
    [[nodiscard]] inline Half half() const;
    [[nodiscard]] inline Row3d basisX() const;
    [[nodiscard]] inline Row3d basisY() const;
    [[nodiscard]] inline Row3d normal() const;
    [[nodiscard]] inline Row3d center() const;
    [[nodiscard]] inline double area() const;
    AdjIter<AdjFH> adjHalfs(bool ccw = true);
    AdjIter<AdjFH> adjHalfs(Half h, bool ccw = true);
};

struct Edge : Elem {
    [[nodiscard]] inline Half half() const;
    [[nodiscard]] inline Vert v0() const;
    [[nodiscard]] inline Vert v1() const;
    [[nodiscard]] inline Face f0() const;
    [[nodiscard]] inline Face f1() const;
    [[nodiscard]] inline double len() const;
    [[nodiscard]] inline double cot() const;
    [[nodiscard]] inline bool isBoundary() const;
};

struct Vert : Elem {
    [[nodiscard]] inline Half half() const;
    [[nodiscard]] inline Row3d pos() const;
    [[nodiscard]] inline Row3d basisX() const;
    [[nodiscard]] inline Row3d basisY() const;
    [[nodiscard]] inline Row3d normal() const;
    [[nodiscard]] inline bool isBoundary() const;
    [[nodiscard]] inline double baryArea() const;
    [[nodiscard]] inline double circArea() const;
    AdjIter<AdjVH> adjHalfs(bool ccw = true);
    AdjIter<AdjVH> adjHalfs(Half h, bool ccw = true);
};

struct Loop : Elem {
    [[nodiscard]] inline Half half() const;
    AdjIter<AdjLH> adjHalfs(bool ccw = true);
};

struct Half : Elem {
    [[nodiscard]] inline Half next() const;
    [[nodiscard]] inline Half prev() const;
    [[nodiscard]] inline Half twin() const;
    [[nodiscard]] inline Vert tail() const;
    [[nodiscard]] inline Vert head() const;
    [[nodiscard]] inline Edge edge() const;
    [[nodiscard]] inline Face face() const;
    [[nodiscard]] inline double len() const;
    [[nodiscard]] inline double cot() const;
    [[nodiscard]] inline double varg() const;
    [[nodiscard]] inline double farg() const;
    [[nodiscard]] inline double darg() const;
    [[nodiscard]] inline bool isBoundary()  const;
    [[nodiscard]] inline bool isCanonical() const;
    [[nodiscard]] inline Eigen::RowVector3d vec() const;
};

enum class MeshType { GraphOnly, FullFeature };

struct Hmsh {
    int nV;
    int nF;
    int nE;
    int nH;
    int nL;
    int nEularChars;

    std::vector<int> vert2half;
    std::vector<int> edge2half;
    std::vector<int> face2half;
    std::vector<int> loop2half;
    std::vector<int> next;
    std::vector<int> prev;
    std::vector<int> twin;
    std::vector<int> tail;
    std::vector<int> head;
    std::vector<int> edge;
    std::vector<int> face;
    std::vector<bool> isBV;
    std::vector<Half> halfs;
    std::vector<Vert> verts;
    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::vector<Loop> loops;
    MatXd pos;
    MatXi idx;
    MatXi edge2vert;
    MatXi edge2face;
    MatXi face2edge;
    MatXd vertBasisX;
    MatXd vertBasisY;
    MatXd vertNormal;
    MatXd faceBasisX;
    MatXd faceBasisY;
    MatXd faceNormal;
    VecXd faceArea;
    VecXd baryDualArea;
    VecXd circDualArea;
    VecXd halfCotan;
    VecXd edgeCotan;
    MatXd baryCenter;
    VecXd angleDefect;
    VecXd dihedralArg;
    VecXd heArgOnVert;
    VecXd heArgOnFace;
    VecXd scalarMeanCurvature;
    MatXd principalCurvatures;

    inline Hmsh(const MatXd& V, const MatXi& F, MeshType type = MeshType::FullFeature);

    [[nodiscard]] inline SprsD d0()  const;
    [[nodiscard]] inline SprsD d1()  const;
    [[nodiscard]] inline SprsD h0()  const { return static_cast<SprsD>(baryDualArea.asDiagonal()); }
    [[nodiscard]] inline SprsD h0i() const { return static_cast<SprsD>(baryDualArea.asDiagonal().inverse()); }
    [[nodiscard]] inline SprsD h1()  const { return static_cast<SprsD>(edgeCotan.asDiagonal()); }
    [[nodiscard]] inline SprsD h1i() const { return static_cast<SprsD>(edgeCotan.asDiagonal().inverse()); }
    [[nodiscard]] inline SprsD h2()  const;
    [[nodiscard]] inline SprsD h2i() const;
};

Half Half::next() const { return {m->next[id], m}; }
Half Half::prev() const { return {m->prev[id], m}; }
Half Half::twin() const { return {m->twin[id], m}; }
Vert Half::tail() const { return {m->tail[id], m}; }
Vert Half::head() const { return {m->head[id], m}; }
Edge Half::edge() const { return {m->edge[id], m}; }
Face Half::face() const { return {m->face[id], m}; }
Half Vert::half() const { return {m->vert2half[id], m}; }
Half Edge::half() const { return {m->edge2half[id], m}; }
Half Face::half() const { return {m->face2half[id], m}; }
Half Loop::half() const { return {m->loop2half[id], m}; }
Vert Edge::v0() const { return {m->edge2vert(id, 0), m}; }
Vert Edge::v1() const { return {m->edge2vert(id, 1), m}; }
Face Edge::f0() const { return {m->edge2face(id, 0), m}; }
Face Edge::f1() const { return {m->edge2face(id, 1), m}; }
bool Half::isCanonical() const { return edge().half().id == id; }
bool Half::isBoundary()  const { return face().id == -1; }
bool Vert::isBoundary()  const { return m->isBV[id]; }
bool Edge::isBoundary()  const { return half().isBoundary() || half().twin().isBoundary(); }
double Edge::len() const { return Half{m->edge2half[id], m}.len(); }
double Half::len() const { return vec().norm(); }
double Half::cot() const { return m->halfCotan[id]; }
double Edge::cot() const { return m->edgeCotan[id]; }
double Half::varg() const { return m->heArgOnVert[id]; }
double Half::farg() const { return m->heArgOnFace[id]; }
double Half::darg() const { return m->dihedralArg[id]; }
double Face::area() const { return m->faceArea[id]; }
double Vert::baryArea() const { return m->baryDualArea[id]; }
double Vert::circArea() const { return m->circDualArea[id]; }
Row3d Half::vec() const { return head().pos() - tail().pos(); }
Row3d Vert::pos() const { return m->pos.row(id); }
Row3d Vert::basisX() const { return m->vertBasisX.row(id); }
Row3d Vert::basisY() const { return m->vertBasisY.row(id); }
Row3d Vert::normal() const { return m->vertNormal.row(id); }
Row3d Face::basisX() const { return m->faceBasisX.row(id); }
Row3d Face::basisY() const { return m->faceBasisY.row(id); }
Row3d Face::normal() const { return m->faceNormal.row(id); }
Row3d Face::center() const { return m->baryCenter.row(id); }

struct AdjBase {
    Half h;
    bool bgn;
    bool ccw;
    AdjBase(Half h, bool ccw) : h(h), bgn(false), ccw(ccw) { }
    Half operator*() const { return h; }
    bool operator!=(const AdjBase& a) const { return !bgn || h.id != a.h.id; }
};

struct AdjVH : AdjBase {
    AdjVH(Hmsh* m, const int hid, bool ccw): AdjBase(Half{hid, m}, ccw) { }
    AdjVH& operator++() { h = ccw ? h.prev().twin() : h.twin().next(); bgn = true; return *this; }
};

struct AdjFH: AdjBase {
    AdjFH(Hmsh* m, const int hid, bool ccw): AdjBase(Half{hid, m}, ccw) { }
    AdjFH& operator++() { h = ccw ? h.next() : h.prev(); bgn = true; return *this; }
};

struct AdjLH : AdjBase {
    AdjLH(Hmsh* m, const int hid, bool ccw): AdjBase(Half{hid, m}, ccw) { }
    AdjLH& operator++() { h = ccw ? h.next() : h.prev(); bgn = true; return *this; }
};

template<typename N>
struct AdjIter {
    AdjIter(Hmsh* m, int iX, bool ccw) : bgnNav(N(m, iX, ccw)), endNav(N(m, iX, ccw)) { }
    N begin() { return bgnNav; }
    N end()   { return endNav; }
    N bgnNav, endNav;
};

inline AdjIter<AdjVH> Vert::adjHalfs(bool ccw) { return {m, m->vert2half[id], ccw}; }
inline AdjIter<AdjFH> Face::adjHalfs(bool ccw) { return {m, m->face2half[id], ccw}; }
inline AdjIter<AdjLH> Loop::adjHalfs(bool ccw) { return {m, m->loop2half[id], ccw}; }
inline AdjIter<AdjVH> Vert::adjHalfs(Half h, bool ccw) { assert(h.tail().id == id); return {m, h.id, ccw}; }
inline AdjIter<AdjFH> Face::adjHalfs(Half h, bool ccw) { assert(h.face().id == id); return {m, h.id, ccw}; }


Hmsh::Hmsh(const MatXd& V, const MatXi& F, MeshType type) {
    etopo(F, edge2vert, face2edge, edge2face);
    pos = V;
    idx = F;
    nV = static_cast<int>(pos.rows());
    nF = static_cast<int>(idx.rows());
    nE = static_cast<int>(edge2vert.rows());
    nH = static_cast<int>(edge2vert.rows() * 2);
    nEularChars = nV - nE + nF;
    vert2half.assign(nV, -1);
    edge2half.assign(nE, -1);
    face2half.assign(nF, -1);
    head.assign(nH, -1);
    tail.assign(nH, -1);
    edge.assign(nH, -1);
    face.assign(nH, -1);
    next.assign(nH, -1);
    prev.assign(nH, -1);
    twin.assign(nH, -1);

    int nP = static_cast<int>(idx.cols());
    auto pair = [&](int a, int b) { twin[a] = b; twin[b] = a; };

    /*--- setup topology ---*/
    MatXi EFi, EH, FH;
    VecXi VH, HV, HE, HF, nextH, prevH, twinH;
    VecXi D = VecXi::Constant(idx.rows(), 3);
    EFi = MatXi::Constant(edge2face.rows(), 2, -1); // number of an edge inside the face
    for (int i = 0; i < edge2face.rows(); i++) {
    for (int k = 0; k < 2; k++) {
        if (edge2face(i, k) == -1) continue;
        for (int j = 0; j < 3; j++) if (face2edge(edge2face(i, k), j) == i) EFi(i, k) = j;
    }}
    dcel(D, idx, edge2vert, edge2face, EFi, VH, EH, FH, HV, HE, HF, nextH, prevH, twinH);

    for (int iV = 0; iV < nV; ++iV) vert2half[iV] = VH[iV];
    for (int iE = 0; iE < nE; ++iE) edge2half[iE] = EH.row(iE)[0];
    for (int iF = 0; iF < nF; ++iF) face2half[iF] = FH.row(iF)[0];
    for (int iH = 0; iH < HV.rows(); ++iH) {
        head[iH] = HV[nextH[iH]];
        tail[iH] = HV[iH];
        edge[iH] = HE[iH];
        face[iH] = HF[iH];
        next[iH] = nextH[iH];
        prev[iH] = prevH[iH];
        if (twinH[iH] != -1) twin[iH] = twinH[iH];
    }

    /*--- setup boundaries ---*/
    for (int iH = 0, iB = nF * nP; iH < nF * nP; ++iH) {
        if (twin[iH] != -1)
            continue;
        int jH = iH;
        int jB = iB;
        while (true) {
            tail[iB] = head[jH];
            head[iB] = tail[jH];
            edge[iB] = edge[jH];
            prev[iB] = iB - 1;
            next[iB] = iB + 1;
            pair(jH, iB);

            jH = prev[jH];
            while (twin[jH] != -1) {
                if (jH == iH) goto loop_done;
                jH = prev[twin[jH]];
            }
            iB++;
        }
        loop_done:
        prev[jB] = iB;
        next[iB] = jB;
        loop2half.push_back(jB);
        iB++;
    }
    nL = static_cast<int>(loop2half.size());

    /*--- setup edge->half on boundary ---*/
    for (int iE = 0; iE < nE; ++iE) if(edge2half[iE] == -1) edge2half[iE] = twin[EH.row(iE)[1]];

    for (int iL = 0; iL < nL; ++iL) loops.emplace_back(Loop{iL, this});
    for (int iF = 0; iF < nF; ++iF) faces.emplace_back(Face{iF, this});
    for (int iV = 0; iV < nV; ++iV) verts.emplace_back(Vert{iV, this});
    for (int iE = 0; iE < nE; ++iE) edges.emplace_back(Edge{iE, this});
    for (int iH = 0; iH < nH; ++iH) halfs.emplace_back(Half{iH, this});

    if (type == MeshType::FullFeature) {
        vertNormal.resize(nV, 3);
        vertBasisX.resize(nV, 3);
        vertBasisY.resize(nV, 3);
        faceBasisX.resize(nF, 3);
        faceBasisY.resize(nF, 3);
        faceNormal.resize(nF, 3);
        faceArea.resize(nF);
        halfCotan.resize(nH);
        edgeCotan.resize(nE);
        angleDefect.resize(nV);
        dihedralArg.resize(nH);
        heArgOnVert.resize(nH);
        heArgOnFace.resize(nH);
        baryCenter.resize(nF, 3);
        baryDualArea.resize(nV);
        circDualArea.resize(nV);
        scalarMeanCurvature.resize(nV);
        principalCurvatures.resize(nV, 2);

        for (Face f: faces) {
            Row3d x = f.half().vec();
            Row3d t = f.half().prev().vec() * -1;
            Row3d n = x.cross(t);
            Row3d p = Row3d::Zero();

            /*--- setup  basis of face ---*/
            faceBasisX.row(f.id) = x.normalized();
            faceBasisY.row(f.id) = -x.cross(n).normalized();
            faceNormal.row(f.id) = n.normalized();

            /*--- setup barycenter of face ---*/
            for (Half h: f.adjHalfs()) p += h.tail().pos();
            baryCenter.row(f.id) = p / nP;

            /*--- setup area of face ---*/
            faceArea[f.id] = n.norm() * 0.5;
        }

        /*--- setup halfedge angle on face coordinate ---*/
        for (Face f: faces) {
            double sum = 0;
            for (Half h: f.adjHalfs()) {
                Row3d v1 = h.vec().normalized();
                Row3d v2 = h.prev().vec().normalized();
                if ((v1 - f.basisX()).norm() > 1e-10) sum += acos(v1.dot(v2));
                heArgOnFace[h.id] = sum;
            }
        }

        /*--- setup vertex orthogonal coordinate ---*/
        for (int iF = 0; iF < nF; ++iF) {
            for (int iP = 0; iP < nP; iP++) {
                vertNormal.row(idx(iF, iP)).array()
                        += faceNormal.row(iF).array() * faceArea[iF];
            }}
        vertNormal.rowwise().normalize();

        for (Vert v: verts) {
            Row3d v_ = pos.row(v.half().head().id) - pos.row(v.id);
            Row3d n = vertNormal.row(v.id);
            Row3d x = (v_ - v_.dot(n) * n).normalized();
            vertBasisX.row(v.id) = x;
            vertBasisY.row(v.id) = n.cross(x);
        }

        /*--- setup dihedral angle for each halfedge ---*/
        for (Half h: halfs) {
            if (h.edge().isBoundary()) { dihedralArg[h.id] = 0; continue; }
            Row3d n1 = faceNormal.row(h.face().id);
            Row3d n2 = faceNormal.row(h.twin().face().id);
            Row3d v = h.vec() / h.len();
            Row3d c = n1.cross(n2);
            dihedralArg[h.id] = atan2(v.dot(c), n1.dot(n2));
        }

        /*--- setup halfedge cotan and edge cotan ---*/
        for (Half h: halfs) {
            if (h.isBoundary()) { halfCotan[h.id] = 0.; continue; };
            Row3d vn = h.next().vec();
            Row3d vp = h.prev().vec() * -1;
            halfCotan[h.id] = vp.dot(vn) / vp.cross(vn).norm();
        }
        for (Edge e: edges) {
            edgeCotan[e.id] = (e.half().cot() + e.half().twin().cot()) * 0.5;
        }

        /*--- setup dual area of vertex ---*/
        for (Vert v: verts) {
            double a1 = 0.;
            double a2 = 0.;
            for (Half h: v.adjHalfs()) {
                if (h.face().id != -1) a1 += faceArea[h.face().id];
                a2 += h.cot() * h.vec().squaredNorm();
                a2 += h.prev().cot() * h.prev().vec().squaredNorm();
            }
            baryDualArea[v.id] = a1 / 3.;
            circDualArea[v.id] = a2 * 0.125;
        }

        /*--- setup a filter to check whether a vertex is on boundary or not ---*/
        isBV.assign(nV, false);
        for (Edge e: edges) {
            if (e.isBoundary()) {
                isBV[e.half().tail().id] = true;
                isBV[e.half().head().id] = true;
            }
        }

        /*--- setup halfedge angle on vertex coordinate, and angle defect on vertex ---*/
        for (Vert v: verts) {
            double sum = 0;
            for (Half h: v.adjHalfs()) {
                heArgOnVert[h.id] = sum;
                Row3d v1 = h.vec().normalized();
                Row3d v2 = h.prev().twin().vec().normalized();
                sum += acos(v1.dot(v2));
            }
            angleDefect[v.id] = v.isBoundary() ? PI - sum : TwoPI - sum;
            for (Half h: v.adjHalfs()) heArgOnVert[h.id] *= TwoPI / sum;
        }

        for (Vert v: verts) {
            /*--- setup scalar mean curvature ---*/
            double sum = 0;
            for (Half h: v.adjHalfs()) {
                sum += h.darg() * h.len();
            }
            scalarMeanCurvature(v.id) = sum * 0.5;

            /*--- setup principal curvatures ---*/
            double a = circDualArea(v.id);
            double h = scalarMeanCurvature(v.id) / a;
            double k = angleDefect(v.id) / a;
            double d = sqrt(h * h - k);
            principalCurvatures.row(v.id) << h - d, h + d;
        }
    }
}

[[nodiscard]] SprsD Hmsh::d0() const {
    SprsD S(nE, nV);
    std::vector<TripD> T(nE);
    for (Edge e: edges) {
        T.emplace_back(e.id, e.half().tail().id, -1.);
        T.emplace_back(e.id, e.half().head().id, 1.);
    }
    S.setFromTriplets(T.begin(), T.end());
    return S;
}

[[nodiscard]] SprsD Hmsh::d1() const {
    SprsD S(nF, nE);
    std::vector<TripD> T(nF);
    for (Face f: faces) {
        for (Half h: f.adjHalfs()) {
            Edge e = h.edge();
            T.emplace_back(f.id, e.id, h.isCanonical() ? 1. : -1.);
        }}
    S.setFromTriplets(T.begin(), T.end());
    return S;
}

[[nodiscard]] SprsD Hmsh::h2()  const {
    VecXd tmp(nF);
    for (Face f: faces) {
        double l1 = f.half().len();
        double l2 = f.half().next().len();
        double l3 = f.half().prev().len();
        double s = (l1 + l2 + l3) * 0.5;
        tmp[f.id] = 1. / sqrt(s * (s - l1) * (s - l2) * l3);
    }
    return static_cast<SprsD>(tmp.asDiagonal());
}

[[nodiscard]] SprsD Hmsh::h2i() const {
    VecXd tmp(nF);
    for (Face f: faces) {
        double l1 = f.half().len();
        double l2 = f.half().next().len();
        double l3 = f.half().prev().len();
        double s = (l1 + l2 + l3) * 0.5;
        tmp[f.id] = 1. / sqrt(s * (s - l1) * (s - l2) * l3);
    }
    return static_cast<SprsD>(tmp.asDiagonal().inverse());
}
}

#endif
