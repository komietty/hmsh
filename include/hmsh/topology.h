// Copyright (C)
// 2013 Alec Jacobson <alecjacobson@gmail.com>,
// 2016 Amir Vaxman <avaxman@gmail.com>,
// 2024 Saki Komikado <komietty@gmail.com>,
// This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef HMSH_TOPOLOGY_H
#define HMSH_TOPOLOGY_H
#include "typedefs.h"

namespace pddg {
inline void etopo(
        const MatXi& F,
        MatXi& EV,
        MatXi& FE,
        MatXi& EF
) {
    std::vector<std::vector<int>> ETT;
    for(int f = 0; f < F.rows(); ++f)
        for (int i = 0; i < 3; ++i) {
            int v1 = F(f, i);
            int v2 = F(f, (i + 1) % 3);
            if (v1 > v2) std::swap(v1, v2);
            std::vector<int> r(4);
            r[0] = v1; r[1] = v2;
            r[2] = f;  r[3] = i;
            ETT.push_back(r);
        }
    std::sort(ETT.begin(), ETT.end());

    int En = 1;
    for(int i = 0; i < int(ETT.size())-1; ++i)
        if (!((ETT[i][0] == ETT[i + 1][0]) && (ETT[i][1] == ETT[i + 1][1])))++En;

    EV = MatXi::Constant((int)(En), 2, -1);
    FE = MatXi::Constant((int)(F.rows()), 3, -1);
    EF = MatXi::Constant((int)(En), 2, -1);
    En = 0;

    for(unsigned i = 0; i < ETT.size(); ++i) {
        if (i == ETT.size() - 1 || !((ETT[i][0] == ETT[i + 1][0]) && (ETT[i][1] == ETT[i + 1][1]))) {
            // Border edge
            std::vector<int>& r1 = ETT[i];
            EV(En, 0) = r1[0];
            EV(En, 1) = r1[1];
            EF(En, 0) = r1[2];
            FE(r1[2], r1[3]) = En;
        } else {
            std::vector<int>& r1 = ETT[i];
            std::vector<int>& r2 = ETT[i+1];
            EV(En, 0) = r1[0];
            EV(En, 1) = r1[1];
            EF(En, 0) = r1[2];
            EF(En, 1) = r2[2];
            FE(r1[2], r1[3]) = En;
            FE(r2[2], r2[3]) = En;
            ++i; // skip the next one
        }
        ++En;
    }

    // Sort the relation EF, accordingly to EV
    // the first one is the face on the left of the edge
    for(unsigned i = 0; i < EF.rows(); ++i) {
        int fid = EF(i, 0);
        bool flip = true;
        // search for edge EV.row(i)
        for (unsigned j = 0; j < 3; ++j)
            if ((F(fid, j) == EV(i, 0)) && (F(fid, (j + 1) % 3) == EV(i, 1))) flip = false;

        if (flip) {
            int tmp = EF(i,0);
            EF(i,0) = EF(i,1);
            EF(i,1) = tmp;
        }
    }
}

// Created a Double-Connected Edge-List (a.k.a. "halfedge structure") from the usual
// libhedra mesh representation. This data structure is very convenient for mesh editing
// and traversing, and the data structure is again only Eigen vectors and matrices.

//input:
//  D           #F by 1 - face degrees
//  F           #F by max(D) - vertex indices in face
//  EV          #E by 2 - edge vertex indices
//  EF          #E by 2 - edge face indices (EF(i,0) is left face, EF(i,1)=-1 if boundary
//  EFi         #E by 2 - position of edge in face by EF
//  innerEdges  vector of inner edges into EV

// Output:
// the number of halfedges can be determined by H=|HV/HE/HF|. It is 2*[Inner edges]+[Boundary Edges]
// VH   #V by 1 - Vertex to outgoing halfedge (into HE)
// EH   #E by 2 - edge to halfedge, where EH(i,0) halfedge is positively oriented, and EH(i,1)=-1 when boundary.
// FH   #F by max(D) - face to (correctly oriented) halfedge s.t. the origin vertex of FH(i,j) is F(i,j)
// HV   #H by 1 - origin vertex of the halfedge
// HE   #H by 1 - edge carrying this halfedge. It does not say which direction.
// HF   #F by 1 - face containing halfedge
// nextH, prevH, twinH - #H by 1 DCEL traversing operations. twinH(i)=-1 for boundary edges.

inline void dcel(
        const VecXi& D,
        const MatXi& F,
        const MatXi& EV,
        const MatXi& EF,
        const MatXi& EFi,
        VecXi& VH,
        MatXi& EH,
        MatXi& FH,
        VecXi& HV,
        VecXi& HE,
        VecXi& HF,
        VecXi& nextH,
        VecXi& prevH,
        VecXi& twinH)
{
    //doing a local halfedge structure for polygonal meshes
    EH = MatXi::Constant(EV.rows(), 2, -1);
    int numH = 0;

    for (int i = 0; i < EF.rows(); i++){
        if (EF(i, 0) != -1)
            EH(i, 0) = numH++;
        if (EF(i, 1) != -1)
            EH(i, 1) = numH++;
    }


    //halfedges to edge
    HE.conservativeResize(numH);
    for (int i = 0; i < EH.rows(); i++){
        if (EH(i, 0) != -1)
            HE(EH(i, 0)) = i;
        if (EH(i, 1) != -1)
            HE(EH(i, 1)) = i;
    }

    //halfedge to vertex and vice versa
    HV.conservativeResize(numH);
    VH.conservativeResize(EV.maxCoeff() + 1);
    for (int i = 0; i < EV.rows(); i++){
        if (EH(i, 0) != -1) {
            HV(EH(i, 0)) = EV(i, 0);
            VH(EV(i, 0)) = EH(i, 0);
        }
        if (EH(i, 1) != -1) {
            HV(EH(i, 1)) = EV(i, 1);
            VH(EV(i, 1)) = EH(i, 1);
        }
    }

    //halfedge to twin
    twinH = Eigen::VectorXi::Constant(numH, -1);
    for (int i = 0; i < EH.rows(); i++)
        if ((EH(i, 0) != -1) && (EH(i, 1) != -1)) {
            twinH(EH(i, 0)) = EH(i, 1);
            twinH(EH(i, 1)) = EH(i, 0);
        }

    //faces to halfedges and vice versa
    FH.resize(F.rows(), F.cols());
    HF.resize(numH);
    for (int i = 0; i < EF.rows(); i++){
        if (EF(i, 0) != -1) {
            FH(EF(i, 0), EFi(i, 0)) = EH(i, 0);
            HF(EH(i, 0)) = EF(i, 0);
        }
        if (EF(i, 1) != -1) {
            FH(EF(i, 1), EFi(i, 1)) = EH(i, 1);
            HF(EH(i, 1)) = EF(i, 1);
        }
    }

    //halfedge to next and prev
    nextH.conservativeResize(HE.rows());
    prevH.conservativeResize(HE.rows());
    for (int i = 0; i < D.rows(); i++){
        for (int j = 0; j < D(i); j++){
            nextH(FH(i, j)) = FH(i, (j + 1) % D(i));
            prevH(FH(i, (j + 1) % D(i))) = FH(i, j);
        }
    }
}
}
#endif
