#include <gtest/gtest.h>
#include <igl/readOBJ.h>
#include <igl/doublearea.h>
#include <igl/barycenter.h>
#include <igl/edge_topology.h>
#include <igl/local_basis.h>
#include "hmsh.h"
#include "hgen.h"

using namespace pddg;

TEST(hmsh, check_graph_structure) {
    //TODO: Use testcase submodule repo & replace to relative path...
    std::string directory = "/Users/saki/dev/models/";
    std::vector<std::string> pathes = {
            "tetrahedron.obj",
            "beetle-alt.obj",
            "lilium.obj",
            "bunny.obj",
    };
    for (auto& path : pathes) {
        MatXd V;
        MatXi F;
        igl::readOBJ(directory + path, V, F);
        auto hmsh = std::make_unique<Hmsh>(V, F);
        //TODO: Better have some tests for homology group...
        auto hgen = std::make_unique<Hgen>(*hmsh);
        hgen->calcHomologyGens();

        MatXi EV;
        MatXi FE;
        MatXi EF;
        igl::edge_topology(V, F, EV, FE, EF);

        //--- check corner structure ---//
        for (Face f: hmsh->faces) {
        for (Half h: f.adjHalfs()) {
            ASSERT_EQ(h.crnr().face().id, f.id);
        }}

        for (Vert v: hmsh->verts) {
        for (Half h: v.adjHalfs()) {
            if (h.isBoundary()) continue;
            ASSERT_EQ(h.next().crnr().vert().id, v.id);
        }}

        ASSERT_LT((EV - hmsh->edge2vert).norm(), 1e-10);
        ASSERT_LT((EF - hmsh->edge2face).norm(), 1e-10);
        ASSERT_LT((FE - hmsh->face2edge).norm(), 1e-10);
    }
}

TEST(hmsh, check_properties) {
    //TODO: Use testcase submodule repo & replace to relative path...
    std::string directory = "/Users/saki/dev/models/";
    std::vector<std::string> pathes = {
            "tetrahedron.obj",
            "beetle-alt.obj",
            "lilium.obj",
            "bunny.obj",
    };
    for (auto& path : pathes) {
        MatXd V;
        MatXi F;
        igl::readOBJ(directory + path, V, F);
        auto hmsh = std::make_unique<Hmsh>(V, F);
        MatXd baryCenter;
        VecXd faceArea;
        MatXd fbx, fby, fbz;
        igl::doublearea(V, F, faceArea);
        igl::barycenter(V, F, baryCenter);
        igl::local_basis(V, F, fbx, fby, fbz);
        ASSERT_LT((baryCenter - hmsh->baryCenter).norm(), 1e-10);
        ASSERT_LT((faceArea * 0.5 - hmsh->faceArea).norm(), 1e-10);
        ASSERT_LT((fbx - hmsh->faceBasisX).norm(), 1e-10);
        ASSERT_LT((fby - hmsh->faceBasisY).norm(), 1e-10);
        ASSERT_LT((fbz - hmsh->faceNormal).norm(), 1e-10);
    }
}
