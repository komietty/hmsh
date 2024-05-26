#include <memory>
#include <igl/readOBJ.h>
#include <polyscope/polyscope.h>
#include <polyscope/curve_network.h>
#include <polyscope/surface_mesh.h>
#include "hmsh.h"
#include "hgen.h"


using namespace pddg;

int main(int argc, char *argv[]) {
    polyscope::options::autocenterStructures = true;

    MatXd V;
    MatXi F;
    igl::readOBJ(argv[1], V, F);
    auto hmsh = std::make_unique<Hmsh>(V, F);
    auto hgen = std::make_unique<Hgen>(*hmsh);
    hgen->calcHomologyGens(false);

    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    auto ps = polyscope::registerSurfaceMesh("mesh", hmsh->pos, hmsh->idx);

    /*--- visuailize curvature ---*/
    ps->setSurfaceColor({0, 10./ 255., 27./ 255.});
    ps->addVertexScalarQuantity("Gaussian Curvature", hmsh->angleDefect);
    ps->addVertexScalarQuantity("Mean Curvature", hmsh->scalarMeanCurvature);
    ps->addVertexScalarQuantity("Principal Curvature 1", hmsh->principalCurvatures.col(0));
    ps->addVertexScalarQuantity("Principal Curvature 2", hmsh->principalCurvatures.col(1));
    ps->resetTransform();
    ps->setSmoothShade(true);

    /*--- visuailize face basis ---*/

    /*--- visuailize vert basis ---*/

    /*--- visuailize boundary ---*/
    for (Loop l: hmsh->loops) {
        std::vector<glm::vec3> nodes;
        std::vector<double> value;
        std::vector<std::array<size_t, 2>> edges;
        size_t n = 0;
        for (Half h: l.adjHalfs()) {
            Row3d p1 = h.tail().pos();
            Row3d p2 = h.head().pos();
            nodes.emplace_back(p1.x(), p1.y(), p1.z());
            nodes.emplace_back(p2.x(), p2.y(), p2.z());
            edges.emplace_back(std::array{ n, n + 1 });
            value.emplace_back(n);
            n += 2;
        }
        auto pn = polyscope::registerCurveNetwork("boundaryloop_" + std::to_string(l.id), nodes, edges);
        pn->addEdgeScalarQuantity("value", value);
        pn->resetTransform();
        pn->setRadius(0.001);
    }

    /*--- visuailize generators ---*/
    {
        std::vector<glm::vec3> nodes;
        std::vector<double> value;
        std::vector<std::array<size_t, 2>> edges;
        size_t n = 0;
        for (Half h: hmsh->halfs) {
            int val = hgen->generators[h.id];
            if (val > 0) {
                Row3d p1 = h.face().center();
                Row3d p2 = h.twin().face().center();
                nodes.emplace_back(p1.x(), p1.y(), p1.z());
                nodes.emplace_back(p2.x(), p2.y(), p2.z());
                edges.emplace_back(std::array{ n, n + 1 });
                value.emplace_back(val);
                n += 2;
            }
        }
        auto pn = polyscope::registerCurveNetwork("generator", nodes, edges);
        pn->addEdgeScalarQuantity("value", value);
        pn->resetTransform();
        pn->setRadius(0.001);
    }

    polyscope::show();
}
