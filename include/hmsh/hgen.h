#ifndef HMSH_HGEN_H
#define HMSH_HGEN_H
#include "hmsh.h"

/// Homology generator class with tree-cotree algorithm

namespace pddg {
class Hgen {
public:
    std::vector<int> generators;
    explicit Hgen(const Hmsh& mesh): mesh(mesh) {}

    void calcHomologyGens(bool includeBoundary = true) {
        auto shared = [](Face a, Face b) {
            for (Half h: a.adjHalfs()) if (h.twin().face().id == b.id) return h;
            std::cerr << "no halfedge is shared";
        };

        auto isInner = [&](const Edge e) {
            return !e.v0().isBoundary() && !e.v1().isBoundary();
        };

        calcPrimalTree(mesh.verts[0]);
        calcDualCotree(mesh.faces[0]);
        generators.resize(mesh.nH);

        int igen = 0;

        for(Edge e: mesh.edges) {
            if (e.isBoundary()) continue;
            if (!includeBoundary && !isInner(e)) continue;
            Half h = e.half();
            if (!isInPrimalTree(h) && !isInDualCotree(h)) {
                igen++;
                std::vector<Half> tmp1;
                std::vector<Half> tmp2;
                Face iF = h.face();
                Face jF = h.twin().face();

                while(fp[iF.id] != iF.id) {
                    auto kF = mesh.faces[fp[iF.id]];
                    tmp1.push_back(shared(iF, kF));
                    iF = kF;
                }

                while(fp[jF.id] != jF.id) {
                    auto kF = mesh.faces[fp[jF.id]];
                    tmp2.push_back(shared(jF, kF));
                    jF = kF;
                }

                size_t m = tmp1.size() - 1;
                size_t n = tmp2.size() - 1;
                while (tmp1[m].id == tmp2[n].id) { m--; n--; }
                std::vector<Half> gen;
                gen.push_back(h);
                for (size_t i = 0; i <= m; ++i) gen.push_back(tmp1[i].twin());
                for (size_t i = 0; i <= n; ++i) gen.push_back(tmp2[i]);
                for (Half h_ : gen) generators[h_.id] = igen;
            }
        }
    }

private:
    const Hmsh& mesh;
    std::vector<int> vp; // vert parent tree
    std::vector<int> fp; // face parent tree

    bool isInPrimalTree(const Half h)  {
        auto a = h.tail();
        auto b = h.head();
        return vp[a.id] == b.id || vp[b.id] == a.id;
    }

    bool isInDualCotree(const Half h) {
        auto a = h.face();
        auto b = h.twin().face();
        return fp[a.id] == b.id || fp[b.id] == a.id;
    }

    void calcPrimalTree(Vert bgn) {
        vp.resize(mesh.nV);
        for (auto v: mesh.verts) vp[v.id] = v.id;
        std::queue<Vert> que;
        que.push(bgn);
        while (!que.empty()) {
            Vert jV = que.front();
            que.pop();
            for (Half h : jV.adjHalfs()) {
                Vert v = h.head();
                if (!v.isBoundary() && vp[v.id] == v.id && v.id != bgn.id) {
                    vp[v.id] = jV.id;
                    que.push(v);
                }
            }
        }
    }

    void calcDualCotree(Face bgn) {
        fp.resize(mesh.nF);
        for(auto f: mesh.faces) fp[f.id] = f.id;
        std::queue<Face> que;
        que.push(bgn);
        while (!que.empty()) {
            Face jF = que.front();
            que.pop();
            for (Half h : jF.adjHalfs()) {
                if (isInPrimalTree(h)) continue;
                Face f = h.twin().face();
                if (!h.twin().isBoundary() && fp[f.id] == f.id && f.id != bgn.id) {
                    fp[f.id] = jF.id;
                    que.push(f);
                }
            }
        }
    }
};
}

#endif
