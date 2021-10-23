#include <cinolib/meshes/meshes.h>
#include <cinolib/winding_number.h>
#include <cinolib/grid_projector.h>
#include <cinolib/feature_network.h>
#include <cinolib/feature_mapping.h>
#include <cinolib/export_surface.h>
#include <cinolib/padding.h>
#include <cinolib/smoother.h>

struct vert_compare{
    bool operator()(const cinolib::vec3d& a, const cinolib::vec3d& b) const {

       double eps = 1e-6;
       if(a.x()-b.x() < 0.0 && abs(a.x()-b.x()) > eps){
           return true;
       }
       else if(abs(a.x()-b.x()) < eps){
           if(a.y()-b.y() < 0.0 && abs(a.y()-b.y()) > eps){
               return true;
           }
           else if(abs(a.y()-b.y()) < eps){
               if(a.z()-b.z() < 0.0 && abs(a.z()-b.z()) > eps){
                   return true;
               }
           }
       }

       return false;
    }
};

void remove_external_polys(cinolib::Trimesh<> &m, cinolib::Hexmesh<> &hm){

    std::cout<<"Removing external polyhedra...."<<std::endl;
    std::vector<bool> polys_to_delete_bool(hm.num_polys());
    std::vector<uint> polys_to_delete;
    cinolib::PARALLEL_FOR(0, hm.num_polys(),1000, [&polys_to_delete_bool, &hm, &m](uint pid){
        cinolib::vec3d centroid = hm.poly_centroid(pid);
        if(cinolib::winding_number(m,centroid) == 0){
            polys_to_delete_bool[pid] = true;

        }

    });

    for(uint pid=0; pid<hm.num_polys(); pid++)
        if(polys_to_delete_bool[pid]) polys_to_delete.push_back(pid);

    cinolib::Hexmesh<> tmp_hm = hm;
    tmp_hm.polys_remove(polys_to_delete);

    std::map<cinolib::vec3d, uint, vert_compare> vmap;

    for(uint vid=0; vid<hm.num_verts(); vid++) vmap[hm.vert(vid)] = vid;

    for(uint eid=0; eid<tmp_hm.num_edges(); eid++){
        if(!tmp_hm.edge_is_manifold(eid)){

            uint eid_ = hm.edge_id(vmap[tmp_hm.edge_vert(eid, 0)], vmap[tmp_hm.edge_vert(eid, 1)]);
            for(uint pid : hm.adj_e2p(eid_)){
                polys_to_delete_bool[pid] = false;
            }

        }
    }

    for(uint vid=0; vid<tmp_hm.num_verts(); vid++){
        if(!tmp_hm.vert_is_manifold(vid)){

            uint vid_ = vmap[tmp_hm.vert(vid)];
            for(uint pid : hm.adj_v2p(vid_)){
                polys_to_delete_bool[pid] = false;
            }
        }
    }
    polys_to_delete.clear();
    for(uint pid=0; pid<hm.num_polys(); pid++){
        if(polys_to_delete_bool[pid]){
            polys_to_delete.push_back(pid);
            hm.poly_data(pid).label = 1;
        }
        else hm.poly_data(pid).label = 0;
    }


    hm.polys_remove(polys_to_delete);
}

void project_on_surface_smoother(cinolib::Hexmesh<> &hm, cinolib::Trimesh<> &target, cinolib::Quadmesh<> &peel){


    std::unordered_map<uint, uint> m_map;
    std::unordered_map<uint, uint> v_map;
    cinolib::export_surface(hm, peel, m_map, v_map);
    target.edge_set_flag(cinolib::MARKED, false);
    peel.edge_set_flag(cinolib::MARKED, false);
    cinolib::Quadmesh tmp_peel = peel;

    target.edge_mark_sharp_creases();
    cinolib::Octree tree;
    tree.build_from_mesh_edges(peel);

    uint num_samples = 10;
    for(uint eid=0; eid<target.num_edges(); eid++){
        break;
        if(target.edge_data(eid).flags[cinolib::MARKED]){

            std::vector<cinolib::vec3d> samples;
            for(uint i=0; i<=num_samples; i++){
                samples.push_back(target.edge_sample_at(eid, static_cast<double>(i)/static_cast<double>(num_samples)));
            }


            for(uint i=0; i<samples.size()-1; i++){
                uint id1, id2;
                double dist1, dist2;
                cinolib::vec3d pos1, pos2;

                tree.closest_point(samples[i], id1, pos1, dist1);
                tree.closest_point(samples[i+1], id2, pos2, dist2);

                uint source = 0;
                if(pos1.dist(peel.edge_vert(id1, 0)) < pos1.dist(peel.edge_vert(id1, 1)))
                    source = peel.edge_vert_id(id1, 0);
                else
                    source = peel.edge_vert_id(id1, 1);
                uint dest = 0;
                if(pos2.dist(peel.edge_vert(id2, 0)) < pos2.dist(peel.edge_vert(id2, 1)))
                    dest = peel.edge_vert_id(id2, 0);
                else
                    dest = peel.edge_vert_id(id2, 1);

                std::vector<uint> path;
                cinolib::dijkstra(peel, source, dest, path);
                if(path.size() == 0) continue;
                for(uint j=0; j<path.size()-1; j++){
                    int edge_to_mark = peel.edge_id(path[j], path[j+1]);
                    if(edge_to_mark < 0) std::cerr<<"This edge doen't exist"<<std::endl;
                    else{
                        peel.edge_set_flag(cinolib::MARKED, true, {static_cast<uint>(edge_to_mark)});
                    }
                }

            }


        }
    }
    cinolib::SmootherOptions options;
    options.n_iters = 3;
    options.laplacian_mode = cinolib::UNIFORM;
    cinolib::mesh_smoother(peel, target, options);


    std::cout<<"Projecting..."<<std::endl;
    for(uint vid=0; vid<peel.num_verts(); vid++){
        hm.vert(v_map[vid]) = peel.vert(vid);
    }


}

void project_on_surface(cinolib::Hexmesh<> &hm, cinolib::Trimesh<> &target){

    std::cout<<"Starting projection"<<std::endl;
    std::vector<std::vector<uint>> fn_tri, fn_quad;
    cinolib::feature_network(target, fn_tri);

    cinolib::Quadmesh<> peel;
    std::unordered_map<uint, uint> h2q, q2h;
    cinolib::export_surface(hm, peel, h2q, q2h);

    cinolib::feature_mapping(target, fn_tri, peel, fn_quad);

    for(auto &f : fn_quad){

        for(uint i=1; i<f.size(); ++i){
            uint v0 = q2h.at(f.at(i));
            uint v1 = q2h.at(f.at(i-1));
            int eid = hm.edge_id(v0, v1);

            assert(eid >= 0 );

            hm.edge_data(eid).flags[cinolib::CREASE] = true;
            hm.edge_data(eid).flags[cinolib::MARKED] = true;
        }
    }

    cinolib::GridProjectorOptions options;
    //options.SJ_thresh = -1;
    cinolib::grid_projector(hm, target, options);

    std::cout<<"Projection completed"<<std::endl;

}

int main(int argc, char *argv[])
{

    std::string mesh_name;
    std::string grid_name;
    std::string output_name;
    std::string projection_method = "mesh_smoother";

    if(argc > 2){
        mesh_name = std::string(argv[1]);
        grid_name = std::string(argv[2]);
        output_name = std::string(argv[3]);
        projection_method = std::string(argv[4]);
    }
    else{
        std::cerr<<"Wrong number of arguments"<<std::endl;
        return -1;
    }

    cinolib::Trimesh<> mesh(mesh_name.c_str());
    cinolib::Hexmesh<> grid(grid_name.c_str());

    remove_external_polys(mesh, grid);
    cinolib::padding(grid, false);


    cinolib::Quadmesh<> peel;

    if(projection_method == "mesh_smoother") project_on_surface_smoother(grid, mesh, peel);
    else project_on_surface(grid, mesh);

    grid.save(output_name.c_str());


    return  0;
}