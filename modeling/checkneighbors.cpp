#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<float, 3, bg::cs::cartesian> point;

struct Atom {
    int id;
    int type;
    point position;

    Atom(int id, int type, float x, float y, float z) : id(id), type(type), position(point(x, y, z)) {}
};

// Define how to get the 'point' from 'Atom'
struct atom_indexable {
    typedef point result_type;
    point operator()(Atom const& a) const { return a.position; }
};

// Define how to compare two 'Atom' objects
struct atom_equal {
    bool operator()(Atom const& a1, Atom const& a2) const {
        return a1.id == a2.id;
    }
};

typedef bgi::rtree<Atom, bgi::quadratic<16>, atom_indexable, atom_equal> rtree;

int main() {
    std::string input_filename = "cu_equilibrated_3nm.data";
    std::string output_filename = "output.txt";
    std::string newlammpsdata_file = "newlammpsdata.lmp";
    
    std::ifstream infile(input_filename);
    std::ofstream outfile(output_filename);
    std::ofstream newlammpsdata(newlammpsdata_file);
    
    if (!infile.is_open()) {
        std::cerr << "Failed to open file " << input_filename << "\n";
        return 1;
    }

    std::string line;
    bool in_atoms_section = false;
    std::vector<Atom> atoms;

    // Read and store coordinates and types
    while (getline(infile, line)) {
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        if (line.find("Atoms") != std::string::npos) {
            in_atoms_section = true;
            getline(infile, line); // Skip the next line which could be a comment or empty
            continue;
        }

        if (line.empty() || line[0] == '#') continue;

        if (in_atoms_section) {
            if (line.find("Bonds") != std::string::npos || line.find("Velocities") != std::string::npos) {
                break;
            } else {
                std::istringstream iss(line);
                int id, atype;
                float x, y, z;
                if (!(iss >> id >> atype >> x >> y >> z)) continue;
                atoms.emplace_back(id, atype, x, y, z);
            }
        }
    }
    infile.close();

    // Build the R-tree from the atom vector
    rtree tree;
    for (const auto& atom : atoms) {
        tree.insert(atom);
    }

    // Perform a radius search and update types
    if (outfile.is_open()) {
        outfile << "Number of atoms: " << atoms.size() << "\n";
        for (auto& atom : atoms) {
            std::vector<Atom> results;
            const float radius = 2.7f;
            bg::model::box<point> query_box(
                point(bg::get<0>(atom.position) - radius, bg::get<1>(atom.position) - radius, bg::get<2>(atom.position) - radius),
                point(bg::get<0>(atom.position) + radius, bg::get<1>(atom.position) + radius, bg::get<2>(atom.position) + radius)
            );
            tree.query(bgi::within(query_box), std::back_inserter(results));
            
            int neighbor_count = results.size() - 1; // Subtract 1 to not count itself
            if (neighbor_count == 12) {
                atom.type = 3;  // Update type to 3
            }
            
            outfile << "Atom " << atom.id << " at (" << bg::get<0>(atom.position) << ", " << bg::get<1>(atom.position) << ", " << bg::get<2>(atom.position)
                    << ") has " << neighbor_count << " neighbors within 2.7 units and is type " << atom.type << "." << std::endl;
        }
        outfile.close();
    } else {
        std::cerr << "Failed to open output file " << output_filename << "\n";
        return 1;
    }

    // Write to new LAMMPS data file
    if (newlammpsdata.is_open()) {
        newlammpsdata << "LAMMPS data file via C++ conversion\n\n";
        newlammpsdata << atoms.size() << " atoms\n\n";
        newlammpsdata << "Atoms\n\n";
        for (const auto& atom : atoms) {
            newlammpsdata << atom.id << " " << atom.type << " " << bg::get<0>(atom.position) << " " << bg::get<1>(atom.position) << " " << bg::get<2>(atom.position) << "\n";
        }
        newlammpsdata.close();
    } else {
        std::cerr << "Failed to open new LAMMPS data file " << newlammpsdata_file << "\n";
        return 1;
    }

    return 0;
}

