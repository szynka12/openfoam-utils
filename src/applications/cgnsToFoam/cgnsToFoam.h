#ifndef CGNSHELPERS_H
#define CGNSHELPERS_H

#include <string>
#include <vector>

// Foam headers
#include "cgnslib.h"
#include "error.H"

namespace cgh
{
void inline cgns_check_error(cgerr_t err)
{
  if (err) cg_error_exit();
}

int open_cgns(const char *file_name)
{
  Foam::Info << "Opening cgns file: " << file_name << Foam::endl;
  int file;
  cgns_check_error(cg_open(file_name, CG_MODE_READ, &file));
  return file;
}

struct Base
{
  int         n;
  std::string name;
  int         cell_dim;
  int         phys_dim;

  Base(int file, int n_) : n(n_)
  {
    char base_name_[256];
    cgns_check_error(cg_base_read(file, n, base_name_, &cell_dim, &phys_dim));
    name = base_name_;
  }

  void info()
  {
    Foam::Info << "Base " << n << " with name: " << name << Foam::endl;
    Foam::Info << "\tCell dim: " << cell_dim << Foam::endl;
    Foam::Info << "\tPhys. dim: " << phys_dim << Foam::endl;
  }
};

struct Zone
{
  int         n;
  int         n_nodes;
  int         n_cells;
  int         n_bc_nodes;
  std::string name;

  Zone(int file, int base, int n_) : n(n_)
  {
    cgsize_t isize[3];
    char     zone_name_[256];
    cgns_check_error(cg_zone_read(file, base, n, zone_name_, isize));
    n_nodes    = isize[0];
    n_cells    = isize[1];
    n_bc_nodes = isize[2];
    name       = zone_name_;
  }

  void info()
  {
    Foam::Info << "Zone " << n << " with name:" << name << Foam::endl;
    Foam::Info << "\t# cells: " << n_cells << Foam::endl;
    Foam::Info << "\t# points: " << n_nodes << Foam::endl;

    // TODO Is this information useful
    Foam::Info << "\t# boundary points: " << n_bc_nodes << Foam::endl;
  }
};

struct Boundary
{
  std::string           name;
  std::vector<cgsize_t> faces;
  bool                  is_contagious;
  BCType_t              type;

  Boundary(int file, int base, int zone, int bc_id)
  {
    char name_[33];

    PointSetType_t ptset_type;

    int        NormalIndex, ndataset;
    DataType_t NormalDataType;
    cgsize_t   NormalListSize;

    cgsize_t n_data;
    cgns_check_error(cg_boco_info(file, base, zone, bc_id, name_, &type,
                                  &ptset_type, &n_data, &NormalIndex,
                                  &NormalListSize, &NormalDataType, &ndataset));

    name = name_;
    faces.resize(n_data);

    cgns_check_error(
        cg_boco_read(file, base, zone, bc_id, faces.data(), nullptr));

    // Renumber since cgns numbers from 1
    for (auto &x : faces) { x -= 1; }

    is_contagious = (faces.back() - faces[0] + 1) == faces.size();
  }

  int size() { return faces.size(); }

  void info()
  {
    Foam::Info << "Boundary " << name << "name" << Foam::endl;
    Foam::Info << "\tSize: " << size() << Foam::endl;
    Foam::Info << "\tRange: [" << faces[0] << ", " << faces.back() << "]"
               << Foam::endl;
    Foam::Info << "\tcontiguous: " << (is_contagious ? "Yes" : "No")
               << Foam::endl;
  }

  void range()
  {
    Foam::Info << "\tRange: [" << faces[0] << ", " << faces.back() << "]"
               << Foam::endl;
  }

  bool operator<(const Boundary &other) const
  {
    return faces[0] < other.faces[0];
  }

  // This is currently unused, consider deleting
  bool has_face(int face_no)
  {
    if (!is_contagious)
    {
      Foam::FatalError
          << "You can only use 'has_face()' on contagious boundaries!"
          << Foam::endl;
    }
    return (face_no > faces[0]) && (face_no > faces.back());
  }
};

template <typename T>
std::vector<T> read_coord(int file, Base base, Zone zone, int direction)
{
  static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value,
                "Only 'double' or 'float' types supported.");

  DataType_t data_type;
  if (std::is_same<T, double>::value) { data_type = RealDouble; }
  else if (std::is_same<T, float>::value)
  {
    data_type = RealSingle;
  }

  std::string coord_name;
  switch (direction)
  {
    case 1:
    {
      coord_name = "CoordinateX";
      break;
    }
    case 2:
    {
      coord_name = "CoordinateY";
      break;
    }
    case 3:
    {
      coord_name = "CoordinateZ";
      break;
    }
    default:
      Foam::FatalError << "Direction number can be from range [1,3]"
                       << Foam::endl;
  }

  cgsize_t imin = 1;
  cgsize_t imax = zone.n_nodes;

  std::vector<T> coord(imax);

  cgns_check_error(cg_coord_read(file, base.n, zone.n, coord_name.c_str(),
                                 data_type, &imin, &imax, coord.data()));
  return coord;
}

std::pair<int, int> count_faces_and_cells(int file, int base, int zone)
{
  int n_sections;
  cgns_check_error(cg_nsections(file, base, zone, &n_sections));

  int n_faces = 0;
  int n_cells = 0;

  for (int sec_i = 1; sec_i <= n_sections; sec_i++)
  {
    char          sectionname[33];
    ElementType_t type;
    cgsize_t      start, end;
    int           nbndry, parent_flag;
    cgns_check_error(cg_section_read(file, base, zone, sec_i, sectionname,
                                     &type, &start, &end, &nbndry,
                                     &parent_flag));

    if (type == CGNS_ENUMV(NGON_n)) { n_faces += end - start + 1; }
    else if (type == CGNS_ENUMV(NFACE_n))
    {
      n_cells += end - start + 1;
    }
  }
  return std::make_pair(n_faces, n_cells);
}

// Read N_NFACE or N_NGON to prepared Foam::List<T>
template <class T, int OFFSET = 1>
void read_ngon(int file, int base, int zone, int n_elements, int section_id,
               typename Foam::List<T>::iterator &data_iterator)
{
  cgsize_t c_size;
  cgns_check_error(cg_ElementDataSize(file, base, zone, section_id, &c_size));

  std::vector<cgsize_t>           offsets_raw(n_elements + 1);
  std::vector<cgsize_t>           connectivity_raw(c_size);
  std::vector<cgsize_t>::iterator conn_it = connectivity_raw.begin();

  cgns_check_error(cg_poly_elements_read(file, base, zone, section_id,
                                         connectivity_raw.data(),
                                         offsets_raw.data(), nullptr));

  T tmp_ngon;

  for (int it = 1; it < offsets_raw.size(); ++it)
  {
    int n_labels = offsets_raw[it] - offsets_raw[it - 1];
    tmp_ngon.resize(n_labels);
    for (int i = 0; i < tmp_ngon.size(); i++)
    {
      // change numbering to start from 0 if the default value of OFFSET is
      // used. OFFSET = 0 is useful when reading cells when negative values are
      // also present
      tmp_ngon[i] = (*conn_it) - OFFSET;
      conn_it++;
    }
    (*data_iterator) = tmp_ngon;
    data_iterator++;
  }

  // Sanity check (we should be at the end of connectivity vector)
  if (conn_it != connectivity_raw.end())
  {
    Foam::FatalError << "Error while reading data from cgns file. "
                     << "Possibly bad ordering or range of elements. "
                     << Foam::endl;
  }
}
}  // namespace cgh

#endif
