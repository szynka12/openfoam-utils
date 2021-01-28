#ifndef CGNS_MESHTOOLS_H
#define CGNS_MESHTOOLS_H

#include "ListOps.H"
#include "face.H"
#include "faceList.H"
#include "label.H"
#include "labelList.H"

namespace mt
{
using cellNeigboursList = Foam::List<Foam::List<Foam::label>>;
using badCellPairsList  = Foam::List<std::pair<Foam::label, Foam::label>>;
using cellFaces         = Foam::List<Foam::List<Foam::label>>;
using cellInfoPair      = std::pair<cellFaces, cellNeigboursList>;

std::pair<Foam::label, Foam::label> find_matching_vertex(const Foam::face& f1,
                                                         const Foam::face& f2)
{
  forAll(f1, f1i)
  {
    forAll(f1, f2i)
    {
      if (f1[f1i] == f2[f2i]) { return std::make_pair(f1i, f2i); }
    }
  }
  Foam::FatalError << "Faces not connected by an edge!" << Foam::endl;
  return std::make_pair(-1, -1);
}

inline Foam::label bound_index(const Foam::label index, const Foam::label size)
{
  if (index >= size) { return index % size; }
  else if (index < 0)
  {
    return size + index;
  }
  else
  {
    return index;
  }
}

Foam::face combine_faces(const Foam::faceList&  faces,
                         const Foam::labelList& merged_faces)
{
  Foam::face recipient_face = faces[merged_faces[0]];
  for (Foam::label facei = 1; facei < merged_faces.size(); facei++)
  {
    Foam::face donor_face = faces[merged_faces[facei]];
    // This is dirty hack used to exit two nested loops

    auto indices = find_matching_vertex(recipient_face, donor_face);

    auto next_i_r = bound_index(indices.first + 1, recipient_face.size());
    auto prev_i_d = bound_index(indices.second - 1, donor_face.size());

    if (recipient_face[next_i_r] == donor_face[prev_i_d])
    {
      indices.first  = next_i_r;
      indices.second = prev_i_d;
    }

    Foam::inplaceRotateList<Foam::List, Foam::label>(recipient_face,
                                                     -indices.first);
    Foam::inplaceRotateList<Foam::List, Foam::label>(donor_face,
                                                     -indices.second - 1);

    Foam::boolList tmp(donor_face.size(), true);
    tmp.first() = false;
    tmp.last()  = false;

    recipient_face.append(Foam::subset(tmp, donor_face));
  }
  return recipient_face;
}

Foam::boolList repair_multiply_connected_cells(
    const cellInfoPair& cell_info, const badCellPairsList& bad_cells,
    const Foam::labelList& owner, const Foam::labelList& neighbour,
    Foam::faceList& faces)
{
  const auto& cell_list = cell_info.first;

  Foam::boolList bad_faces(faces.size(), true);
  forAll(bad_cells, bi)
  {
    //// grab the cells in question
    const auto& c1 = cell_list[bad_cells[bi].first];
    const auto& c2 = cell_list[bad_cells[bi].second];

    Foam::List<Foam::label> multiple_faces;

    forAll(c1, facei)
    {
      if (owner[c1[facei]] == bad_cells[bi].first
          && neighbour[c1[facei]] == bad_cells[bi].second)
      {
        multiple_faces.append(c1[facei]);
      }
    }
    forAll(c2, facei)
    {
      if (owner[c2[facei]] == bad_cells[bi].second
          && neighbour[c2[facei]] == bad_cells[bi].first)
      {
        multiple_faces.append(c2[facei]);
      }
    }
    Foam::Info << "\t"
               << "Correcting faces: " << multiple_faces << Foam::endl;

    faces[multiple_faces[0]] = combine_faces(faces, multiple_faces);
    Foam::Info << faces[multiple_faces[0]] << Foam::endl;

    for (Foam::label i = 1; i < multiple_faces.size(); i++)
    {
      bad_faces[multiple_faces[i]] = false;
    }
  }
  return bad_faces;
}

// Retururns list of faces per cell (only internal ones)  and list of cell
// neighbours (neighbours are sorted)
cellInfoPair cell_neighbours(const Foam::labelList& owners,
                             const Foam::labelList& neighbours)
{
  Foam::label n_cells = std::max(Foam::max(owners), Foam::max(neighbours)) + 1;

  Foam::List<Foam::List<Foam::label>> cell_nei(n_cells);
  cellFaces                           cell_faces(n_cells);

  forAll(neighbours, i)
  {
    Foam::label c1 = neighbours[i];
    Foam::label c2 = owners[i];

    cell_nei[c1].append(c2);
    cell_nei[c2].append(c1);

    cell_faces[c1].append(i);
    cell_faces[c2].append(i);
  }

  // Sort every list, so that finding duplicates would be easy
  std::for_each(cell_nei.begin(), cell_nei.end(),
                [](Foam::labelList& l) { std::sort(l.begin(), l.end()); });

  return std::make_pair(cell_faces, cell_nei);
}

badCellPairsList find_multiply_connected_cells(
    const cellNeigboursList& cell_nei)
{
  // Find multiply connected cells (as in the drawing). Lists are sorted so now
  // we simple compare the consecutive elements. We store here ordered pairs
  // (lower cell is first).
  badCellPairsList bad_cells;
  forAll(cell_nei, celli)
  {
    const auto& curr_nei = cell_nei[celli];
    for (Foam::label neii = 1; neii < curr_nei.size(); neii++)
    {
      if (curr_nei[neii - 1] == curr_nei[neii])
      {
        auto bcell = celli < curr_nei[neii]
                         ? std::make_pair(celli, curr_nei[neii])
                         : std::make_pair(curr_nei[neii], celli);

        bad_cells.append(bcell);
      }
    }
  }

  Foam::inplaceUniqueSort(bad_cells);

  return bad_cells;
}

}  // namespace mt

#endif
