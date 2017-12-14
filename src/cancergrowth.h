#ifndef CANCERGROWTH_H_
#define CANCERGROWTH_H_

#include "biodynamo.h"
#include "interface.h"

namespace bdm {

// Define growth behaviour
struct GrowthModule {
  int curr_celltype;
  double rand_double;

  template <typename T>
  void Run(T* cell) {
    if (cell->GetDiameter() <= 40) {
      cell->ChangeVolume(400);
    } else {
      Divide(*cell);
    }
  }

  bool IsCopied(BmEvent event) const { return true; }
  ClassDefNV(GrowthModule, 1);
};

/// Default Cell implementation does not contain a data member type id
/// Therefore, we extend cell and add this data member along with a Getter and
/// Setter. Furthermore, we overwrite the DivideImpl function which is called
/// for cell division, to decide in which way type id is copied from the mother
/// to the daughter cell.
BDM_SIM_OBJECT(MyCell, Cell) {
  BDM_SIM_OBJECT_HEADER(MyCellExt, 1, type_id_);

 public:
  MyCellExt() {}
  MyCellExt(const std::array<double, 3>& position) : Base(position) {}

  int GetCellTypeID() const { return type_id_[kIdx]; }
  void SetCellTypeID(int type_id) { type_id_[kIdx] = type_id; }

  void DivideImpl(void* daughter_vptr, double volume_ratio, double phi,
                  double theta) {
    auto daughter = static_cast<Self<Scalar>*>(daughter_vptr);
    double rand_double = gRandom.NextDouble();
    if (rand_double < 0.7) {
      daughter->type_id_[0] = 1;
      type_id_[kIdx] = 1;
    } else {
      daughter->type_id_[0] = 2;
      type_id_[kIdx] = 2;
    }
    Base::DivideImpl(daughter, volume_ratio, phi, theta);
  }

 private:
  vec<int> type_id_;
};

// Define compile time parameter
template <typename Backend>
struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
  using BiologyModules = Variant<GrowthModule>;
  using AtomicTypes = VariadicTypedef<MyCell>;
};

inline int Simulate(int argc, const char** argv) {
  InitializeBioDynamo(argc, argv);
  size_t cells_per_dim = 16;

  int curr_celltype;
  double rand_double;

  int type1_counter = 0;
  int type2_counter = 0;
  double mass1_sum = 0.0;
  double mass2_sum = 0.0;

  // Define initial model - in this case 3D grid of cells
  double space = 20;
  auto construct = [&](const std::array<double, 3>& position) {
    MyCell cell(position);
    cell.SetDiameter(30);
    cell.AddBiologyModule(GrowthModule());

    rand_double = gRandom.NextDouble();
    if (rand_double < 0.7) {
      curr_celltype = 1;
      cell.SetMass(1.2);
      // std::cout << "should imprint type 1 for normoxic cells" << std::endl;
    } else {
      curr_celltype = 2;
      cell.SetMass(0.8);
      // std::cout << "should imprint type hypoxic cells" << std::endl;
    }
    cell.SetCellTypeID(curr_celltype);
    return cell;
  };
  ModelInitializer::Grid3D(cells_per_dim, space, construct);

  Scheduler<> scheduler;
  auto cells = ResourceManager<>::Get()->Get<MyCell>();

  for (int curr_time = 0; curr_time < 600; curr_time++) {
    type1_counter = 0;
    type2_counter = 0;
    mass1_sum = 0.0;
    double mass2_sum = 0.0;

    scheduler.Simulate(1);

    // Readout relevant information

    for (size_t i = 0; i < cells->size(); i++) {
      auto&& cell = (*cells)[i];

      curr_celltype = cell.GetCellTypeID();
      if (curr_celltype < 2) {
        type1_counter++;
        mass1_sum += cell.GetMass();
      } else {
        type2_counter++;
        mass2_sum += cell.GetMass();
      }

      // a_rand = gRandom.NextDouble();
      // if (a_rand<0.01) {
      //  std::cout << "a random cell diameter and type: " << cell.GetDiameter()
      //  << " / " << curr_celltype << std::endl;
      //}
    }

    std::cout << "Current time step: " << curr_time
              << ". Number of type 0 and 1 cells: " << type1_counter << " / "
              << type2_counter << std::endl;
    std::cout << "Mass of type 0 and 1 cells: " << mass1_sum << " / "
              << mass2_sum << std::endl;
  }

  DiscontinuousInterfaceData myDisc_fd(1);
  myDisc_fd.normoxic_cells_mass.assign(1, mass1_sum);
  myDisc_fd.hypoxic_cells_mass.assign(1, mass2_sum);
  myDisc_fd.normoxic_cells_population.assign(1, type1_counter);
  myDisc_fd.hypoxic_cells_population.assign(1, type2_counter);

  return 0;
}

}  // namespace bdm

#endif  // CANCERGROWTH_H_
