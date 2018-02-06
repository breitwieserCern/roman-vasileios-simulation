#ifndef CANCERGROWTH_H_
#define CANCERGROWTH_H_

#include "biodynamo.h"
#include "interface.h"

namespace bdm {

  // 0. Define my custom cell, which extends Cell by adding an extra
// data member cell_type.
  //oxygen_level should be between 0 and 1
  BDM_SIM_OBJECT(MyCell, Cell) {
    BDM_SIM_OBJECT_HEADER(MyCellExt, 1, can_divide_, type_id_, oxygen_level_, hypo_division_); // create the header with our new data member

  public:
    MyCellExt() {}
  MyCellExt(const std::array<double, 3>& position) : Base(position) {} // our creator

    // getter and setter for our new data member
    void SetCanDivide(bool d) { can_divide_[kIdx] = d; }
    bool GetCanDivide() { return can_divide_[kIdx]; }
    bool* GetCanDividePtr() { return can_divide_.data(); }

    void SetCellTypeID(int type_id) { type_id_[kIdx] = type_id; }
    int GetCellTypeID() { return type_id_[kIdx]; }
    int* GetCellTypeIDPtr() { return type_id_.data(); }

    void SetOxygenLevel(double oxygenLevel) {oxygen_level_[kIdx] = oxygenLevel; }
    double GetOxygenLevel() { return oxygen_level_[kIdx]; }
    double* GetOxygenLevelPtr() { return oxygen_level_.data(); }

    void SetHypoDiv(bool div) { hypo_division_[kIdx] = div; }
    bool GetHypoDiv() { return hypo_division_[kIdx]; }
    bool* GetHypoDivPtr() { return hypo_division_.data(); }

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
    // declare new data member and define their type
    // private data can only be accessed by public function and not directly
    vec<bool> can_divide_;
    vec<int> type_id_;
    vec<double> oxygen_level_;
    vec<bool> hypo_division_;
  };

// 1. Define growth behaviour
  struct GrowthModule : public BaseBiologyModule {
    
  GrowthModule() : BaseBiologyModule(gAllBmEvents) {}
    
    template <typename T>
      void Run(T* cell) {

      int growthSpeed;
      array<double, 3> cell_movements;
      double divideProba;
      double currentOxygenLevel=cell->GetOxygenLevel();
      // if normoxy: high division rate but low migration
      if (currentOxygenLevel > 0.7) {
        cell->SetHypoDiv(true);
        growthSpeed=100;
        cell_movements={gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1)}; // create an array of 3 ramdom numbers between -1 and 1
        divideProba=0.8;
      }
      // if hypoxy: low division rate but high migration
      else if (currentOxygenLevel > 0.3) {
        growthSpeed = 40;
        //TODO: not a random migration. depending on oxygen gradient. Can get oxygen gradien from feb3
        cell_movements={gTRandom.Uniform(-4, 4), gTRandom.Uniform(-4, 4), gTRandom.Uniform(-4, 4)}; // create an array of 3 ramdom numbers between -4 and 4
        divideProba=0.4;
      }
      // necrose
      else {
        return; // quit bio module
        // don't do anything, but can have a behaviour
//        growthSpeed=0;
//        cell_movements={0, 0, 0};
//        divideProba=0;
//        cell->SetHypoDiv(false);
      }
        
      //cell grow until it reach a diam of 40
      if (cell->GetDiameter() < 40) {
        cell->ChangeVolume(growthSpeed);
        cell->UpdateMassLocation(cell_movements);
        cell->SetPosition(cell->GetMassLocation());
        //Reset biological movement to 0.
        cell->SetTractorForce({0, 0, 0});
      }
      // when diam of 8, it has a chance of dividing (if it can divide)
      else if (cell->GetCanDivide()==true && cell->GetHypoDiv()==true) {
        double aNewRandomDouble=gTRandom.Uniform(0, 1);
        if (aNewRandomDouble <= divideProba) { //0.55 //0.65
          auto&& daughter = Divide(*cell);
          daughter.SetCellTypeID(cell->GetCellTypeID()); // daughter takes the type_id_ value of her mother
          daughter.SetCanDivide(true); // daughter will be able to divide
          daughter.SetHypoDiv(true); // daughter will be able to divide in hypoxy
          daughter.SetOxygenLevel(cell->GetOxygenLevel()); // daughter takes the oxygen_level_ value of her mother
        }
        // if it doesn't divide
        else {
          // it won't divide in hypoxy any more
          if (currentOxygenLevel <= 0.7) {
            cell->SetHypoDiv(false);
          }
          // it will never be able to divide any more
          if (currentOxygenLevel > 0.7) {
            cell->SetCanDivide(false);
          }
        }
      }
      
    } // end of Run()
    
//    bool IsCopied(BmEvent event) const { return true; }
    ClassDefNV(GrowthModule, 1);
    
  }; // end GrowthModule

// Define compile time parameter
  template <typename Backend>
    struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
    using BiologyModules = Variant<GrowthModule>;
    using AtomicTypes = VariadicTypedef<MyCell>;
  };

  template <typename TResourceManager = ResourceManager<>> 
  inline int Simulate(int argc, const char** argv) {
    InitializeBioDynamo(argc, argv);
    size_t cells_per_dim = 16;

    Param::bound_space_ = true;
    Param::min_bound_ = 0;
    Param::max_bound_ = 320; // cells_per_dim * space

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
      cell.SetCanDivide(true);
      cell.SetOxygenLevel(1);
      cell.SetHypoDiv(true);
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
    ModelInitializer::Grid3D(cells_per_dim, space, construct); // call Grid3D to create cells following "construct" model

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
        auto&& cell = (cells)[i];

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

    // OxygenLevel will be updated depending on feb3 simulations
    for (size_t i = 0; i < cells->size(); i++) {
      auto&& cell = (cells)[i];
      cell.SetOxygenLevel(0.5);
    }
    
    for (int i=0; i<400; i++) {
      scheduler.Simulate(1);
    }

    for (size_t i = 0; i < cells->size(); i++) {
      auto&& cell = (cells)[i];
      cell.SetOxygenLevel(0.1);
    }
    
    for (int i=0; i<400; i++) {
      scheduler.Simulate(1);
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
