// Serial Geant4 adapter for controlled interface experiments.
// Output doubles retain IEEE-754 bits; no histogram or scientific acceptance.
#include "G4Event.hh"
#include "G4GDMLParser.hh"
#include "G4GeneralParticleSource.hh"
#include "G4GammaGeneralProcess.hh"
#include "G4HadronicProcess.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4NistManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4VProcess.hh"
#include "G4PhysListFactory.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Region.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UImanager.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4Version.hh"
#include "G4VTouchable.hh"
#include "Randomize.hh"
#include "CLHEP/Random/JamesRandom.h"
#include <openssl/evp.h>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>

namespace fs = std::filesystem;

std::string bits(double value) {
  static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
  std::uint64_t raw;
  std::memcpy(&raw, &value, sizeof(raw));
  std::ostringstream stream;
  stream << std::hex << std::setfill('0') << std::setw(16) << raw;
  return stream.str();
}

#include "hp_readback.hh"
#ifdef XCODE_ENABLE_INCL42_EXPERIMENT
#include "incl_experiment.hh"
#endif

struct Output {
  std::ofstream events, exits, primaries;
  int event = -1, completed = 0, steps = 0, crossings = 0, neutrons = 0;
  double deposit = 0;
  explicit Output(const fs::path& root)
      : events(root/"events.csv"), exits(root/"exits.csv"), primaries(root/"primaries.csv") {
    if (!events || !exits || !primaries) throw std::runtime_error("Cannot open outputs");
    events << "event,steps,exits,neutron_exits,edep_MeV_bits\n";
    exits << "event,track,parent,pdg,energy_MeV_bits,x_mm_bits,y_mm_bits,z_mm_bits,dx_bits,dy_bits,dz_bits,time_ns_bits\n";
    primaries << "event,pdg,x_mm_bits,y_mm_bits,z_mm_bits,px_MeV_bits,py_MeV_bits,pz_MeV_bits\n";
  }
};

class Detector final : public G4VUserDetectorConstruction {
 public:
  explicit Detector(std::string path) : path_(std::move(path)) {}
  G4VPhysicalVolume* Construct() override {
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
    parser_.Read(path_, false);
    return parser_.GetWorldVolume();
  }
  G4LogicalVolume* Target(const std::string& name) { return parser_.GetVolume(name); }
  G4VPhysicalVolume* World() { return parser_.GetWorldVolume(); }
 private:
  std::string path_;
  G4GDMLParser parser_;
};

class Primary final : public G4VUserPrimaryGeneratorAction {
 public:
  explicit Primary(Output& out) : out_(out) {}
  void GeneratePrimaries(G4Event* event) override {
    source_.GeneratePrimaryVertex(event);
    if (event->GetNumberOfPrimaryVertex() != 1) throw std::runtime_error("Expected one primary vertex");
    auto* vertex = event->GetPrimaryVertex();
    if (vertex->GetNumberOfParticle() != 1) throw std::runtime_error("Expected one primary");
    auto* particle = vertex->GetPrimary();
    auto pos = vertex->GetPosition();
    auto momentum = particle->GetMomentum();
    out_.primaries << event->GetEventID() << ',' << particle->GetPDGcode()
      << ',' << bits(pos.x()/mm) << ',' << bits(pos.y()/mm) << ',' << bits(pos.z()/mm)
      << ',' << bits(momentum.x()/MeV) << ',' << bits(momentum.y()/MeV) << ',' << bits(momentum.z()/MeV) << '\n';
  }
 private:
  Output& out_;
  G4GeneralParticleSource source_;
};

class Event final : public G4UserEventAction {
 public:
  explicit Event(Output& out) : out_(out) {}
  void BeginOfEventAction(const G4Event* event) override {
    out_.event = event->GetEventID();
    out_.steps = out_.crossings = out_.neutrons = 0;
    out_.deposit = 0;
  }
  void EndOfEventAction(const G4Event* event) override {
    if (event->IsAborted()) throw std::runtime_error("Aborted event");
    out_.events << out_.event << ',' << out_.steps << ',' << out_.crossings
      << ',' << out_.neutrons << ',' << bits(out_.deposit/MeV) << '\n';
    ++out_.completed;
  }
 private:
  Output& out_;
};

class Step final : public G4UserSteppingAction {
 public:
  Step(Output& out, std::string target) : out_(out), target_(std::move(target)) {}
  void UserSteppingAction(const G4Step* step) override {
    auto* pre = step->GetPreStepPoint();
    auto* volume = pre->GetTouchable()->GetVolume();
    if (!volume || volume->GetLogicalVolume()->GetName() != target_) return;
    ++out_.steps;
    out_.deposit += step->GetTotalEnergyDeposit();
    auto* post = step->GetPostStepPoint();
    if (post->GetStepStatus() != fGeomBoundary) return;
    auto* next = post->GetTouchable()->GetVolume();
    if (next && next->GetLogicalVolume()->GetName() == target_) return;
    auto* track = step->GetTrack();
    auto pos = post->GetPosition();
    auto direction = post->GetMomentumDirection();
    const int pdg = track->GetParticleDefinition()->GetPDGEncoding();
    ++out_.crossings;
    if (pdg == 2112) ++out_.neutrons;
    out_.exits << out_.event << ',' << track->GetTrackID() << ',' << track->GetParentID() << ',' << pdg
      << ',' << bits(post->GetKineticEnergy()/MeV)
      << ',' << bits(pos.x()/mm) << ',' << bits(pos.y()/mm) << ',' << bits(pos.z()/mm)
      << ',' << bits(direction.x()) << ',' << bits(direction.y()) << ',' << bits(direction.z())
      << ',' << bits(post->GetGlobalTime()/ns) << '\n';
  }
 private:
  Output& out_;
  std::string target_;
};

std::string read_source(const fs::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("Cannot read source macro");
  std::ostringstream data;
  data << input.rdbuf();
  if (input.bad()) throw std::runtime_error("Cannot finish reading source macro");
  return data.str();
}

std::string sha256(const std::string& data) {
  std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> context(
      EVP_MD_CTX_new(), EVP_MD_CTX_free);
  unsigned char digest[EVP_MAX_MD_SIZE];
  unsigned int size = 0;
  if (!context || EVP_DigestInit_ex(context.get(), EVP_sha256(), nullptr) != 1 ||
      EVP_DigestUpdate(context.get(), data.data(), data.size()) != 1 ||
      EVP_DigestFinal_ex(context.get(), digest, &size) != 1 || size != 32)
    throw std::runtime_error("SHA-256 failed");
  std::ostringstream text;
  text << std::hex << std::setfill('0');
  for (unsigned int i = 0; i < size; ++i)
    text << std::setw(2) << static_cast<unsigned int>(digest[i]);
  return text.str();
}

void apply_source(const std::string& source) {
  std::istringstream input(source);
  std::string line;
  while (std::getline(input, line)) {
    line = line.substr(0, line.find('#'));
    const auto start = line.find_first_not_of(" \t\r");
    if (start == std::string::npos) continue;
    line = line.substr(start);
    if (line.rfind("/gps/", 0) != 0) throw std::runtime_error("Source macro allows GPS commands only");
    if (G4UImanager::GetUIpointer()->ApplyCommand(line) != 0)
      throw std::runtime_error("Failed GPS command: " + line);
  }
}

int main(int argc, char** argv) {
  try {
#ifdef XCODE_ENABLE_INCL42_EXPERIMENT
    const std::string incl_mode = incl_experiment_mode(argc, argv);
#else
    if (argc != 8) throw std::runtime_error("usage: g4_transport geometry macro physics target events seed output");
#endif
    const std::string physics_name = argv[3], target_name = argv[4];
    const int count = std::stoi(argv[5]);
    const long seed = std::stol(argv[6]);
    const fs::path root = argv[7];
    if (count < 0 || count > 10000 || seed <= 0) throw std::runtime_error("Invalid pilot budget or seed");
    if (!fs::create_directory(root)) throw std::runtime_error("Output directory already exists");
    Output out(root);
    G4Random::setTheEngine(new CLHEP::HepJamesRandom);
    G4Random::setTheSeed(seed);
    G4Random::saveEngineStatus((root/"rng_seed.txt").string().c_str());
    auto manager = std::make_unique<G4RunManager>();
    if (typeid(*manager) != typeid(G4RunManager) || static_cast<int>(manager->GetRunManagerType()) != 0)
      throw std::runtime_error("Expected explicit serial G4RunManager");
    auto* detector = new Detector(argv[1]);
    manager->SetUserInitialization(detector);
    G4PhysListFactory factory;
    if (!factory.IsReferencePhysList(physics_name)) throw std::runtime_error("Unknown physics list");
    auto* physics = factory.GetReferencePhysList(physics_name);
    physics->SetVerboseLevel(0);
    // Sensitivity runs may request another range cut; unset keeps the established 0.1 mm.
    const char* cut_text = std::getenv("XCODE_G4_CUT_MM");
    const double requested_cut_mm = cut_text ? std::stod(cut_text) : 0.1;
    if (!std::isfinite(requested_cut_mm) || requested_cut_mm <= 0) throw std::runtime_error("Invalid requested cut");
    physics->SetDefaultCutValue(requested_cut_mm*mm);
    // Readback probes may switch gamma-nuclear off to show the readback can tell on from off;
    // unset keeps the physics list default. Applied before initialisation, as the messenger requires.
    const char* gamma_nuclear_text = std::getenv("XCODE_G4_GAMMA_NUCLEAR");
    const std::string gamma_nuclear_request = gamma_nuclear_text ? gamma_nuclear_text : "default";
    if (gamma_nuclear_request == "off") {
      if (G4UImanager::GetUIpointer()->ApplyCommand("/physics_lists/em/GammaNuclear false") != 0)
        throw std::runtime_error("GammaNuclear command failed");
    } else if (gamma_nuclear_request != "default") {
      throw std::runtime_error("Invalid XCODE_G4_GAMMA_NUCLEAR");
    }
    manager->SetUserInitialization(physics);
    manager->SetUserAction(new Primary(out));
    manager->SetUserAction(new Event(out));
    manager->SetUserAction(new Step(out, target_name));
    manager->Initialize();
    const auto macro_bytes = read_source(argv[2]);
    apply_source(macro_bytes);
#ifdef XCODE_ENABLE_INCL42_EXPERIMENT
    if (physics_name != "QGSP_INCLXX_HP") throw std::runtime_error("Q062 experiment requires QGSP_INCLXX_HP");
    incl_experiment_configure(incl_mode, root);
#endif
    auto* target = detector->Target(target_name);
    if (!target) throw std::runtime_error("Target volume missing");
    const auto* region = target->GetRegion();
    const auto* cuts = region ? region->GetProductionCuts() : nullptr;
    if (!cuts) throw std::runtime_error("Target production cuts missing");
    const double effective[] = {
      manager->GetUserPhysicsList()->GetDefaultCutValue()/mm,
      cuts->GetProductionCut("gamma")/mm, cuts->GetProductionCut("e-")/mm,
      cuts->GetProductionCut("e+")/mm, cuts->GetProductionCut("proton")/mm
    };
    for (double value : effective)
      if (!std::isfinite(value) || value < 0) throw std::runtime_error("Invalid effective cut");
    std::ofstream identity(root/"identity.tsv");
    identity << "schema_version\t3\nengine\t" << G4Random::getTheEngine()->name()
      << "\nseed\t" << seed << "\nrequested_events\t" << count
      << "\ngeant4_version\t" << G4Version
      << "\nrun_manager_class\tG4RunManager"
      << "\nphysics\t" << physics_name << "\ntarget\t" << target_name
      << "\nrequested_cut_mm_bits\t" << bits(requested_cut_mm)
      << "\neffective_default_cut_mm_bits\t" << bits(effective[0])
      << "\ntarget_gamma_cut_mm_bits\t" << bits(effective[1])
      << "\ntarget_electron_cut_mm_bits\t" << bits(effective[2])
      << "\ntarget_positron_cut_mm_bits\t" << bits(effective[3])
      << "\ntarget_proton_cut_mm_bits\t" << bits(effective[4])
      << "\nmacro_sha256\t" << sha256(macro_bytes)
      << "\nmaterial\t" << target->GetMaterial()->GetName()
      << "\ndensity_g_cm3_bits\t" << bits(target->GetMaterial()->GetDensity()/(g/cm3))
      << "\nvolume_mm3_bits\t" << bits(target->GetSolid()->GetCubicVolume()/mm3) << '\n';
    identity.close();
    if (!identity) throw std::runtime_error("Identity write failed");
    // Keep the established L1 schema byte-for-byte stable. M1 has its own evidence channel.
    const auto core_record = read_source(root/"identity.tsv");
    std::ofstream details(root/"m1-readback.tsv");
    details << "schema_version\tm1-3\n" << core_record.substr(core_record.find('\n')+1)
            << "material_temperature_K_bits\t" << bits(target->GetMaterial()->GetTemperature()/kelvin) << '\n'
            << "neutron_data_path\t" << (std::getenv("G4NEUTRONHPDATA") ? std::getenv("G4NEUTRONHPDATA") : "") << '\n';
    for (int i = 0; const auto* constructor = physics->GetPhysics(i); ++i) {
      details << "physics_constructor_" << i << '\t' << constructor->GetPhysicsName() << '\n';
      details << "physics_constructor_type_" << i << '\t' << typeid(*constructor).name() << '\n';
    }
    auto* neutron = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
    if (!neutron || !neutron->GetProcessManager()) throw std::runtime_error("Neutron process manager missing");
    hp_readback(details, neutron, target->GetMaterial());
    const auto* processes = neutron->GetProcessManager()->GetProcessList();
    for (int i = 0; i < neutron->GetProcessManager()->GetProcessListLength(); ++i)
      details << "neutron_process_" << i << '\t' << (*processes)[i]->GetProcessName() << '\n';
    // Count + every process name proves that a missing named process really is absent.
    for (const auto& particle_name : {"gamma", "e-", "e+"}) {
      auto* particle = G4ParticleTable::GetParticleTable()->FindParticle(particle_name);
      if (!particle || !particle->GetProcessManager()) throw std::runtime_error("EM process manager missing");
      auto* process_manager = particle->GetProcessManager();
      const auto* list = process_manager->GetProcessList();
      details << particle_name << "_process_count\t" << process_manager->GetProcessListLength() << '\n';
      for (int i = 0; i < process_manager->GetProcessListLength(); ++i)
        details << particle_name << "_process_" << i << '\t' << (*list)[i]->GetProcessName() << '\n';
      if (std::string(particle_name) != "gamma") continue;
      // With the general process active, photonNuclear is added inside GammaGeneralProc
      // (11.4.3 G4EmExtraPhysics) and never appears in the list above.
      int general = 0;
      for (int i = 0; i < process_manager->GetProcessListLength(); ++i) {
        auto* combined = dynamic_cast<G4GammaGeneralProcess*>((*list)[i]);
        if (!combined) continue;
        const G4HadronicProcess* nuclear = combined->GetGammaNuclear();
        details << "gamma_general_" << general << "_name\t" << combined->GetProcessName() << '\n'
                << "gamma_general_" << general << "_active\t" << process_manager->GetProcessActivation(combined) << '\n'
                << "gamma_general_" << general << "_nuclear\t" << (nuclear ? nuclear->GetProcessName() : G4String("none")) << '\n';
        ++general;
      }
      details << "gamma_general_count\t" << general << '\n'
              << "gamma_nuclear_request\t" << gamma_nuclear_request << '\n';
    }
    if (auto* tube = dynamic_cast<G4Tubs*>(target->GetSolid())) {
      details << "target_shape\tcylinder\n"
               << "target_radius_cm_bits\t" << bits(tube->GetOuterRadius()/cm) << '\n'
               << "target_length_cm_bits\t" << bits(2*tube->GetZHalfLength()/cm) << '\n';
    }
    auto* world_lv = detector->World()->GetLogicalVolume();
    if (auto* orb = dynamic_cast<G4Orb*>(world_lv->GetSolid()))
      details << "world_radius_cm_bits\t" << bits(orb->GetRadius()/cm) << '\n';
    details << "world_material\t" << world_lv->GetMaterial()->GetName() << '\n';
    for (std::size_t i = 0; i < world_lv->GetNoDaughters(); ++i) {
      auto* volume = world_lv->GetDaughter(i);
      if (volume->GetLogicalVolume() != target) continue;
      const auto position = volume->GetObjectTranslation();
      const auto axis = volume->GetObjectRotationValue() * G4ThreeVector(0,0,1);
      details << "target_x_cm_bits\t" << bits(position.x()/cm) << '\n'
               << "target_y_cm_bits\t" << bits(position.y()/cm) << '\n'
               << "target_z_cm_bits\t" << bits(position.z()/cm) << '\n'
               << "target_axis_x_bits\t" << bits(axis.x()) << '\n'
               << "target_axis_y_bits\t" << bits(axis.y()) << '\n'
               << "target_axis_z_bits\t" << bits(axis.z()) << '\n';
    }
    const auto* material = target->GetMaterial();
    details << "target_element_count\t" << material->GetNumberOfElements() << '\n';
    // m1-2: the atom density transport uses. Density and abundances alone miss declared isotope
    // masses, and the GDML reader may bind an element to preloaded NIST isotopes instead (T031).
    details << "target_atoms_per_cm3_bits\t" << bits(material->GetTotNbOfAtomsPerVolume()*cm3) << '\n';
    if (material->GetNumberOfElements() == 1) {
      const auto* element = (*material->GetElementVector())[0];
      for (std::size_t i = 0; i < element->GetNumberOfIsotopes(); ++i) {
        const auto* isotope = (*element->GetIsotopeVector())[i];
        details << "isotope_" << isotope->GetZ() << '_' << isotope->GetN() << "_fraction_bits\t"
                 << bits(element->GetRelativeAbundanceVector()[i]) << '\n';
        details << "molar_mass_" << isotope->GetZ() << '_' << isotope->GetN() << "_g_mole_bits\t"
                 << bits(isotope->GetA()/(g/mole)) << '\n';
      }
    }
    details.close();
    if (!details) throw std::runtime_error("M1 readback write failed");
    G4Random::saveEngineStatus((root/"rng_before.txt").string().c_str());
    if (count > 0) manager->BeamOn(count);
    G4Random::saveEngineStatus((root/"rng_after.txt").string().c_str());
    if (out.completed != count) throw std::runtime_error("Wrong completed event count");
    // Range cuts become per-material energies only at run start; record what transport used.
    if (count > 0) {
      const auto* table = G4ProductionCutsTable::GetProductionCutsTable();
      std::ofstream energies(root/"cuts-energy.tsv");
      energies << "material\tgamma_keV\telectron_keV\tpositron_keV\tproton_keV\n";
      for (std::size_t i = 0; i < table->GetTableSize(); ++i) {
        energies << table->GetMaterialCutsCouple(i)->GetMaterial()->GetName();
        for (std::size_t particle = 0; particle < 4; ++particle)
          energies << '\t' << std::setprecision(17) << (*table->GetEnergyCutsVector(particle))[i]/keV;
        energies << '\n';
      }
      if (!energies) throw std::runtime_error("Energy cut write failed");
    }
    out.events.flush(); out.exits.flush(); out.primaries.flush();
    if (!out.events || !out.exits || !out.primaries) throw std::runtime_error("Output write failed");
    std::cout << "TRANSPORT_COMPLETED " << out.completed << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "TRANSPORT_FAILED " << error.what() << '\n';
    return 2;
  }
}
