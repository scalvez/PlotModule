// control_plot_module.cc

// Ourselves:
#include </home/calvez/nemo/work_dir/Falaise/trunk/modules/PlotModule/source/falaise/snemo/analysis/control_plot_module.h>

// Standard library:
#include <stdexcept>
#include <sstream>
#include <algorithm>

// Third party:
// - Boost:
#include <boost/fusion/iterator/next.hpp>

// - Bayeux/mctools
#include <mctools/utils.h>
#include <mctools/simulated_data.h>

// - Bayeux/datatools:
#include <datatools/clhep_units.h>
#include <datatools/service_manager.h>
// - Bayeux/mygsl
#include <mygsl/histogram_pool.h>
// - Bayeux/dpp
#include <dpp/histogram_service.h>
// - Bayeux/geomtools:
#include <geomtools/geometry_service.h>
#include <geomtools/manager.h>

// - Falaise
#include <falaise/snemo/processing/services.h>
#include <falaise/snemo/datamodels/data_model.h>
#include <falaise/snemo/datamodels/calibrated_data.h>
#include <falaise/snemo/datamodels/event_header.h>
#include <falaise/snemo/datamodels/particle_track.h>
#include <falaise/snemo/datamodels/particle_track_data.h>

#include <falaise/snemo/datamodels/topology_data.h>
#include <falaise/snemo/datamodels/topology_2e_pattern.h>

namespace analysis {

  // Registration instantiation macro :
  DPP_MODULE_REGISTRATION_IMPLEMENT(control_plot_module,
                                    "analysis::control_plot_module");

  // Character separator between key for histogram dict.
  const char KEY_FIELD_SEPARATOR = '_';

  // Set the histogram pool used by the module :
  void control_plot_module::set_histogram_pool(mygsl::histogram_pool & pool_)
  {
    DT_THROW_IF(is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is already initialized !");
    _histogram_pool_ = &pool_;
    return;
  }

  // Grab the histogram pool used by the module :
  mygsl::histogram_pool & control_plot_module::grab_histogram_pool()
  {
    DT_THROW_IF(! is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");
    return *_histogram_pool_;
  }

  void control_plot_module::_set_defaults()
  {
    _histogram_pool_ = 0;

    return;
  }

  // Initialization :
  void control_plot_module::initialize(const datatools::properties  & config_,
                                       datatools::service_manager   & service_manager_,
                                       dpp::module_handle_dict_type & module_dict_)
  {
    DT_THROW_IF(is_initialized(),
                std::logic_error,
                "Module '" << get_name() << "' is already initialized ! ");

    dpp::base_module::_common_initialize(config_);


    // Service label
    std::string histogram_label;
    if (config_.has_key("Histo_label"))
      {
        histogram_label = config_.fetch_string("Histo_label");
      }
    if (! _histogram_pool_)
      {
        DT_THROW_IF(histogram_label.empty(), std::logic_error,
                    "Module '" << get_name() << "' has no valid 'Histo_label' property !");

        DT_THROW_IF(! service_manager_.has(histogram_label) ||
                    ! service_manager_.is_a<dpp::histogram_service>(histogram_label),
                    std::logic_error,
                    "Module '" << get_name() << "' has no '" << histogram_label << "' service !");
        dpp::histogram_service & Histo
          = service_manager_.get<dpp::histogram_service>(histogram_label);
        set_histogram_pool(Histo.grab_pool());
        if (config_.has_key("Histo_output_files"))
          {
            std::vector<std::string> output_files;
            config_.fetch("Histo_output_files", output_files);
            for (size_t i = 0; i < output_files.size(); i++) {
              Histo.add_output_file(output_files[i]);
            }
          }
        if (config_.has_key("Histo_input_file"))
          {
            const std::string input_file = config_.fetch_string("Histo_input_file");
            Histo.load_from_boost_file(input_file);
          }
        if (config_.has_key("Histo_template_files"))
          {
            std::vector<std::string> template_files;
            config_.fetch("Histo_template_files", template_files);
            for (size_t i = 0; i < template_files.size(); i++) {
              Histo.grab_pool().load(template_files[i]);
            }
          }

        // Tag the module as initialized :
        _set_initialized(true);
        return;
      }
  }

  // Reset :
  void control_plot_module::reset()
  {
    DT_THROW_IF(! is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");

    // Tag the module as un-initialized :
    _set_initialized(false);
    _set_defaults();
    return;
  }

  // Constructor :
  control_plot_module::control_plot_module(datatools::logger::priority logging_priority_)
    : dpp::base_module(logging_priority_)
  {
    _set_defaults();
    return;
  }

  // Destructor :
  control_plot_module::~control_plot_module()
  {
    if (is_initialized()) control_plot_module::reset();
    return;
  }

  // Processing :
  dpp::base_module::process_status control_plot_module::process(datatools::things & data_record_)
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    DT_THROW_IF(! is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");

    // const process_status status = _compute_gamma_track_length(data_record_);
    // if (status != dpp::base_module::PROCESS_OK) {
    //   DT_LOG_ERROR(get_logging_priority(), "Computing the gamma track length fails !");
    //   return status;
    // }

    // Check if some 'topology_data' are available in the data model:
    const std::string td_label = snemo::datamodel::data_info::default_topology_data_label();
    if (! data_record_.has(td_label)) {
      DT_LOG_ERROR(get_logging_priority(), "Missing topology data to be processed !");
      return dpp::base_module::PROCESS_ERROR;
    }

  // Get the 'topology_data' entry from the data model :
    const snemo::datamodel::topology_data & td
      = data_record_.get<snemo::datamodel::topology_data>(td_label);

    DT_LOG_DEBUG(get_logging_priority(), "Topology data : ");
    if (get_logging_priority() >= datatools::logger::PRIO_DEBUG) td.tree_dump();

    const snemo::datamodel::base_topology_pattern & a_pattern = td.get_pattern();
    const std::string & a_pattern_id = a_pattern.get_pattern_id();

    if (a_pattern_id != "2e") {
      DT_LOG_WARNING(get_logging_priority(), "PlotModule only works for '2e' topology for now !");
      return dpp::base_module::PROCESS_ERROR;
    }

    const snemo::datamodel::topology_2e_pattern * ptr_2e_pattern
      = dynamic_cast<const snemo::datamodel::topology_2e_pattern *>(&a_pattern);

    double proba_int;
    if (ptr_2e_pattern->has_internal_probability()) {
      proba_int = ptr_2e_pattern->get_internal_probability();
    }

    double proba_ext;
    if (ptr_2e_pattern->has_external_probability()) {
      proba_ext = ptr_2e_pattern->get_external_probability();
    }

    double delta_y;
    if (ptr_2e_pattern->has_delta_vertices_y()) {
      delta_y = ptr_2e_pattern->get_delta_vertices_y();
    }

    double delta_z;
    if (ptr_2e_pattern->has_delta_vertices_z()) {
      delta_z = ptr_2e_pattern->get_delta_vertices_z();
    }

    double angle;
    if (ptr_2e_pattern->has_angle()) {
    }
    angle = ptr_2e_pattern->get_angle();
    // std::cout << "Angle : " << angle << std::endl;

    std::ostringstream key_pint;
    key_pint << "Internal TOF probability";

    // Getting histogram pool
    mygsl::histogram_pool & a_pool = grab_histogram_pool();

    if (! a_pool.has(key_pint.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key_pint.str(), "", "tof_probability");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d", "tof_probability_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_pint = a_pool.grab_1d(key_pint.str ());

    a_histo_pint.fill(proba_int);

    std::ostringstream key_pext;
    key_pext << "External TOF probability";

    if (! a_pool.has(key_pext.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key_pext.str(), "", "tof_probability");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d", "tof_probability_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_pext = a_pool.grab_1d(key_pext.str ());

    a_histo_pext.fill(proba_ext);

    std::ostringstream key_deltay;
    key_deltay << "#Delta_{y} vertices";

    if (! a_pool.has(key_deltay.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key_deltay.str(), "", "delta_y");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d", "delta_vertices_Y_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_deltay = a_pool.grab_1d(key_deltay.str ());

    a_histo_deltay.fill(delta_y);

    std::ostringstream key_deltaz;
    key_deltaz << "#Delta_{z} vertices";

    if (! a_pool.has(key_deltaz.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key_deltaz.str(), "", "delta_z");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d", "delta_vertices_Z_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_deltaz = a_pool.grab_1d(key_deltaz.str ());

    a_histo_deltaz.fill(delta_z);

    std::ostringstream key_angle;
    key_angle << "Cos(#theta)";

    if (! a_pool.has(key_angle.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key_angle.str(), "", "cos_angle");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d", "cos_angle_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_angle = a_pool.grab_1d(key_angle.str ());

    a_histo_angle.fill(cos(angle));

//-----Energy plots

// Check if the 'particle track' record bank is available :
    const std::string ptd_label = snemo::datamodel::data_info::default_particle_track_data_label();
    if (! data_record_.has(ptd_label))
      {
        DT_LOG_ERROR(get_logging_priority (), "Could not find any bank with label '"
                     << ptd_label << "' !");
        return dpp::base_module::PROCESS_STOP;
      }
    const snemo::datamodel::particle_track_data & ptd
      = data_record_.get<snemo::datamodel::particle_track_data>(ptd_label);

    const snemo::datamodel::particle_track_data::particle_collection_type & the_particles = ptd.get_particles();

//Under the 2e hypothesis
         const snemo::datamodel::particle_track electron_1 = the_particles.at(0).get();
         const snemo::datamodel::particle_track electron_2 = the_particles.at(1).get();

        const snemo::datamodel::calibrated_calorimeter_hit::collection_type &
          the_calorimeters_1 = electron_1.get_associated_calorimeter_hits ();

        const snemo::datamodel::calibrated_calorimeter_hit::collection_type &
          the_calorimeters_2 = electron_2.get_associated_calorimeter_hits ();

        if (the_calorimeters_1.size() > 1 || the_calorimeters_2.size() > 1)
          {
            DT_LOG_WARNING(get_logging_priority(),
                         "A particle is associated to more than 1 calorimeter !");
          }

 double energy_1 = the_calorimeters_1.at(0).get().get_energy();
 double energy_2 = the_calorimeters_2.at(0).get().get_energy();

// std::cout << "E1 " << energy_1 <<  std::endl;
// std::cout << "E2 " << energy_2 <<  std::endl;

 std::ostringstream key_Etot;
 key_Etot << "Etot";

 if (! a_pool.has(key_Etot.str()))
      {
    mygsl::histogram_1d & h = a_pool.add_1d(key_Etot.str(), "", "Etot");
    datatools::properties hconfig;
    hconfig.store_string("mode", "mimic");
    hconfig.store_string("mimic.histogram_1d", "energy_template");
    mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
  }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_Etot = a_pool.grab_1d(key_Etot.str ());

    a_histo_Etot.fill(energy_1 + energy_2);

    std::ostringstream key_Emin;
    key_Emin << "Emin";

    if (! a_pool.has(key_Emin.str()))
      {
    mygsl::histogram_1d & h = a_pool.add_1d(key_Emin.str(), "", "Emin");
    datatools::properties hconfig;
    hconfig.store_string("mode", "mimic");
    hconfig.store_string("mimic.histogram_1d", "energy_template");
    mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
  }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_Emin = a_pool.grab_1d(key_Emin.str ());

    a_histo_Emin.fill(std::min(energy_1,energy_2));


    std::ostringstream key_Emax;
    key_Emax << "Emax";

    if (! a_pool.has(key_Emax.str()))
      {
    mygsl::histogram_1d & h = a_pool.add_1d(key_Emax.str(), "", "Emax");
    datatools::properties hconfig;
    hconfig.store_string("mode", "mimic");
    hconfig.store_string("mimic.histogram_1d", "energy_template");
    mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
  }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_Emax = a_pool.grab_1d(key_Emax.str ());

    a_histo_Emax.fill(std::max(energy_1,energy_2));

    std::ostringstream key_EminEmax;
    key_EminEmax << "EminEmax";

    if (! a_pool.has(key_EminEmax.str()))
      {
    mygsl::histogram_2d & h = a_pool.add_2d(key_EminEmax.str(), "", "EminEmax");
    datatools::properties hconfig;
    hconfig.store_string("mode", "mimic");
    hconfig.store_string("mimic.histogram_2d", "Emin_Emax_template");
    mygsl::histogram_pool::init_histo_2d(h, hconfig, &a_pool);
  }

    // Getting the current histogram
    mygsl::histogram_2d & a_histo_EminEmax = a_pool.grab_2d(key_EminEmax.str ());

    a_histo_EminEmax.fill(std::min(energy_1,energy_2), std::max(energy_1,energy_2));

    DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return dpp::base_module::PROCESS_SUCCESS;
  }

} // namespace analysis

  // end of control_plot_module.cc
  /*
  ** Local Variables: --
  ** mode: c++ --
  ** c-file-style: "gnu" --
  ** tab-width: 2 --
  ** End: --
  */
