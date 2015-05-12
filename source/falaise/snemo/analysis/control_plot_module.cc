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

    // DT_LOG_DEBUG(get_logging_priority(), "Topology data : ");
    // if (get_logging_priority() >= datatools::logger::PRIO_DEBUG) td.tree_dump();

    // const snemo::datamodel::base_topology_pattern & a_pattern = TD.get_pattern();
    // const std::string & a_pattern_id = a_pattern.get_pattern_id();

    // if (a_pattern_id != "2e") {
    //   DT_LOG_WARNING(get_logging_priority(), "PlotModule only works for '2e' topology for now !");
    //   return dpp::base_module::PROCESS_ERROR;
    // }

    // double proba_int;
    // const snemo::datamodel::topology_2e_pattern * ptr_2e_pattern
    //   = dynamic_cast<const snemo::datamodel::topology_2e_pattern *>(&a_pattern);
    // if (ptr_2e_pattern->has_internal_probability()) {
    //   proba_int = ptr_2e_pattern->get_internal_probability();
    // }

    // std::cout << "debug proba int : " << proba_int << std::endl;


    // std::ostringstream key;
    // key << "tof";

    // // Getting histogram pool
    // mygsl::histogram_pool & a_pool = grab_histogram_pool();

    // if (! a_pool.has(key.str()))
    //   {
    //     mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "tof_probability");
    //     datatools::properties hconfig;
    //     hconfig.store_string("mode", "mimic");
    //     hconfig.store_string("mimic.histogram_1d", "tof_probability_template");
    //     mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
    //   }

    // // Getting the current histogram
    // mygsl::histogram_1d & a_histo = a_pool.grab_1d(key.str ());

    // a_histo.fill(proba_int);

    // DT_LOG_TRACE(get_logging_priority(), "Exiting.");
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
