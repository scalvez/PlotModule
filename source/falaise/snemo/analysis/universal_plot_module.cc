// universal_plot_module.cc

// Ourselves:
#include <snemo/analysis/universal_plot_module.h>

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
#include <falaise/snemo/datamodels/topology_1e_pattern.h>
#include <falaise/snemo/datamodels/vertex_measurement.h>

namespace analysis {

  // Registration instantiation macro :
  DPP_MODULE_REGISTRATION_IMPLEMENT(universal_plot_module,
                                    "analysis::universal_plot_module");

  // Character separator between key for histogram dict.
  const char KEY_FIELD_SEPARATOR = '_';

  // Set the histogram pool used by the module :
  void universal_plot_module::set_histogram_pool(mygsl::histogram_pool & pool_)
  {
    DT_THROW_IF(is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is already initialized !");
    _histogram_pool_ = &pool_;
    return;
  }

  // Grab the histogram pool used by the module :
  mygsl::histogram_pool & universal_plot_module::grab_histogram_pool()
  {
    DT_THROW_IF(! is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");
    return *_histogram_pool_;
  }

  void universal_plot_module::_set_defaults()
  {
    _key_fields_.clear ();

    _histogram_pool_ = 0;

    return;
  }

  // Initialization :
  void universal_plot_module::initialize(const datatools::properties  & config_,
                                       datatools::service_manager   & service_manager_,
                                        dpp::module_handle_dict_type & module_dict_)
  {
    DT_THROW_IF(is_initialized(),
                std::logic_error,
                "Module '" << get_name() << "' is already initialized ! ");

    dpp::base_module::_common_initialize(config_);

    // Get the keys from 'Event Header' bank
    if (config_.has_key("key_fields"))
      {
        config_.fetch("key_fields", _key_fields_);
      }

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
        dpp::histogram_service & Histo = service_manager_.grab<dpp::histogram_service>(histogram_label);
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
  void universal_plot_module::reset()
  {
    DT_THROW_IF(! is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");

    // Tag the module as un-initialized :
    _set_initialized(false);
    _set_defaults();
    return;
  }

  // Constructor :
  universal_plot_module::universal_plot_module(datatools::logger::priority logging_priority_)
    : dpp::base_module(logging_priority_)
  {
    _set_defaults();
    return;
  }

  // Destructor :
  universal_plot_module::~universal_plot_module()
  {
    if (is_initialized()) universal_plot_module::reset();
    return;
  }

  // Processing :
  dpp::base_module::process_status universal_plot_module::process(datatools::things & data_record_)
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    DT_THROW_IF(! is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");

    // Check if the 'event header' record bank is available :
    const std::string eh_label = snemo::datamodel::data_info::default_event_header_label();
    if (! data_record_.has(eh_label))
      {
        DT_LOG_ERROR(get_logging_priority(), "Could not find any bank with label '"
                     << eh_label << "' !");
        return dpp::base_module::PROCESS_STOP;
      }
    const snemo::datamodel::event_header & eh
      = data_record_.get<snemo::datamodel::event_header>(eh_label);

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

    // Check if some 'topology_data' are available in the data model:
    // const std::string td_label = snemo::datamodel::data_info::default_topology_data_label();
    const std::string td_label = "TD";
    if (! data_record_.has(td_label)) {
      DT_LOG_ERROR(get_logging_priority(), "Missing topology data to be processed !");
      return dpp::base_module::PROCESS_ERROR;
    }

    // Get the 'topology_data' entry from the data model :
    const snemo::datamodel::topology_data & td
      = data_record_.get<snemo::datamodel::topology_data>(td_label);

    DT_LOG_DEBUG(get_logging_priority(), "Topology data : ");
    if (get_logging_priority() >= datatools::logger::PRIO_DEBUG) td.tree_dump();

    if (! td.has_pattern()) {
      // DT_LOG_ERROR(get_logging_priority(), "Missing pattern !");
      return dpp::base_module::PROCESS_ERROR;
    }
    const snemo::datamodel::base_topology_pattern & a_pattern = td.get_pattern();

    const std::string & a_pattern_id = a_pattern.get_pattern_id();

    if (a_pattern_id != "1e") {
      DT_LOG_WARNING(get_logging_priority(), "PlotModule only works for '2e' topology for now !");
      return dpp::base_module::PROCESS_CONTINUE;
    }

    const snemo::datamodel::topology_1e_pattern & a_1e_pattern
      = dynamic_cast<const snemo::datamodel::topology_1e_pattern &>(a_pattern);

    double energy = a_1e_pattern.get_electron_energy();

    // Getting histogram pool
    mygsl::histogram_pool & a_pool = grab_histogram_pool();


    // Build unique key for histogram map:
    std::ostringstream key;
    // Retrieving info from header bank:
    const datatools::properties & eh_properties = eh.get_properties();
    for (std::vector<std::string>::const_iterator
           ifield = _key_fields_.begin();
         ifield != _key_fields_.end(); ++ifield)
      {
        const std::string & a_field = *ifield;
        if (! eh_properties.has_key(a_field))
          {
            DT_LOG_WARNING(get_logging_priority(),
                           "No properties with key '" << a_field << "' "
                           << "has been found in event header !");
            continue;
          }

        if (eh_properties.is_vector(a_field))
          {
            DT_LOG_WARNING(get_logging_priority (),
                           "Stored properties '" << a_field << "' " << "must be scalar !");
            continue;
          }
        if (eh_properties.is_boolean(a_field))      key << eh_properties.fetch_boolean(a_field);
        else if (eh_properties.is_integer(a_field)) key << eh_properties.fetch_integer(a_field);
        else if (eh_properties.is_real(a_field))    key << eh_properties.fetch_real(a_field);
        else if (eh_properties.is_string(a_field))  key << eh_properties.fetch_string(a_field);
        // Add a underscore separator between fields
        key << KEY_FIELD_SEPARATOR;
      }


    key << "energy";

    if (! a_pool.has(key.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "energy_distrib");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d","energy_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo = a_pool.grab_1d(key.str());

    if(datatools::is_valid(energy))
      a_histo.fill(energy);

    double weight = 1.0;
    if (eh_properties.has_key("event.genbb_label")) {
      if (eh_properties.fetch_string("event.genbb_label").find("0nubb") != std::string::npos)
        // weight /= 1e5;
        weight /= 1e5*(48./50.);
      // if (eh_properties.fetch_string("event.genbb_label").find("2nubb-2MeV") != std::string::npos)
      //   weight /= 1e3;
      if (eh_properties.fetch_string("event.genbb_label").find("2nubb") != std::string::npos)
        // weight /= 1e5;
        weight /= 1e7;
      if (eh_properties.fetch_string("event.genbb_label").find("Tl208") != std::string::npos) {
        weight /= 1e7;
      }
      if (eh_properties.fetch_string("event.genbb_label").find("Bi214_Po214") != std::string::npos &&
          eh_properties.fetch_string("analysis.vertex_origin").find("source") != std::string::npos)
        weight /= 1e7;
      if (eh_properties.fetch_string("event.genbb_label").find("Bi214_Po214") != std::string::npos &&
          eh_properties.fetch_string("analysis.vertex_origin").find("wire_surface") != std::string::npos)
        weight /= 1e7;
    }

    if (eh_properties.has_key(mctools::event_utils::EVENT_GENBB_WEIGHT))
      {
        weight *= eh_properties.fetch_real(mctools::event_utils::EVENT_GENBB_WEIGHT);
      }

    // Store the weight into histogram properties
    if (! a_histo.get_auxiliaries().has_key("weight"))
      {
        a_histo.grab_auxiliaries().update("weight", weight);
      }

    DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return dpp::base_module::PROCESS_SUCCESS;
  }

} // namespace analysis

  // end of universal_plot_module.cc
  /*
  ** Local Variables: --
  ** mode: c++ --
  ** c-file-style: "gnu" --
  ** tab-width: 2 --
  ** End: --
  */
