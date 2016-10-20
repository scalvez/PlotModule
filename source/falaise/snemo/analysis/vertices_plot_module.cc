// vertices_plot_module.cc

// Ourselves:
#include <snemo/analysis/vertices_plot_module.h>

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
#include <falaise/snemo/datamodels/vertex_measurement.h>

namespace analysis {

  // Registration instantiation macro :
  DPP_MODULE_REGISTRATION_IMPLEMENT(vertices_plot_module,
                                    "analysis::vertices_plot_module");

  // Character separator between key for histogram dict.
  const char KEY_FIELD_SEPARATOR = '_';

  // Set the histogram pool used by the module :
  void vertices_plot_module::set_histogram_pool(mygsl::histogram_pool & pool_)
  {
    DT_THROW_IF(is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is already initialized !");
    _histogram_pool_ = &pool_;
    return;
  }

  // Grab the histogram pool used by the module :
  mygsl::histogram_pool & vertices_plot_module::grab_histogram_pool()
  {
    DT_THROW_IF(! is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");
    return *_histogram_pool_;
  }

  void vertices_plot_module::_set_defaults()
  {
    _histogram_pool_ = 0;

    return;
  }

  // Initialization :
  void vertices_plot_module::initialize(const datatools::properties  & config_,
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
  void vertices_plot_module::reset()
  {
    DT_THROW_IF(! is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");

    // Tag the module as un-initialized :
    _set_initialized(false);
    _set_defaults();
    return;
  }

  // Constructor :
  vertices_plot_module::vertices_plot_module(datatools::logger::priority logging_priority_)
    : dpp::base_module(logging_priority_)
  {
    _set_defaults();
    return;
  }

  // Destructor :
  vertices_plot_module::~vertices_plot_module()
  {
    if (is_initialized()) vertices_plot_module::reset();
    return;
  }

  // Processing :
  dpp::base_module::process_status vertices_plot_module::process(datatools::things & data_record_)
  {
    DT_LOG_TRACE(get_logging_priority(), "Entering...");
    DT_THROW_IF(! is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");

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

    if (a_pattern_id != "2e") {
      DT_LOG_WARNING(get_logging_priority(), "PlotModule only works for '2e' topology for now !");
      return dpp::base_module::PROCESS_ERROR;
    }

    // std::ostringstream key;
    // key << "vertices_probability";


    const snemo::datamodel::particle_track_data::particle_collection_type & the_particles = ptd.get_particles();

    // if(the_particles.size() != 1) {
    //   DT_LOG_WARNING(get_logging_priority(), "More than one particle !");
    //   return dpp::base_module::PROCESS_ERROR;
    // }

    // Getting histogram pool
    mygsl::histogram_pool & a_pool = grab_histogram_pool();

    std::ostringstream key;
    key << "vertex_distrib";

    if (! a_pool.has(key.str()))
      {
        mygsl::histogram_2d & h = a_pool.add_2d(key.str(), "", "vertices");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_2d", "vertex_distribution_template");
        mygsl::histogram_pool::init_histo_2d(h, hconfig, &a_pool);
      }
    mygsl::histogram_2d & a_histo = a_pool.grab_2d(key.str ());
    // geomtools::blur_spot tmp = dynamic_cast<const snemo::datamodel::vertex_measurement&> (a_pattern.get_measurement("vertex_e1_e2")).get_vertex();
    double vertex_y = dynamic_cast<const snemo::datamodel::vertex_measurement&> (a_pattern.get_measurement("vertex_e1_e2")).get_vertex().get_position().y();
    double vertex_z = dynamic_cast<const snemo::datamodel::vertex_measurement&> (a_pattern.get_measurement("vertex_e1_e2")).get_vertex().get_position().z();

    if(datatools::is_valid(vertex_y) && datatools::is_valid(vertex_z))
      a_histo.fill(vertex_y,vertex_z);

    DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return dpp::base_module::PROCESS_SUCCESS;

    /*
    std::ostringstream key_tot;
    key_tot << "track_length_distrib";

    if (! a_pool.has(key_tot.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key_tot.str(), "", "track_length_distrib");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d","track_length_efficiency_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_distrib = a_pool.grab_1d(key_tot.str());

    const snemo::datamodel::base_trajectory_pattern & a_track_pattern = the_particles.front().get().get_trajectory().get_pattern();
    double track_length = a_track_pattern.get_shape().get_length();

    if(datatools::is_valid(track_length))
      a_histo_distrib.fill(track_length);

    // if (a_pattern_id != "1e") {
    //   DT_LOG_WARNING(get_logging_priority(), "PlotModule only works for '1e' topology for now !");
    //   return dpp::base_module::PROCESS_SUCCESS;
    //   // return dpp::base_module::PROCESS_ERROR;
    // }

    std::ostringstream key_efficiency;
    key_efficiency << "track_length_efficiency";

    if (! a_pool.has(key_efficiency.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key_efficiency.str(), "", "track_length_efficiency");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d","track_length_efficiency_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_efficiency = a_pool.grab_1d(key_efficiency.str ());

    if(datatools::is_valid(track_length))
      a_histo_efficiency.fill(track_length);
*/

    // if (! a_pool.has(key.str()))
    //   {
    //     mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "vertices");
    //     datatools::properties hconfig;
    //     hconfig.store_string("mode", "mimic");
    //     hconfig.store_string("mimic.histogram_1d", "tof_probability_template");
    //     mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
    //   }

    // // Getting the current histogram
    // mygsl::histogram_1d & a_histo = a_pool.grab_1d(key.str ());
    // // a_histo.fill(a_pattern.get_measurement("vertex_e1_e2").get_probability());
    // // std:: cout << "vertices plot " << dynamic_cast<const snemo::datamodel::vertex_measurement&> (a_pattern.get_measurement("vertex_e1_e2")).get_probability() << std::endl;
    // double vertices_proba = dynamic_cast<const snemo::datamodel::vertex_measurement&> (a_pattern.get_measurement("vertex_e1_e2")).get_probability();
    // if(datatools::is_valid(vertices_proba))
    //    a_histo.fill(vertices_proba);

    // std::ostringstream key_delta_x;
    // key_delta_x << "vertices_distance_x";

    // if (! a_pool.has(key_delta_x.str()))
    //   {
    //     mygsl::histogram_1d & h = a_pool.add_1d(key_delta_x.str(), "", "vertices");
    //     datatools::properties hconfig;
    //     hconfig.store_string("mode", "mimic");
    //     hconfig.store_string("mimic.histogram_1d", "delta_vertices_Y_template");
    //     mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
    //   }

    // // Getting the current histogram
    // mygsl::histogram_1d & a_histo_delta_x = a_pool.grab_1d(key_delta_x.str ());
    // double delta_vertices_x = dynamic_cast<const snemo::datamodel::vertex_measurement&> (a_pattern.get_measurement("vertex_e1_e2")).get_vertices_distance_x();
    // if(datatools::is_valid(delta_vertices_x))
    //    a_histo_delta_x.fill(delta_vertices_x);

    // std::ostringstream key_delta_y;
    // key_delta_y << "vertices_distance_y";

    // if (! a_pool.has(key_delta_y.str()))
    //   {
    //     mygsl::histogram_1d & h = a_pool.add_1d(key_delta_y.str(), "", "vertices");
    //     datatools::properties hconfig;
    //     hconfig.store_string("mode", "mimic");
    //     hconfig.store_string("mimic.histogram_1d", "delta_vertices_Y_template");
    //     mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
    //   }

    // // Getting the current histogram
    // mygsl::histogram_1d & a_histo_delta_y = a_pool.grab_1d(key_delta_y.str ());
    // double delta_vertices_y = dynamic_cast<const snemo::datamodel::vertex_measurement&> (a_pattern.get_measurement("vertex_e1_e2")).get_vertices_distance_y();
    // if(datatools::is_valid(delta_vertices_y))
    //    a_histo_delta_y.fill(delta_vertices_y);

    // std::ostringstream key_delta_z;
    // key_delta_z << "vertices_distance_z";

    // if (! a_pool.has(key_delta_z.str()))
    //   {
    //     mygsl::histogram_1d & h = a_pool.add_1d(key_delta_z.str(), "", "vertices");
    //     datatools::properties hconfig;
    //     hconfig.store_string("mode", "mimic");
    //     hconfig.store_string("mimic.histogram_1d", "delta_vertices_Y_template");
    //     mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
    //   }

    // // Getting the current histogram
    // mygsl::histogram_1d & a_histo_delta_z = a_pool.grab_1d(key_delta_z.str ());
    // double delta_vertices_z = dynamic_cast<const snemo::datamodel::vertex_measurement&> (a_pattern.get_measurement("vertex_e1_e2")).get_vertices_distance_z();
    // if(datatools::is_valid(delta_vertices_z))
    //    a_histo_delta_z.fill(delta_vertices_z);

    /*
    //==========
    std::ostringstream key_vertex;
    key_vertex << "vertex_position";

    double min_track_length = 1000000000;

    for (snemo::datamodel::particle_track_data::particle_collection_type::const_iterator
           iparticle = ptd.get_particles().begin();
         iparticle != ptd.get_particles().end();
         ++iparticle)
      {
        const snemo::datamodel::particle_track & a_particle = iparticle->get();

        if (a_particle.has_trajectory()) {
          const snemo::datamodel::tracker_trajectory & a_trajectory = a_particle.get_trajectory();
          const snemo::datamodel::base_trajectory_pattern & a_track_pattern = a_trajectory.get_pattern();
          if(a_track_pattern.get_shape().get_length() < min_track_length)
          min_track_length = a_track_pattern.get_shape().get_length();
        } else {
          DT_LOG_WARNING(get_logging_priority(), "Particle has no attached trajectory !");
        }
      }
    // std::cout << "particle track length " << min_track_length << std::endl;
    if(min_track_length < 400)
      key_vertex << "_short";
    else if(min_track_length < 600)
      key_vertex << "_medium";
    else
      key_vertex << "_long";

    if (! a_pool.has(key_vertex.str()))
      {
        mygsl::histogram_2d & h = a_pool.add_2d(key_vertex.str(), "", "vertices");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_2d", "vertex_distribution_template");
        mygsl::histogram_pool::init_histo_2d(h, hconfig, &a_pool);
      }

    // // Getting the current histogram
    // mygsl::histogram_2d & a_histo = a_pool.grab_2d(key_vertex.str ());
    // double vertex_y = dynamic_cast<const snemo::datamodel::vertex_measurement&> (a_pattern.get_measurement("vertex_e1_e2")).get_vertex().get_position().y();
    // double vertex_z = dynamic_cast<const snemo::datamodel::vertex_measurement&> (a_pattern.get_measurement("vertex_e1_e2")).get_vertex().get_position().z();
    // if(datatools::is_valid(vertex_y) && datatools::is_valid(vertex_z))
    //   a_histo.fill(vertex_y,vertex_z);

    // Getting the current histogram
    mygsl::histogram_2d & a_histo = a_pool.grab_2d(key_vertex.str ());
    double vertex_y = dynamic_cast<const snemo::datamodel::vertex_measurement&> (a_pattern.get_measurement("vertex_e1")).get_vertex().get_position().y();
    double vertex_z = dynamic_cast<const snemo::datamodel::vertex_measurement&> (a_pattern.get_measurement("vertex_e1")).get_vertex().get_position().z();
    if(datatools::is_valid(vertex_y) && datatools::is_valid(vertex_z))
      a_histo.fill(vertex_y,vertex_z);
    std::cout << "vertex " << vertex_y << "  " << vertex_z << std::endl;
    //=======
*/
    DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return dpp::base_module::PROCESS_SUCCESS;
  }

} // namespace analysis

  // end of vertices_plot_module.cc
  /*
  ** Local Variables: --
  ** mode: c++ --
  ** c-file-style: "gnu" --
  ** tab-width: 2 --
  ** End: --
  */
