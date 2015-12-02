/* halflife_limit_module.h
 * Author(s)     : Steven Calvez <calvez@lal.in2p3.fr>
 * Creation date : 2012-06-13
 * Last modified : 2012-08-27
 *
 * Copyright (C) 2012-2014 Steven Calvez <calvez@lal.in2p3.fr>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 *
 * Description:
 *
 * A dedicated module to compute the limit on the neutrinoless double
 * beta decay halflife.
 *
 * History:
 *
 */

#ifndef ANALYSIS_HALFLIFE_LIMIT_MODULE_H
#define ANALYSIS_HALFLIFE_LIMIT_MODULE_H 1

// Data processing module abstract base class
#include <dpp/base_module.h>

#include <map>
#include <string>
#include <vector>

namespace mygsl {
  class histogram_pool;
}

namespace analysis {

  class halflife_limit_module : public dpp::base_module
  {
  public:

    struct experiment_entry_type
    {
      typedef std::map<std::string, double> background_dict_type;

      double isotope_mass_number;
      double isotope_mass;
      double isotope_bb2nu_halflife;
      double exposure_time;
      background_dict_type background_activities;

      void initialize(const datatools::properties & config_);
    };

  public:

    /// Return the label of for the boolean property associated to 'signal'
    static const std::string & signal_flag();

    /// Return the label of for the boolean property associated to 'background'
    static const std::string & background_flag();

    /// Setting histogram pool
    void set_histogram_pool(mygsl::histogram_pool & pool_);

    /// Grabbing histogram pool
    mygsl::histogram_pool & grab_histogram_pool();

    /// Constructor
    halflife_limit_module(datatools::logger::priority = datatools::logger::PRIO_FATAL);

    /// Destructor
    virtual ~halflife_limit_module();

    /// Initialization
    virtual void initialize(const datatools::properties  & config_,
                            datatools::service_manager   & service_manager_,
                            dpp::module_handle_dict_type & module_dict_);

    /// Reset
    virtual void reset();

    /// Data record processing
    virtual process_status process(datatools::things & data_);

    /// Dump
    void dump_result (std::ostream      & out_    = std::clog,
                      const std::string & title_  = "",
                      const std::string & indent_ = "",
                      bool inherit_               = false) const;

  protected:

    /// Give default values to specific class members.
    void _set_defaults();

    /// Compute topology channel efficiencies.
    void _compute_efficiency();

    /// Compute topology channel efficiencies.
    void _compute_halflife();

  private:

    // The key fields from 'event header' bank to build the histogram key:
    std::vector<std::string> _key_fields_;

    // The histogram pool :
    mygsl::histogram_pool * _histogram_pool_;

    // The experiment running condition
    experiment_entry_type _experiment_conditions_;

    // Macro to automate the registration of the module :
    DPP_MODULE_REGISTRATION_INTERFACE (halflife_limit_module);

  };

} // namespace analysis

#endif // ANALYSIS_HALFLIFE_LIMIT_MODULE_H

// end of halflife_limit_module.h
/*
** Local Variables: --
** mode: c++ --
** c-file-style: "gnu" --
** tab-width: 2 --
** End: --
*/
