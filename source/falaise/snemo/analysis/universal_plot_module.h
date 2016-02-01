/* universal_plot_module.h
 * Author(s)     : Steven Calvez <calvez@lal.in2p3.fr>
 * Creation date : 2015-05-12
 * Last modified : 2014-05-12
 *
 * Copyright (C) 2014 Steven Calvez <calvez@lal.in2p3.fr>
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
 * A module which generates some vertices plots.
 *
 * History:
 *
 */

#ifndef ANALYSIS_UNIVERSAL_PLOT_MODULE_H_
#define ANALYSIS_UNIVERSAL_PLOT_MODULE_H_ 1

// Standard libraires:
#include <set>
#include <map>

// Data processing module abstract base class
#include <dpp/base_module.h>

namespace mygsl {
  class histogram_pool;
}

namespace analysis {

  class universal_plot_module : public dpp::base_module
  {
  public:

    /// Setting histogram pool
    void set_histogram_pool(mygsl::histogram_pool & pool_);

    /// Grabbing histogram pool
    mygsl::histogram_pool & grab_histogram_pool();

    /// Constructor
    universal_plot_module(datatools::logger::priority = datatools::logger::PRIO_FATAL);

    /// Destructor
    virtual ~universal_plot_module();

    /// Initialization
    virtual void initialize(const datatools::properties  & setup_,
                            datatools::service_manager   & service_manager_,
                            dpp::module_handle_dict_type & module_dict_);

    /// Reset
    virtual void reset();

    /// Data record processing
    virtual process_status process(datatools::things & data_);

  protected:

    /// Give default values to specific class members.
    void _set_defaults();

  private:

    // The key fields from 'event header' bank to build the histogram key:
    std::vector<std::string> _key_fields_;

    // The histogram pool :
    mygsl::histogram_pool * _histogram_pool_;

    // Macro to automate the registration of the module :
    DPP_MODULE_REGISTRATION_INTERFACE(universal_plot_module);
  };

} // namespace analysis

#endif // ANALYSIS_UNIVERSAL_PLOT_MODULE_H_

// end of universal_plot_module.h
/*
** Local Variables: --
** mode: c++ --
** c-file-style: "gnu" --
** tab-width: 2 --
** End: --
*/
