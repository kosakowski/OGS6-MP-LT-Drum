/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include <memory>
#include "TwoPhaseComponentialFlowMaterialProperties.h"
namespace BaseLib
{
    class ConfigTree;
}

namespace ProcessLib
{
    namespace TwoPhaseComponentialFlow
    {
        std::unique_ptr<TwoPhaseComponentialFlowMaterialProperties>
            createTwoPhaseComponentialFlowMaterialProperties(
                BaseLib::ConfigTree const& config,
                boost::optional<MeshLib::PropertyVector<int> const&> material_ids,
                std::vector<std::unique_ptr<ParameterBase>> const& parameters);

    }  // end namespace
}  // end namespace
