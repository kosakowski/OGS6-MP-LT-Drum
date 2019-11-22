/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include "TwoPhaseComponentialFlowMaterialProperties.h"
namespace MeshLib
{
    class Element;
}

namespace ProcessLib
{
    template <typename T>
    struct Parameter;

    namespace TwoPhaseComponentialFlow
    {
        struct TwoPhaseComponentialFlowProcessData
        {
            Eigen::VectorXd const _specific_body_force;

            bool const _has_gravity;
            bool const _has_mass_lumping;
            ParameterLib::Parameter<double> const& _diffusion_coeff_component_b;
            ParameterLib::Parameter<double> const& _diffusion_coeff_component_a;
            ParameterLib::Parameter<double> const& _diffusion_coeff_component_c;//diffusion coefficient of CO2 in gas
            ParameterLib::Parameter<double> const& _temperature;
            std::unique_ptr<TwoPhaseComponentialFlowMaterialProperties> _material;
            MathLib::PiecewiseLinearInterpolation const& _interpolated_Q_slow;
            MathLib::PiecewiseLinearInterpolation const& _interpolated_Q_fast;
            MathLib::PiecewiseLinearInterpolation const& _interpolated_kinetic_rate;
            double _dt = 0;

            // mesh properties for scalar output
            MeshLib::PropertyVector<double>* mesh_prop_saturation = nullptr;
            MeshLib::PropertyVector<double>* mesh_prop_bazant_power = nullptr;
            MeshLib::PropertyVector<double>* mesh_prop_mol_density_gas = nullptr;
            MeshLib::PropertyVector<double>* mesh_prop_mol_density_liquid = nullptr;
            MeshLib::PropertyVector<double>* mesh_prop_porosity = nullptr;
            MeshLib::PropertyVector<double>* mesh_prop_pHvalue = nullptr;
            MeshLib::PropertyVector<double>* mesh_prop_co2_cumulate_consume = nullptr;
            MeshLib::PropertyVector<double>* mesh_prop_sio2_cumulate_consume = nullptr;
            MeshLib::PropertyVector<double>* mesh_prop_mol_frac_h2o_vapor = nullptr;
            MeshLib::PropertyVector<double>* mesh_prop_mol_frac_n2 = nullptr;

            //mesh properties for velocity output
            MeshLib::PropertyVector<double>*  mesh_prop_overall_liquid_darcy_volumetric_flux = nullptr;
            MeshLib::PropertyVector<double>*  mesh_prop_overall_gas_darcy_volumetric_flux = nullptr;
            MeshLib::PropertyVector<double>*  mesh_prop_gas_co2_darcy_volumetric_flux = nullptr;
            MeshLib::PropertyVector<double>*  mesh_prop_gas_co2_diffusive_volumetric_flux = nullptr;
            MeshLib::PropertyVector<double>*  mesh_prop_gas_hydrogen_darcy_volumetric_flux = nullptr;
            MeshLib::PropertyVector<double>*  mesh_prop_gas_hydrogen_diffusive_volumetric_flux = nullptr;
            MeshLib::PropertyVector<double>*  mesh_prop_gas_methane_darcy_volumetric_flux = nullptr;
            MeshLib::PropertyVector<double>*  mesh_prop_gas_methane_diffusive_volumetric_flux = nullptr;
            MeshLib::PropertyVector<double>*  mesh_prop_gas_water_vapor_darcy_volumetric_flux = nullptr;
            MeshLib::PropertyVector<double>*  mesh_prop_gas_water_vapor_diffusive_volumetric_flux = nullptr;
            MeshLib::PropertyVector<double>*  mesh_prop_gas_nitrogen_darcy_volumetric_flux = nullptr;
            MeshLib::PropertyVector<double>*  mesh_prop_gas_nitrogen_diffusive_volumetric_flux = nullptr;

        };

    }  // namespace TwoPhaseComponentialFlow
}  // namespace ProcessLib
