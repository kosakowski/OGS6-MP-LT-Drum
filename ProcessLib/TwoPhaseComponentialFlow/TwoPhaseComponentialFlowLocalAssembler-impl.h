/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include "TwoPhaseComponentialFlowLocalAssembler.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "TwoPhaseComponentialFlowProcessData.h"

using MaterialLib::PhysicalConstant::MolarMass::Water;
using MaterialLib::PhysicalConstant::MolarMass::H2;
using MaterialLib::PhysicalConstant::MolarMass::Air;
using MaterialLib::PhysicalConstant::MolarMass::CH4;
using MaterialLib::PhysicalConstant::MolarMass::CO2;

namespace ProcessLib
{
    namespace TwoPhaseComponentialFlow
    {
        template <typename ShapeFunction, typename IntegrationMethod,
            unsigned GlobalDim>
            void TwoPhaseComponentialFlowLocalAssembler<
            ShapeFunction, IntegrationMethod,
            GlobalDim>::assemble(double const t, std::vector<double> const& local_x,
                std::vector<double>& local_M_data,
                std::vector<double>& local_K_data,
                std::vector<double>& local_b_data)
        {
            auto const local_matrix_size = local_x.size();
            auto const n_nodes = ShapeFunction::NPOINTS;
            assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

            auto const p_nodal_values = Eigen::Map<const NodalVectorType>(
                local_x.data() + nonwet_pressure_matrix_index, ShapeFunction::NPOINTS);
            auto const pc_nodal_values = Eigen::Map<const NodalVectorType>(
                local_x.data() + cap_pressure_matrix_index, ShapeFunction::NPOINTS);
            auto const x1_nodal_values = Eigen::Map<const NodalVectorType>(
                local_x.data() + mol_fraction_h2_matrix_index, ShapeFunction::NPOINTS);
            auto const x2_nodal_values = Eigen::Map<const NodalVectorType>(
                local_x.data() + mol_fraction_ch4_matrix_index, ShapeFunction::NPOINTS);
            auto const x3_nodal_values = Eigen::Map<const NodalVectorType>(
                local_x.data() + mol_fraction_co2_matrix_index, ShapeFunction::NPOINTS);

            auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
                local_M_data, local_matrix_size, local_matrix_size);
            auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
                local_K_data, local_matrix_size, local_matrix_size);
            auto local_b = MathLib::createZeroedVector<LocalVectorType>(
                local_b_data, local_matrix_size);

            Eigen::MatrixXd mass_mat_coeff =
                Eigen::MatrixXd::Zero(NUM_NODAL_DOF, NUM_NODAL_DOF);
            Eigen::MatrixXd K_mat_coeff =
                Eigen::MatrixXd::Zero(NUM_NODAL_DOF, NUM_NODAL_DOF);
            Eigen::VectorXd H_vec_coeff = Eigen::VectorXd::Zero(NUM_NODAL_DOF);
            //Eigen::VectorXd F_vec_coeff = Eigen::VectorXd::Zero(NUM_NODAL_DOF);

            NodalMatrixType localMass_tmp;
            localMass_tmp.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
            NodalMatrixType localDispersion_tmp;
            localDispersion_tmp.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
            Eigen::VectorXd localGravity_tmp =
                Eigen::VectorXd::Zero(ShapeFunction::NPOINTS);
            Eigen::VectorXd localSource_tmp =
                Eigen::VectorXd::Zero(ShapeFunction::NPOINTS);

            NodalVectorType localNeumann_tmp;
            localNeumann_tmp.setZero(ShapeFunction::NPOINTS);

            NodalMatrixType mass_operator;
            mass_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

            NodalMatrixType laplace_operator;
            laplace_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);
            NodalMatrixType diffusive_operator;
            diffusive_operator.setZero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

            unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();

            SpatialPosition pos;
            pos.setElementID(_element.getID());

            auto nodes = _element.getNodes();
            auto rx0 = (*nodes[0])[0];
            auto rx1 = (*nodes[1])[0];
            auto rx2 = (*nodes[2])[0];
            auto rx3 = (*nodes[3])[0];
            auto ry0 = (*nodes[0])[1];
            auto ry1 = (*nodes[1])[1];
            auto ry2 = (*nodes[2])[1];
            auto ry3 = (*nodes[3])[1];

            auto area = _element.getContent(); // returns length, area or volume..depending on element dimension...for the drum we have a 2D mesh -> area
             // radial symmetry 2 * 3.1415926*rx ...calculate volume of element according to Guidinsche Regel V = 2 * pi * R * area R: distance of center of gravity of area to rotation axis
            double Rdummy = _element.getCenterOfGravity()[0]; // for us the x coordinate is the one we need
            double node_volume_radial = 2.0 * 3.141592654 * Rdummy * area / _element.getNumberOfEdges(); // aprox the node volume covered by a node..this is not exact due to radial symmetry (I guess)

            const int material_id =
                _process_data._material->getMaterialID(pos.getElementID().get());

            const Eigen::MatrixXd& perm = _process_data._material->getPermeability(
                material_id, t, pos, _element.getDimension());
            assert(perm.rows() == _element.getDimension() || perm.rows() == 1);
            GlobalDimMatrixType permeability = GlobalDimMatrixType::Zero(
                _element.getDimension(), _element.getDimension());
            if (perm.rows() == _element.getDimension())
                permeability = perm;
            else if (perm.rows() == 1)
                permeability.diagonal().setConstant(perm(0, 0));
            MathLib::PiecewiseLinearInterpolation const& interpolated_Q_slow =
                _process_data._interpolated_Q_slow;
            MathLib::PiecewiseLinearInterpolation const& interpolated_Q_fast =
                _process_data._interpolated_Q_fast;
            MathLib::PiecewiseLinearInterpolation const& interpolated_kinetic_rate =
                _process_data._interpolated_kinetic_rate;

            _overall_velocity_gas.clear();
            _overall_velocity_liquid.clear();
            _gas_co2_velocity.clear();
            _gas_hydrogen_velocity.clear();
            _gas_methane_velocity.clear();
            _gas_nitrogen_velocity.clear();
            _gas_vapor_velocity.clear();

            auto cache_mat_gas_vel = MathLib::createZeroedMatrix<
                Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
                    _overall_velocity_gas, GlobalDim, n_integration_points);

            auto cache_mat_liquid_vel = MathLib::createZeroedMatrix<
                Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
                    _overall_velocity_liquid, GlobalDim, n_integration_points);

            auto cache_mat_gas_co2_vel = MathLib::createZeroedMatrix<
                Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
                    _gas_co2_velocity, GlobalDim, n_integration_points);

            auto cache_mat_gas_hydrogen_vel = MathLib::createZeroedMatrix<
                Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
                    _gas_hydrogen_velocity, GlobalDim, n_integration_points);

            auto cache_mat_gas_methane_vel = MathLib::createZeroedMatrix<
                Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
                    _gas_methane_velocity, GlobalDim, n_integration_points);

            auto cache_mat_gas_nitrogen_vel = MathLib::createZeroedMatrix<
                Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
                    _gas_nitrogen_velocity, GlobalDim, n_integration_points);

            auto cache_mat_gas_water_vapor_vel = MathLib::createZeroedMatrix<
                Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
                    _gas_vapor_velocity, GlobalDim, n_integration_points);

            Eigen::VectorXd neumann_vector = Eigen::VectorXd::Zero(local_x.size());
            auto neumann_vec = neumann_vector.segment<ShapeFunction::NPOINTS>(
                0);

            Eigen::VectorXd _neumann_vector_output = Eigen::VectorXd::Zero(local_x.size());
            auto _neumann_vec_output = _neumann_vector_output.segment<ShapeFunction::NPOINTS>(
                0);
            Eigen::VectorXd F_vec_coeff = Eigen::VectorXd::Zero(NUM_NODAL_DOF);

            accelerate_flag = false;
            int gp_carb_neutral_count = 3;
            bool atm_flag = false;
            double deriv_flag = 1;
            double accelerate_factor = 1;
            if (_process_data._material->getMaterialID(pos.getElementID().get()) >=
                2)//backfill
            {
                atm_flag = true;
            }
            if (pos.getElementID().get() == 148 || pos.getElementID().get() == 149 ||
                pos.getElementID().get() == 152 || pos.getElementID().get() == 153 ||
                pos.getElementID().get() == 160 || pos.getElementID().get() == 161 ||
                pos.getElementID().get() == 186 || pos.getElementID().get() == 189 ||
                pos.getElementID().get() == 190 || pos.getElementID().get() == 191 ||
                pos.getElementID().get() == 246 || pos.getElementID().get() == 236)
            {
                accelerate_factor = 2.0;
            }
            if (t > 65)
            {
                accelerate_factor = 1;
            }
            for (unsigned ip = 0; ip < n_integration_points; ip++)
            {
                F_vec_coeff.setZero(NUM_NODAL_DOF);
                auto const& sm = _shape_matrices[ip];

                double pg_int_pt = 0.0;
                double X1_int_pt = 0.0;
                double X2_int_pt = 0.0;
                double X3_int_pt = 0.0;
                double PC_int_pt = 0.0;

                double gas_h2_generation_rate = 0.0;
                double gas_ch4_generation_rate = 0.0;

                NumLib::shapeFunctionInterpolate(local_x, sm.N, pg_int_pt, X1_int_pt,
                    X2_int_pt, X3_int_pt, PC_int_pt);

                const auto _interpolateGaussNode_coord = interpolateNodeCoordinates(
                    _element, sm.N);

                _pressure_wetting[ip] = pg_int_pt - PC_int_pt;
                const double dt = _process_data._dt;
                if (atm_flag&&t > 40&& t<50)
                    X3_int_pt = 1e-4;
                else if(atm_flag && t>50 && t<60)
                    X3_int_pt = 1e-5;
                else if (atm_flag && t>60)
                    X3_int_pt = 1e-6;
                auto const& wp = _integration_method.getWeightedPoint(ip);
                double integration_factor =
                    sm.integralMeasure * sm.detJ * wp.getWeight();
                if (atm_flag)
                    integration_factor = 1e+60;
                const double temperature = _process_data._temperature(t, pos)[0];

                //store the integration value for pressure and molar fraction
                double& pg_ip_data = _ip_data[ip].pressure_cur;
                pg_ip_data = pg_int_pt;
                double& mol_frac_h2_gas_ip_data = _ip_data[ip].mol_frac_h2_cur;
                mol_frac_h2_gas_ip_data = X1_int_pt;
                double const gas_generation_rate = (pg_int_pt - _ip_data[ip].pressure_pre) / dt / R / temperature;
                _gas_generation_rate[ip] = gas_generation_rate;
                double const partial_h2_gas_generation_rate
                    = (pg_int_pt - _ip_data[ip].pressure_pre)*X1_int_pt / dt / R / temperature
                    + pg_int_pt*(X1_int_pt - _ip_data[ip].mol_frac_h2_pre) / dt / R / temperature;
                double X_L_h_gp = pg_int_pt * X1_int_pt / Hen_L_h;  // Henry law
                double X_L_c_gp = pg_int_pt * X2_int_pt / Hen_L_c;
                double X_L_co2_gp = pg_int_pt * X3_int_pt / Hen_L_co2;

                double P_sat_gp = get_P_sat(temperature);

                double K_G_w = pg_int_pt / P_sat_gp;  // henry law ratio
                double K_G_air = pg_int_pt / Hen_L_air;
                double L = 1 - (X_L_h_gp + X_L_c_gp + X_L_co2_gp);
                double G = 1 - X1_int_pt - X2_int_pt - X3_int_pt;
                const double rho_mol_water = rho_l_std / Water;
                const double kelvin_term =
                    exp(PC_int_pt / rho_mol_water / R / temperature);
                const double d_kelvin_term_d_pc =
                    kelvin_term / rho_mol_water / R / temperature;

                double x_nonwet_h2o = get_x_nonwet_h2o(
                    pg_int_pt, X1_int_pt, X2_int_pt, X3_int_pt, P_sat_gp, kelvin_term);
                if (atm_flag)
                    x_nonwet_h2o = P_sat_gp*0.8 / 101325;
                //store the secondary variable of nonwet vapor molar fraction
                _mol_fraction_nonwet_vapor[ip] = x_nonwet_h2o;
                double const x_nonwet_air = (G - x_nonwet_h2o);
                //store the secondary variable of nonwet air molar fraction
                _mol_fraction_nonwet_air[ip] = x_nonwet_air;
                double const x_wet_air = pg_int_pt * x_nonwet_air / Hen_L_air;
                double const x_wet_h2o = 1 - X_L_co2_gp - X_L_c_gp - X_L_h_gp - x_wet_air;
                //pg_int_pt * x_nonwet_h2o * kelvin_term / P_sat_gp;
                /*double const rho_gas =
                _process_data._material->getGasDensity(pg_int_pt, temperature);*/
                double const rho_w = _process_data._material->getLiquidDensity(
                    _pressure_wetting[ip], temperature);

                double Sw = _process_data._material->getSaturation(
                    material_id, t, pos, pg_int_pt, temperature, PC_int_pt);

                _saturation[ip] = Sw;//store the secondary variable
                double const S_G_gp = 1 - Sw;

                double dSwdPc = _process_data._material->getDerivSaturation(
                    material_id, t, pos, pg_int_pt, temperature, Sw);

                const double dSgdPC = -dSwdPc;

                const double rho_mol_nonwet = pg_int_pt / R / temperature;

                const double rho_mass_G_gp =
                    rho_mol_nonwet *
                    (X1_int_pt * H2 + X2_int_pt * CH4 + X3_int_pt * CO2 +
                        x_nonwet_air * Air + x_nonwet_h2o * Water);

                double dLdPG =
                    -X1_int_pt / Hen_L_h - X2_int_pt / Hen_L_c - X3_int_pt / Hen_L_co2;

                double d_x_nonwet_h2o_d_pg = deriv_flag*get_derivative_x_nonwet_h2o_d_pg(
                    pg_int_pt, X1_int_pt, X2_int_pt, X3_int_pt, P_sat_gp, kelvin_term);
                /*double d_x_nonwet_h2o_d_pg_test = (get_x_nonwet_h2o(
                pg_int_pt + 1e-6, X1_int_pt, X2_int_pt, X3_int_pt, P_sat_gp,
                kelvin_term) - get_x_nonwet_h2o(
                pg_int_pt - 1e-6, X1_int_pt, X2_int_pt, X3_int_pt, P_sat_gp,
                kelvin_term)) / 2 / 1e-6;*/

                double const d_x_nonwet_h2o_d_x1 = deriv_flag* get_derivative_x_nonwet_h2o_d_x1(
                    pg_int_pt, X1_int_pt, X2_int_pt, X3_int_pt, P_sat_gp, kelvin_term);
                /*double d_x_nonwet_h2o_d_x1_test = (get_x_nonwet_h2o(
                pg_int_pt, X1_int_pt + 1e-6, X2_int_pt, X3_int_pt, P_sat_gp,
                kelvin_term) - get_x_nonwet_h2o(
                pg_int_pt, X1_int_pt - 1e-6, X2_int_pt, X3_int_pt, P_sat_gp,
                kelvin_term)) / 2 / 1e-6;*/

                double d_x_nonwet_h2o_d_x2 = deriv_flag * get_derivative_x_nonwet_h2o_d_x2(
                    pg_int_pt, X1_int_pt, X2_int_pt, X3_int_pt, P_sat_gp, kelvin_term);
                /*double d_x_nonwet_h2o_d_x2_test = (get_x_nonwet_h2o(
                pg_int_pt, X1_int_pt, X2_int_pt + 1e-6, X3_int_pt, P_sat_gp,
                kelvin_term) - get_x_nonwet_h2o(
                pg_int_pt, X1_int_pt, X2_int_pt - 1e-6, X3_int_pt, P_sat_gp,
                kelvin_term)) / 2 / 1e-6;*/

                double d_x_nonwet_h2o_d_x3 = deriv_flag*get_derivative_x_nonwet_h2o_d_x3(
                    pg_int_pt, X1_int_pt, X2_int_pt, X3_int_pt, P_sat_gp, kelvin_term);
                /*double d_x_nonwet_h2o_d_x3_test = (get_x_nonwet_h2o(
                pg_int_pt, X1_int_pt, X2_int_pt, X3_int_pt + 1e-6, P_sat_gp,
                kelvin_term) - get_x_nonwet_h2o(
                pg_int_pt, X1_int_pt, X2_int_pt, X3_int_pt - 1e-6, P_sat_gp,
                kelvin_term)) / 2 / 1e-6;*/
                double const d_x_nonwet_h2o_d_kelvin =
                    get_derivative_x_nonwet_h2o_d_kelvin(pg_int_pt, X1_int_pt,
                        X2_int_pt, X3_int_pt, P_sat_gp,
                        kelvin_term);
                /*double const d_x_nonwet_h2o_d_kelvin_test = (get_x_nonwet_h2o(
                pg_int_pt, X1_int_pt, X2_int_pt, X3_int_pt, P_sat_gp, kelvin_term +
                1e-6) - get_x_nonwet_h2o(
                pg_int_pt, X1_int_pt, X2_int_pt, X3_int_pt, P_sat_gp,
                kelvin_term - 1e-6)) / 2 / 1e-6;*/
                double const d_x_nonwet_h2o_d_pc = deriv_flag*
                    d_x_nonwet_h2o_d_kelvin * d_kelvin_term_d_pc;

                double const d_x_nonwet_air_d_pg = -d_x_nonwet_h2o_d_pg;
                double const d_x_nonwet_air_d_x1 = -1 - d_x_nonwet_h2o_d_x1;
                double const d_x_nonwet_air_d_x2 = -1 - d_x_nonwet_h2o_d_x2;
                double const d_x_nonwet_air_d_x3 = -1 - d_x_nonwet_h2o_d_x3;
                double const d_x_nonwet_air_d_pc = -d_x_nonwet_h2o_d_pc;

                /// double const d_x_nonwet_h2o_d_x3 =
                /// get_derivative_x_nonwet_h2o_d_x3(pg_int_pt, X1_int_pt, X2_int_pt,
                /// X3_int_pt, P_sat_gp);
                double porosity = _process_data._material->getPorosity(
                    material_id, t, pos, pg_int_pt, temperature, 0);
                if (atm_flag)
                    porosity = 0.5;

                // Assemble M matrix
                // nonwetting
                double const d_rho_mol_nonwet_d_pg = 1 / R / temperature;

                double const rho_mol_wet = rho_mol_water / x_wet_h2o;
                const double rho_mass_wet =
                    rho_mol_wet * (X_L_h_gp * H2 + X_L_c_gp * CH4 + X_L_co2_gp * CO2 +
                        x_wet_air * Air + x_wet_h2o * Water);
                double const d_rho_mol_wet_d_pg =
                    -rho_mol_water * kelvin_term *
                    (x_nonwet_h2o / P_sat_gp +
                        pg_int_pt * d_x_nonwet_h2o_d_pg / P_sat_gp) /
                    x_wet_h2o / x_wet_h2o;
                double const d_rho_mol_wet_d_x1 =
                    -rho_mol_water * kelvin_term *
                    (pg_int_pt * d_x_nonwet_h2o_d_x1 / P_sat_gp) / x_wet_h2o /
                    x_wet_h2o;
                double const d_rho_mol_wet_d_x2 =
                    -rho_mol_water * kelvin_term *
                    (pg_int_pt * d_x_nonwet_h2o_d_x2 / P_sat_gp) / x_wet_h2o /
                    x_wet_h2o;
                double const d_rho_mol_wet_d_x3 =
                    -rho_mol_water * kelvin_term *
                    (pg_int_pt * d_x_nonwet_h2o_d_x3 / P_sat_gp) / x_wet_h2o /
                    x_wet_h2o;
                double const d_rho_mol_wet_d_pc =
                    -rho_mol_water *
                    (pg_int_pt * d_x_nonwet_h2o_d_pc * kelvin_term / P_sat_gp +
                        pg_int_pt * x_nonwet_h2o * d_kelvin_term_d_pc / P_sat_gp) /
                    x_wet_h2o / x_wet_h2o;

                // calculate the carbonation and ASR source/sink term
                // For the backfill part
                double& rho_mol_sio2_wet_backfill = _ip_data[ip].rho_mol_sio2_backfill;
                double& rho_mol_co2_cumul_total_backfill =
                    _ip_data[ip].rho_mol_co2_cumul_total_backfill;  // get cumulative co2
                double& porosity2 = _ip_data[ip].porosity_backfill;
                double rho_mol_total_co2_backfill = 0.;
                double rho_mol_co2_kinetic_rate_backfill = 0.;
                const double rho_co2_max_consume = 8500;
                // for the waste part
                double& porosity3 = _ip_data[ip].porosity_waste;
                double& rho_mol_co2_cumul_total_waste = _ip_data[ip].rho_mol_co2_cumul_total_waste;
                double rho_mol_total_co2_waste = 0.;
//                double M_organic_fast_co2_ini =  _amount_organic_waste_cellulose[ip]; // amount from which gas is produced
//                double M_organic_slow_co2_ini =  _amount_organic_waste_polystyrene[ip]; // amount from which gas is produced
                double& M_organic_fast_co2_ini = _ip_data[ip].amount_organic_waste_prev_cellulose; // amount from which gas is produced
                double& M_organic_slow_co2_ini = _ip_data[ip].amount_organic_waste_prev_polystyrene; // amount from which gas is produced


                //saturation dependent chemical reactivity
                //double const rel_humidity = std::exp(-PC_int_pt*0.018 / rho_mass_wet / 8.314 / temperature);
                double pc_origion = _process_data._material->getCapillaryPressure(material_id, t, pos, pg_int_pt, temperature, Sw);
                double rel_humidity = std::exp(-pc_origion*0.018 / rho_mass_wet / 8.314 / temperature);
                //rel_humidity = pg_int_pt*x_nonwet_h2o / P_sat_gp;
                if (atm_flag)
                    rel_humidity = 0.8;
                _rel_humidity[ip] = rel_humidity;
                double bazant_power = 0.0;

                if (_process_data._material->getMaterialID(pos.getElementID().get()) ==
                    0)//backfill
                {
                    //bazant_power = pow((1 + pow((7.5 - 7.5*rel_humidity), 4)), -1);
                    bazant_power = (rel_humidity - 0.6) / (1 - 0.6);
                    //bazant_power= 5 * rel_humidity - 4;
                    if (rel_humidity<0.6 )
                        bazant_power = 1e-7;
                    if (bazant_power<0.0|| Sw<0.2)
                        bazant_power = 1e-7;
                    // should be only valid for material 0
                    porosity = porosity2;  // use look-up table value
                                           // calculate the current ammount of co2
                                           // calculate the total amount of co2 in this element(gp), which should
                                           // be consumed at this time step
                    rho_mol_total_co2_backfill = _ip_data[ip].porosity_prev_backfill *
                        (rho_mol_nonwet * X3_int_pt * (1 - Sw) +
                            rho_mol_wet * X_L_co2_gp * 0.0);//assume only gas co2 will carbonate
                    double rho_mol_co2_consume_rate_backfill = bazant_power*rho_mol_total_co2_backfill / dt;

                    rho_mol_co2_consume_rate_backfill =
                        (rho_mol_co2_consume_rate_backfill < 0) ? 0.0 : rho_mol_co2_consume_rate_backfill;
                    if (rho_mol_co2_consume_rate_backfill < 0.0)
                        rho_mol_co2_consume_rate_backfill = 0.0;
                    //double const dcarb_rate=
                    //interpolated_kinetic_rate.getValue(rho_mol_co2_cumul_total_backfill*100 / rho_co2_max_consume);
                    //double const dcarb_rate_analytic
                    //= 0.04*(94.32 - rho_mol_co2_cumul_total_backfill*100 / rho_co2_max_consume);
                    rho_mol_co2_kinetic_rate_backfill = rho_mol_co2_consume_rate_backfill;
                    //bazant_power*rho_co2_max_consume*dcarb_rate/100;
                    // impose a max rate bound
                    /*if (rho_mol_co2_kinetic_rate_backfill > rho_mol_co2_consume_rate_backfill)
                    //rho_mol_co2_kinetic_rate_backfill = rho_mol_co2_consume_rate_backfill;
                    rho_mol_co2_kinetic_rate_backfill =
                    (rho_mol_co2_consume_rate_backfill < 0) ? 0.0 : rho_mol_co2_consume_rate_backfill;
                    else
                    rho_mol_co2_kinetic_rate_backfill = rho_mol_co2_kinetic_rate_backfill;*/
                    // get the pH value at this iteration based cumulated dissovled quatz and
                    // cumulated co2
                    double const pH = bi_interpolation(
                        _ip_data[ip].rho_mol_sio2_prev_backfill,
                        _ip_data[ip].rho_mol_co2_cumul_total_prev_backfill, _pH_at_supp_pnt_backfill);
                    if (pH < 10.5)
                        accelerate_flag = true;
                    _pH_value[ip] = pH;//update the secondary variable
                }
                else if (_process_data._material->getMaterialID(pos.getElementID().get()) ==
                    1)//waste matrix
                {
                    // calculate the current ammount of co2
                    // calculate the total amount of co2 in this element(gp), which should
                    // be consumed at this time step
                    bazant_power = (rel_humidity - 0.6) / (1 - 0.6);
                    //bazant_power = 5 * rel_humidity - 4;//for the inner part
                    if (rel_humidity<0.6)
                        bazant_power = 1e-7;
                    if (bazant_power<0.0|| Sw<0.2)
                        bazant_power = 1e-7;
                    rho_mol_total_co2_waste = _ip_data[ip].porosity_prev_waste *
                        (rho_mol_nonwet * X3_int_pt * (1 - Sw) +
                            rho_mol_wet * X_L_co2_gp * Sw);
                    porosity = porosity3;
                    double const pH = 10.0653;// piecewiselinear_interpolation(
                        //_ip_data[ip].rho_mol_co2_cumul_total_prev_waste, _pH_at_supp_pnt_waste);
                    if (pH < 10.5 )
                        accelerate_flag = true;
                    _pH_value[ip] = pH;
                }
                if (bazant_power > 1)
                    bazant_power = 1;
                _reactivity_bazant_power[ip] = bazant_power;
                //consumed CO2 for current step
                _co2_consumed_current_step[ip] = rho_mol_co2_kinetic_rate_backfill*dt;
                //calculate the co2 concentration
                _co2_concentration[ip] = rho_mol_nonwet*X3_int_pt + rho_mol_wet*X_L_co2_gp;
                // store the molar density of gas phase
                _rho_mol_gas_phase[ip] = rho_mol_nonwet;
                // store the molar density of liquid phase
                _rho_mol_liquid_phase[ip] = rho_mol_wet;
                //Assembly
                // H2
                mass_mat_coeff(0, 0) =
                    porosity * ((1 - Sw) * X1_int_pt * d_rho_mol_nonwet_d_pg +
                        Sw * rho_mol_wet * X1_int_pt / Hen_L_h +
                        Sw * X_L_h_gp * d_rho_mol_wet_d_pg);
                mass_mat_coeff(0, 1) =
                    porosity * (rho_mol_nonwet * (1 - Sw) +
                        Sw * rho_mol_wet * pg_int_pt / Hen_L_h +
                        Sw * X_L_h_gp * d_rho_mol_wet_d_x1);
                mass_mat_coeff(0, 2) = porosity * (Sw * X_L_h_gp * d_rho_mol_wet_d_x2);
                mass_mat_coeff(0, 3) = porosity * (Sw * X_L_h_gp * d_rho_mol_wet_d_x3);
                mass_mat_coeff(0, 4) = porosity * (rho_mol_nonwet * X1_int_pt * dSgdPC -
                    rho_mol_wet * X_L_h_gp * dSgdPC +
                    Sw * d_rho_mol_wet_d_pc * X_L_h_gp);
                // CH4
                mass_mat_coeff(1, 0) =
                    porosity * ((1 - Sw) * X2_int_pt * d_rho_mol_nonwet_d_pg +
                        Sw * rho_mol_wet * X2_int_pt / Hen_L_c +
                        Sw * X_L_c_gp * d_rho_mol_wet_d_pg);
                mass_mat_coeff(1, 1) = porosity * (Sw * X_L_c_gp * d_rho_mol_wet_d_x1);
                mass_mat_coeff(1, 2) =
                    porosity * (rho_mol_nonwet * (1 - Sw) +
                        Sw * rho_mol_wet * pg_int_pt / Hen_L_c +
                        Sw * X_L_c_gp * d_rho_mol_wet_d_x2);
                mass_mat_coeff(1, 3) = porosity * (Sw * X_L_c_gp * d_rho_mol_wet_d_x3);

                mass_mat_coeff(1, 4) = porosity * (rho_mol_nonwet * X2_int_pt * dSgdPC -
                    rho_mol_wet * X_L_c_gp * dSgdPC +
                    Sw * d_rho_mol_wet_d_pc * X_L_c_gp);
                // co2
                mass_mat_coeff(2, 0) =
                    porosity * ((1 - Sw) * X3_int_pt * d_rho_mol_nonwet_d_pg +
                        Sw * rho_mol_wet * X3_int_pt / Hen_L_co2 +
                        Sw * X_L_co2_gp * d_rho_mol_wet_d_pg);
                mass_mat_coeff(2, 1) =
                    porosity * (Sw * X_L_co2_gp * d_rho_mol_wet_d_x1);
                mass_mat_coeff(2, 2) =
                    porosity * (Sw * X_L_co2_gp * d_rho_mol_wet_d_x2);
                mass_mat_coeff(2, 3) =
                    porosity * (rho_mol_nonwet * (1 - Sw) +
                        Sw * rho_mol_wet * pg_int_pt / Hen_L_co2 +
                        Sw * X_L_co2_gp * d_rho_mol_wet_d_x3);
                mass_mat_coeff(2, 4) =
                    porosity * (rho_mol_nonwet * X3_int_pt * dSgdPC -
                        rho_mol_wet * X_L_co2_gp * dSgdPC +
                        Sw * d_rho_mol_wet_d_pc * X_L_co2_gp);
                // air
                mass_mat_coeff(3, 0) =
                    porosity * ((1 - Sw) * x_nonwet_air * d_rho_mol_nonwet_d_pg +
                    (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_pg +
                        Sw * rho_mol_wet * d_x_nonwet_air_d_pg * K_G_air +
                        Sw * rho_mol_wet * x_nonwet_air / Hen_L_air +
                        Sw * x_wet_air * d_rho_mol_wet_d_pg);
                mass_mat_coeff(3, 1) =
                    porosity * ((1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x1 +
                        Sw * rho_mol_wet * K_G_air * d_x_nonwet_air_d_x1 +
                        Sw * x_wet_air * d_rho_mol_wet_d_x1);
                mass_mat_coeff(3, 2) =
                    porosity * ((1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x2 +
                        Sw * rho_mol_wet * K_G_air * d_x_nonwet_air_d_x2 +
                        Sw * x_wet_air * d_rho_mol_wet_d_x2);
                mass_mat_coeff(3, 3) =
                    porosity * ((1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x3 +
                        Sw * rho_mol_wet * K_G_air * d_x_nonwet_air_d_x3 +
                        Sw * x_wet_air * d_rho_mol_wet_d_x3);
                mass_mat_coeff(3, 4) =
                    porosity * (rho_mol_nonwet * x_nonwet_air * dSgdPC -
                        rho_mol_wet * K_G_air * x_nonwet_air * dSgdPC +
                        Sw * d_rho_mol_wet_d_pc * x_wet_air +
                        Sw * K_G_air * d_x_nonwet_air_d_pc);

                // h2o
                mass_mat_coeff(4, 0) = porosity * ((1 - Sw) * d_rho_mol_nonwet_d_pg +
                    Sw * d_rho_mol_wet_d_pg);
                mass_mat_coeff(4, 1) = porosity * Sw * d_rho_mol_wet_d_x1;
                mass_mat_coeff(4, 2) = porosity * Sw * d_rho_mol_wet_d_x2;
                mass_mat_coeff(4, 3) = porosity * Sw * d_rho_mol_wet_d_x3;
                mass_mat_coeff(4, 4) =
                    porosity * (rho_mol_nonwet * dSgdPC - rho_mol_wet * dSgdPC +
                        Sw * d_rho_mol_wet_d_pc);
                //-------------debugging------------------------
                // std::cout << "mass_mat_coeff=" << std::endl;
                // std::cout << mass_mat_coeff << std::endl;
                //--------------end debugging-------------------
                for (int ii = 0; ii < NUM_NODAL_DOF; ii++)
                {
                    for (int jj = 0; jj < NUM_NODAL_DOF; jj++)
                    {
                        localMass_tmp.setZero();
                        localMass_tmp.noalias() =
                            mass_mat_coeff(ii, jj) *sm.N.transpose() * sm.N * integration_factor;
                        //_ip_data[ip].massOperator;
                        local_M.block(n_nodes * ii, n_nodes * jj, n_nodes, n_nodes)
                            .noalias() += localMass_tmp;
                    }
                }
                double const k_rel_G =
                    _process_data._material->getNonwetRelativePermeability(
                        t, pos, pg_int_pt, temperature, Sw);
                double const mu_gas =
                    _process_data._material->getGasViscosity(pg_int_pt, temperature);
                double const lambda_G = k_rel_G / mu_gas;
                // diffusion coefficient in water phase
                double const D_L =
                    _process_data._diffusion_coeff_component_b(t, pos)[0];

                // wet
                double k_rel_L =
                    _process_data._material->getWetRelativePermeability(
                        t, pos, _pressure_wetting[ip], temperature, Sw);
                if (atm_flag)
                    k_rel_L = 0;
                double const mu_liquid = _process_data._material->getLiquidViscosity(
                    _pressure_wetting[ip], temperature);
                double const lambda_L = k_rel_L / mu_liquid;
                // diffusion coefficient in gas phase
                double const D_G =
                    _process_data._diffusion_coeff_component_a(t, pos)[0];
                double const D_G_co2 =
                    _process_data._diffusion_coeff_component_c(t, pos)[0];

                K_mat_coeff(0, 0) =
                    (lambda_G * rho_mol_nonwet * X1_int_pt +
                        lambda_L * rho_mol_wet * X_L_h_gp) *
                    permeability(0, 0) +
                    (porosity * D_L * Sw * rho_mol_wet * X1_int_pt / Hen_L_h);
                K_mat_coeff(0, 1) =
                    (porosity * D_G * (1 - Sw) * rho_mol_nonwet +
                        porosity * D_L * Sw * rho_mol_wet * pg_int_pt / Hen_L_h);
                K_mat_coeff(0, 2) = 0.0;
                K_mat_coeff(0, 3) = 0.0;
                K_mat_coeff(0, 4) =
                    (-lambda_L * rho_mol_wet * X_L_h_gp) * permeability(0, 0);
                // ch4
                K_mat_coeff(1, 0) =
                    (lambda_G * rho_mol_nonwet * X2_int_pt +
                        lambda_L * rho_mol_wet * X_L_c_gp) *
                    permeability(0, 0) +
                    (porosity * D_L * Sw * rho_mol_wet * X2_int_pt / Hen_L_c);
                K_mat_coeff(1, 1) = 0.0;
                K_mat_coeff(1, 2) =
                    (porosity * D_G * (1 - Sw) * rho_mol_nonwet +
                        porosity * D_L * Sw * rho_mol_wet * pg_int_pt / Hen_L_c);
                K_mat_coeff(1, 3) = 0.0;
                K_mat_coeff(1, 4) =
                    (-lambda_L * rho_mol_wet * X_L_c_gp) * permeability(0, 0);
                // co2
                K_mat_coeff(2, 0) =
                    (lambda_G * rho_mol_nonwet * X3_int_pt +
                        lambda_L * rho_mol_wet * X_L_co2_gp) *
                    permeability(0, 0) +
                    (porosity * D_L * Sw * rho_mol_wet * X3_int_pt / Hen_L_co2);
                K_mat_coeff(2, 1) = 0.0;
                K_mat_coeff(2, 2) = 0.0;
                K_mat_coeff(2, 3) =
                    (porosity * D_G_co2 * (1 - Sw) * rho_mol_nonwet +
                        porosity * D_L * Sw * rho_mol_wet * pg_int_pt / Hen_L_co2);
                K_mat_coeff(2, 4) =
                    (-lambda_L * rho_mol_wet * X_L_co2_gp) * permeability(0, 0);
                // air
                K_mat_coeff(3, 0) =
                    (lambda_G * rho_mol_nonwet * x_nonwet_air +
                        lambda_L * rho_mol_wet * x_wet_air) *
                    permeability(0, 0) +
                    (porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_pg +
                        porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_pg * K_G_air +
                        porosity * D_L * Sw * rho_mol_wet * x_nonwet_air / Hen_L_air);
                K_mat_coeff(3, 1) =
                    (porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x1 +
                        porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_x1 * K_G_air);
                K_mat_coeff(3, 2) =
                    (porosity * D_G * (1 - Sw) * rho_mol_nonwet * d_x_nonwet_air_d_x2 +
                        porosity * D_L * Sw * rho_mol_wet * d_x_nonwet_air_d_x2 * K_G_air);
                K_mat_coeff(3, 3) =
                    (porosity * D_G * S_G_gp * rho_mol_nonwet * d_x_nonwet_air_d_x3 +
                        porosity * D_L * (1 - S_G_gp) * rho_mol_wet * d_x_nonwet_air_d_x3 *
                        K_G_air);
                K_mat_coeff(3, 4) =
                    (-lambda_L * rho_mol_wet * K_G_air * x_nonwet_air) *
                    permeability(0, 0) +
                    porosity * D_G * S_G_gp * rho_mol_nonwet * d_x_nonwet_air_d_pc +
                    porosity * D_L * (1 - S_G_gp) * rho_mol_wet * K_G_air *
                    d_x_nonwet_air_d_pc;
                // h2o
                K_mat_coeff(4, 0) =
                    (lambda_G * rho_mol_nonwet + lambda_L * rho_mol_wet) *
                    permeability(0, 0);

                K_mat_coeff(4, 1) = 0.0;
                K_mat_coeff(4, 2) = 0.0;
                K_mat_coeff(4, 3) = 0.0;
                K_mat_coeff(4, 4) = (-lambda_L * rho_mol_wet) * permeability(0, 0);

                //-------------debugging------------------------
                // std::cout << "K_mat_coeff=" << std::endl;
                // std::cout << K_mat_coeff << std::endl;
                //--------------end debugging-------------------

                for (int ii = 0; ii < NUM_NODAL_DOF; ii++)
                {
                    for (int jj = 0; jj < NUM_NODAL_DOF; jj++)
                    {
                        localDispersion_tmp.setZero();
                        localDispersion_tmp.noalias() =
                            K_mat_coeff(ii, jj) * sm.dNdx.transpose() * sm.dNdx * integration_factor;
                        //_ip_data[ip].diffusionOperator;
                        local_K.block(n_nodes * ii, n_nodes * jj, n_nodes, n_nodes)
                            .noalias() += localDispersion_tmp;
                    }
                }
                auto const K_mat_coeff_gas = permeability * (k_rel_G / mu_gas);
                auto const K_mat_coeff_liquid = permeability * (k_rel_L / mu_liquid);

                /*for (int nn = 0; nn < ShapeFunction::NPOINTS; nn++) {
                x4_nodal_values[nn] =get_x_nonwet_h2o(
                pg_int_pt, x1_nodal_values[nn], x2_nodal_values[nn], x2_nodal_values[nn], P_sat_gp, kelvin_term);

                x5_nodal_values[nn] = 1 - x4_nodal_values[nn] - x3_nodal_values[nn]
                - x2_nodal_values[nn] - x1_nodal_values[nn];
                }*/
                //calculate the velocity
                GlobalDimVectorType darcy_velocity_gas_phase =
                    -K_mat_coeff_gas * sm.dNdx * p_nodal_values;
                GlobalDimVectorType darcy_velocity_liquid_phase =
                    -K_mat_coeff_liquid * sm.dNdx * (p_nodal_values - pc_nodal_values);
                // for simplify only consider the gaseous specices diffusion velocity
                GlobalDimVectorType diffuse_velocity_h2_gas = -porosity * D_G * (1 - Sw)*sm.dNdx*x1_nodal_values;
                GlobalDimVectorType diffuse_velocity_ch4_gas = -porosity * D_G * (1 - Sw)*sm.dNdx*x2_nodal_values;
                GlobalDimVectorType diffuse_velocity_co2_gas = -porosity * D_G_co2 * (1 - Sw)*sm.dNdx*x3_nodal_values;

                GlobalDimVectorType diffuse_velocity_air_gas = -porosity * D_G * (1 - Sw)*sm.dNdx*(
                    p_nodal_values*d_x_nonwet_air_d_pg
                    + x1_nodal_values*d_x_nonwet_air_d_x1
                    + x2_nodal_values* d_x_nonwet_air_d_x2
                    + (D_G_co2 / D_G)*x3_nodal_values* d_x_nonwet_air_d_x3
                    + pc_nodal_values*d_x_nonwet_air_d_pc);

                GlobalDimVectorType diffuse_velocity_vapor_gas = -porosity * D_G * (1 - Sw)*sm.dNdx*(
                    p_nodal_values*d_x_nonwet_h2o_d_pg
                    + x1_nodal_values*d_x_nonwet_h2o_d_x1
                    + x2_nodal_values* d_x_nonwet_h2o_d_x2
                    + (D_G_co2 / D_G)*x3_nodal_values* d_x_nonwet_h2o_d_x3
                    + pc_nodal_values*d_x_nonwet_h2o_d_pc);

                double co2_degradation_rate = 0;
                if (_process_data._has_gravity)
                {
                    auto const& b = _process_data._specific_body_force;
                    NodalVectorType gravity_operator =
                        sm.dNdx.transpose() * permeability * b * integration_factor;
                    darcy_velocity_gas_phase -= K_mat_coeff_gas * rho_mass_G_gp * b;
                    darcy_velocity_liquid_phase -= K_mat_coeff_liquid*rho_mass_wet*b;
                    H_vec_coeff(0) =
                        (-lambda_G * rho_mol_nonwet * X1_int_pt * rho_mass_G_gp -
                            lambda_L * rho_mol_wet * X1_int_pt * pg_int_pt * rho_mass_wet /
                            Hen_L_h);
                    H_vec_coeff(1) =
                        (-lambda_G * rho_mol_nonwet * X2_int_pt * rho_mass_G_gp -
                            lambda_L * rho_mol_wet * X2_int_pt * pg_int_pt * rho_mass_wet /
                            Hen_L_c);
                    H_vec_coeff(2) =
                        (-lambda_G * rho_mol_nonwet * X3_int_pt * rho_mass_G_gp -
                            lambda_L * rho_mol_wet * X3_int_pt * pg_int_pt * rho_mass_wet /
                            Hen_L_co2);
                    H_vec_coeff(3) =
                        (-lambda_G * rho_mol_nonwet * x_nonwet_air * rho_mass_G_gp -
                            lambda_L * rho_mol_wet * x_nonwet_air * pg_int_pt *
                            rho_mass_wet / Hen_L_air);
                    H_vec_coeff(4) = (-lambda_G * rho_mol_nonwet * rho_mass_G_gp -
                        lambda_L * rho_mol_wet * rho_mass_wet);

                    cache_mat_gas_vel.col(ip).noalias() = darcy_velocity_gas_phase;

                    cache_mat_liquid_vel.col(ip).noalias() = darcy_velocity_liquid_phase;

                    cache_mat_gas_co2_vel.col(ip).noalias() =
                        X3_int_pt * darcy_velocity_gas_phase + diffuse_velocity_co2_gas;

                    cache_mat_gas_hydrogen_vel.col(ip).noalias() =
                        X1_int_pt * darcy_velocity_gas_phase + diffuse_velocity_h2_gas;

                    cache_mat_gas_methane_vel.col(ip).noalias() =
                        X2_int_pt * darcy_velocity_gas_phase + diffuse_velocity_ch4_gas;

                    cache_mat_gas_nitrogen_vel.col(ip).noalias() =
                        x_nonwet_air*darcy_velocity_gas_phase + diffuse_velocity_air_gas;

                    cache_mat_gas_water_vapor_vel.col(ip).noalias() =
                        x_nonwet_h2o*darcy_velocity_gas_phase + diffuse_velocity_vapor_gas;

                    for (unsigned d = 0; d < GlobalDim; ++d)
                    {
                        _total_velocities_gas[d][ip] = darcy_velocity_gas_phase[d];

                        _total_velocities_liquid[d][ip] = darcy_velocity_liquid_phase[d];
                    }

                    for (int idx = 0; idx < NUM_NODAL_DOF; idx++)
                    {
                        // since no primary vairable involved
                        // directly assemble to the Right-Hand-Side
                        // F += dNp^T * K * gz
                        localGravity_tmp.setZero();
                        localGravity_tmp.noalias() =
                            H_vec_coeff(idx) * gravity_operator;
                        local_b.block(n_nodes * idx, 0, n_nodes, 1).noalias() +=
                            localGravity_tmp;
                    }
                }  // end of hasGravityEffect
                   // load the source term
                double porosity_test = bi_interpolation(
                    94.8,
                    1459,
                    _porosity_at_supp_pnts_backfill);  // porosity update
                                                       //store the secondary variable
                if (Sw > -1e-6 && dt > 0)
                {
                    Q_steel_waste_matrix = (9.3682 / 0.13061) * 0.28 * (4.0/3.0) ; // surface area steel in waste / volume waste * reaction rate *4/3 = hydrogen production in mol/(m^3 a)
                                                                              //  Steel in waste is assumed to corrode fast, therefore we can fix rates for H2 here
                                                                              // fortunately water consumption has the same rate as H2 production
                                                                        //interpolated_Q_slow.getValue(0)*
                                                                        /*Eigen::VectorXd F_vec_coeff = Eigen::VectorXd::Zero(NUM_NODAL_DOF);
                                                                        double Q_organic_slow_co2_ini =
                                                                        interpolated_Q_slow.getValue(t);  // read from curves
                                                                        double Q_organic_fast_co2_ini =
                                                                        interpolated_Q_fast.getValue(t);  // read from curves
                                                                        */
                    if (_process_data._material->getMaterialID(
                        pos.getElementID().get()) == 1)//waste matrix
                    {

                                            // here we need for waste matrix the gas production per volume!
                    // surface area of steel [m^2]/ volume of waste compartment [m^3] * gas production [mol / (m^2 * a)]

                       // instead of reading curve, now use analytical formula
                        // quick hack to set start_conditions
//                       double& M_organic_fast_co2_ini = _ip_data[ip].amount_organic_waste_prev_cellulose; // amount from which gas is produced

                       if ((t < 0.00002) &&  (M_organic_fast_co2_ini < m0_cellulose)) M_organic_fast_co2_ini = m0_cellulose;
                       double dummy = (M_organic_fast_co2_ini*k_d_cellulose*bazant_power);  //updated amount for next time step assuming fixed degradation
                       _amount_organic_waste_cellulose[ip] = M_organic_fast_co2_ini - (dummy * dt);
                       _ip_data[ip].amount_organic_waste_cellulose=_amount_organic_waste_cellulose[ip];
                       // read from curvesinterpolated_Q_fast.getValue(0)
  //                   double& M_organic_slow_co2_ini = _ip_data[ip].amount_organic_waste_prev_polystyrene; // amount from which gas is produced
                       if ((t < 0.00002) &&  (M_organic_slow_co2_ini < m0_polystyrene)) M_organic_slow_co2_ini = m0_polystyrene;
                       dummy = M_organic_slow_co2_ini*k_d_polystyrene*bazant_power;  //updated amount for next time step assuming fixed degradation
                       _amount_organic_waste_polystyrene[ip] = M_organic_slow_co2_ini - (dummy * dt);
                       _ip_data[ip].amount_organic_waste_polystyrene=_amount_organic_waste_polystyrene[ip];

                        //calculate the fluid volume change
                        //double& fluid_volume_waste = _ip_data[ip].fluid_volume_waste;
                        //fluid_volume_waste = piecewiselinear_interpolation(
                            //_ip_data[ip].rho_mol_co2_cumul_total_prev_waste, _fluid_volume_suppt_pnt_waste);
                        //double const fluid_volume_rate_waste _ip_data[ip].=
                            //(fluid_volume_waste - _ip_data[ip].fluid_volume_prev_waste) / dt;
                        // steel corrosion rate multiply reactivity
                        Q_steel_waste_matrix *= bazant_power; // multiply with chemical reactivity
                        F_vec_coeff(0) += Q_steel_waste_matrix;  //this is for H2 source/sink
                        //store the gas h2 generation rate
                        _gas_h2_generation_rate[ip] = Q_steel_waste_matrix;
                        const double Q_organic_slow_co2 =
                            M_organic_slow_co2_ini * para_slow*bazant_power; // gas generation rate for CO2 (only) after accounting
                                                                             //for reduction in chemical reactivity due to saturation (bazant_power)
                        const double Q_organic_fast_co2 =
                            M_organic_fast_co2_ini * para_fast*bazant_power;

                        if (_ip_data[ip].rho_mol_co2_cumul_total_prev_waste >=
                            400)  // means carbonation stops, no more co2 will be consumed
                            rho_mol_total_co2_waste = 0.0;
                        //update the cumulated co2 consumption
                        rho_mol_co2_cumul_total_waste =
                            _ip_data[ip].rho_mol_co2_cumul_total_prev_waste;// +rho_mol_total_co2_waste;

                        F_vec_coeff(1) += (Q_organic_slow_co2 * 5/3);  //this is for CH4

                        F_vec_coeff(2) += Q_organic_slow_co2;      //this is for CO2
                        //store the secondary variable
                        co2_degradation_rate += Q_organic_slow_co2;   //adding for degradation (inside the inner pipe it is in current implementation zero)

                        F_vec_coeff(1) += Q_organic_fast_co2;  //this is for CH4

                        F_vec_coeff(2) += Q_organic_fast_co2;  // this is for CO2
                        co2_degradation_rate += Q_organic_fast_co2;

                        //F_vec_coeff(2) -= (rho_mol_total_co2_waste / dt);//consumption of carbonation
                        //F_vec_coeff(4) += Q_steel_waste_matrix;// switch off the water consumption
                        /*F_vec_coeff(4) += (Q_organic_slow_co2 * 8 / 3) +
                            (Q_organic_fast_co2 * 6 / 3);*/
                        F_vec_coeff(4) += (Q_organic_slow_co2 * 2 / 3) +
                            (Q_organic_fast_co2 * 5 / 3);   //KG 44: this adds mol amounts of water and gases water + CO2 + CH4 +H2  ..for H2: water and H2 cancels out!
                            // for cellulose: Q_organic_fast_co2 3/3 + Q_organic_fast_co2 3/3 - 1/3 water = 5/3
                            // for polystyrene: Q_organic_slow_co2 3/3+ Q_organic_slow_co2 5/3 -6/3 = 2/3
                        //F_vec_coeff(4) +=
                        //(fluid_volume_rate_waste) -(rho_mol_total_co2_waste / dt);
                        //update the porosity
                        //= previous porosity + porosity change
                        double const poro = _process_data._material->getPorosity(
                            material_id, t, pos, pg_int_pt, temperature, 0);//initial porosity
                        porosity3 = poro;// +piecewiselinear_interpolation(
                            //_ip_data[ip].rho_mol_co2_cumul_total_prev_waste,
                            //_porosity_change_at_supp_pnt_waste);
                        // porosity should not change in current setup for waste, therefore not calculation of source term
                        // account for change in water content due to porosity change ...we loose mass if porosity gets smaller and saturation is assumed to be constant
                        // fluid_flux [m^3/a] = (porosity_value_old[ip]-porosity_value_new[ip])*saturation_value[i] / dt
                        // positive values (porosity gets smaller): we add the lost water during next time step as source term
                        // negative values (porosity gets bigger): we remove the added water during next time step as sink term
                        // F_vec_coeff(4) += (_porosity_value[ip] -porosity2)*_saturation[ip]/dt;
                        _porosity_value[ip] = porosity3;

                        _rho_mol_co2_cumulated_prev[ip] = rho_mol_co2_cumul_total_waste;
                        //store
                        _h2o_consumed_rate[ip]
                            = -Q_steel_waste_matrix - Q_organic_slow_co2* 2 - Q_organic_fast_co2 / 3 ;
                    }
                    else if (_process_data._material->getMaterialID(
                        pos.getElementID().get()) == 0)//backfill, cement and concrete
                    {
                        double& fluid_volume_backfill = _ip_data[ip].fluid_volume_backfill;
                        fluid_volume_backfill = bi_interpolation(
                            _ip_data[ip].rho_mol_sio2_prev_backfill,
                            _ip_data[ip].rho_mol_co2_cumul_total_prev_backfill, _fluid_volume_suppt_pnt_backfill);
                        double quartz_dissolute_rate_backfill = 31557600 * bi_interpolation(
                            _ip_data[ip].rho_mol_sio2_prev_backfill,
                            _ip_data[ip].rho_mol_co2_cumul_total_prev_backfill, _quartz_rate_suppt_pnt_backfill);
                        double saturation_index_interpolated = bi_interpolation(
                            _ip_data[ip].rho_mol_sio2_prev_backfill,
                            _ip_data[ip].rho_mol_co2_cumul_total_prev_backfill, _saturation_index_suppt_pnts_backfill);
                        // if the "saturation
                        // index" for Quartz is lower than 0.95, indicating that Quartz is supposed
                        // to dissolve(i.e.ASR is active). If the saturation index is close to 1,
                        // the system is in equilibrium (also indicated by the rates which are very
                        // small), for values bigger than 1 Quartz will precipitate...->positive rates.
                        if (saturation_index_interpolated > 0.95)
                            quartz_dissolute_rate_backfill = 0;
                        if (quartz_dissolute_rate_backfill>0)
                            quartz_dissolute_rate_backfill = 0;
                        //quartz_dissolute multiply reactivity
                        quartz_dissolute_rate_backfill = bazant_power*quartz_dissolute_rate_backfill;
                        if (_ip_data[ip].rho_mol_co2_cumul_total_prev_backfill >=
                            7500)  // means carbonation stops, no more co2 will be consumed
                            rho_mol_co2_kinetic_rate_backfill = 0.0;
                        // change of volume due to look-up table Units are m3/a
                        double const fluid_volume_rate =
                            (fluid_volume_backfill - _ip_data[ip].fluid_volume_prev_backfill) / dt;

                        // update the current cumulated co2 consumption
                        rho_mol_co2_cumul_total_backfill =
                            _ip_data[ip].rho_mol_co2_cumul_total_prev_backfill +
                            (rho_mol_co2_kinetic_rate_backfill*dt*accelerate_factor);// *2 * 3.1415926*_interpolateGaussNode_coord[0];
                                                                   // co2 consumption
                        F_vec_coeff(2) -= rho_mol_co2_kinetic_rate_backfill;
                        //(rho_mol_total_co2_backfill / dt);
                        // water source/sink term
                        F_vec_coeff(4) +=
                            (fluid_volume_rate*rho_mol_wet*_element.getContent() / 4)-rho_mol_co2_kinetic_rate_backfill;//
                         //-rho_mol_co2_kinetic_rate_backfill;//switch off the water consumption
                        // update the amount of dissolved sio2
                        rho_mol_sio2_wet_backfill =
                            _ip_data[ip].rho_mol_sio2_prev_backfill -
                            (quartz_dissolute_rate_backfill * dt);// * 2 * 3.1415926*_interpolateGaussNode_coord[0];

                                                                  // porosity update
                        porosity2 = bi_interpolation(
                            _ip_data[ip].rho_mol_sio2_prev_backfill,
                            _ip_data[ip].rho_mol_co2_cumul_total_prev_backfill,
                            _porosity_at_supp_pnts_backfill);  // porosity update

                        // account for change in water content due to porosity change ...we loose mass if porosity gets smaller and saturation is assumed to be constant
                        // fluid_flux [m^3/a] = (porosity_value_old[ip]-porosity_value_new[ip])*saturation_value[i] / dt
                        // positive values (porosity gets smaller): we add the lost water during next time step as source term
                        // negative values (porosity gets bigger): we remove the added water during next time step as sink term
                        // volume change should be multiplied with molar density to get mol/m3/a
                        F_vec_coeff(4) += _element.getContent() / 4
                            *(_porosity_value[ip] -porosity2)*_saturation[ip]/dt*rho_mol_wet;
                        // do a similar correction for the gases
                        // add the whole gas phase (includes all components and humidity)
                        F_vec_coeff(4) += _element.getContent() / 4
                            *(_porosity_value[ip] -porosity2)*(1-_saturation[ip])/dt*rho_mol_nonwet;
                        // now add source/sinks for the gases: add volumes multiplied with molar density and with mol fraction
                        // first CH4
                        F_vec_coeff(1) += _element.getContent() / 4
                            *(_porosity_value[ip] -porosity2)*(1-_saturation[ip])/dt*rho_mol_nonwet* X2_int_pt;
                        // now CO2
                        F_vec_coeff(2) += _element.getContent() / 4
                            *(_porosity_value[ip] -porosity2)*(1-_saturation[ip])/dt*rho_mol_nonwet* X3_int_pt;
                        // N2
                        F_vec_coeff(3) += _element.getContent() / 4
                            *(_porosity_value[ip] -porosity2)*(1-_saturation[ip])/dt*rho_mol_nonwet*x_nonwet_air;
                        // then hydrogen
                        F_vec_coeff(0) += _element.getContent() / 4
                            *(_porosity_value[ip] -porosity2)*(1-_saturation[ip])/dt*rho_mol_nonwet* X1_int_pt;

                        //store the secondary variable porosity
                        _porosity_value[ip] = porosity2;

                        _rho_mol_co2_cumulated_prev[ip] = rho_mol_co2_cumul_total_backfill;
                        _rho_mol_sio2_cumulated_prev[ip] = rho_mol_sio2_wet_backfill;

                        //store the gas h2 generation rate
                        _gas_h2_generation_rate[ip] = F_vec_coeff(0);
                        //store the h2o consumption/release rate due to asr&carbonation
                        _h2o_consumed_rate[ip] = fluid_volume_rate*rho_mol_wet*_element.getContent() / 4
                            + _element.getContent() / 4
                            * (_porosity_value[ip] - porosity2)*_saturation[ip] / dt * rho_mol_wet;
                    }
                    //store the source term for each component
                    // thanks to the facts that the source terms are all for gas phase
                    //gas_h2_generation_rate = F_vec_coeff(0);
                    _gas_ch4_generation_rate[ip] = F_vec_coeff(1);
                    _gas_co2_generation_rate[ip] = F_vec_coeff(2);
                    _gas_co2_degradation_rate[ip] = co2_degradation_rate;
                    for (int idx = 0; idx < NUM_NODAL_DOF; idx++)
                    {
                        // since no primary vairable involved
                        // directly assemble to the Right-Hand-Side
                        // F += dNp^T * K * gz
                        localSource_tmp.setZero();
                        localSource_tmp.noalias() =
                            sm.N.transpose() * F_vec_coeff(idx) * integration_factor;
                        local_b.block(n_nodes * idx, 0, n_nodes, 1).noalias() +=
                            localSource_tmp;
                    }
                } // end of loading the source and sink term
                  //_gas_h2_generation_rate[ip] = gas_h2_generation_rate;
            }  // end of GP asm
               //for the outer drum
               //apply the neumann boundary condition on the line element
               //first search the nodes
               //axissymmetric
            double ele_bazant_power = 0;
            for (auto const& ip : _reactivity_bazant_power)
            {
                ele_bazant_power += ip;
            }
            ele_bazant_power /= static_cast<double>(n_integration_points);
            if (ele_bazant_power > 1)
                ele_bazant_power = 1;
            else if (ele_bazant_power < 0)
                ele_bazant_power = 0;

            if (abs( rx0-0.303) < eps  && abs(rx1 - 0.303) < eps)
            {
                //indicates edge 0-1 located on the boundary
                length = std::sqrt(std::pow(rx0 - rx1, 2) + std::pow(ry0 - ry1, 2));
                neumann_vec[2] = 0.0;
                neumann_vec[3] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx0;
                neumn_h2 = 0.003733333*ele_bazant_power;
                if (accelerate_flag) {
                    if (gp_carb_neutral_count <= 2) {
                        if (ry0 > ry1)
                        {
                            neumann_vec[0] = 100 * neumn_h2;
                            neumann_vec[1] = neumn_h2;
                        }
                        else
                        {
                            neumann_vec[0] = neumn_h2;
                            neumann_vec[1] = 100 * neumn_h2;
                        }
                    }
                    else {
                        neumann_vec[0] = 100 * neumn_h2;
                        neumann_vec[1] = 100 * neumn_h2;
                    }
                }
                else {
                    neumann_vec[0] = neumn_h2;
                    neumann_vec[1] = neumn_h2;
                }
                localNeumann_tmp = neumann_vec* radial_sym_fac* length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2 / node_volume_radial;
            }
            else if (abs(rx1- 0.303) <eps  && abs(rx2 - 0.303) <eps)
            {
                length = std::sqrt(std::pow(rx1 - rx2, 2) + std::pow(ry1 - ry2, 2));
                neumann_vec[0] = 0.0;
                neumann_vec[3] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx1;
                neumn_h2 = 0.003733333*ele_bazant_power;
                if (accelerate_flag) {
                    if (gp_carb_neutral_count <= 2) {
                        if (ry1 > ry2)
                        {
                            neumann_vec[1] = 100 * neumn_h2;
                            neumann_vec[2] = neumn_h2;
                        }
                        else
                        {
                            neumann_vec[1] = neumn_h2;
                            neumann_vec[2] = 100 * neumn_h2;
                        }
                    }
                    else {
                        neumann_vec[1] = 100 * neumn_h2;
                        neumann_vec[2] = 100 * neumn_h2;
                    }
                }
                else {
                    neumann_vec[1] = neumn_h2;
                    neumann_vec[2] = neumn_h2;
                }
                localNeumann_tmp = neumann_vec * radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (abs(rx2 - 0.303) < eps && abs(rx0 - 0.303) < eps)
            {
                length = std::sqrt(std::pow(rx0 - rx2, 2) + std::pow(ry0 - ry2, 2));
                neumann_vec[1] = 0.0;
                neumann_vec[3] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx2;
                neumn_h2 = 0.003733333*ele_bazant_power;
                if (accelerate_flag) {
                    if (gp_carb_neutral_count <= 2) {
                        if (ry0 > ry2)
                        {
                            neumann_vec[0] = 100 * neumn_h2;
                            neumann_vec[2] = neumn_h2;
                        }
                        else
                        {
                            neumann_vec[0] = neumn_h2;
                            neumann_vec[2] = 100 * neumn_h2;
                        }
                    }
                    else {
                        neumann_vec[0] = 100 * neumn_h2;
                        neumann_vec[2] = 100 * neumn_h2;
                    }
                }
                else {
                    neumann_vec[0] = neumn_h2;
                    neumann_vec[2] = neumn_h2;
                }
                localNeumann_tmp = neumann_vec* radial_sym_fac*length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (abs(rx3 - 0.303) < eps && abs(rx0 - 0.303) < eps)
            {
                length = std::sqrt(std::pow(rx0 - rx3, 2) + std::pow(ry0 - ry3, 2));
                neumann_vec[1] = 0.0;
                neumann_vec[2] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx2;
                neumn_h2 = 0.003733333*ele_bazant_power;
                if (accelerate_flag) {
                    if (gp_carb_neutral_count <= 2) {
                        if (ry0 > ry3)
                        {
                            neumann_vec[0] = 100 * neumn_h2;
                            neumann_vec[3] = neumn_h2;
                        }
                        else
                        {
                            neumann_vec[0] = neumn_h2;
                            neumann_vec[3] = 100 * neumn_h2;
                        }
                    }
                    else {
                        neumann_vec[0] = 100 * neumn_h2;
                        neumann_vec[3] = 100 * neumn_h2;
                    }
                }
                else {
                    neumann_vec[0] = neumn_h2;
                    neumann_vec[3] = neumn_h2;
                }
                localNeumann_tmp = neumann_vec * radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (abs(rx3 - 0.303) < eps && abs(rx1 - 0.303) < eps)
            {
                length = std::sqrt(std::pow(rx1 - rx3, 2) + std::pow(ry1 - ry3, 2));
                neumann_vec[0] = 0.0;
                neumann_vec[2] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx1;
                neumn_h2 = 0.003733333*ele_bazant_power;
                if (accelerate_flag) {
                    if (gp_carb_neutral_count <= 2) {
                        if (ry1 > ry3)
                        {
                            neumann_vec[1] = 100 * neumn_h2;
                            neumann_vec[3] = neumn_h2;
                        }
                        else
                        {
                            neumann_vec[1] = neumn_h2;
                            neumann_vec[3] = 100 * neumn_h2;
                        }
                    }
                    else {
                        neumann_vec[1] = 100 * neumn_h2;
                        neumann_vec[3] = 100 * neumn_h2;
                    }
                }
                else {
                    neumann_vec[1] = neumn_h2;
                    neumann_vec[3] = neumn_h2;
                }
                localNeumann_tmp = neumann_vec * radial_sym_fac* length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (abs(rx3 - 0.303) < eps && abs(rx2 - 0.303) < eps)
            {
                length = std::sqrt(std::pow(rx2 - rx3, 2) + std::pow(ry2 - ry3, 2));
                neumann_vec[0] = 0.0;
                neumann_vec[1] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx2;
                neumn_h2 = 0.003733333*ele_bazant_power;
                if (accelerate_flag) {
                    if (gp_carb_neutral_count <= 2) {
                        if (ry2 > ry3)
                        {
                            neumann_vec[2] = 100 * neumn_h2;
                            neumann_vec[3] = neumn_h2;
                        }
                        else
                        {
                            neumann_vec[2] = neumn_h2;
                            neumann_vec[3] = 100 * neumn_h2;
                        }
                    }
                    else {
                        neumann_vec[2] = 100 * neumn_h2;
                        neumann_vec[3] = 100 * neumn_h2;
                    }
                }
                else {
                    neumann_vec[2] = neumn_h2;
                    neumann_vec[3] = neumn_h2;
                }
                localNeumann_tmp = neumann_vec * radial_sym_fac* length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            // for the second Neumann boundary condition
            if (std::abs(rx0 - 0.245)<0.0025 + eps && std::abs(rx1 - 0.245)<0.0025 + eps &&
                ry1<0.795 + eps && ry1>0.088 - eps && ry0<0.795 + eps && ry0>0.088 - eps)
            {
                //indicates edge 0-1 located on the boundary
                length = std::sqrt(std::pow(rx0 - rx1, 2) + std::pow(ry0 - ry1, 2));
                neumann_vec[2] = 0;
                radial_sym_fac = 2 * 3.1415926*rx0;
                neumn_h2 = 0.003733333*ele_bazant_power;
                if (accelerate_flag) {
                    if (gp_carb_neutral_count <= 2) {
                        if (ry0 > ry1)
                        {
                            neumann_vec[0] = 100 * neumn_h2;
                            neumann_vec[1] = neumn_h2;
                        }
                        else
                        {
                            neumann_vec[0] = neumn_h2;
                            neumann_vec[1] = 100 * neumn_h2;
                        }
                    }
                    else {
                        neumann_vec[0] = 100 * neumn_h2;
                        neumann_vec[1] = 100 * neumn_h2;
                    }
                }
                else {
                    neumann_vec[0] = neumn_h2;
                    neumann_vec[1] = neumn_h2;
                }
                localNeumann_tmp = neumann_vec* radial_sym_fac*length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(rx1 - 0.245)<0.0025 + eps && std::abs(rx2 - 0.245)<0.0025 + eps &&
                ry1<0.795 + eps && ry1>0.088 - eps && ry2<0.795 + eps && ry2>0.088 - eps)
            {
                length = std::sqrt(std::pow(rx1 - rx2, 2) + std::pow(ry1 - ry2, 2));
                neumann_vec[0] = 0.0;
                radial_sym_fac = 2 * 3.14159*rx1;
                neumn_h2 = 0.003733333*ele_bazant_power;
                if (accelerate_flag) {
                    if (gp_carb_neutral_count <= 2) {
                        if (ry1 > ry2)
                        {
                            neumann_vec[1] = 100 * neumn_h2;
                            neumann_vec[2] = neumn_h2;
                        }
                        else
                        {
                            neumann_vec[1] = neumn_h2;
                            neumann_vec[2] = 100 * neumn_h2;
                        }
                    }
                    else {
                        neumann_vec[1] = 100 * neumn_h2;
                        neumann_vec[2] = 100 * neumn_h2;
                    }
                }
                else {
                    neumann_vec[1] = neumn_h2;
                    neumann_vec[2] = neumn_h2;
                }
                localNeumann_tmp = neumann_vec* radial_sym_fac*length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(rx2 - 0.245)<0.0025 + eps && std::abs(rx0 - 0.245)<0.0025 + eps &&
                ry2<0.795 + eps && ry2>0.088 - eps && ry0<0.795 + eps && ry0>0.088 - eps)
            {
                length = std::sqrt(std::pow(rx0 - rx2, 2) + std::pow(ry0 - ry2, 2));
                neumann_vec[1] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx2;
                neumn_h2 = 0.003733333*ele_bazant_power;
                if (accelerate_flag) {
                    if (gp_carb_neutral_count <= 2) {
                        if (ry0 > ry2)
                        {
                            neumann_vec[0] = 100 * neumn_h2;
                            neumann_vec[2] = neumn_h2;
                        }
                        else
                        {
                            neumann_vec[0] = neumn_h2;
                            neumann_vec[2] = 100 * neumn_h2;
                        }
                    }
                    else {
                        neumann_vec[0] = 100 * neumn_h2;
                        neumann_vec[2] = 100 * neumn_h2;
                    }
                }
                else {
                    neumann_vec[0] = neumn_h2;
                    neumann_vec[2] = neumn_h2;
                }
                localNeumann_tmp = neumann_vec* radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(rx0 - 0.245)<0.0025 + eps && std::abs(rx3 - 0.245)<0.0025 + eps &&
                ry0<0.795 + eps && ry3>0.088 - eps && ry0<0.795 + eps && ry3>0.088 - eps)
            {
                length = std::sqrt(std::pow(rx0 - rx3, 2) + std::pow(ry0 - ry3, 2));
                neumann_vec[1] = 0.0;
                neumann_vec[2] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx3;
                neumn_h2 = 0.003733333*ele_bazant_power;
                if (accelerate_flag) {
                    if (gp_carb_neutral_count <= 2) {
                        if (ry0 > ry3)
                        {
                            neumann_vec[0] = 100 * neumn_h2;
                            neumann_vec[3] = neumn_h2;
                        }
                        else
                        {
                            neumann_vec[0] = neumn_h2;
                            neumann_vec[3] = 100 * neumn_h2;
                        }
                    }
                    else {
                        neumann_vec[0] = 100 * neumn_h2;
                        neumann_vec[3] = 100 * neumn_h2;
                    }
                }
                else {
                    neumann_vec[0] = neumn_h2;
                    neumann_vec[3] = neumn_h2;
                }
                localNeumann_tmp = neumann_vec* radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(rx2 - 0.245)<0.0025 + eps && std::abs(rx3 - 0.245)<0.0025 + eps &&
                ry2<0.795 + eps && ry3>0.088 - eps && ry2<0.795 + eps && ry3>0.088 - eps)
            {
                length = std::sqrt(std::pow(rx2 - rx3, 2) + std::pow(ry2 - ry3, 2));
                neumann_vec[0] = 0.0;
                neumann_vec[1] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx2;
                neumn_h2 = 0.003733333*ele_bazant_power;
                if (accelerate_flag) {
                    if (gp_carb_neutral_count <= 2) {
                        if (ry2 > ry3)
                        {
                            neumann_vec[2] = 100 * neumn_h2;
                            neumann_vec[3] = neumn_h2;
                        }
                        else
                        {
                            neumann_vec[2] = neumn_h2;
                            neumann_vec[3] = 100 * neumn_h2;
                        }
                    }
                    else {
                        neumann_vec[2] = 100 * neumn_h2;
                        neumann_vec[3] = 100 * neumn_h2;
                    }
                }
                else {
                    neumann_vec[2] = neumn_h2;
                    neumann_vec[3] = neumn_h2;
                }
                localNeumann_tmp = neumann_vec* radial_sym_fac *length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(rx1 - 0.245)<0.0025 + eps && std::abs(rx3 - 0.245)<0.0025 + eps &&
                ry1<0.795 + eps && ry3>0.088 - eps && ry1<0.795 + eps && ry3>0.088 - eps)
            {
                length = std::sqrt(std::pow(rx1 - rx3, 2) + std::pow(ry1 - ry3, 2));
                neumann_vec[0] = 0.0;
                neumann_vec[2] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx1;
                neumn_h2 = 0.003733333*ele_bazant_power;
                if (accelerate_flag) {
                    if (gp_carb_neutral_count <= 2) {
                        if (ry1 > ry3)
                        {
                            neumann_vec[1] = 100 * neumn_h2;
                            neumann_vec[3] = neumn_h2;
                        }
                        else
                        {
                            neumann_vec[1] = neumn_h2;
                            neumann_vec[3] = 100 * neumn_h2;
                        }
                    }
                    else {
                        neumann_vec[1] = 100 * neumn_h2;
                        neumann_vec[3] = 100 * neumn_h2;
                    }
                }
                else {
                    neumann_vec[1] = neumn_h2;
                    neumann_vec[3] = neumn_h2;
                }
                localNeumann_tmp = neumann_vec * radial_sym_fac* length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }

            // for the third boundary condition
            if (std::abs(rx0 - 0.24)<0.0025 + eps && std::abs(rx1 - 0.24)<0.0025 + eps &&
                ry1<0.795 + eps && ry1>0.088 - eps && ry0<0.795 + eps && ry0>0.088 - eps)
            {
                //indicates edge 0-1 located on the boundary
                length = std::sqrt(std::pow(rx0 - rx1, 2) + std::pow(ry0 - ry1, 2));
                neumann_vec[2] = 0;
                neumann_vec[3] = 0;
                radial_sym_fac = 2 * 3.1415926*rx0;
                neumn_h2 = 0.3733333*ele_bazant_power;
                neumann_vec[0] = neumn_h2;
                neumann_vec[1] = neumn_h2;
                localNeumann_tmp = neumann_vec* radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(rx1 - 0.24)<0.0025 + eps && std::abs(rx2 - 0.24)<0.0025 + eps &&
                ry1<0.795 + eps && ry1>0.088 - eps && ry2<0.795 + eps && ry2>0.088 - eps)
            {
                length = std::sqrt(std::pow(rx1 - rx2, 2) + std::pow(ry1 - ry2, 2));
                neumann_vec[0] = 0.0;
                neumann_vec[3] = 0;
                radial_sym_fac = 2 * 3.14159*rx1;
                neumn_h2 = 0.3733333* ele_bazant_power;
                neumann_vec[1] = neumn_h2;
                neumann_vec[2] = neumn_h2;
                localNeumann_tmp = neumann_vec* radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec * radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(rx2 - 0.24)<0.0025 + eps && std::abs(rx0 - 0.24)<0.0025 + eps &&
                ry2<0.795 + eps && ry2>0.088 - eps && ry0<0.795 + eps && ry0>0.088 - eps)
            {
                length = std::sqrt(std::pow(rx0 - rx2, 2) + std::pow(ry0 - ry2, 2));
                neumann_vec[1] = 0.0;
                neumann_vec[3] = 0;
                radial_sym_fac = 2 * 3.1415926*rx2;
                neumn_h2 = 0.3733333* ele_bazant_power;
                neumann_vec[0] = neumn_h2;
                neumann_vec[2] = neumn_h2;
                localNeumann_tmp = neumann_vec* radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(rx0 - 0.24)<0.0025 + eps && std::abs(rx3 - 0.24)<0.0025 + eps &&
                ry0<0.795 + eps && ry3>0.088 - eps && ry0<0.795 + eps && ry3>0.088 - eps)
            {
                length = std::sqrt(std::pow(rx0 - rx3, 2) + std::pow(ry0 - ry3, 2));
                neumann_vec[1] = 0.0;
                neumann_vec[2] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx0;
                neumn_h2 = 0.3733333* ele_bazant_power;
                neumann_vec[0] = neumn_h2;
                neumann_vec[3] = neumn_h2;
                localNeumann_tmp = neumann_vec* radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(rx2 - 0.24)<0.0025 + eps && std::abs(rx3 - 0.24)<0.0025 + eps &&
                ry2<0.795 + eps && ry3>0.088 - eps && ry2<0.795 + eps && ry3>0.088 - eps)
            {
                length = std::sqrt(std::pow(rx2 - rx3, 2) + std::pow(ry2 - ry3, 2));
                neumann_vec[0] = 0.0;
                neumann_vec[1] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx2;
                neumn_h2 = 0.3733333* ele_bazant_power;
                neumann_vec[2] = neumn_h2;
                neumann_vec[3] = neumn_h2;
                localNeumann_tmp = neumann_vec* radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(rx1 - 0.24)<0.0025 + eps && std::abs(rx3 - 0.24)<0.0025 + eps &&
                ry1<0.795 + eps && ry3>0.088 - eps && ry1<0.795 + eps && ry3>0.088 - eps)
            {
                length = std::sqrt(std::pow(rx1 - rx3, 2) + std::pow(ry1 - ry3, 2));
                neumann_vec[0] = 0.0;
                neumann_vec[2] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx1;
                neumn_h2 = 0.3733333* ele_bazant_power;

                neumann_vec[1] = neumn_h2;
                neumann_vec[3] = neumn_h2;

                localNeumann_tmp = neumann_vec* radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }

            //for the bottom boundary condition
            if (std::abs(ry0)<0.0 + eps && std::abs(ry1)<0.00 + eps )
            {
                //indicates edge 0-1 located on the boundary
                length = std::sqrt(std::pow(rx0 - rx1, 2) + std::pow(ry0 - ry1, 2));
                neumann_vec[2] = 0;
                neumann_vec[3] = 0;
                radial_sym_fac = 2 * 3.1415926*rx0;
                neumn_h2 = 0.003733333*ele_bazant_power;
                neumann_vec[0] = neumn_h2;
                neumann_vec[1] = neumn_h2;
                localNeumann_tmp = neumann_vec* radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec*  radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(ry1)<0.00 + eps && std::abs(ry2)<0.00 + eps)
            {
                //indicates edge 1-2 located on the boundary
                length = std::sqrt(std::pow(rx1 - rx2, 2) + std::pow(ry1 - ry2, 2));
                neumann_vec[0] = 0.0;
                neumann_vec[3] = 0;
                radial_sym_fac = 2 * 3.14159*rx1;
                neumn_h2 = 0.003733333* ele_bazant_power;
                neumann_vec[1] = neumn_h2;
                neumann_vec[2] = neumn_h2;
                localNeumann_tmp = neumann_vec *radial_sym_fac *length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(ry2)<0.00 + eps && std::abs(ry0)<0.00 + eps)
            {
                //indicates edge 2-0 located on the boundary
                length = std::sqrt(std::pow(rx0 - rx2, 2) + std::pow(ry0 - ry2, 2));
                neumann_vec[1] = 0.0;
                neumann_vec[3] = 0;
                radial_sym_fac = 2 * 3.1415926*rx2;
                neumn_h2 = 0.003733333* ele_bazant_power;
                neumann_vec[0] = neumn_h2;
                neumann_vec[2] = neumn_h2;
                localNeumann_tmp = neumann_vec* radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(ry0 )<0.00 + eps && std::abs(ry3)<0.00 + eps)
            {
                //indicates edge 0-3 located on the boundary
                length = std::sqrt(std::pow(rx0 - rx3, 2) + std::pow(ry0 - ry3, 2));
                neumann_vec[1] = 0.0;
                neumann_vec[2] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx0;
                neumn_h2 = 0.003733333* ele_bazant_power;
                neumann_vec[0] = neumn_h2;
                neumann_vec[3] = neumn_h2;
                localNeumann_tmp = neumann_vec * radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec*radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(ry2)<0.00 + eps && std::abs(ry3)<0.00 + eps)
            {
                //indicates edge 2-3 located on the boundary
                length = std::sqrt(std::pow(rx2 - rx3, 2) + std::pow(ry2 - ry3, 2));
                neumann_vec[0] = 0.0;
                neumann_vec[1] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx2;
                neumn_h2 = 0.003733333* ele_bazant_power;
                neumann_vec[2] = neumn_h2;
                neumann_vec[3] = neumn_h2;
                localNeumann_tmp = neumann_vec * radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec* radial_sym_fac * length / 2/ node_volume_radial;
            }
            else if (std::abs(ry1)<0.00 + eps && std::abs(ry3)<0.00 + eps)
            {
                //indicates edge 1-3 located on the boundary
                length = std::sqrt(std::pow(rx1 - rx3, 2) + std::pow(ry1 - ry3, 2));
                neumann_vec[0] = 0.0;
                neumann_vec[2] = 0.0;
                radial_sym_fac = 2 * 3.1415926*rx1;
                neumn_h2 = 0.003733333* ele_bazant_power;

                neumann_vec[1] = neumn_h2;
                neumann_vec[3] = neumn_h2;

                localNeumann_tmp = neumann_vec * radial_sym_fac * length / 2;
                _neumann_vec_output = neumann_vec * radial_sym_fac * length / 2/ node_volume_radial;
            }
            local_b.block(n_nodes * 0, 0, n_nodes, 1).noalias() += localNeumann_tmp; // This is for hydrogen-> which is created
            local_b.block(n_nodes * 4, 0, n_nodes, 1).noalias() -= localNeumann_tmp; // This is for water-> which is consumed

            //output secondary variable
            for (unsigned ip = 0; ip < n_integration_points; ip++)
            {
                // interpolate the node value on the element
                auto const& sm = _shape_matrices[ip];

                double h2_flux = 0.0;
                double h2_flux2 = 0.0;
                double h2_flux3 = 0.0;
                double h2_flux4 = 0.0;
                double h2_flux5 = 0.0;

                NumLib::shapeFunctionInterpolate(_neumann_vector_output, sm.N, h2_flux, h2_flux2,
                    h2_flux3, h2_flux4, h2_flux5);
                _gas_h2_boundary_generation_rate[ip] = h2_flux; //KG: This is the flux across a surface (boundary)/ divided by element volume
                _gas_h2_overall_generation_rate[ip] =
                    _gas_h2_boundary_generation_rate[ip] + _gas_h2_generation_rate[ip];
                _h2o_consumed_rate[ip] -= _gas_h2_boundary_generation_rate[ip]; // this is only for output, or...where in the model the water sink term for boundary fluxes is set?
            }

            int n = n_integration_points;
            //Eigen::VectorXd ele_saturation = Eigen::VectorXd::Zero(n);
            double ele_saturation = 0;
            for (auto const& ip : _saturation)
            {
                ele_saturation += ip;
            }
            ele_saturation /= static_cast<double>(n_integration_points);
            auto const element_id = _element.getID();
            (*_process_data.mesh_prop_saturation)[element_id] = ele_saturation;

            double ele_porosity = 0;
            for (auto const& ip : _porosity_value)
            {
                ele_porosity += ip;
            }
            ele_porosity /=static_cast<double>(n_integration_points);
            (*_process_data.mesh_prop_porosity)[element_id] = ele_porosity;

            double ele_pH = 0;
            for (auto const& ip : _pH_value)
            {
                ele_pH += ip;
            }
            ele_pH /= static_cast<double>(n_integration_points);
            (*_process_data.mesh_prop_pHvalue)[element_id] = ele_pH;

            double ele_bazant_power_2 = 0;
            for (auto const& ip : _reactivity_bazant_power)
            {
                ele_bazant_power_2 += ip;
            }
            ele_bazant_power_2 /= static_cast<double>(n);
            (*_process_data.mesh_prop_bazant_power)[element_id] = ele_bazant_power_2;

            double ele_mol_density_gas = 0;
            for (auto const& ip : _rho_mol_gas_phase)
            {
                ele_mol_density_gas += ip;
            }
            ele_mol_density_gas /= static_cast<double>(n);
            (*_process_data.mesh_prop_mol_density_gas)[element_id] = ele_mol_density_gas;

            double ele_mol_density_liquid = 0;
            for (auto const& ip : _rho_mol_liquid_phase)
            {
                ele_mol_density_liquid += ip;
            }
            ele_mol_density_liquid /= static_cast<double>(n);
            (*_process_data.mesh_prop_mol_density_liquid)[element_id]
                = ele_mol_density_liquid;

            double ele_co2_cumulate_consume = 0;
            for (auto const& ip : _rho_mol_co2_cumulated_prev)
            {
                ele_co2_cumulate_consume += ip;
            }
            ele_co2_cumulate_consume/= static_cast<double>(n);
            (*_process_data.mesh_prop_co2_cumulate_consume)[element_id]
                = ele_co2_cumulate_consume;

            double ele_sio2_cumulate_consume = 0;
            for (auto const& ip : _rho_mol_sio2_cumulated_prev)
            {
                ele_sio2_cumulate_consume += ip;
            }
            ele_sio2_cumulate_consume/= static_cast<double>(n);
            (*_process_data.mesh_prop_sio2_cumulate_consume)[element_id]
                = ele_sio2_cumulate_consume;

            double ele_mol_frac_h2o_vapor = 0;
            for (auto const& ip : _mol_fraction_nonwet_vapor)
            {
                ele_mol_frac_h2o_vapor += ip;
            }
            ele_mol_frac_h2o_vapor /= static_cast<double>(n);
            (*_process_data.mesh_prop_mol_frac_h2o_vapor)[element_id]
                = ele_mol_frac_h2o_vapor;

            double ele_mol_frac_n2_air = 0.0;
            for (auto const& ip : _mol_fraction_nonwet_air)
            {
                ele_mol_frac_n2_air += ip;
            }
            ele_mol_frac_n2_air /= static_cast<double>(n);
            (*_process_data.mesh_prop_mol_frac_n2)[element_id]
                = ele_mol_frac_n2_air;
            // store the velocity value on each cell
            //Eigen::Vector3d ele_liquid_velocity = Eigen::Vector3d::Zero();
            Eigen::VectorXd ele_liquid_velocity = Eigen::VectorXd::Zero(GlobalDim);
            Eigen::VectorXd ele_gas_velocity = Eigen::VectorXd::Zero(GlobalDim);
            Eigen::VectorXd ele_co2_gas_velocity = Eigen::VectorXd::Zero(GlobalDim);
            Eigen::VectorXd ele_h2_gas_velocity = Eigen::VectorXd::Zero(GlobalDim);
            Eigen::VectorXd ele_CH4_gas_velocity = Eigen::VectorXd::Zero(GlobalDim);
            Eigen::VectorXd ele_h2o_vapor_gas_velocity = Eigen::VectorXd::Zero(GlobalDim);
            for (unsigned ip = 0; ip < n_integration_points; ip++)
            {
                ele_liquid_velocity += cache_mat_liquid_vel.col(ip);
                ele_gas_velocity += cache_mat_gas_vel.col(ip);
                ele_co2_gas_velocity += cache_mat_gas_co2_vel.col(ip);
                ele_h2_gas_velocity += cache_mat_gas_hydrogen_vel.col(ip);
                ele_CH4_gas_velocity += cache_mat_gas_methane_vel.col(ip);
                ele_h2o_vapor_gas_velocity += cache_mat_gas_water_vapor_vel.col(ip);
            }
            for (unsigned i = 0; i < GlobalDim; i++) {
                (*_process_data.mesh_prop_overall_liquid_vel)[element_id * 3 + i] =
                    ele_liquid_velocity[i];
                (*_process_data.mesh_prop_overall_gas_vel)[element_id * 3 + i] =
                    ele_gas_velocity[i];
                (*_process_data.mesh_prop_gas_co2_vel)[element_id * 3 + i] =
                    ele_co2_gas_velocity[i];
                (*_process_data.mesh_prop_gas_hydrogen_vel)[element_id * 3 + i] =
                    ele_h2_gas_velocity[i];
                (*_process_data.mesh_prop_gas_methane_vel)[element_id * 3 + i] =
                    ele_CH4_gas_velocity[i];
                (*_process_data.mesh_prop_gas_water_vapor_vel)[element_id * 3 + i] =
                    ele_h2o_vapor_gas_velocity[i];
            }
            (*_process_data.mesh_prop_overall_liquid_vel)[element_id * 3 + 2] =
                0.0;
            (*_process_data.mesh_prop_overall_gas_vel)[element_id * 3 + 2] =
                0.0;
            (*_process_data.mesh_prop_gas_co2_vel)[element_id * 3 + 2] =
                0.0;
            (*_process_data.mesh_prop_gas_hydrogen_vel)[element_id * 3 + 2] =
                0.0;
            (*_process_data.mesh_prop_gas_methane_vel)[element_id * 3 + 2] =
                0.0;
            (*_process_data.mesh_prop_gas_water_vapor_vel)[element_id * 3 + 2] =
                0.0;
            if (_process_data._has_mass_lumping)
            {
                auto Mhpg =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);

                auto Mhmolh2 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        nonwet_pressure_matrix_index,
                        nonwet_pressure_size * mol_fraction_h_coeff_index);

                auto Mhmolch4 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        nonwet_pressure_matrix_index,
                        nonwet_pressure_size * mol_fraction_ch4_coeff_index);
                auto Mhmolco2 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        nonwet_pressure_matrix_index,
                        nonwet_pressure_size * mol_fraction_co2_coeff_index);

                auto Mhpc =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        nonwet_pressure_matrix_index,
                        nonwet_pressure_size * cap_pressure_coeff_index);

                auto Mch4pg =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        nonwet_pressure_size, nonwet_pressure_matrix_index);

                auto Mch4molh2 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        nonwet_pressure_size,
                        nonwet_pressure_size * mol_fraction_h_coeff_index);

                auto Mch4molch4 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        nonwet_pressure_size,
                        nonwet_pressure_size * mol_fraction_ch4_coeff_index);
                auto Mch4molco2 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        nonwet_pressure_size,
                        nonwet_pressure_size * mol_fraction_co2_coeff_index);

                auto Mch4pc =
                    local_M.template block<nonwet_pressure_size, cap_pressure_size>(
                        nonwet_pressure_size,
                        nonwet_pressure_size * cap_pressure_coeff_index);

                auto Mco2pg =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        2 * nonwet_pressure_size, nonwet_pressure_matrix_index);

                auto Mco2molh2 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        2 * nonwet_pressure_size,
                        nonwet_pressure_size * mol_fraction_h_coeff_index);

                auto Mco2molch4 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        2 * nonwet_pressure_size,
                        nonwet_pressure_size * mol_fraction_ch4_coeff_index);
                auto Mco2molco2 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        2 * nonwet_pressure_size,
                        nonwet_pressure_size * mol_fraction_co2_coeff_index);

                auto Mco2pc =
                    local_M.template block<nonwet_pressure_size, cap_pressure_size>(
                        2 * nonwet_pressure_size,
                        nonwet_pressure_size * cap_pressure_coeff_index);

                auto Mairpg =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        3 * nonwet_pressure_size, nonwet_pressure_matrix_index);

                auto Mairmolh2 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        3 * nonwet_pressure_size,
                        nonwet_pressure_size * mol_fraction_h_coeff_index);

                auto Mairmolch4 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        3 * nonwet_pressure_size,
                        nonwet_pressure_size * mol_fraction_ch4_coeff_index);
                auto Mairmolco2 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        3 * nonwet_pressure_size,
                        nonwet_pressure_size * mol_fraction_co2_coeff_index);

                auto Mairpc =
                    local_M.template block<nonwet_pressure_size, cap_pressure_size>(
                        3 * nonwet_pressure_size,
                        nonwet_pressure_size * cap_pressure_coeff_index);

                auto Mh2opg =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        4 * nonwet_pressure_size, nonwet_pressure_matrix_index);

                auto Mh2omolh2 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        4 * nonwet_pressure_size,
                        nonwet_pressure_size * mol_fraction_h_coeff_index);

                auto Mh2omolch4 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        4 * nonwet_pressure_size,
                        nonwet_pressure_size * mol_fraction_ch4_coeff_index);
                auto Mh2omolco2 =
                    local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
                        4 * nonwet_pressure_size,
                        nonwet_pressure_size * mol_fraction_co2_coeff_index);

                auto Mh2opc =
                    local_M.template block<nonwet_pressure_size, cap_pressure_size>(
                        4 * nonwet_pressure_size,
                        nonwet_pressure_size * cap_pressure_coeff_index);
                for (unsigned row = 0; row < Mhpg.cols(); row++)
                {
                    for (unsigned column = 0; column < Mhpg.cols(); column++)
                    {
                        if (row != column)
                        {
                            Mhpg(row, row) += Mhpg(row, column);
                            Mhpg(row, column) = 0.0;
                            Mhmolh2(row, row) += Mhmolh2(row, column);
                            Mhmolh2(row, column) = 0.0;
                            Mhmolch4(row, row) += Mhmolch4(row, column);
                            Mhmolch4(row, column) = 0.0;
                            Mhmolco2(row, row) += Mhmolco2(row, column);
                            Mhmolco2(row, column) = 0.0;
                            Mhpc(row, row) += Mhpc(row, column);
                            Mhpc(row, column) = 0.0;

                            Mch4pg(row, row) += Mch4pg(row, column);
                            Mch4pg(row, column) = 0.0;
                            Mch4molh2(row, row) += Mch4molh2(row, column);
                            Mch4molh2(row, column) = 0.0;
                            Mch4molch4(row, row) += Mch4molch4(row, column);
                            Mch4molch4(row, column) = 0.0;
                            Mch4molco2(row, row) += Mch4molco2(row, column);
                            Mch4molco2(row, column) = 0.0;
                            Mch4pc(row, row) += Mch4pc(row, column);
                            Mch4pc(row, column) = 0.0;

                            Mco2pg(row, row) += Mco2pg(row, column);
                            Mco2pg(row, column) = 0.0;
                            Mco2molh2(row, row) += Mco2molh2(row, column);
                            Mco2molh2(row, column) = 0.0;
                            Mco2molch4(row, row) += Mco2molch4(row, column);
                            Mco2molch4(row, column) = 0.0;
                            Mco2molco2(row, row) += Mco2molco2(row, column);
                            Mco2molco2(row, column) = 0.0;
                            Mco2pc(row, row) += Mco2pc(row, column);
                            Mco2pc(row, column) = 0.0;

                            Mairpg(row, row) += Mairpg(row, column);
                            Mairpg(row, column) = 0.0;
                            Mairmolh2(row, row) += Mairmolh2(row, column);
                            Mairmolh2(row, column) = 0.0;
                            Mairmolch4(row, row) += Mairmolch4(row, column);
                            Mairmolch4(row, column) = 0.0;
                            Mairmolco2(row, row) += Mairmolco2(row, column);
                            Mairmolco2(row, column) = 0.0;
                            Mairpc(row, row) += Mairpc(row, column);
                            Mairpc(row, column) = 0.0;

                            Mh2opg(row, row) += Mh2opg(row, column);
                            Mh2opg(row, column) = 0.0;
                            Mh2omolh2(row, row) += Mh2omolh2(row, column);
                            Mh2omolh2(row, column) = 0.0;
                            Mh2omolch4(row, row) += Mh2omolch4(row, column);
                            Mh2omolch4(row, column) = 0.0;
                            Mh2omolco2(row, row) += Mh2omolco2(row, column);
                            Mh2omolco2(row, column) = 0.0;
                            Mh2opc(row, row) += Mh2opc(row, column);
                            Mh2opc(row, column) = 0.0;
                        }
                    }
                }
            }  // end of mass-lumping
        }

    }  // end of namespace
}  // end of namespace
