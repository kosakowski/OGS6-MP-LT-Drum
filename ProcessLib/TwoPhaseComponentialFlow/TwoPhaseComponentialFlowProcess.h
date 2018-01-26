/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Process.h"
#include "TwoPhaseComponentialFlowLocalAssembler.h"

namespace MeshLib
{
    class Element;
    class Mesh;
    template <typename PROP_VAL_TYPE>
    class PropertyVector;
}

namespace ProcessLib
{
    namespace TwoPhaseComponentialFlow
    {
        /**
        * \brief A class to simulate the two-phase flow process with P-rho model in
        * porous media
        */
        class TwoPhaseComponentialFlowProcess final : public Process
        {
        public:
            TwoPhaseComponentialFlowProcess(
                MeshLib::Mesh& mesh,
                std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
                std::vector<std::unique_ptr<ParameterBase>> const& parameters,
                unsigned const integration_order,
                std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
                process_variables,
                TwoPhaseComponentialFlowProcessData&& process_data,
                SecondaryVariableCollection&& secondary_variables,
                NumLib::NamedFunctionCaller&& named_function_caller,
                BaseLib::ConfigTree const& config,
                std::map<std::string,
                std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
                curves);

            bool isLinear() const override { return false; }
        private:
            void initializeConcreteProcess(
                NumLib::LocalToGlobalIndexMap const& dof_table,
                MeshLib::Mesh const& mesh, unsigned const integration_order) override;

            void assembleConcreteProcess(const double t, GlobalVector const& x,
                GlobalMatrix& M, GlobalMatrix& K,
                GlobalVector& b) override;

            void assembleWithJacobianConcreteProcess(
                const double t, GlobalVector const& x, GlobalVector const& xdot,
                const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
                GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

            void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                double const dt,
                const int /*process_id*/) override
            {
                DBUG("PreTimestep TwoPhaseCarbonation.");

                _process_data._dt = dt;

                GlobalExecutor::executeMemberOnDereferenced(
                    &TwoPhaseComponentialFlowLocalAssemblerInterface::preTimestep,
                    _local_assemblers, *_local_to_global_index_map, x, t, dt);
            }

            TwoPhaseComponentialFlowProcessData _process_data;

            std::vector<std::unique_ptr<TwoPhaseComponentialFlowLocalAssemblerInterface>>
                _local_assemblers;
        };

    }  // end of namespace
}  // end of namespace
