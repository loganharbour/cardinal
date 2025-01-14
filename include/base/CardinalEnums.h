/********************************************************************/
/*                  SOFTWARE COPYRIGHT NOTIFICATION                 */
/*                             Cardinal                             */
/*                                                                  */
/*                  (c) 2021 UChicago Argonne, LLC                  */
/*                        ALL RIGHTS RESERVED                       */
/*                                                                  */
/*                 Prepared by UChicago Argonne, LLC                */
/*               Under Contract No. DE-AC02-06CH11357               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*             Prepared by Battelle Energy Alliance, LLC            */
/*               Under Contract No. DE-AC07-05ID14517               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*                 See LICENSE for full restrictions                */
/********************************************************************/

#pragma once

#include "MooseEnum.h"

MooseEnum getNekOrderEnum();
MooseEnum getBinnedVelocityComponentEnum();
MooseEnum getNekFieldEnum();
MooseEnum getOperationEnum();
MooseEnum getTallyTypeEnum();
MooseEnum getEigenvalueEnum();
MooseEnum getChannelTypeEnum();
MooseEnum getRelaxationEnum();

namespace order
{
  /// Enumeration of possible surface order reconstructions for nekRS solution transfer
  enum NekOrderEnum
  {
    first,
    second
  };
}

namespace component
{
  /// Directions in which to evaluate velocity
  enum BinnedVelocityComponentEnum
  {
    normal,
    user
  };
}

namespace field
{
  /// Enumeration of possible fields to read from nekRS
  enum NekFieldEnum
  {
    velocity_component,
    velocity_x,
    velocity_y,
    velocity_z,
    velocity,
    temperature,
    pressure,
    unity
  };

  /// Enumeration of possible fields to write in nekRS
  enum NekWriteEnum
  {
    flux,
    heat_source,
    x_displacement,
    y_displacement,
    z_displacement
  };
}

namespace operation
{
  /// Enumeration of possible operations to perform in global postprocessors
  enum OperationEnum
  {
    max,
    min
  };
}

namespace tally
{
  /// Type of tally to construct for the OpenMC model
  enum TallyTypeEnum
  {
    cell,
    mesh
  };
}

namespace coupling
{
  /// Type of feedback in Monte Carlo simulation
  enum CouplingFields
  {
    temperature,
    density_and_temperature,
    none
  };
}

namespace eigenvalue
{
  /// Type of OpenMC k-eigenvalue global tally
  enum EigenvalueEnum
  {
    collision,
    absorption,
    tracklength,
    combined
  };
}

namespace channel_type
{
  /// Type of subchannel
  enum ChannelTypeEnum
  {
    interior,
    edge,
    corner
  };
}

namespace relaxation
{
  /// Type of relaxation
  enum RelaxationEnum
  {
    constant,
    robbins_monro,
    dufek_gudowski,
    none
  };
}
