#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#***********************************************************************

set (sources
     backtrans_k.F90
     covarmatrix.F90
     covarvector.F90
     dispersion_kriging.F90
     dispersion_kriging_layer.F90
     energy_kriging.F90
     energy_kriging_layer.F90
     finishing_kriging.F90
     gradient_kriging.F90
     gradient_kriging_layer.F90
     hessian_kriging.F90
     hessian_kriging_layer.F90
     kriging_model.F90
     kriging_mod.F90
     matderiv.F90
     matern.F90
     pgek.F90
     predict.F90
     set_l_kriging.F90
     start_kriging.F90
     trans_k.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.txt)
