************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Logical Function LDF_IntegralPrescreeningInfoIsSet()
      Implicit None
#include "ldf_integral_prescreening_info.fh"
      LDF_IntegralPrescreeningInfoIsSet=
     &    l_GDiag_1C.gt.0 .or. l_GDiag_1C_Mx.gt.0 .or.
     &    l_GDiag_2C.gt.0 .or. l_GDiag_2C_Mx.gt.0 .or.
     &    l_IDiag.gt.0    .or. l_IDiag_Mx.gt.0
      End
