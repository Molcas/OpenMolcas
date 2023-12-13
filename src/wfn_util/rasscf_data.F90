!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

! This is just an encapsulation of the common block in
! src/Include/rasdim.fh
! src/Include/rasscf.fh
! src/Include/output_ras.fh
! into a data module

module rasscf_data

implicit none

#include "rasdim.fh"
#include "rasscf.fh"
#include "output_ras.fh"

end module rasscf_data
