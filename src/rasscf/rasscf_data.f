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

! This is just an encapsulation of the common block in
! src/Include/rasscf.fh
! src/Include/rasdim.fh
! into a data module

      module rasscf_data
      implicit none
! Order of inclusion matters!
#include "rasdim.fh"
#include "rasscf.fh"
      save
      end module rasscf_data
