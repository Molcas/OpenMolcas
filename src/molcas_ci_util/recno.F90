!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

function RecNo(itype,iRoot)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Compute the Record number of a vector                            *
!                                                                      *
!     calling arguments:                                               *
!     itype   : integer                                                *
!               vector type: 1 = H_diag                                *
!                            2 = CI_vec                                *
!                            3 = Sig_vec                               *
!     iRoot   : integer                                                *
!               root number                                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use davctl_mod, only: n_Roots, nkeep
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: RecNo
integer(kind=iwp), intent(in) :: itype, iRoot
integer(kind=iwp), external :: PageNo
#include "rasdim.fh"

select case (itype)
  case (1)
    !H_diag
    RecNo = 1
  case (2)
    !CI_vec
    RecNo = 1+PageNo(iRoot)
  case (3)
    !Sig_vec
    RecNo = 1+nkeep+PageNo(iRoot)
  case (4)
    !tmp_CI_vec
    RecNo = 1+2*nKeep+iRoot
  case (5)
    !tmp_Sig_vec
    RecNo = 1+2*nKeep+n_Roots+iRoot
  case default
    RecNo = 0
    write(u6,*) 'RecNo: itype does not match'
    write(u6,*) 'itype = ',itype
    call Abend()
end select

return

end function RecNo
