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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine depend_cvb(chr1,chr2)

use casvb_global, only: i_dep_on_j, ioffs, iprint, j_dep_on_i, joffs, ndep_ij, ndep_ji, nobj
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: chr1, chr2
integer(kind=iwp) :: ii

call mkafter_cvb(chr1,chr2)
call touchdepend_cvb(chr1,chr2)

if (iprint >= 10) then
  write(u6,*) ' IOFFS :',(ioffs(ii),ii=1,nobj+1)
  write(u6,*) ' JOFFS :',(joffs(ii),ii=1,nobj+1)
  write(u6,*) ' I_DEP_ON_J :',(i_dep_on_j(ii),ii=1,ndep_ij)
  write(u6,*) ' J_DEP_ON_I :',(j_dep_on_i(ii),ii=1,ndep_ji)
end if

return

end subroutine depend_cvb
