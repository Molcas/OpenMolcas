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

subroutine SIGMADET_CVB(C,HC,IREFSM,PERMS2,NCI)

implicit real*8(A-H,O-Z)
logical PERMS2
#include "WrkSpc.fh"
#include "rasscf_lucia.fh"
dimension C(NCI), HC(NCI)
dimension DUMMY(1)

! Export arguments to be used in sigma_master_cvb

C_POINTER = KCI_POINTER
call DCOPY_(NCI,C,1,WORK(C_POINTER),1)
! Call the sigma routine
call LUCIA_UTIL('SIGMA_CVB',IREFSM,IDUMMY,DUMMY)
call DCOPY_(NCI,WORK(KSIGMA_POINTER),1,HC,1)

return
! Avoid unused argument warnings
if (.false.) call Unused_logical(PERMS2)

end subroutine SIGMADET_CVB
