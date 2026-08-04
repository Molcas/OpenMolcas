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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------
! 1998  PER-AAKE MALMQUIST
! DEPARTMENT OF THEORETICAL CHEMISTRY
! UNIVERSITY OF LUND
! SWEDEN
!--------------------------------------------

subroutine MKRHS_SAVE(ICASE,ISYM,IVEC,LW)
!SVC: special routine to save the RHS array. MKRHS works in serial, so
! in case of a true parallel run we need to put the local array in a
! global array and then save that to disk in a distributed fashion.

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use fake_GA, only: GA_Arrays
#endif
use caspt2_module, only: NASUP, NISUP
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ICASE, ISYM, IVEC, LW
integer(kind=iwp) :: lg_w, NAS, NIS

NAS = NASUP(ISYM,ICASE)
NIS = NISUP(ISYM,ICASE)

#ifdef _MOLCAS_MPP_
if (IS_REAL_PAR()) then
  call RHS_ALLO(NAS,NIS,lg_W)
  call RHS_PUT(NAS,NIS,lg_W,GA_Arrays(LW)%A)
else
#endif
  lg_W = LW
#ifdef _MOLCAS_MPP_
end if
#endif

call RHS_SAVE(NAS,NIS,lg_W,iCASE,iSYM,iVEC)

#ifdef _MOLCAS_MPP_
if (IS_REAL_PAR()) call RHS_FREE(lg_W)
#endif

end subroutine MKRHS_SAVE
