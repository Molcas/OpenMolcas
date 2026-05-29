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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine CASPT2_ResD(Mode,NIN,NIS,lg_W1,lg_W2,DIN,DIS)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use fake_GA, only: GA_Arrays
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Mode, NIN, NIS, lg_W1, lg_W2
real(kind=wp), intent(in) :: DIN(NIN), DIS(NIS)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: iHi1, iHi2, iLo1, iLo2, jHi1, jHi2, jLo1, jLo2, LDW1, LDW2, mW1, mW2, myRank, NCOL, NROW
#include "global.fh"
#include "mafdecls.fh"
#endif

! Apply the resolvent of the diagonal part of H0 to an RHS array

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  call GA_Sync()
  myRank = GA_NodeID()
  !-SVC: get the local vertical stripes of the lg_W vector
  call GA_Distribution(lg_W1,myRank,iLo1,iHi1,jLo1,jHi1)
  call GA_Distribution(lg_W2,myRank,iLo2,iHi2,jLo2,jHi2)
  !! Well, assume the same dimension
  if ((iLo1 > 0) .and. (jLo1 > 0) .and. (iLo2 > 0) .and. (jLo2 > 0)) then
    NROW = iHi1-iLo1+1
    NCOL = jHi1-jLo1+1
    call GA_Access(lg_W1,iLo1,iHi1,jLo1,jHi1,mW1,LDW1)
    call GA_Access(lg_W2,iLo2,iHi2,jLo2,jHi2,mW2,LDW2)
    call CASPT2_ResD2(MODE,NROW,NCOL,DBL_MB(mW1),DBL_MB(mW2),LDW1,DIN(iLo1),DIS(jLo1))
    call GA_Release_Update(lg_W1,iLo1,iHi1,jLo1,jHi1)
    call GA_Release_Update(lg_W2,iLo2,iHi2,jLo2,jHi2)
  end if
  call GA_Sync()
else
#endif
  call CASPT2_ResD2(MODE,NIN,NIS,GA_Arrays(lg_W1)%A,GA_Arrays(lg_W2)%A,NIN,DIN,DIS)
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine CASPT2_ResD
