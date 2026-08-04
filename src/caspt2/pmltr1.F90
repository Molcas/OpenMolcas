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

subroutine PMLTR1(KOD,IMLTOP,LST1,X,nX,NAS1,NIS1,JXOFF,F,NFI,NFJ,lg_Y,NAS2,NIS2)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use GA_Wrapper, only: DBL_MB, GA_NodeId
#endif
use Sigma_data, only: NLST1
use fake_GA, only: GA_Arrays
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: KOD, IMLTOP, LST1(4,NLST1), nX, NAS1, NIS1, JXOFF, NFI, NFJ, lg_Y, NAS2, NIS2
real(kind=wp), intent(inout) :: X(nX), F(NFI,NFJ)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: iYHi, iYLo, jYHi, jYLo, LDY, mY, myRank

! SVC: Determine the index ranges of the local chunks of lg_X and lg_Y.
! The boundaries and leading dimension are stored in a common block for
! access inside the lower-level routines.
! For now, only case H is handled as a distributed array, which is
! always the Y array.
if (Is_Real_Par()) then
  call GA_Sync()
  myRank = GA_NodeID()
  !call GA_Distribution(lg_X,myRank,iXLo,iXHi,jXLo,jXHi)
  !if ((iXLo /= 0) .and. (jXLo /= 0)) call GA_Access(lg_X,iXLo,iXHi,jXLo,jXHi,mX,LDX)
  call GA_Distribution(lg_Y,myRank,iYLo,iYHi,jYLo,jYHi)
  if ((iYLo /= 0) .and. (jYLo /= 0)) then
    call GA_Access(lg_Y,iYLo,iYHi,jYLo,jYHi,mY,LDY)
    if ((KOD == 17) .or. (KOD == 18)) then
      call MLTR1_EH(IMLTOP,LST1,X,NAS1,NIS1,JXOFF,F,NFI,NFJ,DBL_MB(mY),NAS2,jYLo,jYHi)
    else if ((KOD == 21) .or. (KOD == 22)) then
      call MLTR1_GH(IMLTOP,LST1,X,NAS1,NIS1,JXOFF,F,NFI,NFJ,DBL_MB(mY),NAS2,jYLo,jYHi)
    end if
    call GA_Release_Update(lg_Y,iYLo,iYHi,jYLo,jYHi)
  end if
  call GA_Sync()
else
#endif
  if ((KOD == 17) .or. (KOD == 18)) then
    call MLTR1_EH(IMLTOP,LST1,X,NAS1,NIS1,JXOFF,F,NFI,NFJ,GA_Arrays(lg_Y)%A,NAS2,1,NIS2)
  else if ((KOD == 21) .or. (KOD == 22)) then
    call MLTR1_GH(IMLTOP,LST1,X,NAS1,NIS1,JXOFF,F,NFI,NFJ,GA_Arrays(lg_Y)%A,NAS2,1,NIS2)
  end if
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine PMLTR1
