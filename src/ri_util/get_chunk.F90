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

subroutine Get_Chunk(LenVec,NumVec_,iChoVec,iSym,iVec_Global)

use RI_glob, only: Chunk
#ifdef _MOLCAS_MPP_
use RI_glob, only: iMap, ip_Chunk
use ChoSwp, only: InfVec
use Para_Info, only: MyRank, Is_Real_Par
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: LenVec, NumVec_, iSym, iVec_Global
integer(kind=iwp), intent(inout) :: iChoVec
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"
integer(kind=iwp) :: Indx, iVec, J_e, J_s, ld, MuNu_e, MuNu_s, nJ
#endif

#ifndef _MOLCAS_MPP_
#include "macros.fh"
unused_var(iVec_Global)
#endif

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  call GA_Sync()
  call Cho_RI_SwapVecUnit(iSym)

  ! Get the subrange and check that we are within range.

  J_s = iMap(1+MyRank)
  if (J_s <= NumVec_) then
    J_e = iMap(1+MyRank+1)-1
    J_e = min(J_e,NumVec_)
    nJ = J_e-J_s+1
    if (nJ > 0) then
      MuNu_s = 1
      MuNu_e = LenVec
      call GA_Access(ip_Chunk,MuNu_s,MuNu_e,J_s,J_e,Indx,ld)
      !call RecPrt('Dbl_Mb(Indx)',' ',Dbl_Mb(Indx),LenVec,nJ)
      call Cho_PutVec(DBL_MB(Indx),LenVec,nJ,iChoVec+1,iSym)
      ! Set mapping from local to global vector index
      do iVec=1,nJ
        InfVec(iChoVec+iVec,5,iSym) = iVec_Global+J_s-2+iVec
      end do
      call GA_Release(ip_Chunk,MuNu_s,MuNu_e,J_s,J_e)
      iChoVec = iChoVec+nJ
    end if
  end if
  call Cho_RI_SwapVecUnit(iSym)
else
#endif
  call Cho_PutVec(Chunk,LenVec,NumVec_,iChoVec+1,iSym)
  iChoVec = iChoVec+NumVec_
#ifdef _MOLCAS_MPP_
end if
#endif

return

end subroutine Get_Chunk
