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

subroutine Put_Chunk(MuNu_s,MuNu_e,j_s,j_e,Rv,nMuNu,LenVec)

use RI_glob, only: Chunk
#ifdef _MOLCAS_MPP_
use RI_glob, only: ip_Chunk
use Para_Info, only: Is_Real_Par
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: MuNu_s, MuNu_e, j_s, j_e, nMuNu, LenVec
real(kind=wp), intent(inout) :: Rv(nMuNu,j_e-j_s+1)
integer(kind=iwp) :: jp_ChoVec, jVec, mMuNu, NumVec_

#ifndef _MOLCAS_MPP_
#include "macros.fh"
unused_var(MuNu_e)
#endif

!                                                                      *
!***********************************************************************
!                                                                      *

NumVec_ = j_e-j_s+1

#ifdef _MOLCAS_MPP_
if (NumVec_ > 0) then
  if (Is_Real_Par()) then
    call GA_Put(ip_Chunk,MuNu_s,MuNu_e,j_s,j_e,Rv,nMuNu)
  else
#endif
    mMuNu = MuNu_s-1
    jp_ChoVec = 1+mMuNu
    do jVec=1,NumVec_
      Chunk(jp_ChoVec:jp_ChoVec+nMuNu-1) = Rv(:,jVec)
      jp_ChoVec = jp_ChoVec+LenVec
    end do
#ifdef _MOLCAS_MPP_
  end if
end if
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Put_Chunk
