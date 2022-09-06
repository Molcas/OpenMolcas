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

use Chunk_Mod
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif

implicit real*8(A-H,O-Z)
real*8 Rv(nMuNu,(j_e-j_s+1))

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _MOLCAS_MPP_

NumVec_ = j_e-j_s+1
if (NumVec_ > 0) then
  if (Is_Real_Par()) then
    call GA_Put(ip_Chunk,MuNu_s,MuNu_e,j_s,j_e,Rv,nMuNu)
  else
    mMuNu = MuNu_s-1
    jp_ChoVec = 1+mMuNu
    do jVec=1,NumVec_
      call dcopy_(nMuNu,Rv(1,jVec),1,Chunk(jp_ChoVec),1)
      jp_ChoVec = jp_ChoVec+LenVec
    end do
  end if
end if

#else

mMuNu = MuNu_s-1
NumVec_ = j_e-j_s+1

jp_ChoVec = 1+mMuNu
do jVec=1,NumVec_
  call dcopy_(nMuNu,Rv(1,jVec),1,Chunk(jp_ChoVec),1)
  jp_ChoVec = jp_ChoVec+LenVec
end do

! Avoid unused argument warnings
if (.false.) call Unused_integer(MuNu_e)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Put_Chunk
