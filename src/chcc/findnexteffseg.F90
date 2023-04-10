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

subroutine findNextEffSeg(NvGrp,eff,Nprocs,eff_thrs,maxGrp,printkey)

implicit none
integer Nprocs, NvGrp, printkey
real*8 eff_thrs, eff, tmp0, tmp1, tmp2
real*8 dceil
external dceil
integer maxGrp

11 continue
! calculate theoretical efficiency for o2v4 step
tmp0 = NvGrp*NvGrp/2.0d0
tmp1 = dceil(tmp0) ! otazne ... zisti
tmp0 = 1.0d0*tmp1/Nprocs
tmp2 = dceil(tmp0)
eff = (NvGrp*NvGrp/2.0d0)/(1.0d0*tmp2*Nprocs) ! mozno daj tu tmp

! report Np and efficiency
if (printkey >= 10) write(6,'(A,i4,A,f6.2)') 'Efficiency check: ',NvGrp,', efficiency: ',eff*100

if (((eff*100) < eff_thrs) .and. (NvGrp < maxGrp)) then
  NvGrp = NvGrp+1
  goto 11
end if

return

end subroutine findNextEffSeg
