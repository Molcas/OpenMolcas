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

use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: NvGrp
real(kind=wp), intent(out) :: eff
integer(kind=iwp), intent(in) :: Nprocs, maxGrp, printkey
real(kind=wp), intent(in) :: eff_thrs
integer(kind=iwp) :: tmp1, tmp2
real(kind=wp) :: tmp0

do
  ! calculate theoretical efficiency for o2v4 step
  tmp0 = Half*real(NvGrp**2,kind=wp)
  tmp1 = ceiling(tmp0) ! otazne ... zisti
  tmp0 = real(tmp1,kind=wp)/real(Nprocs,kind=wp)
  tmp2 = ceiling(tmp0)
  eff = Half*real(NvGrp**2,kind=wp)/real(tmp2*Nprocs,kind=wp) ! mozno daj tu tmp

  ! report Np and efficiency
  if (printkey >= 10) write(u6,'(A,i4,A,f6.2)') 'Efficiency check: ',NvGrp,', efficiency: ',eff*100

  if (((eff*100) >= eff_thrs) .or. (NvGrp >= maxGrp)) exit
  NvGrp = NvGrp+1
end do

return

end subroutine findNextEffSeg
