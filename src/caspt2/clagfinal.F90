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

subroutine CLagFinal(nConf,nRoots,nState,CLag,SLag)

use PrintLevel, only: VERBOSE
use caspt2_global, only: iPrGlb
use caspt2_module, only: ISCF, REFENE
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nConf, nRoots, nState
real(kind=wp), intent(inout) :: CLag(nConf,nRoots), SLag(nState**2)
integer(kind=iwp) :: ijst, ilStat, jlStat
real(kind=wp) :: Ovl, Scal
real(kind=wp), allocatable :: CI1(:), CI2(:)
real(kind=wp), external :: DDOT_

call mma_allocate(CI1,nConf,Label='CI1')
call mma_allocate(CI2,nConf,Label='CI2')

!! Construct SLag
ijst = 0
do ilStat=1,nState
  if (ISCF == 0) then
    call LoadCI(CI1,ilStat)
  else
    CI1(1) = One
  end if
  do jlStat=1,ilStat !! -1
    ijst = ilStat+nState*(jlStat-1)
    if (ilStat == jlStat) cycle
    if (ISCF == 0) then
      call LoadCI(CI2,jlStat)
    else
      CI2(1) = One
    end if
    Scal = DDOT_(nConf,CI1,1,CLag(1,jlStat),1)-DDOT_(nConf,CI2,1,CLag(1,ilStat),1)
    Scal = Scal/(REFENE(jlStat)-REFENE(ilStat))
    SLag(ijst) = SLag(ijst)+Scal
    if (IPRGLB >= VERBOSE) then
      write(u6,*)
      write(u6,'(1x,"SLag for State ",i1,"-",i1," = ",f20.10)') ilstat,jlstat,slag(ijst)
      write(u6,*)
    end if
  end do
end do

!! This projection is required to get convergence in MCLR.
!! also the linear equation for non-invariant CASPT2
do ilStat=1,nState
  CI1(1:nConf) = CLag(1:nConf,ilStat)
  !do i=1,nconf
  !  write(u6,'(i3,f20.10)') i,clag(i,ilstat)
  !end do
  do jlStat=1,nState
    if (ISCF == 0) then
      call LoadCI(CI2,jlStat)
    else
      CI2(1) = One
    end if
    Ovl = DDot_(nConf,CI1,1,CI2,1)
    !write(u6,*) 'projection coeff = ',ovl
    CLag(1:nConf,ilStat) = CLag(1:nConf,ilStat)-Ovl*CI2(1:nConf)
  end do
  !write(u6,*) 'clag after projection'
  !write(u6,*) 'state = ',ilstat
  !do i=1,nconf
  !  write(u6,'(i3,f20.10)') i,clag(i,ilstat)
  !end do
end do

call mma_deallocate(CI1)
call mma_deallocate(CI2)

return

end subroutine CLagFinal
