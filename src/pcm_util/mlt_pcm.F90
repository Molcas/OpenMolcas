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

subroutine Mlt_PCM(nAt,nTs,nOrdOp,Tessera,AtmC,V,EF_n,EF_e)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAt, nTs, nOrdOp
real(kind=wp), intent(in) :: Tessera(4,nTs), AtmC(3,nAt)
real(kind=wp), intent(out) :: V(nTs), EF_n(3,nTs)
real(kind=wp), intent(inout) :: EF_e(3,nTs)
integer(kind=iwp) :: iTs, nDens
real(kind=wp) :: Temp(3)
logical(kind=iwp) :: Found
integer(kind=iwp), allocatable :: lOper(:)
real(kind=wp), allocatable :: Chrg(:), D1ao(:), FactOp(:)

call mma_allocate(Chrg,nat)
call Get_dArray('Nuclear charge',Chrg,nAt)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the desired multipole on the tiles

! 1) The nuclear contribution

do iTs=1,nTs
  call EFNuc(Tessera(1,iTs),Chrg,AtmC,nAt,Temp,nOrdOp)
  if (nOrdOp == 0) then
    V(iTs) = Temp(1)
  else if (nOrdOp == 1) then
    EF_n(:,iTs) = Temp(:)
  end if
end do

call mma_deallocate(Chrg)

! 2) The electronic contribution

! Get the total 1st order AO density matrix

call Qpg_dArray('D1ao',Found,nDens)
if (Found .and. (nDens /= 0)) then
  call mma_allocate(D1ao,nDens,Label='D1ao')
else
  write(u6,*) 'Mlt_pcm: D1ao not found.'
  call Abend()
end if
call Get_dArray_chk('D1ao',D1ao,nDens)

call mma_allocate(FactOp,nTs,label='FactOp')
call mma_allocate(lOper,nTs,label='lOper')
FactOp(:) = One
lOper(:) = 255

call drv_ef_PCM(FactOp,nTs,D1ao,nDens,Tessera,lOper,EF_e,nOrdOp)
if (nOrdOp == 0) then
  do iTs=1,nTs
    V(iTs) = Ef_e(1,iTs)
  end do
end if

call mma_deallocate(D1ao)
call mma_deallocate(FactOp)
call mma_deallocate(lOper)

return

end subroutine Mlt_PCM
