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

subroutine NatOrb_MCLR(Dens,CMOO,CMON,OCCN)

use MCLR_Data, only: ipCM, ipMat, nDens
use input_mclr, only: kPrint, nBas, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: Dens(*), CMOO(*)
real(kind=wp), intent(_OUT_) :: CMON(*), OCCN(*)
integer(kind=iwp) :: i, iEnd, ii, ij, iO, iS, iSt, j
real(kind=wp), allocatable :: EVal(:), EVec(:)

call mma_allocate(EVec,nDens,Label='EVec')
call mma_allocate(EVal,nDens,Label='EVal')

! Diagonalize the density matrix and transform orbitals

if (btest(kprint,3)) then
  write(u6,*)
  write(u6,*) '           Effective natural population'
  write(u6,*) '           ============================'
  write(u6,*)
end if
io = 0
do is=1,nsym
  ij = 0
  do i=0,nbas(is)-1
    do j=0,i
      ij = ij+1
      Eval(ij) = Dens(ipMat(is,is)+i+j*nbas(is))
    end do
  end do
  call unitmat(EVec,nbas(is))
  call JACOB(EVal,EVec,nbas(is),nbas(is))
  ii = 0
  do i=1,nbas(is)
    ii = ii+i
    OCCN(io+i) = Eval(ii)
  end do
  IST = IO+1
  IEND = IO+NBAS(is)
  if (btest(kprint,1)) write(u6,'(6X,A3,I2,A1,10F11.6,/,(12X,10F11.6))') 'sym',iS,':',(OCCN(I),I=IST,IEND)
  if (nBas(is) >= 1) &
    call DGEMM_('N','N',NBAS(is),NBAS(is),NBAS(is),One,CMOO(ipCM(is)),NBAS(is),EVec,NBAS(is),Zero,CMON(ipCM(is)),NBAS(is))
  io = io+nbas(is)
end do

call mma_deallocate(EVec)
call mma_deallocate(Eval)

end subroutine NatOrb_MCLR
