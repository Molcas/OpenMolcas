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

! Not properly worked through. Do not use!
subroutine Mbpt2Corr(nBas,Cmo)

use qmstat_global, only: DenCorrD, iOrb, iPrint, Trace_MP2
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBas
real(kind=wp), intent(in) :: Cmo(nBas,nBas)
integer(kind=iwp) :: i, iB1, iT, j, kaunt1
real(kind=wp) :: Det
real(kind=wp), allocatable :: Diff(:), Inv(:,:), RedSq(:,:), SqD(:,:), SqE(:), TEMP(:,:)
#include "warnings.h"

write(u6,*)
write(u6,*) 'MP2 density correction is requested.'
write(u6,*) ' -- perturbative correlation correction to the solute density.'

! No-no zone!

write(u6,*)
write(u6,*) 'THIS OPTION IS NOT PROPERLY WORKED THROUGH! SHOULD NOT BE USED!'
call Quit(_RC_GENERAL_ERROR_)
! Check that the density difference is sound.
iT = nTri_Elem(nBas)
call mma_allocate(Diff,iT,Label='Diff')
call Get_dArray_chk('D1ao',Diff,iT)
if (iPrint >= 10) call TriPrt('Non-reduced difference density matrix',' ',Diff,nBas)
! Transform density difference to orbital basis.
call mma_allocate(SqD,nBas,nBas,label='SqDenA')
call mma_allocate(SqE,nBas**2,label='SqDenM')
call mma_allocate(TEMP,nBas,nBas,label='TEMP')
call mma_allocate(Inv,nBas,nBas,label='Inv')
call mma_allocate(RedSq,nBas,nBas,label='RedSq')
! Do not forget the density matrix convention in Molcas.
call Dsq(Diff,SqD,1,nBas,nBas)
! Inverse of orbital file and transformation.
call Minv(Cmo,Inv,Det,nBas)
call Dgemm_('N','N',nBas,nBas,nBas,One,Inv,nBas,SqD,nBas,Zero,TEMP,nBas)
call Dgemm_('N','T',nBas,nBas,nBas,One,TEMP,nBas,Inv,nBas,Zero,SqE,nBas)
! Remove all except the suck-out orbitals.
do i=1,nBas
  do j=1,nBas
    if ((i <= iOrb(1)) .and. (j <= iOrb(1))) then
      RedSq(j,i) = SqE(j+(i-1)*nBas)
    else
      RedSq(j,i) = Zero
    end if
  end do
end do
! Make a check of the trace. Should be small.
Trace_MP2 = Zero
do iB1=1,nBas
  Trace_MP2 = Trace_MP2+RedSq(iB1,iB1)
end do
if (iPrint >= 10) then
  write(u6,*) 'Trace: ',Trace_MP2
end if
! Make things a bit more tidy.
kaunt1 = 0
do i=1,iOrb(1)
  do j=1,iOrb(1)
    kaunt1 = kaunt1+1
    SqE(kaunt1) = RedSq(j,i)
  end do
end do
call mma_allocate(DenCorrD,nTri_Elem(iOrb(1)),label='DenCorrD')
call SqToTri_q(SqE,DenCorrD,iOrb(1))

! Transform back if we want to keep things in AO-basis. Not
! used in QMSTAT at the present. If you wish, comment away the
! code below 'make things a bit more tidy' and you are in
! ready to rumble.
!call Dgemm_('N','N',nBas,nBas,nBas,One,Cmo,nBas,RedSq,nBas,Zero,TEMP,nBas)
!call Dgemm_('N','T',nBas,nBas,nBas,One,TEMP,nBas,Cmo,nBas,Zero,SqE,nBas)
!do i=1,nBas
!  do j=1,nBas
!    if (i /= j) SqE(j+(i-1)*nBas) = SqE(j+(i-1)*nBas)*Two
!  end do
!end do
!call SqToTri_q(SqE,TrDiffD,nBas)
!if (iPrint >= 10) call TriPrt('Reduced difference density matrix',' ',TrDiffD,nBas)

call mma_deallocate(Diff)
call mma_deallocate(SqD)
call mma_deallocate(SqE)
call mma_deallocate(TEMP)
call mma_deallocate(Inv)
call mma_deallocate(RedSq)

return

end subroutine Mbpt2Corr
