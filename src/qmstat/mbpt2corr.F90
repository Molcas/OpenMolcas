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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "maxi.fh"
#include "qm1.fh"
#include "qminp.fh"
#include "WrkSpc.fh"
#include "warnings.h"
integer(kind=iwp) :: nBas
real(kind=wp) :: Cmo(MxBas**2)
integer(kind=iwp) :: i, iB1, iI, ipSqD, ipSqE, ipTEMP, iRedSq, Ising, iT, j, jjj, kaunt1, kaunt2, kaunter
real(kind=wp) :: Det
real(kind=wp), allocatable :: Diff(:)

write(u6,*)
write(u6,*) 'MP2 density correction is requested.'
write(u6,*) ' -- perturbative correlation correction to the solute density.'

! No-no zone!

write(u6,*)
write(u6,*) 'THIS OPTION IS NOT PROPERLY WORKED THROUGH! SHOULD NOT BE USED!'
call Quit(_RC_GENERAL_ERROR_)
! Check that the density difference is sound.
iT = nBas*(nBas+1)/2
call mma_allocate(Diff,iT,Label='Diff')
call Get_D1ao(Diff,iT)
if (iPrint >= 10) then
  call TriPrt('Non-reduced difference density matrix',' ',Diff,nBas)
end if
! Transform density difference to orbital basis.
call GetMem('SqDenA','Allo','Real',ipSqD,nBas**2)
call GetMem('SqDenM','Allo','Real',ipSqE,nBas**2)
call GetMem('TEMP','Allo','Real',ipTEMP,nBas**2)
call GetMem('Inv','Allo','Real',iI,nBas**2)
call GetMem('RedSq','Allo','Real',iRedSq,nBas**2)
call dcopy_(nBas**2,[Zero],0,Work(ipSqD),1)
call dcopy_(iOrb(1)**2,[Zero],0,Work(ipSqE),1)
call dcopy_(nBas*iOrb(1),[Zero],0,Work(ipTEMP),1)
! Do not forget the density matrix convention in Molcas.
call Dsq(Diff,Work(ipSqD),1,nBas,nBas)
! Inverse of orbital file and transformation.
call Minv(Cmo,Work(iI),Ising,Det,nBas)
call Dgemm_('N','N',nBas,nBas,nBas,One,Work(iI),nBas,Work(ipSqD),nBas,Zero,Work(ipTEMP),nBas)
call Dgemm_('N','T',nBas,nBas,nBas,One,Work(ipTEMP),nBas,Work(iI),nBas,Zero,Work(ipSqE),nBas)
! Remove all except the suck-out orbitals.
kaunt1 = 0
do i=1,nBas
  do j=1,nBas
    if ((i <= iOrb(1)) .and. (j <= iOrb(1))) then
      Work(iRedSq+kaunt1) = Work(ipSqE+kaunt1)
    else
      Work(iRedSq+kaunt1) = Zero
    end if
    kaunt1 = kaunt1+1
  end do
end do
! Make a check of the trace. Should be small.
kaunter = 0
Trace_MP2 = Zero
do iB1=1,nBas
  do jjj=1,nBas
    if (iB1 == jjj) Trace_MP2 = Trace_MP2+Work(iRedSq+kaunter)
    kaunter = kaunter+1
  end do
end do
if (iPrint >= 10) then
  write(u6,*) 'Trace: ',Trace_MP2
end if
! Make things a bit more tidy.
kaunt1 = 0
kaunt2 = 0
do i=1,iOrb(1)
  do j=1,nBas
    if (j <= iOrb(1)) then
      Work(ipSqE+kaunt1) = Work(iRedSq+kaunt2)
      kaunt1 = kaunt1+1
    end if
    kaunt2 = kaunt2+1
  end do
end do
call SqToTri_q(Work(ipSqE),DenCorrD,iOrb(1))

! Transform back if we want to keep things in AO-basis. Not
! used in QMSTAT at the present. If you wish, comment away the
! code below 'make things a bit more tidy' and you are in
! ready to rumble.
!call Dgemm_('N','N',nBas,nBas,nBas,One,Cmo,nBas,Work(iRedSq),nBas,Zero,Work(ipTEMP),nBas)
!call Dgemm_('N','T',nBas,nBas,nBas,One,Work(ipTEMP),nBas,Cmo,nBas,Zero,Work(ipSqE),nBas)
!k = 0
!do i=1,nBas
!  do j=1,nBas
!    if (i /= j) Work(ipSqE+k) = Work(ipSqE+k)*2
!      k = k+1
!  end do
!end do
!call SqToTri_q(Work(ipSqE),Work(ipTrDiffD),nBas)
!if (iPrint >= 10) then
!  call TriPrt('Reduced difference density matrix',' ',Work(ipTrDiffD),nBas)
!end if

call mma_deallocate(Diff)
call GetMem('SqDenA','Free','Real',ipSqD,nBas**2)
call GetMem('SqDenM','Free','Real',ipSqE,nBas**2)
call GetMem('TEMP','Free','Real',ipTEMP,nBas**2)
call GetMem('Inv','Free','Real',iI,nBas**2)
call GetMem('RedSq','Free','Real',iRedSq,nBas**2)

return

end subroutine Mbpt2Corr
