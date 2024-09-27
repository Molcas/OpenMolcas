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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

!#define _DEBUGPRINT_
subroutine OneEl_Integrals(Kernel,KrnlMm,Label,ip,lOper,nComp,CCoor,nOrdOp,rHrmt,iChO,Integrals)

use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
procedure(int_kernel) :: Kernel
procedure(int_mem) :: KrnlMm
character(len=8), intent(in) :: Label
integer(kind=iwp), intent(in) :: nComp, lOper(nComp), nOrdOp, iChO(nComp)
integer(kind=iwp), intent(out) :: ip(nComp)
real(kind=wp), intent(in) :: CCoor(3,nComp), rHrmt
real(kind=wp), allocatable, intent(out) :: Integrals(:)
integer(kind=iwp) :: iComp, iIrrep, iStabO(0:7), LenInt, LenTot, llOper, nIC, nStabO
real(kind=wp) :: dum(1)
integer(kind=iwp), external :: n2Tri

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) ' In OneEl: Label',Label
write(u6,*) ' In OneEl: nComp'
write(u6,'(1X,8I5)') nComp
write(u6,*) ' In OneEl: lOper'
write(u6,'(1X,8I5)') lOper
write(u6,*) ' In OneEl: n2Tri'
do iComp=1,nComp
  ip(iComp) = n2Tri(lOper(iComp))
end do
write(u6,'(1X,8I5)') (ip(iComp),iComp=1,nComp)
call RecPrt(' CCoor',' ',CCoor,3,nComp)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the number of blocks from each component of the operator
! and the irreps it will span.

nIC = 0
llOper = 0
do iComp=1,nComp
  llOper = ior(llOper,lOper(iComp))
  do iIrrep=0,nIrrep-1
    if (btest(lOper(iComp),iIrrep)) nIC = nIC+1
  end do
end do
#ifdef _DEBUGPRINT_
write(u6,*) ' nIC =',nIC
#endif
if (nIC == 0) then
  call WarningMessage(2,'OneEl_Integrals: nIC == 0')
  call Abend()
end if
call SOS(iStabO,nStabO,llOper)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate memory for symmetry adapted one electron integrals.
! Will just store the unique elements, i.e. low triangular blocks
! and lower triangular elements in the diagonal blocks.

ip(:) = -1
LenTot = 0
do iComp=1,nComp
  ip(iComp) = 1+LenTot
  LenInt = n2Tri(lOper(iComp))
  LenTot = LenTot+LenInt+4
end do
call mma_allocate(Integrals,LenTot,Label='Integrals')
Integrals(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute all SO integrals for all components of the operator.

call OneEl_Inner(Kernel,KrnlMm,Label,ip,lOper,nComp,CCoor,nOrdOp,rHrmt,iChO,iStabO,nStabO,nIC,Dum,1,0,Integrals,LenTot)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine OneEl_Integrals
