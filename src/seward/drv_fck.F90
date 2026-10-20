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
! Copyright (C) 1990,1991,1993,1999, Roland Lindh                      *
!               1990, IBM                                              *
!***********************************************************************

subroutine Drv_Fck(Label,ip,lOper,nComp,CCoor,nOrdOp,rNuc,rHrmt,iChO,opmol,ipad,opnuc,iopadr,idirect,isyop,PtChrg,nGrid,iAddPot)

use PAM2, only: iPAMcount
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
character(len=8), intent(inout) :: Label
integer(kind=iwp), intent(in) :: nComp, lOper(nComp), nOrdOp, iChO(nComp), ipad, iopadr(*), idirect, isyop, nGrid, iAddPot
integer(kind=iwp), intent(out) :: ip(nComp)
real(kind=wp), intent(in) :: CCoor(3,nComp), rNuc(nComp), rHrmt, opmol(*), opnuc(*), PtChrg(nGrid)
integer(kind=iwp) :: iadr, iComp, iIrrep, iOpt, iRC, iSmLbl, iStabO(0:7), LenInt, LenTot, llOper, nIC, nStabO
real(kind=wp), allocatable :: Int1El(:)
integer(kind=iwp), external :: n2Tri

#include "warnings.h"
#include "macros.fh"
unused_var(CCoor)
unused_var(nOrdOp)
unused_var(iChO)
unused_var(opmol(1))
unused_var(opnuc(1))
unused_var(ipad)
unused_var(iopadr(1))
unused_var(idirect)
unused_var(isyop)
unused_var(PtChrg)
unused_var(iAddPot)

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
!-----Compute the number of blocks from each component of the operator
!     and the irreps it will span.

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
if (nIC == 0) return
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
  LenInt = n2Tri(lOper(iComp))
  LenTot = LenTot+LenInt+4
end do
call mma_allocate(Int1El,LenTot)
Int1El(:) = Zero
ip(1) = 1
iadr = ip(1)
do iComp=1,nComp
  LenInt = n2Tri(lOper(iComp))
  ip(icomp) = iadr
  iadr = iadr+LenInt+4
  ! Copy center of operator to work area.
  call DCopy_(3,Ccoor(1,iComp),1,Int1El(ip(iComp)+LenInt),1)
  ! Copy nuclear contribution to work area.
  Int1El(ip(iComp)+LenInt+3) = rNuc(iComp)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Compute all SO integrals for all components of the operator.

call Drv_Fck_Inner(ip,Int1El,LenTot,lOper,nComp,rHrmt,iStabO,nStabO,nIC)
!                                                                      *
!***********************************************************************
!                                                                      *
!                    P O S T P R O C E S S I N G                       *
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call PrMtrx(Label,lOper,nComp,ip,Int1El)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Write integrals to disc.

do iComp=1,nComp
  iSmLbl = lOper(iComp)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !----- Write integrals to disc

  iOpt = 0
  iRC = -1
  if (Label(1:3) == 'PAM') write(Label,'(A5,I3.3)') 'PAM  ',iPAMcount
  !write(u6,*) ' oneel *',Label,'*'

  call WrOne(iRC,iOpt,Label,iComp,Int1El(ip(iComp)),iSmLbl)

  if (Label(1:3) == 'PAM') call WrOne(iRC,iOpt,Label,1,Int1El(ip(iComp)),iSmLbl)
  iPAMcount = iPAMcount+1

  if (iRC /= 0) then
    write(u6,*) ' *** Error in subroutine ONEEL ***'
    write(u6,*) '     Abend in subroutine WrOne'
    call Quit(_RC_IO_ERROR_WRITE_)
  end if
end do  ! iComp
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Deallocate memory for integral

call mma_deallocate(Int1El)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Drv_Fck
