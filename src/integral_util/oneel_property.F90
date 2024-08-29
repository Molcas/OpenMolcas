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

!#define _DEBUGPRINT_
subroutine OneEl_Property(Kernel,KrnlMm,Label,ip,lOper,nComp,CCoor,nOrdOp,rNuc,rHrmt,iChO,D_tot,nDens,Property,Sig)

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Integral_Interfaces, only: int_kernel, int_mem, OneEl_Integrals
use stdalloc, only: mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
procedure(int_kernel) :: Kernel
procedure(int_mem) :: KrnlMm
character(len=8), intent(in) :: Label
integer(kind=iwp), intent(in) :: nComp, lOper(nComp), nOrdOp, iChO(nComp), nDens
integer(kind=iwp), intent(out) :: ip(nComp)
real(kind=wp), intent(in) :: CCoor(3,nComp), rNuc(nComp), rHrmt, D_tot(nDens), Sig
real(kind=wp), intent(out) :: Property(nComp)
integer(kind=iwp) :: iComp, iSmLbl, LenTot, n_Int
real(kind=wp), allocatable :: Integrals(:)
integer(kind=iwp), external :: n2Tri
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
if (rHrmt /= One) then
  call WarningMessage(2,'OneEl_Property: rHrmt /= One')
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the one-electron integrals

call OneEl_Integrals(Kernel,KrnlMm,Label,ip,lOper,nComp,CCoor,nOrdOp,rHrmt,iChO,Integrals)
!                                                                      *
!***********************************************************************
!                                                                      *
!                    P O S T P R O C E S S I N G                       *
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call PrMtrx(Label,lOper,nComp,ip,Integrals)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute properties

LenTot = 0
do iComp=1,nComp
  iSmLbl = lOper(iComp)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute properties directly from integrals

  n_Int = n2Tri(iSmLbl)
  LenTot = LenTot+n_Int+4
  if (n_Int /= 0) then
    call CmpInt(Integrals(ip(iComp)),n_Int,nBas,nIrrep,iSmLbl)
    if (n_Int /= nDens) then
      call WarningMessage(2,'OneEl_Property: n_Int /= nDens')
      write(u6,*) 'n_Int=',n_Int
      write(u6,*) 'nDens',nDens
      call Abend()
    end if
    Property(iComp) = rNuc(iComp)-Sig*DDot_(nDens,D_Tot,1,Integrals(ip(iComp)),1)
  else
    Property(iComp) = rNuc(iComp)
  end if

end do  ! iComp
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate memory for integral

call mma_deallocate(Integrals)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine OneEl_Property
