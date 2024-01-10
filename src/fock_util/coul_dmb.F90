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

subroutine Coul_DMB(GetFM,nDM,Rep_EN,FM,DMA,DMB,lFDM)

use Cholesky, only: nBas, nSym
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: GetFM
integer(kind=iwp), intent(in) :: nDM, lFDM
real(kind=wp), intent(out) :: Rep_EN
real(kind=wp), intent(inout) :: FM(lFDM)
real(kind=wp), intent(in) :: DMA(lFDM), DMB(lFDM)
integer(kind=iwp) :: irc
type(DSBA_Type) :: DLT(1), FLT(1)
real(kind=wp), external :: ddot_

if ((nDM > 2) .or. (nDM < 1)) then
  write(u6,*) ' In Coul_DMB: wrong value of nDM= ',nDM
  call SysAbendMsg('Coul_DMB ',' nDM must be 1 or 2 ',' ')
end if

if (GetFM) then

  call Allocate_DT(FLT(1),nBas,nBas,nSym,aCase='TRI',Ref=FM)

  call NameRun('AUXRFIL') ! switch RUNFILE name

  call Allocate_DT(DLT(1),nBas,nBas,nSym,aCase='TRI')
  call get_dArray('D1ao',DLT(1)%A0,lFDM)

  FLT(1)%A0(:) = Zero
  call CHO_FOCK_DFT_RED(irc,DLT,FLT)
  if (irc /= 0) then
    call SysAbendMsg('Coul_DMB ',' non-zero rc ',' ')
  end if
  call GADSum(FM,lFDM)

  call Deallocate_DT(DLT(1))
  call Deallocate_DT(FLT(1))

  call NameRun('#Pop')  ! switch back RUNFILE name

end if

Rep_EN = ddot_(lFDM,DMA,1,FM,1)
if (nDM == 2) then
  Rep_EN = Rep_EN+ddot_(lFDM,DMB,1,FM,1)
end if

return

end subroutine Coul_DMB
