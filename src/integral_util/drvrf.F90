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
subroutine DrvRF(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq,iCharge)

use External_Centers, only: iXPolType
use rctfld_module, only: lLangevin, lRF, lRFCav, PCM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nh1
real(kind=wp), intent(inout) :: h1(nh1), TwoHam(nh1), RepNuc
real(kind=wp), intent(in) :: D(nh1)
logical(kind=iwp), intent(in) :: First, Dff, NonEq
integer(kind=iwp), intent(inout) :: iCharge
integer(kind=iwp) :: iComp, iOpt, iRC, iSyLbl
real(kind=wp) :: EEE, ERfSelf, RepNuc_RF, RepNucXX(1)
real(kind=wp), save :: RepNuc_Temp
character(len=8) :: Label
real(kind=wp), allocatable :: h1_RF(:), h1_XX(:), RFld(:,:)
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. lRF) return

call Set_Basis_Mode('Valence')
call Setup_iSD()

call Init_RctFld(NonEq,iCharge)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(RFld,nh1,2,Label='RFld')
RFld(:,2) = Zero
if (First) RepNuc_Temp = RepNuc
!                                                                      *
!***********************************************************************
!                                                                      *
if (lLangevin .or. (iXPolType > 0)) then

  ! Reaction field a la polarizabilities and Langevin dipole
  ! moments on a grid in a cavity in a dielectric medium.

  call Langevin(h1,RFld(:,2),D,RepNuc,nh1,First,Dff,NonEq)

else if (PCM) then

  ! Reaction field a la PCM

  call DrvPCM(h1,RFld(:,2),D,RepNuc,nh1,First,Dff,NonEq)

else if (lRFCav) then

  ! Reaction field a la cavity in dielectric medium.

  call RctFld(h1,RFld(:,2),D,RepNuc,nh1,First,Dff,NonEq)

else
  call WarningMessage(2,'I do not know what reaction field type to use.')
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the original one-electron hamiltonian and add the reaction
! field contribution to it.

Label = 'h1    XX'
call Get_Temp(Label,RFld(:,1),nh1)
RFld(:,1) = h1(:)-RFld(:,1)
! Add contribution to the TwoHam array
TwoHam(:) = TwoHam(:)+RFld(:,2)

! Store away information for perturbative calculations

RFld(:,1)= RFld(:,1)+RFld(:,2)

ERFSelf = RepNuc-RepNuc_Temp
EEE = DDot_(nh1,RFld(1,2),1,D,1)
ERFSelf = ERFSelf-Half*EEE

call Put_dScalar('RF Self Energy',ERFSelf)
call Put_dArray('Reaction field',RFld(1,1),nh1)
call mma_deallocate(RFld)
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out the matrix elements of the reaction field potential
! to be exploited in subsequent codes like, CC, MP2, and CASPT2.

Label = 'PotNucXX'
call Get_Temp(Label,RepNucXX,1)
RepNuc_RF = RepNuc-RepNucXX(1)

call mma_allocate(h1_RF,nh1+4,Label='h1_RF')
call mma_allocate(h1_XX,nh1,Label='h1_XX')

Label = 'h1    XX'
call Get_Temp(Label,h1_XX,nh1)
h1_RF(:) = h1(:)-h1_XX(:)
call mma_deallocate(h1_XX)

h1_RF(nh1+3) = RepNuc_RF

iSyLbl = 1
iRc = -1
iOpt = 0
iComp = 1
Label = 'OneHamRF'
call WrOne(iRc,iOpt,Label,iComp,h1_RF,iSyLbl)

call mma_deallocate(h1_RF)
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_iSD()

return

end subroutine DrvRF
