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

subroutine DrvXV(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq,lRF,KSDFT,ExFac,iCharge,iSpin,D1I,D1A,nD1,DFTFOCK,Do_DFT)

use OFembed, only: Do_OFemb, OFE_KSDFT
use Constants, only: Zero

implicit none
integer nh1, nD1, iCharge, iSpin
real*8 h1(nh1), TwoHam(nh1), D(nh1,2)
real*8 D1I(nD1), D1A(nD1)
real*8 RepNuc, ExFac
logical First, Dff, lRF, NonEq, Do_Grad, Do_DFT
character(len=*) KSDFT
character(len=4) DFTFOCK
real*8 Grad(1), RN(1)
integer nGrad
logical Do_ESPF
#ifdef _EFP_
logical EFP_On
#endif
character(len=8) Label

!                                                                      *
!***********************************************************************
!                                                                      *
! If first iterations save original h1 and RepNuc, the reaction
! field and/or the dft section might later modify these. Any static
! contributions to these are updated internally in the subroutine
! rctfld.

! We also save one set which are not to be modified at all.

if (First) then
  Label = 'PotNuc00'
  call Put_Temp(Label,[RepNuc],1)
  Label = 'h1_raw  '
  call Put_Temp(Label,h1,nh1)

  Label = 'PotNucXX'
  call Put_Temp(Label,[RepNuc],1)
  Label = 'h1    XX'
  call Put_Temp(Label,h1,nh1)
end if

! Start always with the original h1 and RepNuc! Observe that they
! now might contain static updates from the RF or the Langevin
! codes through the action of the subroutine rctfld.

Label = 'PotNuc00'
call Get_Temp(Label,RN,1)
RepNuc = RN(1)
Label = 'h1_raw  '
call Get_Temp(Label,h1,nh1)
!nf
!                                                                      *
!***********************************************************************
!                                                                      *
!                        ESPF section                                  *
!                                                                      *
!***********************************************************************
!                                                                      *
call DecideOnESPF(Do_ESPF)
if (Do_ESPF) call h1_espf(h1,RepNuc,nh1,First,Do_DFT)
!nf
!                                                                      *
!***********************************************************************
!                                                                      *
!              Reaction Field section                                  *
!                                                                      *
!***********************************************************************
!                                                                      *
if (lRF) call DrvRF(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq,iCharge)
!                                                                      *
!***********************************************************************
!                                                                      *
!                         DFT section                                  *
!                                                                      *
!***********************************************************************
!                                                                      *
Do_Grad = .false.
Grad = Zero
nGrad = 1
if ((KSDFT /= 'SCF') .and. Do_DFT) call DrvDFT(h1,nh1,KSDFT,ExFac,Do_Grad,Grad,nGrad,iSpin,DFTFOCK)
!                                                                      *
!***********************************************************************
!                                                                      *
!                         Orbital-Free Embedding section               *
!                                                                      *
!***********************************************************************
!                                                                      *
if (Do_OFemb) call DrvEMB(nh1,OFE_KSDFT,Do_Grad,Grad,nGrad,DFTFOCK)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _EFP_
if (EFP_On()) call DrvEFP(First)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(D1I)
  call Unused_real_array(D1A)
end if

end subroutine DrvXV
