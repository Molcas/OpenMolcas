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
! Copyright (C) 2002, Roland Lindh                                     *
!***********************************************************************

subroutine DrvDFTg(Grad,Temp,nGrad)
!***********************************************************************
!                                                                      *
! Object: driver for computation of gradient with respect to the DFT   *
!         energy.                                                      *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chem. Phys.                       *
!             University of Lund, SWEDEN                               *
!             August 2002                                              *
!***********************************************************************

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Para_Info, only: King

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "disp.fh"
#include "nq_info.fh"
character Label*80, KSDFT*16
real*8 Grad(nGrad), Temp(nGrad)
logical First, Dff, Do_Grad, l_casdft
character*4 DFTFOCK
dimension Dummy1(1), Dummy2(1), Dummy3(1), Dumm0(1), Dumm1(1)

!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu1,TWall1)
!                                                                      *
!***********************************************************************
!                                                                      *
!...  Prologue
DFTFOCK = 'SCF '
iRout = 131
iPrint = nPrint(iRout)
LuWr = 6

nDens = 0
do iIrrep=0,nIrrep-1
  nDens = nDens+nBas(iIrrep)*(nBas(iIrrep)+1)/2
end do

!     D F T - g r a d i e n t                                          *
!***********************************************************************
!8)                                                                    *
!     D F T - g r a d i e n t

!call Get_iOption(iDFT)

call Get_cArray('DFT functional',KSDFT,16)
l_casdft = (KSDFT(1:5) == 'TLSDA') .or. (KSDFT(1:6) == 'TLSDA5') .or. (KSDFT(1:5) == 'TBLYP') .or. (KSDFT(1:5) == 'TOPBE') .or. &
           (KSDFT(1:6) == 'TSSBSW') .or. (KSDFT(1:5) == 'TSSBD') .or. (KSDFT(1:5) == 'TS12G') .or. (KSDFT(1:4) == 'TPBE') .or. &
           (KSDFT(1:5) == 'FTPBE') .or. (KSDFT(1:7) == 'TREVPBE') .or. (KSDFT(1:8) == 'FTREVPBE') .or. &
           (KSDFT(1:6) == 'FTLSDA') .or. (KSDFT(1:6) == 'FTOPBE') .or. (KSDFT(1:6) == 'FTBLYP)')

if (l_casdft) then
  DFTFOCK = 'ROKS'
  call Get_iScalar('System BitSwitch',iOpt)
  iOpt = ior(iOpt,2**6)
  call Put_iScalar('System BitSwitch',iOpt)
end if

call Get_iScalar('System BitSwitch',iDFT)
if (iand(iDFT,2**6) /= 0) then

  call StatusLine(' Alaska:',' Computing DFT gradients')

  First = .true.
  Dff = .false.
  call Get_cArray('DFT functional',KSDFT,16)
  ExFac = Zero ! Set to proper value at retrun!
  Do_Grad = .true.
  call Get_iScalar('Multiplicity',iSpin)
  !write(LuWr,*) 'DrvDFTg: KSDFT=',KSDFT
  !write(LuWr,*) 'DrvDFTg: ExFac=',ExFac
  iDumm = 1
  call DrvDFT(Dummy1,Dummy2,Dummy3,Dummy4,nDens,First,Dff,lRF,KSDFT,ExFac,Do_Grad,Temp,nGrad,iSpin,Dumm0,Dumm1,iDumm,DFTFOCK)

  iEnd = 1
99 continue
  if (KSDFT(iEnd:iEnd) == ' ') then
    iEnd = iEnd-1
  else
    iEnd = iEnd+1
    Go To 99
  end if
  Label = 'The DFT('//KSDFT(1:iEnd)//') contribution'
  jPrint = nPrint(112)
  !AMS
  !jPrint = 15
  if (jPrint >= 15) call PrGrad(Label,Temp,nGrad,ChDisp,5)
  if (king()) call DaXpY_(nGrad,One,Temp,1,Grad,1)
  if (iPrint < 6) Go To 777
  write(LuWr,*)
  if (Grid_Type == Moving_Grid) then
    write(LuWr,*) 'DFT contribution computed for a moving grid.'
  else
    write(LuWr,*) 'DFT contribution computed for a fixed grid.'
  end if
  write(LuWr,*)
777 continue

end if
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu2,TWall2)
call SavTim(5,TCpu2-TCpu1,TWall2-TWall1)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine DrvDFTg
