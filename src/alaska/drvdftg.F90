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
!               2021, Jie J. Bao                                       *
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
!             Jie Bao, Dept. of Chem.                                  *
!             University of Minnesota, US                              *
!             October 2021                                             *
!***********************************************************************

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Para_Info, only: King
use nq_Info, only: nAshT, Grid_Type, Moving_Grid
use Disp, only: ChDisp
use NAC, only: isNAC, NACStates
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Temp(nGrad)
#include "print.fh"
integer(kind=iwp) :: iDFT, iEnd, iI, iJ, iIrrep, IK, iOpt, iPrint, iRout, iSpin, jPrint, LuWr, nAct(nIrrep), nDens, ng1, ng2, nRoots
real(kind=wp) :: Dummy(1), ExFac, TCpu1, TCpu2, TWall1, TWall2
logical(kind=iwp) :: Do_Grad, l_casdft
character(len=80) :: Label
character(len=80) :: KSDFT
character(len=8) :: Method
character(len=4) :: DFTFOCK
real(kind=wp), allocatable :: Temp2(:), R(:), G1qs(:), G2qs(:), G1qt(:), G2qt(:), D1AOMS(:), D1SAOMS(:), D1AOt(:), D1SAOt(:)

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
LuWr = u6

nDens = 0
do iIrrep=0,nIrrep-1
  nDens = nDens+nBas(iIrrep)*(nBas(iIrrep)+1)/2
end do

!     D F T - g r a d i e n t                                          *
!***********************************************************************
!8)                                                                    *
!     D F T - g r a d i e n t

!call Get_iOption(iDFT)

call Get_cArray('DFT functional',KSDFT,80)
l_casdft = (KSDFT(1:2) == 'T:') .or. (KSDFT(1:3) == 'FT:')

if (l_casdft) then
  DFTFOCK = 'ROKS'
  call Get_iScalar('System BitSwitch',iOpt)
  iOpt = ibset(iOpt,6)
  call Put_iScalar('System BitSwitch',iOpt)
end if

call Get_iScalar('System BitSwitch',iDFT)
if (btest(iDFT,6)) then

  call StatusLine(' Alaska:',' Computing DFT gradients')

  call Get_cArray('DFT functional',KSDFT,80)
  ExFac = Zero ! Set to proper value at retrun!
  Do_Grad = .true.
  call Get_iScalar('Multiplicity',iSpin)
  !write(LuWr,*) 'DrvDFTg: KSDFT=',KSDFT
  !write(LuWr,*) 'DrvDFTg: ExFac=',ExFac
  call Get_cArray('Relax Method',Method,8)
  if (Method /= 'MSPDFT') then
    call DrvDFT(Dummy,nDens,KSDFT,ExFac,Do_Grad,Temp,nGrad,iSpin,DFTFOCK)
  else
    ! modifications for MS-PDFT gradient starting here
    call Get_iScalar('Number of roots',nRoots)
    call Get_iArray('nAsh',nAct,nIrrep)
    nasht = 0
    do iIrrep=1,nIrrep
      NASHT = NASHT+nAct(IIrrep)
    end do
    NG1 = NASHT*(NASHT+1)/2
    NG2 = NG1*(NG1+1)/2
    call mma_allocate(R,nRoots**2)
    call mma_allocate(Temp2,nGrad)
    call mma_allocate(G1qt,nG1)
    call mma_allocate(G2qt,nG2)
    call mma_allocate(G1qs,nG1*nRoots)
    call mma_allocate(G2qs,nG2*nRoots)
    call mma_allocate(D1AOMS,nDens*nRoots)
    call mma_allocate(D1AOt,nDens)
    if (iSpin /= 1) then
      call mma_allocate(D1SAOMS,nDens*nRoots)
      call mma_allocate(D1SAOt,nDens)
    end if
    call Get_DArray('MS_FINAL_ROT',R,nRoots**2)
    Temp(:) = Zero
    call Get_dArray_chk('D1mo',G1qt,nG1)
    call Get_dArray_chk('P2mo',G2qt,nG2)
    call Get_DArray('D1INTER',G1qs,ng1*nRoots)
    call Get_DArray('P2INTER',G2qs,ng2*nRoots)
    call Get_DArray('D1AO_MS',D1AOMS,nDens*nRoots)
    call Get_dArray_chk('D1ao',D1AOt,nDens)
    if (iSpin /= 1) then
      call Get_DArray('D1SAO_MS',D1SAOMS,nDens*nRoots)
      call Get_dArray_chk('D1sao',D1SAOt,nDens)
    end if
    do IK=1,nRoots
      call Put_dArray('D1mo',G1qs((IK-1)*nG1+1),nG1)
      call Put_dArray('P2mo',G2qs((IK-1)*nG2+1),nG2)
      call Put_dArray('D1ao',D1AOMS((IK-1)*nDens+1),nDens)
      if (iSpin /= 1) then
        call Put_dArray('D1sao',D1SAOMS((IK-1)*nDens+1),nDens)
      end if
      Temp2(:) = Zero
      call DrvDFT(Dummy,nDens,KSDFT,ExFac,Do_Grad,Temp2,nGrad,iSpin,DFTFOCK)
      jPrint = nPrint(112)
      if (isNAC) then
        iI = NACstates(1)
        iJ = NACstates(2)
        call DAXPY_(nGrad,R((II-1)*nRoots+IK)*R((iJ-1)*nRoots+IK),Temp2,1,Temp,1)
        if (jPrint >= 15) then
          Label = 'DFT Int Contribution'
          write(u6,*) 'state, coeff i, coeff j',IK,R((II-1)*nRoots+IK),R((iJ-1)*nRoots+IK)
          call PrGrad(Label,Temp2,nGrad,ChDisp)
        end if
      else
        call Get_iScalar('Relax CASSCF root',iI)
        call DAXPY_(nGrad,R((II-1)*nRoots+IK)**2,Temp2,1,Temp,1)
        if (jPrint >= 15) then
          Label = 'DFT Int Contribution'
          write(u6,*) 'state, coeff',IK,R((II-1)*nRoots+IK)
          call PrGrad(Label,Temp2,nGrad,ChDisp)
        end if
      end if
    end do
    call Put_dArray('D1mo',G1qt,nG1)
    call Put_dArray('P2mo',G2qt,nG2)
    call Put_dArray('D1ao',D1AOt,nDens)
    if (ISpin /= 1) call Put_dArray('D1sao',D1SAOt,nDens)
    call mma_deallocate(R)
    call mma_deallocate(Temp2)
    call mma_deallocate(G1qt)
    call mma_deallocate(G2qt)
    call mma_deallocate(G1qs)
    call mma_deallocate(G2qs)
    call mma_deallocate(D1AOMS)
    call mma_deallocate(D1AOt)
    if (iSpin /= 1) then
      call mma_deallocate(D1SAOMS)
      call mma_deallocate(D1SAOt)
    end if
    ! End of condition for MS-PDFT gradient
  end if
  iEnd = 1
  do
    if (KSDFT(iEnd:iEnd) == ' ') then
      iEnd = iEnd-1
      exit
    else
      iEnd = iEnd+1
    end if
  end do
  Label = 'The DFT('//KSDFT(1:iEnd)//') contribution'
  jPrint = nPrint(112)
  !AMS
  !jPrint = 15
  if (jPrint >= 15) call PrGrad(Label,Temp,nGrad,ChDisp)
  if (king()) call DaXpY_(nGrad,One,Temp,1,Grad,1)
  if (iPrint >= 6) then
    write(LuWr,*)
    if (Grid_Type == Moving_Grid) then
      write(LuWr,*) 'DFT contribution computed for a moving grid.'
    else
      write(LuWr,*) 'DFT contribution computed for a fixed grid.'
    end if
    write(LuWr,*)
  end if

end if
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu2,TWall2)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine DrvDFTg
