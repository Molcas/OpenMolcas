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
! Copyright (C) 1989-1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine Alaska(LuSpool,ireturn)
!***********************************************************************
!                                                                      *
!  Object: Driver for the one and two electron integral gradient       *
!          program ALASKA.                                             *
!                                                                      *
!          Alaska is a derivative code of Seward 3.1.                  *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA     *
!          July '89 - May '90                                          *
!                                                                      *
!          Roland Lindh, Dept. of Theoretical Chemistry, University of *
!          Lund, SWEDEN. Modified to gradient calculations September   *
!          1991 - February 1992.                                       *
!***********************************************************************

use Alaska_Info, only: Am
use Basis_Info, only: dbsc, nCnttp
use Gateway_global, only: Onenly, Test
use RICD_Info, only: Do_RI, Cholesky
use Para_Info, only: nProcs, King
use OFembed, only: Do_OFemb
use k2_arrays, only: DeDe
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LuSpool
integer(kind=iwp), intent(out) :: ireturn
#include "Molcas.fh"
#include "disp.fh"
#include "print.fh"
#include "rctfld.fh"
#include "columbus_gamma.fh"
#include "nac.fh"
integer(kind=iwp) :: i, iCar, iCnt, iCnttp, iPrint, irlxroot1, irlxroot2, iRout, l1, mdc, nCnttp_Valence, ndc, nDiff, nsAtom
real(kind=wp) :: TCpu1, TCpu2, TWall1, TWall2
logical(kind=iwp) :: DoRys, Found
character(len=180) :: Label
real(kind=wp), allocatable :: Grad(:), Temp(:), Tmp(:), Rlx(:,:), CSFG(:)
integer(kind=iwp), external :: isFreeUnit
real(kind=wp), external :: dnrm2_
logical(kind=iwp), external :: RF_On
!*********** columbus interface ****************************************
integer(kind=iwp) :: Columbus, colgradmode, lcartgrd, iatom, icen, j
real(kind=wp), allocatable :: Cgrad(:,:)
character(len=LenIn5), allocatable :: CNames(:)
character(len=80) :: Lab

!                                                                      *
!***********************************************************************
!                                                                      *
!call Alaska_banner()

call CWTime(TCpu1,TWall1)
!                                                                      *
!***********************************************************************
!                                                                      *
! Prologue
!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 1
!                                                                      *
!***********************************************************************
!                                                                      *
! Print program header
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the input information as Seward dumped on INFO.

nDiff = 1
DoRys = .true.
call IniSew(DoRys,nDiff)
if (RF_On()) then
  if (NonEq_Ref) then
    call WarningMessage(2,'Error in Alaska')
    write(u6,*) 'NonEq=.True., invalid option'
    call Abend()
  end if
  call Init_RctFld(.false.,iCharge_Ref)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Input specific for the gradient calculation.

call Inputg(LuSpool)

!-- Since the input has changed some of the shell information
!   regenerate the tabulated shell information.

iPrint = nPrint(iRout)

call mma_allocate(Grad,lDisp(0),Label='Grad')
call mma_allocate(Temp,lDisp(0),Label='Temp')
Grad(:) = Zero

! remove LuSpool

call Close_LuSpool(LuSpool)

! identify a Columbus calculation
! Columbus=1
! colgradmode=0   standard gradient written to GRAD
! colgradmode=1   standard gradient written to Grad State1
! colgradmode=2   standard gradient written to Grad State1
! colgradmode=3   non-adiabatic coupling vector written to NADC

call Get_iScalar('Columbus',Columbus)
call Get_iScalar('colgradmode',colgradmode)

!-- Start computing the gradients
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute nuclear contributions.

if (king() .or. HF_Force) then

  ! per default NADC must not have nuclear contributions added

  if (NO_NUC .or. ((Columbus == 1) .and. (colgradmode == 3))) then
    write(u6,*) 'Skipping Nuclear Charge Contribution'
  else
    call DrvN1(Grad,Temp,lDisp(0))
    if (iPrint >= 15) then
      Lab = ' Total Nuclear Contribution'
      call PrGrad(Lab,Grad,lDisp(0),ChDisp)
    end if
  end if
end if

!iPrint = 16
!                                                                      *
!***********************************************************************
!                                                                      *
if (Do_OFemb) then
  ! RepNuc term to the Orbital-Free Embedding gradient
  call DrvN1_EMB(Grad,Temp,lDisp(0))
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. Test) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !-- Compute contribution due to the derivative of the one-
  !   electron hamiltonian and the contribution due the re-
  !   normalization.

  if ((nProcs <= 1) .or. Onenly) then
    call Drvh1(Grad,Temp,lDisp(0))
  end if
  if (Do_OFemb) then
    ! NucAtt term to the Orbital-Free Embedding gradient
    call Drvh1_EMB(Grad,Temp,lDisp(0))
  end if
  !Lab = 'Nuc + One-electron Contribution'
  !call PrGrad(Lab,Grad,lDisp(0),ChDisp)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the DFT contribution to the gradient

  call DrvDFTg(Grad,Temp,lDisp(0))
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (.not. Onenly) then
    call mma_allocate(DeDe,[-1,-1],label='DeDe') ! Dummy allocation
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! DFT-type Orbital-Free Embedding term to the gradient

    if (Do_OFemb) call DrvEMBg(Grad,Temp,lDisp(0))
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    !-- Compute contribution due to 2-electron integrals.

    if (Cholesky .or. Do_RI) then
      if (Cholesky) then
        if (iPrint >= 6) write(u6,*) 'Cholesky-ERI gradients!'
      else
        if (iPrint >= 6) write(u6,*) 'RI-ERI gradients!'
      end if
      call Drvg1_RI(Grad,Temp,lDisp(0))
    else
      if (iPrint >= 6) write(u6,*) 'Conventional ERI gradients!'
      call Drvg1(Grad,Temp,lDisp(0))
    end if

    call DScal_(lDisp(0),Half,Temp,1)
    if (iPrint >= 15) then
      Lab = ' Two-electron Contribution'
      call PrGrad(Lab,Temp,lDisp(0),ChDisp)
    end if

    !-- Accumulate contribution to the gradient

    call GR_DArray(Grad,lDisp(0))
    call DaXpY_(lDisp(0),One,Temp,1,Grad,1)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    call mma_deallocate(DeDe)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call CWTime(TCpu2,TWall2)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !-- Apply the translational and rotational invariance of the energy.

  if (TRSymm) then
    if (iPrint >= 99) then
      call PrGrad(' Molecular gradients (no TR) ',Grad,lDisp(0),ChDisp)
      call RecPrt(' The A matrix',' ',Am,lDisp(0),lDisp(0))
    end if
    Temp(1:lDisp(0)) = Grad(1:lDisp(0))

    call dGeMV_('N',lDisp(0),lDisp(0),One,Am,lDisp(0),Temp,1,Zero,Grad,1)
    call mma_deallocate(Am)
  end if ! TRSymm
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !-- Equivalence option

  if (lEq) then
    do i=1,lDisp(0)
      if (IndxEq(i) /= i) Grad(i) = Grad(IndxEq(i))
    end do
  end if ! lEq
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
nCnttp_Valence = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%Aux) exit
  nCnttp_Valence = nCnttp_Valence+1
end do

! f^AB is the "total derivative coupling"
! h^AB is the "CI derivative coupling"
! f^AB = <B|dA/dR> = h^AB/(E_A-E_B) + f_CSF^AB = -f^BA
! h^AB = <B|dH/dR|A> = h^BA
! f_CSF^AB = -f_CSF^BA
!
! Note that we store h^AB + f_CSF^AB*(E_A-E-B), or just h^AB if
! NOCSF was given, to avoid division by (nearly) zero

if (isNAC) then
  call PrGrad('CI derivative coupling ',Grad,lDisp(0),ChDisp)
  if (DoCSF) then
    call mma_Allocate(CSFG,lDisp(0),Label='CSFG')
    call CSFGrad(CSFG,lDisp(0))
    call PrGrad('CSF derivative coupling ',CSFG,lDisp(0),ChDisp)
    call daxpy_(lDisp(0),EDiff,CSFG,1,Grad,1)
    call mma_deallocate(CSFG)
  end if
  write(u6,'(15X,A,ES13.6)') 'Energy difference: ',EDiff
  Label = ''
  Label = 'Total derivative coupling'//trim(Label)
  call mma_allocate(Tmp,lDisp(0),Label='Tmp')
  Tmp(:) = Grad(:)/EDiff
  call PrGrad(trim(Label),Tmp,lDisp(0),ChDisp)
  write(u6,'(15X,A,F12.4)') 'norm: ',dnrm2_(lDisp(0),Tmp,1)
  call mma_deallocate(Tmp)
else if (iPrint >= 4) then
  if (HF_Force) then
    call PrGrad('Hellmann-Feynman Forces ',Grad,lDisp(0),ChDisp)
  else
    call PrGrad(' Molecular gradients',Grad,lDisp(0),ChDisp)
  end if
end if
if (isNAC) then
  ! For NAC, the sign is undefined (because the wave functions can
  ! change sign), check only absolute values
  call mma_allocate(Tmp,lDisp(0),Label='Tmp')
  Tmp(:) = abs(Grad(:))
  call Add_Info('Grad',Tmp,lDisp(0),6)
  call mma_deallocate(Tmp)
else
  call Add_Info('Grad',Grad,lDisp(0),6)
end if

!-- Molcas format

!-- Write gradient to runfile.
!
! Save the gradient

call Get_iScalar('Unique atoms',nsAtom)
l1 = 3*nsAtom
call mma_allocate(Rlx,3,nsAtom,Label='Rlx')
mdc = 0
ndc = 0
do iCnttp=1,nCnttp_Valence

  ! Skip gradients for pseudo atoms

  if (dbsc(iCnttp)%pChrg .or. (dbsc(iCnttp)%nFragType > 0) .or. dbsc(iCnttp)%Frag) then
    mdc = mdc+dbsc(iCnttp)%nCntr
  else
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      ndc = ndc+1
      do iCar=1,3
        if (InxDsp(mdc,iCar) /= 0) then
          Rlx(iCar,ndc) = Grad(InxDsp(mdc,iCar))
        else

          ! Put in explicit zero if gradient is zero by symmetry.

          Rlx(iCar,ndc) = Zero
        end if
      end do
    end do
  end if
end do
!
if (HF_Force) then
  call Put_dArray('HF-forces',Rlx,l1)
else if (Columbus == 1) then
  call Put_nadc(colgradmode,Rlx,l1)
else
  call Put_dArray('GRAD',Rlx,l1)
end if
call mma_deallocate(Rlx)

!*********** columbus interface ****************************************
! print full cartesian gradient in Columbus format

if (Columbus == 1) then
  ! real*8 Cgrad(3,mxatom)
  ! character CNames(MxAtom)*9
  ! integer lcartgrd, iatom,icen,j
  call mma_allocate(CGrad,3,MxAtom,label='CGrad')
  call mma_allocate(CNames,MxAtom,label='CNames')
  call TrGrd_Alaska_(CGrad,CNames,Grad,lDisp(0),iCen)
  lcartgrd = 60
  lcartgrd = isFreeUnit(lcartgrd)
  call Molcas_Open(lcartgrd,'cartgrd')
  do IATOM=1,iCen
    write(60,1010) (CGrad(j,iatom),j=1,3)
  end do
  close(lcartgrd)
  call mma_deallocate(CGrad)
  call mma_deallocate(CNames)
end if

!-- At the end of the calculation free all memory to check for
!   corruption of the memory.

call mma_deallocate(Temp)
call mma_deallocate(Grad)

! Restore iRlxRoot if changed as set by the RASSCF module.

call qpg_iScalar('Relax CASSCF root',Found)
if (Found) then
  call Get_iScalar('Relax CASSCF root',irlxroot1)
  call qpg_iScalar('Relax Original root',Found)
  if (Found) then
    call Get_iScalar('Relax Original root',irlxroot2)
    if (iRlxRoot1 /= iRlxRoot2) then
      call Put_iScalar('Relax CASSCF root',irlxroot2)
      call Put_iScalar('NumGradRoot',irlxroot2)
    end if
  end if
end if

! Epilogue

call ClsSew()

if (iPrint >= 6) then
  call FastIO('STATUS')
end if

if (Test) then
  ireturn = 20
else
  ireturn = 0
end if

return

1010 format(3d15.6)

end subroutine Alaska
