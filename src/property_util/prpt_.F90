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
! Copyright (C) 1991, Roland Lindh                                     *
!               2018, Sijia S. Dong                                    *
!***********************************************************************

subroutine Prpt_(nIrrep,nBas,nDim,Occ,n2Tot,Vec,var,Short,iUHF,ifallorb)
!***********************************************************************
!                                                                      *
! Purpose: calculation of expectation values of different              *
!          operators as available on the 'ONEINT' file                 *
!                                                                      *
! Caution: before calling this subroutine one needs to                 *
!          open the ONEINT file                                        *
!                                                                      *
! Calling parameters:                                                  *
!   nIrRep            number of irreducible representations            *
!   nBas(0:nIrRep)    number of basis functions in each repre-         *
!                     sentation                                        *
!   ndim              total number of basis functions                  *
!   occ(ndim,*)       occupation number for all eigenvectors           *
!   ifallorb          logical option for whether the property of       *
!                     all orbitals are printed (and not weighted by    *
!                     occupation number)in property calculation        *
!                     (S.S.Dong, 2018)                                 *
!                                                                      *
! 1991 R. Lindh, Dept. of Theor. Chem. Univ. of Lund, Sweden.          *
! Modified by S.S.Dong, 2018, Univ. of Minnesota                       *
! - Enable properties to be printed for all orbitals                   *
! (including virtuals) and not weighted by occupation numbers          *
!***********************************************************************

use OneDat, only: sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nIrrep, nBas(0:nIrrep-1), nDim, n2Tot, iUHF
real(kind=wp), intent(in) :: Occ(nDim,*), Vec(n2Tot,*)
logical(kind=iwp), intent(in) :: var, Short, ifallorb
integer(kind=iwp) :: i, iCmp, iComp, idum(1), iEF, iOcc, iopt, iPL, iS, irc, iSmLbl, iTol, j, jRC, lpole, maxCen, maxGG, &
                     mBas(0:nIrrep-1), mDim, mInt, nbast, nblock, nCen, nComp, nfblock, nS, tNUC
real(kind=wp) :: C1(3), C2(3)
logical(kind=iwp) :: NxtOpr
character(len=8) :: label
real(kind=wp), allocatable :: Den(:), El(:,:,:), ElSum(:), Nuc(:), NucSum(:), Opr(:)
real(kind=wp), parameter :: Thrs = 1.0e-6_wp
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt
#include "hfc_logical.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0

if (iPL >= 2) then
  write(u6,*)
  call CollapseOutput(1,'   Molecular properties:')
  write(u6,'(3X,A)') '   ---------------------'
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
nS = 1
if (Short) then
  mDim = 1
else
  mDim = nDim
  if (iUHF == 1) nS = 2
end if
mBas(:) = nS*nBas(:)
!                                                                      *
!***********************************************************************
!                                                                      *
nblock = 0
do i=0,nirrep-1
  nblock = nblock+nbas(i)*(nbas(i)+1)/2
end do
nfblock = ndim*(ndim+1)/2

!                                                                      *
!***********************************************************************
!                                                                      *
if (Short) then
  ! calculate the density matrix with all off-diagonal elements
  ! multipled by 2

  !if (iUHF ==  0) then
  call mma_allocate(Den,nBlock,label='Den')
  if (var) then
    call Get_D1ao_Var(Den,nBlock)
  else
    call Get_dArray_chk('D1ao',Den,nBlock)
  end if
  !end if

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Scan the ONEINT file for multipole moment operators

call mma_allocate(Opr,nfblock+4,label='Opr')

!write(u6,*) ' Starting scan of ONEINT for multipole moments'
do i=0,99
  NxtOpr = .false.
  nComp = (i+1)*(i+2)/2
  if (allocated(Nuc)) call mma_deallocate(Nuc)
  call mma_allocate(Nuc,nComp,label='Nuc')
  if (allocated(El)) call mma_deallocate(El)
  call mma_allocate(El,mDim,nS,nComp,label='El')

  Nuc(:) = Zero
  El(:,:,:) = Zero
  write(label,'(a,i2)') 'MLTPL ',i
  do iComp=1,nComp
    iCmp = iComp
    irc = -1
    iopt = ibset(0,sOpSiz)
    call iRdOne(irc,iopt,label,iCmp,idum,iSmLbl)
    if (irc == 0) mInt = idum(1)
    if (irc /= 0) cycle
    NxtOpr = .true.
    irc = -1
    iopt = 0
    call RdOne(irc,iopt,label,iCmp,Opr,iSmLbl)
    if (irc /= 0) cycle
    if (mInt /= 0) call CmpInt(Opr,mInt,nBas,nIrrep,iSmLbl)
    Nuc(iComp) = Opr(mInt+4)
    if (iComp == 1) C1(:) = Opr(mInt+1:mInt+3)
    if (mInt == 0) cycle
    if (Short) then
      call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Den,nDim,Occ(:,1),nblock,Opr,El(:,1,iComp))
    else
      call Xprop(short,ifallorb,nIrrep,nBas,n2Tot,Vec(:,1),nDim,Occ(:,1),nblock,Opr,El(:,1,iComp))
      if (iUHF == 1) call Xprop(short,ifallorb,nIrrep,nBas,n2Tot,Vec(:,2),nDim,Occ(:,2),nblock,Opr,El(:,2,iComp))
    end if
  end do
  if (.not. NxtOpr) exit

  call prop(short,label,C1,C1,nirrep,mBas,nS*mDim,Occ,Thrs,El,Nuc,i,ifallorb)
end do
if (allocated(Nuc)) call mma_deallocate(Nuc)
if (allocated(El)) call mma_deallocate(El)

! Scan 'ONEINT' for magnetic hyperfine integrals

!                                                                      *
!***********************************************************************
!                                                                      *
! This check if the magnetic hypferfine integrals are available
! and only calculates hyperfine tensor matrix for UHF wavefunction
! in C1 symmetry. Since for UHF the spin density matrix is easily
! obtained.

MAG_X2C = .false.
irc = -1
iopt = ibset(0,sOpSiz)
label = 'MAGXP  1'
iComp = 1
call iRdOne(irc,iopt,label,iComp,idum,iSmLbl)
if (irc == 0) then
  MAG_X2C = .true.
  !if (Method == 'UHF-SCF ') then
  if (iUHF == 1) then
    if (UHF_HFC) then
      if (nIrrep == 1) then
        tNUC = 0
        nbast = nbas(0)
        call Get_iScalar('LP_nCenter',tNUC)
        call cmp_hfc(nbast,tNUC)
      else
        write(u6,'(/,/,6X,A)') 'Skipping Hyperfine tensor matrix for UHF with symmetry'
      end if
    else
      write(u6,'(/,/,6X,A)') 'Skipping Hyperfine tensor matrix'
    end if
  end if
end if

! Scan 'ONEINT' for electric field integrals

!                                                                      *
!***********************************************************************
!                                                                      *

!write(u6,*) ' Starting scan of ONEINT for various elec. field integrals'

do iEF=0,2
  nComp = (iEF+1)*(iEF+2)/2
  if (allocated(Nuc)) call mma_deallocate(Nuc)
  call mma_allocate(Nuc,nComp,label='Nuc')
  if (allocated(El)) call mma_deallocate(El)
  call mma_allocate(El,mDim,nS,nComp,label='El')

  ! create vectors to store the sums of electronic and nuclear components over all centers
  call mma_allocate(ElSum,nComp,label='ElSum')
  call mma_allocate(NucSum,nComp,label='NucSum')
  ElSum(:) = Zero
  NucSum(:) = Zero

  ! loop over different operator origins (max. 99999)

  maxCen = 99999
  nCen = 0
  do i=1,maxCen
    Nuc(:) = Zero
    El(:,:,:) = Zero
    write(label,'(a,i1,i5)') 'EF',iEF,i
    NxtOpr = .false.
    do iComp=1,nComp
      iCmp = iComp
      irc = -1
      iopt = ibset(0,sOpSiz)
      call iRdOne(irc,iopt,label,iCmp,idum,iSmLbl)
      if (irc == 0) mInt = idum(1)
      if (irc /= 0) cycle
      NxtOpr = .true.
      irc = -1
      iopt = 0
      call RdOne(irc,iopt,label,iCmp,Opr,iSmLbl)
      if (irc /= 0) cycle
      if (mInt /= 0) call CmpInt(Opr,mInt,nBas,nIrrep,iSmLbl)
      Nuc(iComp) = Opr(mInt+4)
      if (iComp == 1) C1(:) = Opr(mInt+1:mInt+3)
      if (mInt == 0) cycle
      if (Short) then
        call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Den,nDim,Occ(:,1),nblock,Opr,El(:,1,iComp))
      else
        call Xprop(short,ifallorb,nIrrep,nBas,n2Tot,Vec(:,1),nDim,Occ(:,1),nblock,Opr,El(:,1,iComp))
        if (iUHF == 1) call Xprop(short,ifallorb,nIrrep,nBas,n2Tot,Vec(:,2),nDim,Occ(:,2),nblock,Opr,El(:,2,iComp))
      end if
    end do
    if (.not. NxtOpr) exit

    call Prop(short,label,C1,C1,nirrep,mBas,nS*mDim,Occ,Thrs,El,Nuc,iEF,ifallorb)
    ! add the components to the sums, and update the total number of centers
    do iComp=1,nComp
      do iS=1,nS
        do iOcc=1,mDim
          ElSum(iComp) = ElSum(iComp)+El(iOcc,iS,iComp)
        end do
      end do
    end do
    NucSum(:) = NucSum(:)+Nuc(:)
    nCen = i
  end do

  if (nCen > 0) then
    ! set the tolerance according to the total number of centers
    ! (assuming error scales with sqrt(ncen))
    if (Label(1:3) == 'EF2') then
      iTol = 4
    else
      iTol = 5
    end if
    iTol = iTol-nint(Half*log10(real(nCen,kind=wp)))
    ! set MAG_X2C to avoid tests of electric field properties when
    ! wavefunction is X2C transformed (there is no way to tell but we
    ! can kind of tell by reading MAG x2c integrals) This is a workaround.
    if (.not. MAG_X2C) then
      write(label,'(a,i1,a)') 'EF',iEF,'   el'
      call Add_Info(label,ElSum,nComp,iTol)
      write(label,'(a,i1,a)') 'EF',iEF,'  nuc'
      call Add_Info(label,NucSum,nComp,iTol)
    end if
  end if
  call mma_deallocate(ElSum)
  call mma_deallocate(NucSum)
end do
if (allocated(Nuc)) call mma_deallocate(Nuc)
if (allocated(El)) call mma_deallocate(El)
!                                                                      *
!***********************************************************************
!                                                                      *

!write(u6,*) ' Starting scan of ONEINT for various contact term integrals'

nComp = 1
call mma_allocate(Nuc,nComp,label='Nuc')
call mma_allocate(El,mDim,nS,nComp,label='El')

! create vectors to store the sums of electronic and nuclear components over all centers
call mma_allocate(ElSum,nComp,label='ElSum')
call mma_allocate(NucSum,nComp,label='NucSum')
ElSum(:) = Zero
NucSum(:) = Zero

! loop over different operator origins (max. 99999)

maxCen = 99999
nCen = 0
do i=1,maxCen
  Nuc(:) = Zero
  El(:,:,:) = Zero
  write(label,'(a,i5)') 'CNT',i
  NxtOpr = .false.

  ! dummy loop, but we keep the structure
  do iComp=1,nComp
    iCmp = iComp
    irc = -1
    iopt = ibset(0,sOpSiz)
    call iRdOne(irc,iopt,label,iCmp,idum,iSmLbl)
    if (irc == 0) mInt = idum(1)
    if (irc /= 0) cycle
    NxtOpr = .true.
    irc = -1
    iopt = 0
    call RdOne(irc,iopt,label,iCmp,Opr,iSmLbl)
    if (irc /= 0) cycle
    if (mInt /= 0) call CmpInt(Opr,mInt,nBas,nIrrep,iSmLbl)
    Nuc(iComp) = Opr(mInt+4)
    C1(:) = Opr(mInt+1:mInt+3)
    if (mInt == 0) cycle
    if (Short) then
      call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Den,nDim,Occ(:,1),nblock,Opr,El(:,1,iComp))
    else
      call Xprop(short,ifallorb,nIrrep,nBas,n2Tot,Vec(:,1),nDim,Occ(:,1),nblock,Opr,El(:,1,iComp))
      if (iUHF == 1) call Xprop(short,ifallorb,nIrrep,nBas,n2Tot,Vec(:,2),nDim,Occ(:,2),nblock,Opr,El(:,2,iComp))
    end if
  end do
  if (.not. NxtOpr) exit

  call Prop(short,label,C1,C1,nirrep,mBas,nS*mDim,Occ,Thrs,El,Nuc,0,ifallorb)
  ! add the components to the sums, and update the total number of centers
  do iComp=1,nComp
    do iS=1,nS
      do iOcc=1,mDim
        ElSum(iComp) = ElSum(iComp)+El(iOcc,iS,iComp)
      end do
    end do
  end do
  NucSum(:) = NucSum(:)+Nuc(:)
  nCen = i
end do
call mma_deallocate(Nuc)
call mma_deallocate(El)

if (nCen > 0) then
  ! set the tolerance according to the total number of centers
  ! (assuming error scales with sqrt(ncen))
  iTol = 5
  iTol = iTol-nint(Half*log10(real(nCen,kind=wp)))
  if (.not. MAG_X2C) then
    write(label,'(a,a)') 'CNT','   el'
    call Add_Info(label,ElSum,nComp,iTol)
    write(label,'(a,a)') 'CNT','  nuc'
    call Add_Info(label,NucSum,nComp,iTol)
  end if
end if
call mma_deallocate(ElSum)
call mma_deallocate(NucSum)
!                                                                      *
!***********************************************************************
!                                                                      *
!write(u6,*) ' Starting scan of ONEINT diamagnetic shielding'

nComp = 9
lpole = 2
call mma_allocate(Nuc,nComp,label='Nuc')
call mma_allocate(El,mDim,nS,nComp,label='El')

maxGG = 99
maxCen = 99
! loop over different gauge origins (max.99)
do j=1,maxGG
  Nuc(:) = Zero
  El(:,:,:) = Zero
  jRC = 0
  ! loop over different operator origins (max.99)
  do i=1,maxCen
    write(label,'(a,i2,i2)') 'DMS ',j,i
    NxtOpr = .false.
    do iComp=1,nComp
      iCmp = iComp
      irc = -1
      iopt = ibset(0,sOpSiz)
      call iRdOne(irc,iopt,label,iCmp,idum,iSmLbl)
      if (irc == 0) mInt = idum(1)
      if (irc /= 0) cycle
      NxtOpr = .true.
      irc = -1
      iopt = 0
      call RdOne(irc,iopt,label,iCmp,Opr,iSmLbl)
      if (irc /= 0) cycle
      if (mInt /= 0) call CmpInt(Opr,mInt,nBas,nIrrep,iSmLbl)
      Nuc(iComp) = Opr(mInt+4)
      if (iComp == 1) then
        C1(:) = Opr(mInt+1:mInt+3)
        C2(:) = C1(:)
      else if (iComp == 2) then
        C2(:) = Opr(mInt+1:mInt+3)
      end if
      if (mInt == 0) cycle
      if (Short) then
        call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Den,nDim,Occ(:,1),nblock,Opr,El(:,1,iComp))
      else
        call Xprop(short,ifallorb,nIrrep,nBas,n2Tot,Vec(:,1),nDim,Occ(:,1),nblock,Opr,El(:,1,iComp))
        if (iUHF == 1) call Xprop(short,ifallorb,nIrrep,nBas,n2Tot,Vec(:,2),nDim,Occ(:,2),nblock,Opr,El(:,2,iComp))
      end if
    end do
    if (.not. NxtOpr) exit

    call prop(short,label,C1,C2,nirrep,mBas,nS*mDim,Occ,Thrs,El,Nuc,lpole,ifallorb)
    jRC = 1
  end do
  if (jRC == 0) exit
end do
call mma_deallocate(Nuc)
call mma_deallocate(El)
call mma_deallocate(Opr)
if (Short) call mma_deallocate(Den)

!                                                                      *
!***********************************************************************
!                                                                      *
if (iPL >= 2) then
  call CollapseOutput(0,'   Molecular properties:')
  write(u6,*)
end if

return

end subroutine Prpt_
