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
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Prpt_old(nirrep,nbas,ndim,n2dim,vec,occ)
!***********************************************************************
!                                                                      *
! Purpose: calculation of expectation values of different              *
!          operators as available on the 'ONEINT' file                 *
!                                                                      *
! Caution: before calling this subroutine one needs to                 *
!          open the ONEINT file                                        *
!                                                                      *
! Calling parameters:                                                  *
!                                                                      *
!   nIrRep            number of irreducible representations            *
!   nBas(0:nIrRep)    number of basis functions in each repre-         *
!                     sentation                                        *
!   ndim              total number of basis functions                  *
!   n2dim             sum(i,i=0,nirrep-1)(nbas(i)**2): elements        *
!                     of all vectors in all representations            *
!   vec(n2dim)        eigenvectors, all for each representation        *
!   occ(ndim)         occupation number for all eigenvectors           *
!                                                                      *
! 1991 R. Lindh, Dept. of Theor. Chem. Univ. of Lund, Sweden.          *
!***********************************************************************

use OneDat, only: sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nirrep, nbas(0:nirrep), ndim, n2dim
real(kind=wp), intent(in) :: vec(n2dim), occ(ndim)
integer(kind=iwp) :: i, iCmp, iComp, iCount, idum(1), iEF, ii, il, iOcc, iopt, ir, irc, iSmLbl, iVec, j, jCount, jRC, maxCen, &
                     maxGG, mDim, mInt, nblock, nComp, nfblock
real(kind=wp) :: C1(3), C2(3), dummy
logical(kind=iwp) :: short, NxtOpr, ifallorb
character(len=8) :: label
real(kind=wp), allocatable :: Den(:), El(:), Nuc(:), Opr(:)

write(u6,*)
call CollapseOutput(1,'   Molecular properties:')
write(u6,'(3X,A)') '   ---------------------'
write(u6,*)
short = .true.
ifallorb = .false.
mDim = 1

nblock = 0
do i=0,nirrep-1
  nblock = nblock+nbas(i)*(nbas(i)+1)/2
end do
nfblock = ndim*(ndim+1)/2

! calculate the density matrix with all off-diagonal elements
! multipled by 2

call mma_allocate(Den,nblock,label='Den')
Den(:) = Zero
iCount = 1
iVec = 0
iOcc = 0
do i=0,nIrrep-1
  do ii=1,nBas(i)
    iOcc = iOcc+1
    jCount = iCount
    do il=1,nBas(i)
      do ir=1,il-1
        Den(jCount) = Den(jCount)+Two*Vec(iVec+il)*Vec(iVec+ir)*Occ(iOcc)
        jCount = jCount+1
      end do
      Den(jCount) = Den(jCount)+Vec(iVec+il)*Vec(iVec+il)*Occ(iOcc)
      jCount = jCount+1
    end do
    iVec = iVec+nBas(i)
  end do
  iCount = iCount+nbas(i)*(nBas(i)+1)/2
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Scan the ONEINT file for multipole moment operators

call mma_allocate(Opr,nfblock+4,label='Opr')

!write(u6,*) ' Starting scan of ONEINT for multipole moments'
do i=1,99
  NxtOpr = .false.
  nComp = (i+1)*(i+2)/2
  if (allocated(Nuc)) call mma_deallocate(Nuc)
  call mma_allocate(Nuc,nComp,label='Nuc')
  if (allocated(El)) call mma_deallocate(El)
  call mma_allocate(El,nComp,label='El')

  Nuc(:) = Zero
  El(:) = Zero
  write(label,'(a,i2)') 'MLTPL ',i
  do iComp=1,nComp
    iCmp = iComp
    irc = -1
    iopt = ibset(0,sOpSiz)
    call iRdOne(irc,iopt,label,iCmp,idum,iSmLbl)
    if (irc /= 0) cycle
    mInt = idum(1)
    NxtOpr = .true.
    irc = -1
    iopt = 0
    call RdOne(irc,iopt,label,iCmp,Opr,iSmLbl)
    if (irc /= 0) cycle
    if (mInt /= 0) call CmpInt(Opr,mInt,nBas,nIrrep,iSmLbl)
    Nuc(iComp) = Opr(mInt+4)
    if (iComp == 1) C1(:) = Opr(mInt+1:mInt+3)
    if (mInt == 0) cycle
    call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Den,nDim,Occ,nblock,Opr,El(iComp))
  end do
  if (.not. NxtOpr) exit

  call prop(short,label,C1,C1,nirrep,nBas,mDim,occ,dummy,El,Nuc,i,ifallorb)
end do
if (allocated(Nuc)) call mma_deallocate(Nuc)
if (allocated(El)) call mma_deallocate(El)
!                                                                      *
!***********************************************************************
!                                                                      *
! Scan 'ONEINT' for electric field integrals

!write(u6,*) ' Starting scan of ONEINT for various elec. field integrals'

do iEF=0,2
  nComp = (iEF+1)*(iEF+2)/2
  if (allocated(Nuc)) call mma_deallocate(Nuc)
  call mma_allocate(Nuc,nComp,label='Nuc')
  if (allocated(El)) call mma_deallocate(El)
  call mma_allocate(El,nComp,label='El')

  ! loop over different operator origins (max. 9999)

  maxCen = 9999
  do i=1,maxCen
    Nuc(:) = Zero
    El(:) = Zero
    write(label,'(a,i1,i5)') 'EF',iEF,i
    NxtOpr = .false.
    do iComp=1,nComp
      iCmp = iComp
      irc = -1
      iopt = ibset(0,sOpSiz)
      call iRdOne(irc,iopt,label,iCmp,idum,iSmLbl)
      if (irc /= 0) cycle
      mInt = idum(1)
      NxtOpr = .true.
      irc = -1
      iopt = 0
      call RdOne(irc,iopt,label,iCmp,Opr,iSmLbl)
      if (irc /= 0) cycle
      if (mInt /= 0) call CmpInt(Opr,mInt,nBas,nIrrep,iSmLbl)
      Nuc(iComp) = Opr(mInt+4)
      if (iComp == 1) C1(:) = Opr(mInt+1:mInt+3)
      if (mInt == 0) cycle
      call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Den,nDim,Occ,nblock,Opr,El(iComp))
    end do
    if (.not. NxtOpr) exit

    call prop(short,label,C1,C1,nirrep,nBas,mDim,occ,dummy,El,Nuc,i,ifallorb)
  end do

end do
if (allocated(Nuc)) call mma_deallocate(Nuc)
if (allocated(El)) call mma_deallocate(El)
!                                                                      *
!***********************************************************************
!                                                                      *
!write(u6,*) ' Starting scan of ONEINT diamagnetic shielding'

nComp = 9
call mma_allocate(Nuc,nComp,label='Nuc')
call mma_allocate(El,nComp,label='El')

maxGG = 99
maxCen = 99
! loop over differnt gauge origins (max.99)
do j=1,maxGG
  Nuc(:) = Zero
  El(:) = Zero
  jRC = 0
  ! loop over differnt operator origins (max.99)
  do i=1,maxCen
    write(label,'(a,i2,i2)') 'DMS ',j,i
    NxtOpr = .false.
    do iComp=1,nComp
      iCmp = iComp
      irc = -1
      iopt = ibset(0,sOpSiz)
      call iRdOne(irc,iopt,label,iCmp,idum,iSmLbl)
      if (irc /= 0) cycle
      mInt = idum(1)
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
      call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Den,nDim,Occ,nblock,Opr,El(iComp))
    end do
    if (.not. NxtOpr) exit

    call prop(short,label,C1,C2,nirrep,nBas,mDim,occ,dummy,El,Nuc,i,ifallorb)
    jRC = 1
  end do
  if (jRC == 0) exit
end do
call mma_deallocate(Nuc)
call mma_deallocate(El)
call mma_deallocate(Opr)
call mma_deallocate(Den)
!                                                                      *
!***********************************************************************
!                                                                      *

call CollapseOutput(0,'   Molecular properties:')
write(u6,*)

end subroutine Prpt_old
