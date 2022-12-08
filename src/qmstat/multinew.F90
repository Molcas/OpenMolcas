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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************
!  MultiNew
!
!> @brief
!>   Perform the MME in contracted AO-basis
!> @author A. Ohrn
!>
!> @details
!> (i) Read in the multipole integrals from Seward. (ii) Construct
!> some data to simplify accessing the computed data. (iii) Make
!> the actual MME.
!>
!> @note
!> Requires numbers taken from ::qfread. We also need some integrals
!> that supposedly have been computed by Seward.
!>
!> @param[in]  nAt      Number of atoms in QM-molecule
!> @param[in]  nBas     Number of contracted basis functions
!> @param[in]  nOcc     Number of basis functions of the \f$ i \f$ -th atom-type
!> @param[in]  natyp    Number of atoms of the \f$ i \f$ -th atom-type
!> @param[in]  nntyp    Number of atom-types
!> @param[out] MME      The multicenter multipole expanded densities of unique pairs of contracted basis functions
!> @param[out] iCenTri  Set of indices that tells to which center the \f$ i \f$ -th unique pair of basis functions in a lower
!>                      triangularly stored matrix belongs
!> @param[out] iCenTriT Just like \p iCenTri, but in square shape
!> @param[out] nMlt     Highest multipole in MME
!> @param[out] outxyz   Expansion centers in molecule
!> @param[in]  lSlater
!***********************************************************************

subroutine MultiNew(nAt,nBas,nOcc,natyp,nntyp,MME,iCenTri,iCenTriT,nMlt,outxyz,lSlater)

use qmstat_global, only: MxMltp
use Index_Functions, only: nTri3_Elem, nTri_Elem
use Data_Structures, only: Alloc1DArray_Type, Allocate_DT, Deallocate_DT
use OneDat, only: sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAt, nBas, nntyp, nOcc(nntyp), natyp(nntyp)
type(Alloc1DArray_Type), intent(out) :: MME(nTri3_Elem(MxMltp))
integer(kind=iwp), intent(out) :: iCenTri(nTri_Elem(nBas)), iCenTriT(nBas,nBas), nMlt
real(kind=wp), intent(out) :: outxyz(3,nTri_Elem(nAt))
logical(kind=iwp), intent(in) :: lSlater
integer(kind=iwp) :: i, iAt, iB1, iB2, iCmp, iComp, iDum(1), iMlt, Ind, Indie, IndiePrev, iOpt, irc, iSmLbl, j, kaunt, kaunter, &
                     LMltSlq, Lu_One, nB1Prev, nB2Prev, nBasA, nComp, nMul, nSize
real(kind=wp) :: CordMul(MxMltp,3), Corr, CorrDip1, CorrDip2, CorrOvl
logical(kind=iwp) :: Changed1, Changed2, Lika
character(len=20) :: MemLab
character(len=8) :: Label
character(len=2) :: ChCo, ChCo2
integer(kind=iwp), allocatable :: nBasAt(:)
real(kind=wp), allocatable :: xyz(:,:,:)
type(Alloc1DArray_Type), allocatable :: Mult(:,:)
integer(kind=iwp), parameter :: iX(6) = [1,1,1,2,2,3], iY(6) = [1,2,3,2,3,3]
character(len=*), parameter :: Integrals(3) = ['MLTPL  0','MLTPL  1','MLTPL  2']
integer(kind=iwp), external :: IsFreeUnit
#include "warnings.h"
!Jose.No Nuclear charges in Slater

!----------------------------------------------------------------------*
! Read the multipole integrals in contracted AO-basis.                 *
!----------------------------------------------------------------------*
irc = -1
Lu_One = IsFreeUnit(49)
iOpt = 0
call OpnOne(irc,iOpt,'ONEINT',Lu_One)
if (irc /= 0) then
  write(u6,*)
  write(u6,*) 'ERROR! Could not open one-electron integral file.'
  call Quit(_RC_IO_ERROR_READ_)
end if

call Allocate_DT(Mult,[1,MxMltp],[1,nTri_Elem(MxMltp)],label='Mult')

! This loop will terminate when no more multipole integrals are
! available, hence there is not a problem that we apparently loop
! over MxMltpl, which is a fixed number.

outer: do iMlt=1,MxMltp
  nComp = nTri_Elem(iMlt)
  do iComp=1,nComp
    iCmp = iComp
    irc = -1
    iOpt = ibset(0,sOpSiz)
    iSmLbl = 1
    Label = integrals(iMlt)
    call iRdOne(irc,iOpt,Label,iCmp,iDum,iSmLbl)
    if (irc == 0) nSize = iDum(1)
    if (irc /= 0) then
      if (iComp /= 1) then
        write(u6,*)
        write(u6,*) 'ERROR! Failed to read number of one-electron integrals.'
        call Quit(_RC_IO_ERROR_READ_)
      else  !Normal exit here.
        nMlt = iMlt-1
        exit outer
      end if
    end if
    if (nSize /= 0) then
      write(ChCo,'(I2.2)') iMlt
      write(ChCo2,'(I2.2)') iComp
      write(MemLab,*) 'MEM'//ChCo//ChCo2
      call mma_allocate(Mult(iMlt,iComp)%A,nSize+4,label=MemLab)
      irc = -1
      iOpt = 0
      iSmLbl = 0
      call RdOne(irc,iOpt,Label,iCmp,Mult(iMlt,iComp)%A,iSmLbl) !Collect integrals
    else
      write(u6,*)
      write(u6,*) 'ERROR! Problem reading ',integrals(iMlt)
      call Quit(_RC_IO_ERROR_READ_)
    end if
  end do
  CordMul(iMlt,:) = Mult(iMlt,1)%A(nSize+1:3)
  nMlt = MxMltp
end do outer

!----------------------------------------------------------------------*
! Collect centers from preceeding MpProp calculation. Compute two      *
! index vectors. First one gives index of atom on which the ith basis  *
! function is centered. The other (iCenTri) gives to which center the  *
! ith unique basis function product belong.                            *
!----------------------------------------------------------------------*
call mma_allocate(xyz,3,nAt,nAt,label='xyz')
call Get_Centers(nAt,xyz)
do i=1,nAt
  outxyz(:,i) = xyz(:,i,i)
end do
kaunt = nAt
do i=1,nAt
  do j=1,i-1
    kaunt = kaunt+1
    outxyz(:,kaunt) = xyz(:,i,j)
  end do
end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
!Jose. Collect data of the Slater representation of the Quantum System *
! Prefactors, Exponents, PointNuclearCharges.                          *
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
if (lSlater) then
  call Get_Slater(LMltSlQ,outxyz,nAt)

  if (LMltSlQ+1 /= nMlt) then
    write(u6,*) 'ERROR! Multipole order',LMltSlQ,' in DiffPr file is different from order',nMlt-1, &
                ' in One-electron file. Check your files.'
    call Quit(_RC_GENERAL_ERROR_)
  end if
end if

call mma_allocate(nBasAt,nBas,label='nBasAt')

kaunter = 0
iAt = 0
do i=1,nntyp
  nBasA = nOcc(i)/natyp(i)
  do j=1,natyp(i)
    iAt = iAt+1
    nBasAt(kaunter+1:kaunter+nBasA) = iAt
    kaunter = kaunter+nBasA
  end do
end do

kaunter = 0
Indie = nAt
IndiePrev = 1
nB1Prev = 1
nB2Prev = 1
do iB1=1,nBas !Count over unique pairs of bas.func.
  do iB2=1,iB1
    kaunter = kaunter+1
    Lika = nBasAt(iB1) == nBasAt(iB2)
    if (Lika) then !If equal indices, then take that number.
      iCenTri(kaunter) = nBasAt(iB1)
      nB1Prev = nBasAt(iB1)
      nB2Prev = nBasAt(iB2)
    else
      Changed1 = nB1Prev /= nBasAt(iB1) !Check if changed atom.
      Changed2 = nB2Prev /= nBasAt(iB2)
      if (Changed1 .and. (.not. Changed2)) then !Case when from center 1 to nAt+1.
        Indie = Indie+1
        nB1Prev = nBasAt(iB1)
        IndiePrev = Indie
      else if (Changed2 .and. (.not. Changed1)) then
        ! Moving to new atom horizontally in lower triangular
        ! matrix. If it is a jump back to left corner,
        ! do not increase index, but get old one.
        if (iB2 == 1) then
          Indie = IndiePrev
        else
          Indie = Indie+1
        end if
        nB2Prev = nBasAt(iB2)
      else if (Changed1 .and. Changed2) then !Changing both atoms.
        Indie = Indie+1
        nB1Prev = nBasAt(iB1)
        nB2Prev = nBasAt(iB2)
        IndiePrev = Indie
      end if
      iCenTri(kaunter) = Indie
    end if
  end do
end do
Ind = 0
do i=1,nBas !Let's be square.
  do j=1,i
    Ind = Ind+1
    iCenTriT(i,j) = iCenTri(ind)
    iCenTriT(j,i) = iCenTri(ind)
  end do
end do

!----------------------------------------------------------------------*
! Start the MME. To get a MME-dipole, we want <psi_i|x-x_o|psi_j> but  *
! we have <psi_i|x-x_M|psi_j> where x_o is the chosen MME-center, while*
! x_M is the center that Molcas uses. We transform in this manner:     *
! <psi_i|x-x_o|psi_j> = <psi_i|x-x_M|psi_j>+(x_M-x_o)*<psi_i|psi_j>.   *
! The quadrupole contains a further complication: not only must we     *
! include more terms, but the dipole correction may not be the MME-    *
! dipole due to that Molcas may not have used the same center for      *
! dipoles and quadrupoles. Let have a look:                            *
! <psi_i|(x-x_o)(y-y_o)|psi_j>=                                        *
! <psi_i|(x-x_M)(y-y_M)|psi_j>+(x_M-x_o)*<psi_i|y-y_M|psi_j>+...       *
! But the last dipole term may need to be transformed further if y_M   *
! for the quadrupoles are not the same as y_M for the dipoles. This is *
! the explanation for the somewhat "sliskiga" expression below for the *
! quadrupoles.                                                         *
!----------------------------------------------------------------------*
if (nMlt > 3) then !This number is connected to for how high order of multipole we have implemented below.
  write(u6,*)
  write(u6,*) 'Too high order of multipole in MME.'
  call Quit(_RC_INTERNAL_ERROR_)
end if
nMul = 0
do i=1,nMlt
  nMul = nMul+nTri_Elem(i)
end do
do iMlt=1,nMul
  write(ChCo,'(I2.2)') iMlt
  call mma_allocate(MME(iMlt)%A,nSize,label='MME'//ChCo)
end do

! The MME.

kaunt = 0
do iB1=1,nBas
  do iB2=1,iB1
    kaunt = kaunt+1

    ! The charge. No translation.

    MME(1)%A(kaunt) = Mult(1,1)%A(kaunt)

    ! The dipole. Translation gives rise to charge.

    do i=1,3
      Corr = (CordMul(2,i)-xyz(i,nBasAt(iB1),nBasAt(iB2)))*Mult(1,1)%A(kaunt)
      MME(i+1)%A(kaunt) = Mult(2,i)%A(kaunt)+Corr
    end do

    ! The quadrupole. Translation gives rise to dipoles and charges.
    ! Also we have to keep track of the centers for the integrals computed
    ! by Seward.

    do i=1,6
      CorrDip1 = (CordMul(3,iX(i))-xyz(iX(i),nBasAt(iB1),nBasAt(iB2)))* &
                 (Mult(2,iY(i))%A(kaunt)+(CordMul(2,iY(i))-CordMul(3,iY(i)))*Mult(1,1)%A(kaunt))
      CorrDip2 = (CordMul(3,iY(i))-xyz(iY(i),nBasAt(iB1),nBasAt(iB2)))* &
                 (Mult(2,iX(i))%A(kaunt)+(CordMul(2,iX(i))-CordMul(3,iX(i)))*Mult(1,1)%A(kaunt))
      CorrOvl = (CordMul(3,iX(i))-xyz(iX(i),nBasAt(iB1),nBasAt(iB2)))* &
                (CordMul(3,iY(i))-xyz(iY(i),nBasAt(iB1),nBasAt(iB2)))*Mult(1,1)%A(kaunt)
      MME(i+4)%A(kaunt) = Mult(3,i)%A(kaunt)+CorrDip1+CorrDip2+CorrOvl
    end do
  end do
end do

! Deallocations.

call mma_deallocate(nBasAt)
call mma_deallocate(xyz)
call Deallocate_DT(Mult)

call ClsOne(irc,Lu_One)

return

end subroutine MultiNew
