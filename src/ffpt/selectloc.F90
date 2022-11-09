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
!  SelectLoc
!
!> @brief
!>   Localize the perturbation LoProp-style. This way a perturbation can be applied selectively on parts of a molecule
!> @author A. Ohrn
!>
!> @details
!> Collect \p H0 as it is, clean the vacuum part so only perturbation
!> is there. Collect LoProp transformation and transform. Put zeros
!> according to user specification and transform back. The localized
!> perturbation is added to the one-electron Hamiltonian and \p H0 is
!> returned.
!>
!> @param[in,out] H0    The one-electron Hamiltonain with perturbations so far.
!>                      On output the localized perturbed one-electron Hamiltonian
!> @param[in]     nSize Size of the triangular \p H0 with the additional origo and nuclear contribution
!***********************************************************************

subroutine SelectLoc(H0,nSize)

use FFPT_Global, only: LCumulate, nSets, nSym, iSelection, Bonds, nBas, Atoms
use OneDat, only: sNoOri, sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSize
real(kind=wp), intent(inout) :: H0(nSize)
character(len=8) :: Label
logical(kind=iwp) :: OneOrNot1, OneOrNot2, OneOrNot3, OneOrNot4, OneOrNot, CrazySet
integer(kind=iwp) :: iComp, idum(1), i, iOpt0, iOpt1, iOpt2, iRc, iSymLbl, j, k, kaunter, l, nInts
real(kind=wp) :: H01, H02, H03, H04, Siff
integer(kind=iwp), allocatable :: oType(:), Center(:)
real(kind=wp), allocatable :: STr(:), SSq(:,:), T(:,:), Tinv(:,:), HVac(:), V(:), VS(:,:), TEMP(:,:), VLoP(:,:)
logical(kind=iwp), parameter :: Debug = .false.

!-- Commence!

write(u6,*)
write(u6,*) ' The perturbation will be localized "LoProp style".'
write(u6,*)
write(u6,*) ' -- Number of basis subsets:',nSets
do k=1,nSets
  write(u6,*) '    ',iSelection(1,k),iSelection(2,k)
end do
write(u6,*) ' -- Atoms and bonds logical flags:'
do i=1,nSets
  write(u6,*) '      Set atom:  ',i,Atoms(i)
  do j=i+1,nSets
    write(u6,*) '      Sets bond: ',i,j,Bonds(i,j)
  end do
end do

!-- No symmetry allowed.

call Get_iScalar('nSym',nSym)
if (nSym /= 1) then
  write(u6,*)
  write(u6,*) ' You have specified symmetry. The keyword "SELEctive" in FFPT is incompatible with this.'
  write(u6,*) ' Aborting....'
  call Abend()
end if

!-- Collect the overlap and some auxiliaries for LoProp.

call mma_allocate(oType,nBas(1),label='Orbital_Type')
call mma_allocate(Center,nBas(1),label='Center_Index')
call Get_iArray('Orbital Type',oType,nBas(1))
call Get_iArray('Center Index',Center,nBas(1))
do i=1,nBas(1)
  if (oType(i) /= 1 .and. oType(i) /= 0) then
    write(u6,*) 'Orbital type vector is corrupted!'
    call Abend()
  end if
end do

iOpt2 = ibset(0,sNoOri)
iOpt1 = ibset(0,sOpSiz)
iOpt0 = 0
iComp = 1
Label = 'MltPl  0'
iRc = -1
iSymLbl = 1
nInts = 0
call iRdOne(iRc,iOpt1,Label,iComp,idum,iSymLbl)
if (iRc == 0) nInts = idum(1)
call mma_allocate(STr,nInts+4,label='SMatTr')
call RdOne(iRc,iOpt0,Label,iComp,STr,iSymLbl)
if (iRc /= 0) then
  write(u6,*) 'Error reading overlap matrix in SELECTLOC!'
  call Abend()
end if
!-- Let's be square.
call mma_allocate(SSq,nBas(1),nBas(1),label='SMatSq')
call Square(STr,SSq,1,nBas(1),nBas(1))

!-- Call the localization utility and get the transformation matrix.

call mma_allocate(T,nBas(1),nBas(1),label='T')
call mma_allocate(Tinv,nBas(1),nBas(1),label='Tinv')
call Localize_LoProp(T,Tinv,nBas(1),SSq,Center,oType)
if (DeBug) then
  call RecPrt('Total transfMat',' ',T,nBas(1),nBas(1))
end if

!-- Transform the perturbation to the LoProp basis. FFPT accumulates the
!   perturbation to H0, but we only want the perturbation V, hence first
!   a subtraction is necessary.

call mma_allocate(HVac,nInts+4,label='VacH0')
if (LCumulate) then
  Label = 'OneHam  '
else
  Label = 'OneHam 0'
end if
iRc = -1
call RdOne(iRc,iOpt2,Label,iComp,HVac,iSymLbl)
if (iRc /= 0) then
  write(u6,*) 'Error reading H0 in SELECTLOC!'
  call Abend()
end if
call mma_allocate(V,nInts,label='Pert')
V(1:nInts) = H0(1:nInts)-HVac(1:nInts)
H01 = H0(nInts+1)
H02 = H0(nInts+2)
H03 = H0(nInts+3)
H04 = H0(nInts+4)
!----But first translate the perturbation origo to the relevant centre
call TransNow(V,STr)
!----You may proceed.
call mma_allocate(VS,nBas(1),nBas(1),label='PertSq')
call Square(V,VS,1,nBas(1),nBas(1))
if (DeBug) then
  call RecPrt('Pert:(Basis:ord)',' ',VS,nBas(1),nBas(1))
end if

call mma_allocate(TEMP,nBas(1),nBas(1),label='TEMP')
call mma_allocate(VLoP,nBas(1),nBas(1),label='PertL')
!----Go to basis where overlap matrix, S, is diagonal.
call DGEMM_('T','N',nBas(1),nBas(1),nBas(1),One,T,nBas(1),VS,nBas(1),Zero,TEMP,nBas(1))
call DGEMM_('N','N',nBas(1),nBas(1),nBas(1),One,TEMP,nBas(1),T,nBas(1),Zero,VLoP,nBas(1))
if (DeBug) then
  call RecPrt('Pert:(Basis:LoP)',' ',VLoP,nBas(1),nBas(1))
end if

!-- Set elements to zero as designated in input. The routine below is
!   far from optimal, but we are not in need of great speed here so....

do i=1,nBas(1)
  do j=1,nBas(1)
    Siff = Zero
    do k=1,nSets
      do l=k,nSets
        if (k /= l) then
          if (.not. Bonds(k,l)) cycle
        else
          if (.not. Atoms(k)) cycle
        end if
        OneOrNot1 = (i >= iSelection(1,k)) .and. (i <= iSelection(2,k))
        OneOrNot2 = (j >= iSelection(1,l)) .and. (j <= iSelection(2,l))
        OneOrNot3 = (i >= iSelection(1,l)) .and. (i <= iSelection(2,l))
        OneOrNot4 = (j >= iSelection(1,k)) .and. (j <= iSelection(2,k))
        OneOrNot = (OneOrNot1 .and. OneOrNot2) .or. (OneOrNot3 .and. OneOrNot4)
        CrazySet = (OneOrNot1 .and. OneOrNot2) .and. (OneOrNot3 .and. OneOrNot4)
        if (CrazySet .and. Bonds(k,l) .and. (.not. Atoms(k))) then
          write(u6,*) 'Your set selection is not exclusive!'
        end if
        if (OneOrNot) Siff = One
        !-- Here we enable to set the weight in the bond-domain to some
        !   other number than one.
        if (OneOrNot .and. (Atoms(k) .and. Bonds(k,l))) then
          ! FIXME
          write(u6,*) 'Bug! SiffBond is uninitialized!'
          call Abend()
          !Siff = SiffBond
        end if
      end do
    end do
    VLoP(j,i) = VLoP(j,i)*Siff
  end do
end do

!-- Transform back. Due to the non-unitarian and non-orthogonal basis
!   the inverse is is contravariant (if the transformation was
!   covariant). See the Book by Lanczos.

call DGEMM_('T','N',nBas(1),nBas(1),nBas(1),One,Tinv,nBas(1),VLoP,nBas(1),Zero,TEMP,nBas(1))
call DGEMM_('N','N',nBas(1),nBas(1),nBas(1),One,TEMP,nBas(1),Tinv,nBas(1),Zero,VS,nBas(1))
if (DeBug) then
  call RecPrt('Pert.Zeroed',' ',VS,nBas(1),nBas(1))
end if

!-- Add this perturbation to the one-electron hamiltonian.

kaunter = 0
do i=1,nBas(1)
  do j=1,i
    kaunter = kaunter+1
    H0(kaunter) = HVac(kaunter)+VS(j,i)
  end do
end do
!----And don't forget the orgio and the nuclear repulsion.
H0(nInts+1) = H01
H0(nInts+2) = H02
H0(nInts+3) = H03
H0(nInts+4) = H04

!-- Deallocations en masse.

call mma_deallocate(oType)
call mma_deallocate(Center)
call mma_deallocate(STr)
call mma_deallocate(SSq)
call mma_deallocate(T)
call mma_deallocate(Tinv)
call mma_deallocate(HVac)
call mma_deallocate(V)
call mma_deallocate(VS)
call mma_deallocate(TEMP)
call mma_deallocate(VLoP)

!-- Exit

write(u6,*)
write(u6,*) '  ....Done!'
write(u6,*)

return

end subroutine SelectLoc
