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
!> @param[in,out] H0    The one-electron Hamiltonain with perturbations so far. On output the localized perturbed one-electron Hamiltonian
!> @param[in]     nSize Size of the triangular \p H0 with the additional origo and nuclear contribution
!***********************************************************************

subroutine SelectLoc(H0,nSize)

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSize
real(kind=wp), intent(inout) :: H0(nSize)
#include "input.fh"
#include "WrkSpc.fh"
character(len=8) :: Label
logical(kind=iwp) :: OneOrNot1, OneOrNot2, OneOrNot3, OneOrNot4, OneOrNot, CrazySet
integer(kind=iwp) :: idum(1), i, iHVac, ind, iOpt0, iOpt1, iOpt2, ipCenter, ipSSq, ipSTr, ipT, ipTinv, ipType, iRc, iSymLbl, &
                     iTEMP, iV, iVLoP, iVS, j, k, kaunter, l, nInts
real(kind=wp) :: H01, H02, H03, H04, Siff
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

call GetMem('Orbital_Type','Allo','Inte',ipType,nBas(1))
call GetMem('Center_Index','Allo','Inte',ipCenter,nBas(1))
call Get_iArray('Orbital Type',iWork(ipType),nBas(1))
call Get_iArray('Center Index',iWork(ipCenter),nBas(1))
do i=ipType,ipType+nBas(1)-1
  if (iWork(i) /= 1 .and. iWork(i) /= 0) then
    write(u6,*) 'Orbital type vector is corrupted!'
    call Abend()
  end if
end do

iOpt2 = 2
iOpt1 = 1
iOpt0 = 0
Label = 'MltPl  0'
iRc = -1
iSymLbl = 1
nInts = 0
call iRdOne(iRc,iOpt1,Label,1,idum,iSymLbl)
if (iRc == 0) nInts = idum(1)
call GetMem('SMatTr','Allo','Real',ipSTr,nInts+4)
call RdOne(iRc,iOpt0,Label,1,Work(ipSTr),iSymLbl)
if (iRc /= 0) then
  write(u6,*) 'Error reading overlap matrix in SELECTLOC!'
  call Abend()
end if
!-- Let's be square.
call GetMem('SMatSq','Allo','Real',ipSSq,nBas(1)**2)
call Square(Work(ipSTr),Work(ipSSq),1,nBas(1),nBas(1))

!-- Call the localization utility and get the transformation matrix.

call GetMem('T','Allo','Real',ipT,nBas(1)**2)
call GetMem('Tinv','Allo','Real',ipTinv,nBas(1)**2)
call Localize_LoProp(Work(ipT),Work(ipTinv),nBas(1),Work(ipSSq),iWork(ipCenter),iWork(ipType))
if (DeBug) then
  call RecPrt('Total transfMat',' ',Work(ipT),nBas(1),nBas(1))
end if

!-- Transform the perturbation to the LoProp basis. FFPT accumulates the
!   perturbation to H0, but we only want the perturbation V, hence first
!   a subtraction is necessary.

call GetMem('VacH0','Allo','Real',iHVac,nInts+4)
if (LCumulate) then
  Label = 'OneHam  '
else
  Label = 'OneHam 0'
end if
iRc = -1
call RdOne(iRc,iOpt2,Label,1,Work(iHVac),iSymLbl)
if (iRc /= 0) then
  write(u6,*) 'Error reading H0 in SELECTLOC!'
  call Abend()
end if
call GetMem('Pert','Allo','Real',iV,nInts)
do i=1,nInts
  Work(iV+i-1) = H0(i)-Work(iHVac+i-1)
end do
H01 = H0(nInts+1)
H02 = H0(nInts+2)
H03 = H0(nInts+3)
H04 = H0(nInts+4)
!----But first translate the perturbation origo to the relevant centre
call TransNow(iV,ipSTr)
!----You may proceed.
call GetMem('PertSq','Allo','Real',iVS,nBas(1)**2)
call Square(Work(iV),Work(iVS),1,nBas(1),nBas(1))
if (DeBug) then
  call RecPrt('Pert:(Basis:ord)',' ',Work(iVS),nBas(1),nBas(1))
end if

call GetMem('TEMP','Allo','Real',iTEMP,nBas(1)**2)
call GetMem('PertL','Allo','Real',iVLoP,nBas(1)**2)
!----Go to basis where overlap matrix, S, is diagonal.
call DGEMM_('T','N',nBas(1),nBas(1),nBas(1),One,Work(ipT),nBas(1),Work(iVS),nBas(1),Zero,Work(iTEMP),nBas(1))
call DGEMM_('N','N',nBas(1),nBas(1),nBas(1),One,Work(iTEMP),nBas(1),Work(ipT),nBas(1),Zero,Work(iVLoP),nBas(1))
if (DeBug) then
  call RecPrt('Pert:(Basis:LoP)',' ',Work(iVLop),nBas(1),nBas(1))
end if

!-- Set elements to zero as designated in input. The routine below is
!   far from optimal, but we are not in need of great speed here so....

kaunter = 0
do i=1,nBas(1)
  do j=1,nBas(1)
    Siff = Zero
    do k=1,nSets
      do l=k,nSets
        if (k /= l) then
          if (.not. Bonds(k,l)) goto 999
        else
          if (.not. Atoms(k)) goto 999
        end if
        OneOrNot1 = i >= iSelection(1,k) .and. i <= iSelection(2,k)
        OneOrNot2 = j >= iSelection(1,l) .and. j <= iSelection(2,l)
        OneOrNot3 = i >= iSelection(1,l) .and. i <= iSelection(2,l)
        OneOrNot4 = j >= iSelection(1,k) .and. j <= iSelection(2,k)
        OneOrNot = (OneOrNot1 .and. OneOrNot2) .or. (OneOrNot3 .and. OneOrNot4)
        CrazySet = (OneOrNot1 .and. OneOrNot2) .and. (OneOrNot3 .and. OneOrNot4)
        if (CrazySet .and. Bonds(k,l) .and. .not. Atoms(k)) then
          write(u6,*) 'Your set selection is not exclusive!'
        end if
        if (OneOrNot) Siff = One
        !-- Here we enable to set the weight in the bond-domain to some
        !   other number than one.
        if (OneOrNot .and. (Atoms(k) .and. Bonds(k,l))) Siff = SiffBond
999     continue
      end do
    end do
    Work(iVLop+kaunter) = Work(iVLop+kaunter)*Siff
    kaunter = kaunter+1
  end do
end do

!-- Transform back. Due to the non-unitarian and non-orthogonal basis
!   the inverse is is contravariant (if the transformation was
!   covariant). See the Book by Lanczos.

call DGEMM_('T','N',nBas(1),nBas(1),nBas(1),One,Work(ipTinv),nBas(1),Work(iVLop),nBas(1),Zero,Work(iTEMP),nBas(1))
call DGEMM_('N','N',nBas(1),nBas(1),nBas(1),One,Work(iTEMP),nBas(1),Work(ipTinv),nBas(1),Zero,Work(iVS),nBas(1))
if (DeBug) then
  call RecPrt('Pert.Zeroed',' ',Work(iVS),nBas(1),nBas(1))
end if

!-- Add this perturbation to the one-electron hamiltonian.

kaunter = 0
do i=1,nBas(1)
  do j=1,i
    kaunter = kaunter+1
    ind = (i-1)*nBas(1)+j
    H0(kaunter) = Work(iHVac+kaunter-1)+Work(iVS+ind-1)
  end do
end do
!----And don't forget the orgio and the nuclear repulsion.
H0(nInts+1) = H01
H0(nInts+2) = H02
H0(nInts+3) = H03
H0(nInts+4) = H04

!-- Deallocations en masse.

call GetMem('Orbital_Type','Free','Inte',ipType,nBas(1))
call GetMem('Center_Index','Free','Inte',ipCenter,nBas(1))
call GetMem('SMatSq','Free','Real',ipSSq,nBas(1)**2)
call GetMem('T','Free','Real',ipT,nBas(1)**2)
call GetMem('Tinv','Free','Real',ipTinv,nBas(1)**2)
call GetMem('VacH0','Free','Real',iHVac,nInts+4)
call GetMem('Pert','Free','Real',iV,nInts)
call GetMem('PertSq','Free','Real',iVS,nBas(1)**2)
call GetMem('TEMP','Free','Real',iTEMP,nBas(1)**2)
call GetMem('PertL','Free','Real',iVLoP,nBas(1)**2)
call GetMem('SMatTr','Free','Real',ipSTr,nInts+4)

!-- Exit

write(u6,*)
write(u6,*) '  ....Done!'
write(u6,*)

return

end subroutine SelectLoc
