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

subroutine ContRASBas(nBas,nStatePrim,iNonH,iNonS,iEig2)

use Constants, only: Zero, one
use Definitions, only: wp, iwp, u6

implicit none
#include "maxi.fh"
#include "qminp.fh"
#include "qm2.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: nBas(MxSym), nStatePrim, iNonH, iNonS, iEig2
integer(kind=iwp) :: i, iEig1, ii, ind, iRedHSq, iRedHTr, iS, iSqH, iState, iT, iTEMP, j, jState, k, kaunt, kaunter, nLvlInd, nTri
real(kind=wp) :: sss, x

! Hi y'all

write(u6,*) '     ----- Constructing CASSI eigenstates.'

! Diagonalize overlap matrix.

call GetMem('EigV1','Allo','Real',iEig1,nStatePrim**2)
kaunter = 0
do i=1,nStatePrim
  do j=1,nStatePrim
    if (i == j) then
      Work(iEig1+kaunter) = One
    else
      Work(iEig1+kaunter) = Zero
    end if
    kaunter = kaunter+1
  end do
end do
call Jacob(Work(iNonS),Work(iEig1),nStatePrim,nStatePrim)
if (iPrint >= 15) then
  call TriPrt('Diagonal RASSCF overlap matrix',' ',Work(iNonS),nStatePrim)
end if

! Construct TS^(-1/2) for canonical orthogonalization.

ii = 0
do i=1,nStatePrim
  ii = ii+i
  x = One/sqrt(max(1.0e-14_wp,Work(iNonS-1+ii)))
  do k=1,nStatePrim
    ind = k+nStatePrim*(i-1)-1
    Work(iEig1+ind) = x*Work(iEig1+ind)
  end do
end do

! Make reductions if requested.

call GetMem('RedEigV1','Allo','Real',iEig2,nStatePrim**2)
iT = 0
if (ContrStateB) then
  do iS=1,nStatePrim
    kaunt = iS*(iS+1)/2-1
    sss = Work(iNonS+kaunt)
    if (sss > ThrsCont) then
      iT = iT+1
      call dcopy_(nStatePrim,Work(iEig1+nStatePrim*(iS-1)),1,Work(iEig2+nStatePrim*(iT-1)),1)
    end if
  end do
  nStateRed = iT
  write(u6,6199) '  ----- Contraction:',nStatePrim,' ---> ',nStateRed
else
  call dcopy_(nStatePrim**2,Work(iEig1),1,Work(iEig2),1)
  nStateRed = nStatePrim
end if

! Transform H and diagonalize in the original basis.

nTri = nStateRed*(nStateRed+1)/2
call GetMem('TEMP','Allo','Real',iTEMP,nStatePrim**2)
call GetMem('SqH','Allo','Real',iSqH,nStatePrim**2)
call GetMem('RedHSq','Allo','Real',iRedHSq,nStateRed**2)
call GetMem('RedHTr','Allo','Real',iRedHTr,nTri)
call Square(Work(iNonH),Work(iSqH),1,nStatePrim,nStatePrim)
call Dgemm_('N','N',nStatePrim,nStateRed,nStatePrim,One,Work(iSqH),nStatePrim,Work(iEig2),nStatePrim,Zero,Work(iTEMP),nStatePrim)
call Dgemm_('T','N',nStateRed,nStateRed,nStatePrim,One,Work(iEig2),nStatePrim,Work(iTEMP),nStatePrim,Zero,Work(iRedHSq),nStateRed)
call SqToTri_Q(Work(iRedHSq),Work(iRedHTr),nStateRed)
call Jacob(Work(iRedHTr),Work(iEig2),nStateRed,nStatePrim)
call JacOrd(Work(iRedHTr),Work(iEig2),nStateRed,nStatePrim)

! At this stage we have eigenvectors to the CASSI states and their
! eigenenergies, hence time to construct the first H_0 and store the
! eigenvectors for subsequent transformations.

kaunter = 0
nLvlInd = 1
do iState=1,nStateRed
  do jState=1,iState
    kaunter = kaunter+1
    HMatState(kaunter) = Zero
  end do
  HMatState(kaunter) = Work(iRedHTr+kaunter-1)
  ! If requested, introduce level-shift of states.
  if (nLvlShift /= 0) then
    if (iState == iLvlShift(nLvlInd)) then
      HMatState(kaunter) = HMatState(kaunter)+dLvlShift(nLvlInd)
      nLvlInd = nLvlInd+1
    end if
  end if
end do

! Print.

if (iPrint >= 10) then
  call TriPrt('RASSI Hamiltonian',' ',HMatState,nStateRed)
  write(u6,*)
  call RecPrt('RASSI eigenvectors',' ',Work(iEig2),nStatePrim,nStateRed)
end if

! Deallocate.

call GetMem('EigV1','Free','Real',iEig1,nStatePrim**2)
call GetMem('TEMP','Free','Real',iTEMP,nStatePrim**2)
call GetMem('SqH','Free','Real',iSqH,nStatePrim**2)
call GetMem('RedHSq','Free','Real',iRedHSq,nStateRed**2)
call GetMem('RedHTr','Free','Real',iRedHTr,nTri)

! OBSERVE! CAUTION! ATTENTION! The variable nState is defined.

nState = nStateRed

! No parasan!

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(nBas)

6199 format(A,I3,A,I3)

end subroutine ContRASBas
