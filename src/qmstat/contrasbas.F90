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

subroutine ContRASBas(nStatePrim,NonH,NonS,Eig2)

use qmstat_global, only: ContrStateB, dLvlShift, HmatSOld, HmatState, iLvlShift, iPrint, nLvlShift, nState, ThrsCont
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, one
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nStatePrim
real(kind=wp), intent(in) :: NonH(nTri_Elem(nStatePrim))
real(kind=wp), intent(inout) :: NonS(nTri_Elem(nStatePrim))
real(kind=wp), intent(out) :: Eig2(nStatePrim,nStatePrim)
integer(kind=iwp) :: i, ii, iS, iState, iT, kaunt, kaunter, nLvlInd, nStateRed, nTri
real(kind=wp) :: sss, x
real(kind=wp), allocatable :: Eig1(:,:), RedHSq(:,:), RedHTr(:), SqH(:,:), TEMP(:,:)

! Hi y'all

write(u6,*) '     ----- Constructing CASSI eigenstates.'

! Diagonalize overlap matrix.

call mma_allocate(Eig1,nStatePrim,nStatePrim,label='EigV1')
call unitmat(Eig1,nStatePrim)
call Jacob(NonS,Eig1,nStatePrim,nStatePrim)
if (iPrint >= 15) call TriPrt('Diagonal RASSCF overlap matrix',' ',NonS,nStatePrim)

! Construct TS^(-1/2) for canonical orthogonalization.

ii = 0
do i=1,nStatePrim
  ii = ii+i
  x = One/sqrt(max(1.0e-14_wp,NonS(ii)))
  Eig1(:,i) = x*Eig1(:,i)
end do

! Make reductions if requested.

iT = 0
if (ContrStateB) then
  do iS=1,nStatePrim
    kaunt = nTri_Elem(iS)+1
    sss = NonS(kaunt)
    if (sss > ThrsCont) then
      iT = iT+1
      Eig2(:,iT) = Eig1(:,iS)
    end if
  end do
  nStateRed = iT
  write(u6,6199) '  ----- Contraction:',nStatePrim,' ---> ',nStateRed
else
  Eig2(:,:) = Eig1
  nStateRed = nStatePrim
end if

! Transform H and diagonalize in the original basis.

nTri = nTri_Elem(nStateRed)
call mma_allocate(TEMP,nStatePrim,nStatePrim,label='TEMP')
call mma_allocate(SqH,nStatePrim,nStatePrim,label='SqH')
call mma_allocate(RedHSq,nStateRed,nStateRed,label='RedHSq')
call mma_allocate(RedHTr,nTri,label='RedHTr')
call Square(NonH,SqH,1,nStatePrim,nStatePrim)
call Dgemm_('N','N',nStatePrim,nStateRed,nStatePrim,One,SqH,nStatePrim,Eig2,nStatePrim,Zero,TEMP,nStatePrim)
call Dgemm_('T','N',nStateRed,nStateRed,nStatePrim,One,Eig2,nStatePrim,TEMP,nStatePrim,Zero,RedHSq,nStateRed)
call SqToTri_Q(RedHSq,RedHTr,nStateRed)
call Jacob(RedHTr,Eig2,nStateRed,nStatePrim)
call JacOrd(RedHTr,Eig2,nStateRed,nStatePrim)

! At this stage we have eigenvectors to the CASSI states and their
! eigenenergies, hence time to construct the first H_0 and store the
! eigenvectors for subsequent transformations.

call mma_allocate(HMatState,nTri_Elem(nStateRed),label='HMatState')
call mma_allocate(HMatSOld,nTri_Elem(nStateRed),label='HMatSOld')

kaunter = 0
nLvlInd = 1
HMatState(:) = Zero
do iState=1,nStateRed
  kaunter = nTri_Elem(iState)
  HMatState(kaunter) = RedHTr(kaunter)
  ! If requested, introduce level-shift of states.
  if (nLvlShift > 0) then
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
  call RecPrt('RASSI eigenvectors',' ',Eig2,nStatePrim,nStateRed)
end if

! Deallocate.

call mma_deallocate(Eig1)
call mma_deallocate(TEMP)
call mma_deallocate(SqH)
call mma_deallocate(RedHSq)
call mma_deallocate(RedHTr)

! OBSERVE! CAUTION! ATTENTION! The variable nState is defined.

nState = nStateRed

! No parasan!

return

6199 format(A,I3,A,I3)

end subroutine ContRASBas
