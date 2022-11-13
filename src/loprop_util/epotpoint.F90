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

subroutine EPotPoint(Potte,nPick,Pick,DPick,T,Ti,NucNr,nB,iAtom,jAtom,Center)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nPick, Pick(nPick), NucNr, nB, iAtom, jAtom, Center(nB)
real(kind=wp), intent(out) :: Potte(nPick)
real(kind=wp), intent(in) :: DPick(nPick), T(nB,nB), Ti(nB,nB)
integer(kind=iwp) :: iB1, iB2, iC1, iC2, iComp, iOpt, iPo, iPoint, irc, iSmLbl, nB2, nDens
real(kind=wp) :: dEx
character(len=10) :: Label
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: D1ao(:), D_Sq(:,:), DTrans(:,:), PP(:), PSq(:,:), PTr(:,:), TEMP(:,:)

! Loop through all points, pick out the relevant ones, obtain the
! partial expectation value from the appropriate basis and return.

nB2 = nB*(nB+1)/2
call mma_allocate(D_Sq,nB,nB,label='DSq')
call Qpg_darray('D1ao',Found,nDens)
if (Found .and. (nDens /= 0)) then
  call mma_allocate(D1ao,nDens,Label='D1ao')
else
  write(u6,*) 'EPotPoint: D1ao not found.'
  call abend()
end if
call Get_dArray_chk('D1ao',D1ao,nDens)
call Dsq(D1ao,D_Sq,1,nB,nB)
call mma_deallocate(D1ao)
call mma_allocate(TEMP,nB,nB,label='TEMP')
call mma_allocate(DTrans,nB,nB,label='DTrans')

! Contravariant transformation of density matrix.

call DGEMM_('N','N',nB,nB,nB,One,Ti,nB,D_Sq,nB,Zero,TEMP,nB)
call DGEMM_('N','T',nB,nB,nB,One,TEMP,nB,Ti,nB,Zero,DTrans,nB)
call mma_allocate(PP,nB2+4,label='Points')
call mma_allocate(PSq,nB,nB,label='PointsSq')
call mma_allocate(PTr,nB,nB,label='PointsTr')
do iPoint=1,nPick
  iPo = Pick(iPoint)
  write(Label,'(A3,I5)') 'EF0',iPo
  irc = -1
  iOpt = 0
  iSmLbl = 0
  iComp = 1
  call RdOne(irc,iOpt,Label,iComp,PP,iSmLbl)
  call Square(PP,PSq,1,nB,nB)

  ! Covariant transformation of the matrix for the potential in this particular point.

  call DGEMM_('T','N',nB,nB,nB,One,T,nB,PSq,nB,Zero,TEMP,nB)
  call DGEMM_('N','N',nB,nB,nB,One,TEMP,nB,T,nB,Zero,PTr,nB)
  dEx = Zero

  ! The usual stuff to get the localized value.

  do iB1=1,nB
    do iB2=1,nB
      iC1 = Center(iB1)
      iC2 = Center(iB2)
      if (((iC1 == iAtom) .and. (iC2 == jAtom)) .or. ((iC1 == jAtom) .and. (iC2 == iAtom))) then
        dEx = dEx+DTrans(iB2,iB1)*PTr(iB2,iB1)
      end if
    end do
  end do

  ! Accumulate in the return vector.

  if (iAtom == jAtom) then
    Potte(iPoint) = -dEx+real(NucNr,kind=wp)/DPick(iPoint)
  else
    Potte(iPoint) = -dEx
  end if
end do

! Deallocate.

call mma_deallocate(D_Sq)
call mma_deallocate(TEMP)
call mma_deallocate(DTrans)
call mma_deallocate(PP)
call mma_deallocate(PSq)
call mma_deallocate(PTr)

return

end subroutine EPotPoint
