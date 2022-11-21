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

subroutine InitB(nMult,natom,nAtQM,nGrdPt,Cord,Grid,T,TT,TTT,Ext,B,IsMM)
! Compute the electrostatic tensor matrix between QM atoms and grid points
! Then form the B = ExtPot[TtT^-1]Tt vector

use espf_global, only: MxExtPotComp
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nMult, natom, nAtQM, nGrdPt, IsMM(natom)
real(kind=wp), intent(in) :: Cord(3,natom), Grid(3,nGrdPt), Ext(MxExtPotComp,natom)
real(kind=wp), intent(out) :: T(nMult,nGrdPt), TT(nMult,nMult), TTT(nGrdPt,nMult), B(nGrdPt)
integer(kind=iwp) :: iMlt, iPL, iPnt, iQM, jAt, jMlt, jPnt, kMlt, kPnt, nOrd
real(kind=wp) :: Det, R, R3, X, Y, Z
real(kind=wp), allocatable :: Scr(:,:)
integer(kind=iwp), external :: iPL_espf

iPL = iPL_espf()
nOrd = nMult/nAtQM

! T

do iPnt=1,nGrdPt
  iQM = 0
  do jAt=1,natom
    if (IsMM(jAt) == 1) cycle
    iQM = iQM+1
    X = Grid(1,iPnt)-Cord(1,jAt)
    Y = Grid(2,iPnt)-Cord(2,jAt)
    Z = Grid(3,iPnt)-Cord(3,jAt)
    R = sqrt(X*X+Y*Y+Z*Z)
    T(nOrd*(iQM-1)+1,iPnt) = One/R
    if (nOrd > 1) then
      R3 = R*R*R
      T(nOrd*(iQM-1)+2,iPnt) = X/R3
      T(nOrd*(iQM-1)+3,iPnt) = Y/R3
      T(nOrd*(iQM-1)+4,iPnt) = Z/R3
    end if
  end do
end do
if (iQM /= nAtQM) then
  write(u6,'(A,I4,A4,I4)') ' Error in espf/initb: iQM != nAtQM ',iQM,' != ',nAtQM
  call Abend()
end if

!call RecPrt('T',' ',T,nMult,nGrdPt)
!nMax = Max(nGrdPt,nMult)
!call mma_allocate(U,nMult,nMax,label='U')
!call mma_allocate(W,nMult,label='W')
!call mma_allocate(V,nMult,nMax,label='V')
!call mma_allocate(Scr,nMult,1,label='Scr')
!call SVD(nMax,nGrdPt,nMult,T,W,.true.,U,.true.,V,iErr,Scr)
!write(u6,*) 'iErr=',iErr
!call RecPrt('U',' ',U,nMult,nMax)
!call RecPrt('w',' ',W,nMult,1)
!call RecPrt('V',' ',V,nMult,nMax)
!call mma_deallocate(U)
!call mma_deallocate(W)
!call mma_deallocate(V)
!call mma_deallocate(Scr)

! TtT

TT(:,:) = Zero
do iMlt=1,nMult
  do jMlt=1,nMult
    do kPnt=1,nGrdPt
      TT(jMlt,iMlt) = TT(jMlt,iMlt)+T(iMlt,kPnt)*T(jMlt,kPnt)
    end do
  end do
end do

! TtT^-1

call mma_allocate(Scr,nMult,nMult,label='Scr')
call minv(TT,Scr,Det,nMult)
TT(:,:) = Scr
call mma_deallocate(Scr)

! [TtT^-1]Tt

TTT(:,:) = Zero
do iMlt=1,nMult
  do jPnt=1,nGrdPt
    do kMlt=1,nMult
      TTT(jPnt,iMlt) = TTT(jPnt,iMlt)+TT(kMlt,iMlt)*T(kMlt,jPnt)
    end do
  end do
end do
if (iPL >= 4) call RecPrt('(TtT)^(-1)Tt matrix in InitB',' ',TTT,nGrdPt,nMult)

! B = ExtPot[TtT^-1]Tt

B(:) = Zero
do iPnt=1,nGrdPt
  iQM = 0
  do jAt=1,natom
    if (IsMM(jAt) == 1) cycle
    iQM = iQM+1
    B(iPnt) = B(iPnt)+Ext(1,jAt)*TTT(iPnt,nOrd*(iQM-1)+1)
    if (nOrd > 1) &
      B(iPnt) = B(iPnt)+Ext(2,jAt)*TTT(iPnt,nOrd*(iQM-1)+2)+Ext(3,jAt)*TTT(iPnt,nOrd*(iQM-1)+3)+Ext(4,jAt)*TTT(iPnt,nOrd*(iQM-1)+4)
  end do
end do
if (iPL >= 4) then
  write(u6,'(A)') ' In InitB (grid coordinates, B value)'
  do iPnt=1,nGrdPt
    write(u6,1234) iPnt,Grid(:,iPnt),B(iPnt)
  end do
end if

return

1234 format(I4,4F12.6)

end subroutine InitB
