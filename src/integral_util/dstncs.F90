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
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Dstncs(Lbls,xyz,mCentr,Angstr,Max_Center,iCols)
!***********************************************************************
!                                                                      *
! Object: to compute distances from a coordinate list                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Molcas, only: LenIn
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Three
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: mCentr, Max_Center, iCols
character(len=LenIn), intent(in) :: Lbls(mCentr)
real(kind=wp), intent(in) :: xyz(3,mCentr), Angstr
integer(kind=iwp) :: i, iC, iCC, ii, iiBst, IsFirst, jC, jCC, MoreToGo
real(kind=wp) :: Fact, r, rr, Thr_D, Thr_R, x1, x2, y1, y2, z1, z2
integer(kind=iwp), allocatable :: iBST(:,:)
real(kind=wp), allocatable :: BST(:)

if (mCentr <= Max_Center) then

  do i=1,2
    write(u6,*)
    if (i == 1) then
      Fact = One
      write(u6,'(19X,A)') ' *************************************** '
      write(u6,'(19X,A)') ' *    InterNuclear Distances / bohr    * '
      write(u6,'(19X,A)') ' *************************************** '
    else
      Fact = Angstr
      write(u6,'(19X,A)') ' ******************************************* '
      write(u6,'(19X,A)') ' *    InterNuclear Distances / angstrom    * '
      write(u6,'(19X,A)') ' ******************************************* '
    end if
    do icc=1,mCentr,iCols
      write(u6,*)
      if (iCols == 6) write(u6,'( 9X,6(5X,I2,1X,A,2X))') (ic,Lbls(ic),ic=icc,min(icc+5,mCentr))
      if (iCols == 5) write(u6,'( 9X,5(5X,I2,1X,A,2X))') (ic,Lbls(ic),ic=icc,min(icc+4,mCentr))

      do jc=icc,mCentr
        x1 = xyz(1,jc)
        y1 = xyz(2,jc)
        z1 = xyz(3,jc)
        if (iCols == 6) &
          write(u6,101) jc,Lbls(jc),(Fact*sqrt((xyz(1,ic)-x1)**2+(xyz(2,ic)-y1)**2+(xyz(3,ic)-z1)**2),ic=icc,min(jc,icc+5,mCentr))
        if (iCols == 5) &
          write(u6,102) jc,Lbls(jc),(Fact*sqrt((xyz(1,ic)-x1)**2+(xyz(2,ic)-y1)**2+(xyz(3,ic)-z1)**2),ic=icc,min(jc,icc+4,mCentr))
      end do
    end do
  end do
  return
else

  write(u6,*)
  write(u6,'(19X,A)') ' ************************************************* '
  write(u6,'(19X,A)') ' **** InterNuclear Distances / bohr, angstrom **** '
  write(u6,'(19X,A)') ' ************************************************* '
  write(u6,*)
  write(u6,'(A)') '     Atom centers         bohr        angstrom'

  !VV Set .false. to get faster printing without sorting.

  if (.true.) then
    Thr_R = (Three/Angstr)**2
    Thr_D = 1.0e-4_wp
    call mma_allocate(BST,mCentr**2,Label='BST')
    call mma_allocate(iBST,2,mCentr**2,Label='iBST')
    iiBST = 0
    do icc=1,mCentr
      x1 = xyz(1,icc)
      y1 = xyz(2,icc)
      z1 = xyz(3,icc)
      do jcc=1,icc-1
        x2 = xyz(1,jcc)
        y2 = xyz(2,jcc)
        z2 = xyz(3,jcc)
        R = (x2-x1)**2+(y2-y1)**2+(z2-z1)**2
        if (R <= Thr_R) then
          iiBST = iiBST+1
          BST(iiBST) = R
          iBST(1,iiBST) = icc
          iBST(2,iiBST) = jcc
        end if
      end do
    end do
#   ifdef _DEBUGPRINT_
    do ii=1,iiBST
      R = sqrt(BST(ii))
      write(u6,'(2(I5,1X,A4),2(F10.6,6X))') iBST(1,ii),Lbls(iBST(1,ii)),iBST(2,ii),Lbls(iBST(2,ii)),R,R*Angstr
    end do
#   endif
    MoreToGo = 1
    do while (MoreToGo == 1)

      ! Find the shortest distance between any atoms

      R = 100.0_wp
      do ii=1,iiBST
        R = min(R,BST(ii))
      end do

      if (R > 90.0_wp) exit

      moretogo = 0
      isfirst = 1
      do ii=1,iiBST

        if (abs(R-BST(ii)) < Thr_D) then
          if (isfirst == 1) then
            RR = sqrt(R)
            write(u6,'(2(I5,1X,A),2(F10.6,6X))') iBST(1,ii),Lbls(iBST(1,ii)),iBST(2,ii),Lbls(iBST(2,ii)),RR,RR*Angstr
            isfirst = 0
          else
            write(u6,'(2(I5,1X,A),2(F10.6,6X))') iBST(1,ii),Lbls(iBST(1,ii)),iBST(2,ii),Lbls(iBST(2,ii))
          end if
          BST(ii) = 100.0_wp ! Effectively remove from the list
          moretogo = 1
        end if
      end do
    end do
    call mma_deallocate(iBST)
    call mma_deallocate(BST)

  else

    Thr_R = Three
    do icc=1,mCentr
      x1 = xyz(1,icc)
      y1 = xyz(2,icc)
      z1 = xyz(3,icc)
      do jcc=1,icc-1
        x2 = xyz(1,jcc)
        y2 = xyz(2,jcc)
        z2 = xyz(3,jcc)
        R = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
        if (R*Angstr <= Thr_R) write(u6,'(2(I5,1X,A),2(F10.6,6X))') icc,Lbls(icc),jcc,Lbls(jcc),R,R*Angstr
      end do
    end do
  end if
end if

return

101 format (I5,1X,A,1X,6(F10.6,6X))
102 format (I5,1X,A,1X,5(F10.6,6X))

end subroutine Dstncs
