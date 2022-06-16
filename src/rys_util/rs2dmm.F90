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
! Copyright (C) 1990,1991, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine Rs2Dmm(xyz2D,nArg,lRys,nabMax,ncdMax,PAWP,QCWQ,B10,B00,B01,la,lb,lc,ld,IfHss,IfGrd)
!***********************************************************************
!                                                                      *
! Object: to compute the 2-dimensional integrals of the Rys            *
!         quadrature. The z components are assumed to be pre-          *
!         conditioned with the weights of the roots of the             *
!         Rys polynomial.                                              *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
! Modified loop structure for RISC 1991 R. Lindh Dept. of Theoretical  *
! Chemistry, University of Lund, Sweden.                               *
!***********************************************************************

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nArg, lRys, nabMax, ncdMax, la, lb, lc, ld
real(kind=wp), intent(inout) :: xyz2D(nArg*lRys,3,0:nabMax,0:ncdMax)
real(kind=wp), intent(in) :: PAWP(nArg*lRys,3), QCWQ(nArg*lRys,3), B10(nArg*lRys,3), B00(nArg*lRys,3), B01(nArg*lRys,3)
logical(kind=iwp), intent(in) :: IfHss(4,3,4,3), IfGrd(3,4)
integer(kind=iwp) :: iab, iCar, icd, llab, llcd, mabMax, mcdMax
real(kind=wp) :: Fac1, Fac2, Fact

!iRout = 15
!iPrint = nPrint(iRout)
!if (iPrint >= 99) then
!  if (nabMax > 0) call RecPrt('PAWP',' ',PAWP,nArg,lRys*3)
!  if (ncdMax > 0) call RecPrt('QCWQ',' ',QCWQ,nArg,lRys*3)
!  call RecPrt(' B10',' ',B10,nArg*lRys,3)
!  call RecPrt(' B00',' ',B00,nArg*lRys,3)
!  call RecPrt(' B01',' ',B01,nArg*lRys,3)
!end if

! Compute 2D integrals with index (0,0). Observe that the z
! component already contains the weight factor.

xyz2D(:,1:2,0,0) = One

! Compute 2D integrals with index (i,0)

do iCar=1,3
  llab = 0
  llcd = 0
  if (IfHss(2,iCar,2,iCar) .or. IfHss(1,iCar,1,iCar)) llab = 2
  if (IfGrd(iCar,2) .or. IfGrd(iCar,1)) llab = max(llab,1)
  if (IfHss(3,iCar,3,iCar) .or. IfHss(4,iCar,4,iCar)) llcd = 2
  if (IfGrd(iCar,3) .or. IfGrd(iCar,4)) llcd = max(llcd,1)

  mabMax = la+lb+llab
  mcdMax = lc+ld+llcd

  if (mabMax /= 0) then
    xyz2D(:,iCar,1,0) = PAWP(:,iCar)*xyz2D(:,iCar,0,0)
  end if
  if (mabMax-1 == 1) then
    xyz2D(:,iCar,2,0) = PAWP(:,iCar)*xyz2D(:,iCar,1,0)+B10(:,iCar)*xyz2D(:,iCar,0,0)
  else if (mabMax-1 > 1) then
    Fact = One
    do iab=1,mabMax-1
      xyz2D(:,iCar,iab+1,0) = PAWP(:,iCar)*xyz2D(:,iCar,iab,0)+Fact*B10(:,iCar)*xyz2D(:,iCar,iab-1,0)
      Fact = Fact+One
    end do
  end if

  ! Compute 2D integrals with index (0,i)

  if (mcdMax /= 0) then
    xyz2D(:,iCar,0,1) = QCWQ(:,iCar)*xyz2D(:,iCar,0,0)
  end if
  if (mcdMax-1 == 1) then
    xyz2D(:,iCar,0,2) = QCWQ(:,iCar)*xyz2D(:,iCar,0,1)+B01(:,iCar)*xyz2D(:,iCar,0,0)
  else if (mcdMax-1 > 1) then
    Fact = One
    do icd=1,mcdMax-1
      xyz2D(:,iCar,0,icd+1) = QCWQ(:,iCar)*xyz2D(:,iCar,0,icd)+Fact*B01(:,iCar)*xyz2D(:,iCar,0,icd-1)
      Fact = Fact+One
    end do
  end if

  ! Compute 2D integrals with index (i,iCar,j)

  if (mcdMax <= mabMax) then
    Fac1 = One
    do icd=1,mcdMax
      xyz2D(:,iCar,1,icd) = PAWP(:,iCar)*xyz2D(:,iCar,0,icd)+Fac1*B00(:,iCar)*xyz2D(:,iCar,0,icd-1)
      if (mabMax-1 == 1) then
        xyz2D(:,iCar,2,icd) = PAWP(:,iCar)*xyz2D(:,iCar,1,icd)+B10(:,iCar)*xyz2D(:,iCar,0,icd)+ &
                              Fac1*B00(:,iCar)*xyz2D(:,iCar,1,icd-1)
      else if (mabMax-1 > 1) then
        Fac2 = One
        do iab=1,mabMax-1
          xyz2D(:,iCar,iab+1,icd) = PAWP(:,iCar)*xyz2D(:,iCar,iab,icd)+Fac2*B10(:,iCar)*xyz2D(:,iCar,iab-1,icd)+ &
                                    Fac1*B00(:,iCar)*xyz2D(:,iCar,iab,icd-1)
          Fac2 = Fac2+One
        end do
      end if
      Fac1 = Fac1+One
    end do
  else
    Fac1 = One
    do iab=1,mabMax
      xyz2D(:,iCar,iab,1) = QCWQ(:,iCar)*xyz2D(:,iCar,iab,0)+Fac1*B00(:,iCar)*xyz2D(:,iCar,iab-1,0)
      if (mcdMax-1 == 1) then
        xyz2D(:,iCar,iab,2) = QCWQ(:,iCar)*xyz2D(:,iCar,iab,1)+B01(:,iCar)*xyz2D(:,iCar,iab,0)+ &
                              Fac1*B00(:,iCar)*xyz2D(:,iCar,iab-1,1)
      else if (mcdMax-1 > 1) then
        Fac2 = One
        do icd=1,mcdMax-1
          xyz2D(:,iCar,iab,icd+1) = QCWQ(:,iCar)*xyz2D(:,iCar,iab,icd)+Fac2*B01(:,iCar)*xyz2D(:,iCar,iab,icd-1)+ &
                                    Fac1*B00(:,iCar)*xyz2D(:,iCar,iab-1,icd)
          Fac2 = Fac2+One
        end do
      end if
      Fac1 = Fac1+One
    end do
  end if
end do

!if (iPrint >= 99) then
!  do iab=0,nabMax
!    do icd=0,ncdMax
!      write(Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(x)'
!      call RecPrt(Label,' ',xyz2D(:,1,iab,icd),nArg,lRys)
!      write(Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(y)'
!      call RecPrt(Label,' ',xyz2D(:,2,iab,icd),nArg,lRys)
!      write(Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(z)'
!      call RecPrt(Label,' ',xyz2D(:,3,iab,icd),nArg,lRys)
!    end do
!  end do
!end if

return

end subroutine Rs2Dmm
