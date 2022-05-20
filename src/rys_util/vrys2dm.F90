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

subroutine vRys2Dm(xyz2D,nArg,lRys,nabMax,ncdMax,PAWP,QCWQ,B10,B00,B01,la,lb,lc,ld,IfGrad)
!***********************************************************************
!                                                                      *
!     Object: to compute the 2-dimensional integrals of the Rys        *
!             quadrature. The z components are assumed to be pre-      *
!             conditioned with the weights of the roots of the         *
!             Rys polynomial.                                          *
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
integer(kind=iwp) :: nArg, lRys, nabMax, ncdMax, la, lb, lc, ld
real(kind=wp) :: xyz2D(nArg*lRys,3,0:nabMax,0:ncdMax), PAWP(nArg*lRys,3), QCWQ(nArg*lRys,3), B10(nArg*lRys), B00(nArg*lRys), &
                 B01(nArg*lRys)
logical(kind=iwp) :: IfGrad(3,4)
integer(kind=iwp) :: i, iab, iCar, icd, llab, llcd, mabMax, mcdMax
real(kind=wp) :: Fac1, Fac2, Fact, temp1, temp2, temp3
#ifdef _DEBUGPRINT_
character(len=30) :: Label
#endif

#ifdef _DEBUGPRINT_
if (nabMax > 0) call RecPrt('PAWP',' ',PAWP,nArg,lRys*3)
if (ncdMax > 0) call RecPrt('QCWQ',' ',QCWQ,nArg,lRys*3)
call RecPrt(' B10',' ',B10,nArg*lRys,3)
call RecPrt(' B00',' ',B00,nArg*lRys,3)
call RecPrt(' B01',' ',B01,nArg*lRys,3)
#endif

! Compute 2D integrals with index (0,0). Observe that the z
! component already contains the weight factor.

call dcopy_(2*nArg*lRys,[One],0,xyz2D(1,1,0,0),1)

! Compute 2D integrals with index (i,0)

do iCar=1,3
  llab = 0
  if (IfGrad(iCar,1) .or. IfGrad(iCar,2)) llab = 1
  llcd = 0
  if (IfGrad(iCar,3) .or. IfGrad(iCar,4)) llcd = 1
  mabMax = la+lb+llab
  mcdMax = lc+ld+llcd

  if (mabMax /= 0) then
    do i=1,nArg*lRys
      xyz2D(i,iCar,1,0) = PAWP(i,iCar)*xyz2D(i,iCar,0,0)
    end do
  end if
  if (mabMax-1 == 1) then
    do i=1,nArg*lRys
      xyz2D(i,iCar,2,0) = PAWP(i,iCar)*xyz2D(i,iCar,1,0)+B10(i)*xyz2D(i,iCar,0,0)
    end do
  else if (mabMax-1 > 1) then
    Fact = One
    do iab=1,mabMax-1
      do i=1,nArg*lRys
        temp1 = PAWP(i,iCar)*xyz2D(i,iCar,iab,0)
        temp2 = Fact*B10(i)*xyz2D(i,iCar,iab-1,0)
        xyz2D(i,iCar,iab+1,0) = temp1+temp2
      end do
      Fact = Fact+One
    end do
  end if

  ! Compute 2D integrals with index (0,i)

  if (mcdMax /= 0) then
    do i=1,nArg*lRys
      xyz2D(i,iCar,0,1) = QCWQ(i,iCar)*xyz2D(i,iCar,0,0)
    end do
  end if
  if (mcdMax-1 == 1) then
    do i=1,nArg*lRys
      xyz2D(i,iCar,0,2) = QCWQ(i,iCar)*xyz2D(i,iCar,0,1)+B01(i)*xyz2D(i,iCar,0,0)
    end do
  else if (mcdMax-1 > 1) then
    Fact = One
    do icd=1,mcdMax-1
      do i=1,nArg*lRys
        temp1 = QCWQ(i,iCar)*xyz2D(i,iCar,0,icd)
        temp2 = Fact*B01(i)*xyz2D(i,iCar,0,icd-1)
        xyz2D(i,iCar,0,icd+1) = temp1+temp2
      end do
      Fact = Fact+One
    end do
  end if

  ! Compute 2D integrals with index (i,iCar,j)

  if (mcdMax <= mabMax) then
    Fac1 = One
    do icd=1,mcdMax
      do i=1,nArg*lRys
        xyz2D(i,iCar,1,icd) = PAWP(i,iCar)*xyz2D(i,iCar,0,icd)+Fac1*B00(i)*xyz2D(i,iCar,0,icd-1)
      end do
      if (mabMax-1 == 1) then
        do i=1,nArg*lRys
          xyz2D(i,iCar,2,icd) = PAWP(i,iCar)*xyz2D(i,iCar,1,icd)+B10(i)*xyz2D(i,iCar,0,icd)+Fac1*B00(i)*xyz2D(i,iCar,1,icd-1)
        end do
      else if (mabMax-1 > 1) then
        Fac2 = One
        do iab=1,mabMax-1
          do i=1,nArg*lRys
            temp1 = PAWP(i,iCar)*xyz2D(i,iCar,iab,icd)
            temp2 = Fac2*B10(i)*xyz2D(i,iCar,iab-1,icd)
            temp3 = Fac1*B00(i)*xyz2D(i,iCar,iab,icd-1)
            xyz2D(i,iCar,iab+1,icd) = temp1+temp2+temp3
          end do
          Fac2 = Fac2+One
        end do
      end if
      Fac1 = Fac1+One
    end do
  else
    Fac1 = One
    do iab=1,mabMax
      do i=1,nArg*lRys
        xyz2D(i,iCar,iab,1) = QCWQ(i,iCar)*xyz2D(i,iCar,iab,0)+Fac1*B00(i)*xyz2D(i,iCar,iab-1,0)
      end do
      if (mcdMax-1 == 1) then
        do i=1,nArg*lRys
          xyz2D(i,iCar,iab,2) = QCWQ(i,iCar)*xyz2D(i,iCar,iab,1)+B01(i)*xyz2D(i,iCar,iab,0)+Fac1*B00(i)*xyz2D(i,iCar,iab-1,1)
        end do
      else if (mcdMax-1 > 1) then
        Fac2 = One
        do icd=1,mcdMax-1
          do i=1,nArg*lRys
            temp1 = QCWQ(i,iCar)*xyz2D(i,iCar,iab,icd)
            temp2 = Fac2*B01(i)*xyz2D(i,iCar,iab,icd-1)
            temp3 = Fac1*B00(i)*xyz2D(i,iCar,iab-1,icd)
            xyz2D(i,iCar,iab,icd+1) = temp1+temp2+temp3
          end do
          Fac2 = Fac2+One
        end do
      end if
      Fac1 = Fac1+One
    end do
  end if
end do

#ifdef _DEBUGPRINT_
do iab=0,nabMax
  do icd=0,ncdMax
    write(Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(x)'
    call RecPrt(Label,' ',xyz2D(:,1,iab,icd),nArg,lRys)
    write(Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(y)'
    call RecPrt(Label,' ',xyz2D(:,2,iab,icd),nArg,lRys)
    write(Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(z)'
    call RecPrt(Label,' ',xyz2D(:,3,iab,icd),nArg,lRys)
  end do
end do
#endif

return

end subroutine vRys2Dm
