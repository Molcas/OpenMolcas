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
! Copyright (C) 1990,1991,1994, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

subroutine vRys2D(xyz2D,nArg,lRys,nabMax,ncdMax,PAWP,QCWQ,B10,B00,B01)
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
! Further modifications in Jan-Feb. 1994.                              *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
real*8 xyz2D(nArg*lRys*3,0:nabMax,0:ncdMax), PAWP(nArg*lRys*3), QCWQ(nArg*lRys*3), B10(nArg*lRys), B00(nArg*lRys), B01(nArg*lRys)
logical lPAWP, lQCWQ
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
character*30 Label

if (nabMax > 0) call RecPrt('PAWP',' ',PAWP,lRys,nArg*3)
if (ncdMax > 0) call RecPrt('QCWQ',' ',QCWQ,lRys,nArg*3)
call RecPrt(' B10',' ',B10,lRys,nArg)
call RecPrt(' B00',' ',B00,lRys,nArg)
call RecPrt(' B01',' ',B01,lRys,nArg)
#endif

iOffy = nArg*lRys
iOffz = 2*nArg*lRys
if ((nabMax /= 0) .or. (ncdMax /= 0)) then

  ! General code

  ! Store away PAWPz and QCPQz

  lPAWP = nabMax > 2
  lQCWQ = ncdMax > 2
  if ((nabMax >= 1) .and. (ncdMax >= 1)) then
    if (ncdMax <= nabMax) then
      lPAWP = .true.
    else
      lQCWQ = .true.
    end if
  end if
  if (lPAWP) call dcopy_(nArg*lRys,PAWP(1+iOffz),1,xyz2D(1,0,0),1)
  if (lQCWQ) call dcopy_(nArg*lRys,QCWQ(1+iOffz),1,xyz2D(1+iOffy,0,0),1)

  ! Compute 2D integrals with index (i,0)

  if (nabMax == 0) then

  else if (nabMax == 1) then
    do i=1,nArg*lRys
      PAWPz = xyz2D(i+iOffz,1,0)
      xyz2D(i+iOffz,1,0) = PAWPz*xyz2D(i+iOffz,0,0)
    end do
  else if (nabMax >= 2) then
    do i=1,nArg*lRys
      PAWPx = xyz2D(i,1,0)
      xyz2D(i,2,0) = PAWPx*xyz2D(i,1,0)+B10(i)

      PAWPy = xyz2D(i+iOffy,1,0)
      xyz2D(i+iOffy,2,0) = PAWPy*xyz2D(i+iOffy,1,0)+B10(i)

      PAWPz = xyz2D(i+iOffz,1,0)
      xyz2D(i+iOffz,1,0) = PAWPz*xyz2D(i+iOffz,0,0)
      xyz2D(i+iOffz,2,0) = PAWPz*xyz2D(i+iOffz,1,0)+B10(i)*xyz2D(i+iOffz,0,0)
    end do
    if (nabMax > 2) then
      Fact = Two
      do iab=2,nabMax-1
        do i=1,nArg*lRys
          PAWPx = xyz2D(i,1,0)
          PAWPy = xyz2D(i+iOffy,1,0)
          PAWPz = xyz2D(i,0,0)
          temp1x = PAWPx*xyz2D(i,iab,0)
          temp1y = PAWPy*xyz2D(i+iOffy,iab,0)
          temp1z = PAWPz*xyz2D(i+iOffz,iab,0)
          temp2x = Fact*B10(i)*xyz2D(i,iab-1,0)
          temp2y = Fact*B10(i)*xyz2D(i+iOffy,iab-1,0)
          temp2z = Fact*B10(i)*xyz2D(i+iOffz,iab-1,0)
          xyz2D(i,iab+1,0) = temp1x+temp2x
          xyz2D(i+iOffy,iab+1,0) = temp1y+temp2y
          xyz2D(i+iOffz,iab+1,0) = temp1z+temp2z
        end do
        Fact = Fact+One
      end do
    end if
  end if

  ! Compute 2D integrals with index (0,i)

  if (ncdMax == 1) then
    do i=1,nArg*lRys

      QCWQz = xyz2D(i+iOffz,0,1)
      xyz2D(i+iOffz,0,1) = QCWQz*xyz2D(i+iOffz,0,0)
    end do
  else if (ncdMax >= 2) then
    do i=1,nArg*lRys
      QCWQx = xyz2D(i,0,1)
      xyz2D(i,0,2) = QCWQx*xyz2D(i,0,1)+B01(i)

      QCWQy = xyz2D(i+iOffy,0,1)
      xyz2D(i+iOffy,0,2) = QCWQy*xyz2D(i+iOffy,0,1)+B01(i)

      QCWQz = xyz2D(i+iOffz,0,1)
      xyz2D(i+iOffz,0,1) = QCWQz*xyz2D(i+iOffz,0,0)
      xyz2D(i+iOffz,0,2) = QCWQz*xyz2D(i+iOffz,0,1)+B01(i)*xyz2D(i+iOffz,0,0)
    end do
    if (ncdMax > 2) then
      Fact = Two
      do icd=2,ncdMax-1
        do i=1,nArg*lRys
          QCWQx = xyz2D(i,0,1)
          QCWQy = xyz2D(i+iOffy,0,1)
          QCWQz = xyz2D(i+iOffy,0,0)
          temp1x = QCWQx*xyz2D(i,0,icd)
          temp1y = QCWQy*xyz2D(i+iOffy,0,icd)
          temp1z = QCWQz*xyz2D(i+iOffz,0,icd)
          temp2x = Fact*B01(i)*xyz2D(i,0,icd-1)
          temp2y = Fact*B01(i)*xyz2D(i+iOffy,0,icd-1)
          temp2z = Fact*B01(i)*xyz2D(i+iOffz,0,icd-1)
          xyz2D(i,0,icd+1) = temp1x+temp2x
          xyz2D(i+iOffy,0,icd+1) = temp1y+temp2y
          xyz2D(i+iOffz,0,icd+1) = temp1z+temp2z
        end do
        Fact = Fact+One
      end do
    end if
  end if

  ! Compute 2D integrals with index (i,j)

  if (ncdMax <= nabMax) then
    Fac1 = One
    do icd=1,ncdMax
      if (icd == 1) then
        do i=1,nArg*lRys
          PAWPx = xyz2D(i,1,0)
          PAWPy = xyz2D(i+iOffy,1,0)
          PAWPz = xyz2D(i,0,0)
          xyz2D(i,1,1) = PAWPx*xyz2D(i,0,1)+B00(i)
          xyz2D(i+iOffy,1,1) = PAWPy*xyz2D(i+iOffy,0,1)+B00(i)
          xyz2D(i+iOffz,1,1) = PAWPz*xyz2D(i+iOffz,0,1)+B00(i)*xyz2D(i+iOffz,0,0)
        end do
      else
        do i=1,nArg*lRys
          PAWPx = xyz2D(i,1,0)
          PAWPy = xyz2D(i+iOffy,1,0)
          PAWPz = xyz2D(i,0,0)
          xyz2D(i,1,icd) = PAWPx*xyz2D(i,0,icd)+Fac1*B00(i)*xyz2D(i,0,icd-1)
          xyz2D(i+iOffy,1,icd) = PAWPy*xyz2D(i+iOffy,0,icd)+Fac1*B00(i)*xyz2D(i+iOffy,0,icd-1)
          xyz2D(i+iOffz,1,icd) = PAWPz*xyz2D(i+iOffz,0,icd)+Fac1*B00(i)*xyz2D(i+iOffz,0,icd-1)
        end do
      end if
      if (nabMax == 2) then
        do i=1,nArg*lRys
          PAWPx = xyz2D(i,1,0)
          PAWPy = xyz2D(i+iOffy,1,0)
          PAWPz = xyz2D(i,0,0)
          xyz2D(i,2,icd) = PAWPx*xyz2D(i,1,icd)+B10(i)*xyz2D(i,0,icd)+Fac1*B00(i)*xyz2D(i,1,icd-1)
          xyz2D(i+iOffy,2,icd) = PAWPy*xyz2D(i+iOffy,1,icd)+B10(i)*xyz2D(i+iOffy,0,icd)+Fac1*B00(i)*xyz2D(i+iOffy,1,icd-1)
          xyz2D(i+iOffz,2,icd) = PAWPz*xyz2D(i+iOffz,1,icd)+B10(i)*xyz2D(i+iOffz,0,icd)+Fac1*B00(i)*xyz2D(i+iOffz,1,icd-1)
        end do
      else if (nabMax > 2) then
        Fac2 = One
        do iab=1,nabMax-1
          do i=1,nArg*lRys
            PAWPx = xyz2D(i,1,0)
            PAWPy = xyz2D(i+iOffy,1,0)
            PAWPz = xyz2D(i,0,0)
            temp1x = PAWPx*xyz2D(i,iab,icd)
            temp1y = PAWPy*xyz2D(i+iOffy,iab,icd)
            temp1z = PAWPz*xyz2D(i+iOffz,iab,icd)
            temp2x = Fac2*B10(i)*xyz2D(i,iab-1,icd)
            temp2y = Fac2*B10(i)*xyz2D(i+iOffy,iab-1,icd)
            temp2z = Fac2*B10(i)*xyz2D(i+iOffz,iab-1,icd)
            temp3x = Fac1*B00(i)*xyz2D(i,iab,icd-1)
            temp3y = Fac1*B00(i)*xyz2D(i+iOffy,iab,icd-1)
            temp3z = Fac1*B00(i)*xyz2D(i+iOffz,iab,icd-1)
            xyz2D(i,iab+1,icd) = temp1x+temp2x+temp3x
            xyz2D(i+iOffy,iab+1,icd) = temp1y+temp2y+temp3y
            xyz2D(i+iOffz,iab+1,icd) = temp1z+temp2z+temp3z
          end do
          Fac2 = Fac2+One
        end do
      end if
      Fac1 = Fac1+One
    end do
  else
    Fac1 = One
    do iab=1,nabMax
      if (iab == 1) then
        do i=1,nArg*lRys
          QCWQx = xyz2D(i,0,1)
          QCWQy = xyz2D(i+iOffy,0,1)
          QCWQz = xyz2D(i+iOffy,0,0)
          xyz2D(i,1,1) = QCWQx*xyz2D(i,1,0)+B00(i)
          xyz2D(i+iOffy,1,1) = QCWQy*xyz2D(i+iOffy,1,0)+B00(i)
          xyz2D(i+iOffz,1,1) = QCWQz*xyz2D(i+iOffz,1,0)+B00(i)*xyz2D(i+iOffz,0,0)
        end do
      else
        do i=1,nArg*lRys
          QCWQx = xyz2D(i,0,1)
          QCWQy = xyz2D(i+iOffy,0,1)
          QCWQz = xyz2D(i+iOffy,0,0)
          xyz2D(i,iab,1) = QCWQx*xyz2D(i,iab,0)+Fac1*B00(i)*xyz2D(i,iab-1,0)
          xyz2D(i+iOffy,iab,1) = QCWQy*xyz2D(i+iOffy,iab,0)+Fac1*B00(i)*xyz2D(i+iOffy,iab-1,0)
          xyz2D(i+iOffz,iab,1) = QCWQz*xyz2D(i+iOffz,iab,0)+Fac1*B00(i)*xyz2D(i+iOffz,iab-1,0)
        end do
      end if
      if (ncdMax == 2) then
        do i=1,nArg*lRys
          QCWQx = xyz2D(i,0,1)
          QCWQy = xyz2D(i+iOffy,0,1)
          QCWQz = xyz2D(i+iOffy,0,0)
          xyz2D(i,iab,2) = QCWQx*xyz2D(i,iab,1)+B01(i)*xyz2D(i,iab,0)+Fac1*B00(i)*xyz2D(i,iab-1,1)
          xyz2D(i+iOffy,iab,2) = QCWQy*xyz2D(i+iOffy,iab,1)+B01(i)*xyz2D(i+iOffy,iab,0)+Fac1*B00(i)*xyz2D(i+iOffy,iab-1,1)
          xyz2D(i+iOffz,iab,2) = QCWQz*xyz2D(i+iOffz,iab,1)+B01(i)*xyz2D(i+iOffz,iab,0)+Fac1*B00(i)*xyz2D(i+iOffz,iab-1,1)
        end do
      else if (ncdMax > 2) then
        Fac2 = One
        do icd=1,ncdmax-1
          do i=1,nArg*lRys
            QCWQx = xyz2D(i,0,1)
            QCWQy = xyz2D(i+iOffy,0,1)
            QCWQz = xyz2D(i+iOffy,0,0)
            temp1x = QCWQx*xyz2D(i,iab,icd)
            temp1y = QCWQy*xyz2D(i+iOffy,iab,icd)
            temp1z = QCWQz*xyz2D(i+iOffz,iab,icd)
            temp2x = Fac2*B01(i)*xyz2D(i,iab,icd-1)
            temp2y = Fac2*B01(i)*xyz2D(i+iOffy,iab,icd-1)
            temp2z = Fac2*B01(i)*xyz2D(i+iOffz,iab,icd-1)
            temp3x = Fac1*B00(i)*xyz2D(i,iab-1,icd)
            temp3y = Fac1*B00(i)*xyz2D(i+iOffy,iab-1,icd)
            temp3z = Fac1*B00(i)*xyz2D(i+iOffz,iab-1,icd)
            xyz2D(i,iab,icd+1) = temp1x+temp2x+temp3x
            xyz2D(i+iOffy,iab,icd+1) = temp1y+temp2y+temp3y
            xyz2D(i+iOffz,iab,icd+1) = temp1z+temp2z+temp3z
          end do
          Fac2 = Fac2+One
        end do
      end if
      Fac1 = Fac1+One
    end do
  end if
end if

#ifdef _DEBUGPRINT_
do iab=0,nabMax
  do icd=0,ncdMax
    write(Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(x)'
    call RecPrt(Label,' ',xyz2D(1,iab,icd),lRys,nArg)
    write(Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(y)'
    call RecPrt(Label,' ',xyz2D(1+nArg*lRys,iab,icd),lRys,nArg)
    write(Label,'(A,I2,A,I2,A)') ' 2D(',iab,',',icd,')(z)'
    call RecPrt(Label,' ',xyz2D(1+2*nArg*lRys,iab,icd),lRys,nArg)
  end do
end do
#endif

return

end subroutine vRys2D
