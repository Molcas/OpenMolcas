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

subroutine Rys2D(xyz2D,nArg,lRys,nabMax,ncdMax,PAWP,QCWQ,B10,B00,B01)
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
! VV: improve loop structure                                           *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
real*8 xyz2D(nArg*lRys*3,0:nabMax,0:ncdMax), PAWP(nArg*lRys*3), QCWQ(nArg*lRys*3), B10(nArg*lRys*3), B00(nArg*lRys*3), &
       B01(nArg*lRys*3)

#ifdef _DEBUGPRINT_
character*30 Label
if (nabMax > 0) call RecPrt('PAWP',' ',PAWP,nArg,lRys*3)
if (ncdMax > 0) call RecPrt('QCWQ',' ',QCWQ,nArg,lRys*3)
call RecPrt(' B10',' ',B10,nArg*lRys,3)
call RecPrt(' B00',' ',B00,nArg*lRys,3)
call RecPrt(' B01',' ',B01,nArg*lRys,3)
#endif

! Compute 2D integrals with index (0,0). Observe that the z
! component already contains the weight factor.

call dcopy_(2*nArg*lRys,[One],0,xyz2D(1,0,0),1)

! Compute 2D integrals with index (i,0)

if (nabMax /= 0) then
  do i=1,nArg*lRys*3
    xyz2D(i,1,0) = PAWP(i)*xyz2D(i,0,0)
  end do
end if
if (nabMax == 2) then
  do i=1,nArg*lRys*3
    xyz2D(i,2,0) = PAWP(i)*xyz2D(i,1,0)+B10(i)*xyz2D(i,0,0)
  end do
else if (nabMax > 2) then
  do iab=1,nabMax-1
    do i=1,nArg*lRys*3
      temp1 = PAWP(i)*xyz2D(i,iab,0)
      temp2 = dble(iab)*B10(i)*xyz2D(i,iab-1,0)
      xyz2D(i,iab+1,0) = temp1+temp2
    end do
  end do
end if

! Compute 2D integrals with index (0,i)

if (ncdMax /= 0) then
  do i=1,nArg*lRys*3
    xyz2D(i,0,1) = QCWQ(i)*xyz2D(i,0,0)
  end do
end if
if (ncdMax == 2) then
  do i=1,nArg*lRys*3
    xyz2D(i,0,2) = QCWQ(i)*xyz2D(i,0,1)+B01(i)*xyz2D(i,0,0)
  end do
else if (ncdMax > 2) then
  do icd=1,ncdMax-1
    do i=1,nArg*lRys*3
      temp1 = QCWQ(i)*xyz2D(i,0,icd)
      temp2 = dble(icd)*B01(i)*xyz2D(i,0,icd-1)
      xyz2D(i,0,icd+1) = temp1+temp2
    end do
  end do
end if

! Compute 2D integrals with index (i,j)

if (ncdMax <= nabMax) then
  do icd=1,ncdMax
    do i=1,nArg*lRys*3
      xyz2D(i,1,icd) = PAWP(i)*xyz2D(i,0,icd)+dble(icd)*B00(i)*xyz2D(i,0,icd-1)
    end do
    if (nabMax == 2) then
      do i=1,nArg*lRys*3
        xyz2D(i,2,icd) = PAWP(i)*xyz2D(i,1,icd)+B10(i)*xyz2D(i,0,icd)+dble(icd)*B00(i)*xyz2D(i,1,icd-1)
      end do
    else if (nabMax > 2) then
      do iab=1,nabMax-1
        do i=1,nArg*lRys*3
          temp1 = PAWP(i)*xyz2D(i,iab,icd)
          temp2 = dble(iab)*B10(i)*xyz2D(i,iab-1,icd)
          temp3 = dble(icd)*B00(i)*xyz2D(i,iab,icd-1)
          xyz2D(i,iab+1,icd) = temp1+temp2+temp3
        end do
      end do
    end if
  end do
else
  do iab=1,nabMax
    do i=1,nArg*lRys*3
      xyz2D(i,iab,1) = QCWQ(i)*xyz2D(i,iab,0)+dble(iab)*B00(i)*xyz2D(i,iab-1,0)
    end do
    if (ncdMax == 2) then
      do i=1,nArg*lRys*3
        xyz2D(i,iab,2) = QCWQ(i)*xyz2D(i,iab,1)+B01(i)*xyz2D(i,iab,0)+dble(iab)*B00(i)*xyz2D(i,iab-1,1)
      end do
    else if (ncdMax > 2) then
      do icd=1,ncdmax-1
        do i=1,nArg*lRys*3
          temp1 = QCWQ(i)*xyz2D(i,iab,icd)
          temp2 = dble(icd)*B01(i)*xyz2D(i,iab,icd-1)
          temp3 = dble(iab)*B00(i)*xyz2D(i,iab-1,icd)
          xyz2D(i,iab,icd+1) = temp1+temp2+temp3
        end do
      end do
    end if
  end do

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

end subroutine Rys2D
