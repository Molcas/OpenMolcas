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

use Constants, only: One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nArg, lRys, nabMax, ncdMax
real(kind=wp), intent(inout) :: xyz2D(nArg*lRys,3,0:nabMax,0:ncdMax)
real(kind=wp), intent(in) :: PAWP(nArg*lRys,3), QCWQ(nArg*lRys,3), B10(nArg*lRys), B00(nArg*lRys), B01(nArg*lRys)
integer(kind=iwp) :: iab, icd
real(kind=wp) :: Fac1, Fac2, Fact
logical(kind=iwp) :: lPAWP, lQCWQ
#ifdef _DEBUGPRINT_
character(len=30) :: Label
#endif

#ifdef _DEBUGPRINT_
if (nabMax > 0) call RecPrt('vRys2D: PAWP',' ',PAWP,lRys,nArg*3)
if (ncdMax > 0) call RecPrt('vRys2D: QCWQ',' ',QCWQ,lRys,nArg*3)
call RecPrt('vRys2D:  B10',' ',B10,lRys,nArg)
call RecPrt('vRys2D:  B00',' ',B00,lRys,nArg)
call RecPrt('vRys2D:  B01',' ',B01,lRys,nArg)
#endif

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
  if (lPAWP) xyz2D(:,1,0,0) = PAWP(:,3)
  if (lQCWQ) xyz2D(:,2,0,0) = QCWQ(:,3)

  ! Compute 2D integrals with index (i,0)

  if (nabMax == 0) then

  else if (nabMax == 1) then
    xyz2D(:,3,1,0) = xyz2D(:,3,1,0)*xyz2D(:,3,0,0)
  else if (nabMax >= 2) then
    xyz2D(:,1,2,0) = xyz2D(:,1,1,0)**2+B10(:)
    xyz2D(:,2,2,0) = xyz2D(:,2,1,0)**2+B10(:)
    xyz2D(:,3,2,0) = (xyz2D(:,3,1,0)**2+B10(:))*xyz2D(:,3,0,0)
    xyz2D(:,3,1,0) = xyz2D(:,3,1,0)*xyz2D(:,3,0,0)
    if (nabMax > 2) then
      Fact = Two
      do iab=2,nabMax-1
        xyz2D(:,1,iab+1,0) = xyz2D(:,1,1,0)*xyz2D(:,1,iab,0)+Fact*B10(:)*xyz2D(:,1,iab-1,0)
        xyz2D(:,2,iab+1,0) = xyz2D(:,2,1,0)*xyz2D(:,2,iab,0)+Fact*B10(:)*xyz2D(:,2,iab-1,0)
        xyz2D(:,3,iab+1,0) = xyz2D(:,1,0,0)*xyz2D(:,3,iab,0)+Fact*B10(:)*xyz2D(:,3,iab-1,0)
        Fact = Fact+One
      end do
    end if
  end if

  ! Compute 2D integrals with index (0,i)

  if (ncdMax == 1) then
    xyz2D(:,3,0,1) = xyz2D(:,3,0,1)*xyz2D(:,3,0,0)
  else if (ncdMax >= 2) then
    xyz2D(:,1,0,2) = xyz2D(:,1,0,1)**2+B01(:)
    xyz2D(:,2,0,2) = xyz2D(:,2,0,1)**2+B01(:)
    xyz2D(:,3,0,2) = (xyz2D(:,3,0,1)**2+B01(:))*xyz2D(:,3,0,0)
    xyz2D(:,3,0,1) = xyz2D(:,3,0,1)*xyz2D(:,3,0,0)
    if (ncdMax > 2) then
      Fact = Two
      do icd=2,ncdMax-1
        xyz2D(:,1,0,icd+1) = xyz2D(:,1,0,1)*xyz2D(:,1,0,icd)+Fact*B01(:)*xyz2D(:,1,0,icd-1)
        xyz2D(:,2,0,icd+1) = xyz2D(:,2,0,1)*xyz2D(:,2,0,icd)+Fact*B01(:)*xyz2D(:,2,0,icd-1)
        xyz2D(:,3,0,icd+1) = xyz2D(:,2,0,0)*xyz2D(:,3,0,icd)+Fact*B01(:)*xyz2D(:,3,0,icd-1)
        Fact = Fact+One
      end do
    end if
  end if

  ! Compute 2D integrals with index (i,j)

  if (ncdMax <= nabMax) then
    Fac1 = One
    do icd=1,ncdMax
      if (icd == 1) then
        xyz2D(:,1,1,1) = xyz2D(:,1,1,0)*xyz2D(:,1,0,1)+B00(:)
        xyz2D(:,2,1,1) = xyz2D(:,2,1,0)*xyz2D(:,2,0,1)+B00(:)
        xyz2D(:,3,1,1) = xyz2D(:,1,0,0)*xyz2D(:,3,0,1)+B00(:)*xyz2D(:,3,0,0)
      else
        xyz2D(:,1,1,icd) = xyz2D(:,1,1,0)*xyz2D(:,1,0,icd)+Fac1*B00(:)*xyz2D(:,1,0,icd-1)
        xyz2D(:,2,1,icd) = xyz2D(:,2,1,0)*xyz2D(:,2,0,icd)+Fac1*B00(:)*xyz2D(:,2,0,icd-1)
        xyz2D(:,3,1,icd) = xyz2D(:,1,0,0)*xyz2D(:,3,0,icd)+Fac1*B00(:)*xyz2D(:,3,0,icd-1)
      end if
      if (nabMax == 2) then
        xyz2D(:,1,2,icd) = xyz2D(:,1,1,0)*xyz2D(:,1,1,icd)+B10(:)*xyz2D(:,1,0,icd)+Fac1*B00(:)*xyz2D(:,1,1,icd-1)
        xyz2D(:,2,2,icd) = xyz2D(:,2,1,0)*xyz2D(:,2,1,icd)+B10(:)*xyz2D(:,2,0,icd)+Fac1*B00(:)*xyz2D(:,2,1,icd-1)
        xyz2D(:,3,2,icd) = xyz2D(:,1,0,0)*xyz2D(:,3,1,icd)+B10(:)*xyz2D(:,3,0,icd)+Fac1*B00(:)*xyz2D(:,3,1,icd-1)
      else if (nabMax > 2) then
        Fac2 = One
        do iab=1,nabMax-1
          xyz2D(:,1,iab+1,icd) = xyz2D(:,1,1,0)*xyz2D(:,1,iab,icd)+Fac2*B10(:)*xyz2D(:,1,iab-1,icd)+Fac1*B00(:)*xyz2D(:,1,iab,icd-1)
          xyz2D(:,2,iab+1,icd) = xyz2D(:,2,1,0)*xyz2D(:,2,iab,icd)+Fac2*B10(:)*xyz2D(:,2,iab-1,icd)+Fac1*B00(:)*xyz2D(:,2,iab,icd-1)
          xyz2D(:,3,iab+1,icd) = xyz2D(:,1,0,0)*xyz2D(:,3,iab,icd)+Fac2*B10(:)*xyz2D(:,3,iab-1,icd)+Fac1*B00(:)*xyz2D(:,3,iab,icd-1)
          Fac2 = Fac2+One
        end do
      end if
      Fac1 = Fac1+One
    end do
  else
    Fac1 = One
    do iab=1,nabMax
      if (iab == 1) then
        xyz2D(:,1,1,1) = xyz2D(:,1,0,1)*xyz2D(:,1,1,0)+B00(:)
        xyz2D(:,2,1,1) = xyz2D(:,2,0,1)*xyz2D(:,2,1,0)+B00(:)
        xyz2D(:,3,1,1) = xyz2D(:,2,0,0)*xyz2D(:,3,1,0)+B00(:)*xyz2D(:,3,0,0)
      else
        xyz2D(:,1,iab,1) = xyz2D(:,1,0,1)*xyz2D(:,1,iab,0)+Fac1*B00(:)*xyz2D(:,1,iab-1,0)
        xyz2D(:,2,iab,1) = xyz2D(:,2,0,1)*xyz2D(:,2,iab,0)+Fac1*B00(:)*xyz2D(:,2,iab-1,0)
        xyz2D(:,3,iab,1) = xyz2D(:,2,0,0)*xyz2D(:,3,iab,0)+Fac1*B00(:)*xyz2D(:,3,iab-1,0)
      end if
      if (ncdMax == 2) then
        xyz2D(:,1,iab,2) = xyz2D(:,1,0,1)*xyz2D(:,1,iab,1)+B01(:)*xyz2D(:,1,iab,0)+Fac1*B00(:)*xyz2D(:,1,iab-1,1)
        xyz2D(:,2,iab,2) = xyz2D(:,2,0,1)*xyz2D(:,2,iab,1)+B01(:)*xyz2D(:,2,iab,0)+Fac1*B00(:)*xyz2D(:,2,iab-1,1)
        xyz2D(:,3,iab,2) = xyz2D(:,2,0,0)*xyz2D(:,3,iab,1)+B01(:)*xyz2D(:,3,iab,0)+Fac1*B00(:)*xyz2D(:,3,iab-1,1)
      else if (ncdMax > 2) then
        Fac2 = One
        do icd=1,ncdmax-1
          xyz2D(:,1,iab,icd+1) = xyz2D(:,1,0,1)*xyz2D(:,1,iab,icd)+Fac2*B01(:)*xyz2D(:,1,iab,icd-1)+Fac1*B00(:)*xyz2D(:,1,iab-1,icd)
          xyz2D(:,2,iab,icd+1) = xyz2D(:,2,0,1)*xyz2D(:,2,iab,icd)+Fac2*B01(:)*xyz2D(:,2,iab,icd-1)+Fac1*B00(:)*xyz2D(:,2,iab-1,icd)
          xyz2D(:,3,iab,icd+1) = xyz2D(:,2,0,0)*xyz2D(:,3,iab,icd)+Fac2*B01(:)*xyz2D(:,3,iab,icd-1)+Fac1*B00(:)*xyz2D(:,3,iab-1,icd)
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
    write(Label,'(A,I2,A,I2,A)') 'vRys2D:  2D(',iab,',',icd,')(x)'
    call RecPrt(Label,' ',xyz2D(:,1,iab,icd),lRys,nArg)
    write(Label,'(A,I2,A,I2,A)') 'vRys2D:  2D(',iab,',',icd,')(y)'
    call RecPrt(Label,' ',xyz2D(:,2,iab,icd),lRys,nArg)
    write(Label,'(A,I2,A,I2,A)') 'vRys2D:  2D(',iab,',',icd,')(z)'
    call RecPrt(Label,' ',xyz2D(:,3,iab,icd),lRys,nArg)
  end do
end do
#endif

xyz2D(:,1:2,0,0) = One

return

end subroutine vRys2D
