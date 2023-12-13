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

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nArg, lRys, nabMax, ncdMax
real(kind=wp), intent(inout) :: xyz2D(nArg*lRys,3,0:nabMax,0:ncdMax)
real(kind=wp), intent(in) :: PAWP(nArg*lRys,3), QCWQ(nArg*lRys,3), B10(nArg*lRys,3), B00(nArg*lRys,3), B01(nArg*lRys,3)
integer(kind=iwp) :: iab, icd
#ifdef _DEBUGPRINT_
character(len=30) :: Label
#endif

#ifdef _DEBUGPRINT_
if (nabMax > 0) call RecPrt('Rys2D: PAWP',' ',PAWP,nArg,lRys*3)
if (ncdMax > 0) call RecPrt('Rys2D: QCWQ',' ',QCWQ,nArg,lRys*3)
call RecPrt('Rys2D:  B10',' ',B10,nArg*lRys,3)
call RecPrt('Rys2D:  B00',' ',B00,nArg*lRys,3)
call RecPrt('Rys2D:  B01',' ',B01,nArg*lRys,3)
#endif

! Compute 2D integrals with index (0,0). Observe that the z
! component already contains the weight factor.

xyz2D(:,1:2,0,0) = One

! Compute 2D integrals with index (i,0)

if (nabMax /= 0) then
  xyz2D(:,:,1,0) = PAWP(:,:)*xyz2D(:,:,0,0)
  if (nabMax == 2) then
    xyz2D(:,:,2,0) = PAWP(:,:)*xyz2D(:,:,1,0)+B10(:,:)*xyz2D(:,:,0,0)
  else if (nabMax > 2) then
    do iab=1,nabMax-1
      xyz2D(:,:,iab+1,0) = PAWP(:,:)*xyz2D(:,:,iab,0)+real(iab,kind=wp)*B10(:,:)*xyz2D(:,:,iab-1,0)
    end do
  end if
end if

! Compute 2D integrals with index (0,i)

if (ncdMax /= 0) then
  xyz2D(:,:,0,1) = QCWQ(:,:)*xyz2D(:,:,0,0)
  if (ncdMax == 2) then
    xyz2D(:,:,0,2) = QCWQ(:,:)*xyz2D(:,:,0,1)+B01(:,:)*xyz2D(:,:,0,0)
  else if (ncdMax > 2) then
    do icd=1,ncdMax-1
      xyz2D(:,:,0,icd+1) = QCWQ(:,:)*xyz2D(:,:,0,icd)+real(icd,kind=wp)*B01(:,:)*xyz2D(:,:,0,icd-1)
    end do
  end if
end if

! Compute 2D integrals with index (i,j)

if (ncdMax <= nabMax) then
  do icd=1,ncdMax
    xyz2D(:,:,1,icd) = PAWP(:,:)*xyz2D(:,:,0,icd)+real(icd,kind=wp)*B00(:,:)*xyz2D(:,:,0,icd-1)
    if (nabMax == 2) then
      xyz2D(:,:,2,icd) = PAWP(:,:)*xyz2D(:,:,1,icd)+B10(:,:)*xyz2D(:,:,0,icd)+real(icd,kind=wp)*B00(:,:)*xyz2D(:,:,1,icd-1)
    else if (nabMax > 2) then
      do iab=1,nabMax-1
        xyz2D(:,:,iab+1,icd) = PAWP(:,:)*xyz2D(:,:,iab,icd)+real(iab,kind=wp)*B10(:,:)*xyz2D(:,:,iab-1,icd)+ &
                               real(icd,kind=wp)*B00(:,:)*xyz2D(:,:,iab,icd-1)
      end do
    end if
  end do
else
  do iab=1,nabMax
    xyz2D(:,:,iab,1) = QCWQ(:,:)*xyz2D(:,:,iab,0)+real(iab,kind=wp)*B00(:,:)*xyz2D(:,:,iab-1,0)
    if (ncdMax == 2) then
      xyz2D(:,:,iab,2) = QCWQ(:,:)*xyz2D(:,:,iab,1)+B01(:,:)*xyz2D(:,:,iab,0)+real(iab,kind=wp)*B00(:,:)*xyz2D(:,:,iab-1,1)
    else if (ncdMax > 2) then
      do icd=1,ncdmax-1
        xyz2D(:,:,iab,icd+1) = QCWQ(:,:)*xyz2D(:,:,iab,icd)+real(icd,kind=wp)*B01(:,:)*xyz2D(:,:,iab,icd-1)+ &
                               real(iab,kind=wp)*B00(:,:)*xyz2D(:,:,iab-1,icd)
      end do
    end if
  end do

end if

#ifdef _DEBUGPRINT_
do iab=0,nabMax
  do icd=0,ncdMax
    write(Label,'(A,I2,A,I2,A)') 'Rys2D:  2D(',iab,',',icd,')(x)'
    call RecPrt(Label,' ',xyz2D(:,1,iab,icd),lRys,nArg)
    write(Label,'(A,I2,A,I2,A)') 'Rys2D:  2D(',iab,',',icd,')(y)'
    call RecPrt(Label,' ',xyz2D(:,2,iab,icd),lRys,nArg)
    write(Label,'(A,I2,A,I2,A)') 'Rys2D:  2D(',iab,',',icd,')(z)'
    call RecPrt(Label,' ',xyz2D(:,3,iab,icd),lRys,nArg)
  end do
end do
#endif

return

end subroutine Rys2D
