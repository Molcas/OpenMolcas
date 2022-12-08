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
! Copyright (C) 1990,1991,1995, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

subroutine XRys2D(xyz2D,nArg,lRys,nabMax,ncdMax,PAWP,QCWQ,B10,B00,B01)
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
! Modified for external field version, Feb '95.                        *
! VV: improve loop structure                                           *
!***********************************************************************

use Constants, only: One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nArg, lRys, nabMax, ncdMax
real(kind=wp), intent(inout) :: xyz2D(nArg*lRys,3,0:nabMax,0:ncdMax)
real(kind=wp), intent(in) :: PAWP(nArg*lRys,3), QCWQ(nArg*lRys,3), B10(nArg*lRys,3), B00(nArg*lRys,3), B01(nArg*lRys,3)
integer(kind=iwp) :: iab, in_
#ifdef _DEBUGPRINT_
character(len=30) :: Label
#endif

#ifdef _DEBUGPRINT_
iRout = 15
iPrint = nPrint(iRout)
if (iPrint >= 59) then
  if (nabMax > 0) call RecPrt('PAWP',' ',PAWP,nArg,lRys*3)
  if (ncdMax > 0) call RecPrt('QCWQ',' ',QCWQ,nArg,lRys*3)
  call RecPrt(' B10',' ',B10,nArg*lRys,3)
  call RecPrt(' B00',' ',B00,nArg*lRys,3)
  call RecPrt(' B01',' ',B01,nArg*lRys,3)
end if
#endif

! Compute 2D integrals with index (0,0). Observe that the z
! component already contains the weight factor.

xyz2D(:,1:2,0,0) = One

! Span first I(i,0)

if (nabMax >= 1) then
  xyz2D(:,:,1,0) = PAWP(:,:)*xyz2D(:,:,0,0)
  do iab=1,nabMax-1
    xyz2D(:,:,iab+1,0) = PAWP(:,:)*xyz2D(:,:,iab,0)+real(iab,kind=wp)*B10(:,:)*xyz2D(:,:,iab-1,0)
  end do
end if

! Now do the rest!

if (ncdMax >= 1) then
  xyz2D(:,:,0,1) = QCWQ(:,:)*xyz2D(:,:,0,0)
  do iab=1,nabMax
    xyz2D(:,:,iab,1) = QCWQ(:,:)*xyz2D(:,:,iab,0)+real(iab,kind=wp)*B00(:,:)*xyz2D(:,:,iab-1,0)
  end do
end if
do in_=1,ncdMax-1
  xyz2D(:,:,0,in_+1) = QCWQ(:,:)*xyz2D(:,:,0,in_)-real(in_,kind=wp)*B01(:,:)*xyz2D(:,:,0,in_-1)
  do iab=1,nabMax
    xyz2D(:,:,iab,in_+1) = QCWQ(:,:)*xyz2D(:,:,iab,in_)+real(iab,kind=wp)*B00(:,:)*xyz2D(:,:,iab-1,in_)- &
                           real(in_,kind=wp)*B01(:,:)*xyz2D(:,:,iab,in_-1)
  end do
end do

#ifdef _DEBUGPRINT_
if (iPrint >= 99) then
  write(u6,*) ' 2D-integral computed in XRys2D'
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
end if
#endif

return

end subroutine XRys2D
