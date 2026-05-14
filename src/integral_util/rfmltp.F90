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
! Copyright (C) 1992, Roland Lindh                                     *
!               1994, Markus P. Fuelscher                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine RFmltp()
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!                                                                      *
!     modified by M. P. Fuelscher, 94/04/28                            *
!***********************************************************************

use Index_Functions, only: nTri3_Elem1, nTri_Elem1
use rctfld_module, only: Eps, EpsInf, lMax, lRF, lRFCav, MM, NonEq_ref, PCM, rds
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, iElem, iM, iOff, jM, l, nComp, nElem, nM
real(kind=wp) :: dESolv, ESolv
real(kind=wp), allocatable :: QTot(:), VTot(:)
real(kind=wp), external :: DDot_

if (.not. lRF) return
nComp = nTri3_Elem1(lMax)
call mma_allocate(VTot,nComp,Label='VTot')
call mma_allocate(QTot,nComp,Label='QTot')

if (lRF .and. (.not. PCM) .and. lRFCav) then
  QTot(:) = MM(:,1)+MM(:,2)
# ifdef _DEBUGPRINT_
  call RecPrt('Total Multipole Moments',' ',QTot,1,nComp)
# endif
  VTot(:) = QTot(:)
  ! Compute the electric field due to the total charge distribution.
  call AppFld(VTot,rds,Eps,lMax,EpsInf,NonEq_ref)
# ifdef _DEBUGPRINT_
  call RecPrt('Total Electric Field',' ',VTot,1,nComp)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  write(u6,*) '     Multipole analysis of the contributions to the dielectric solvation energy'
  write(u6,*)
  write(u6,*) '     --------------------------------------'
  write(u6,*) '        l             dE '
  write(u6,*) '     --------------------------------------'
  Esolv = Zero
  iOff = 1
  do l=0,lMax
    nElem = nTri_Elem1(l)
    dEsolv = -Half*DDot_(nElem,QTot(iOff),1,VTot(iOff),1)
    write(u6,'(8X,I2,10X,F13.10)') l,dEsolv
    iOff = iOff+nElem
    Esolv = Esolv+dEsolv
  end do
  write(u6,*) '     --------------------------------------'
  write(u6,*)
  write(u6,*)
  write(u6,*) '     Total Multipole Moments (cartesian)'
  write(u6,*) '     -----------------------------------'
  iM = 1
  do l=0,lMax
    nElem = nTri_Elem1(l)
    jM = iM
    do iElem=1,nElem,7
      nM = min(7,nElem-iElem+1)
      write(u6,'(8X,7ES14.5)') (QTot(i),i=jM,jM+nM-1)
      jM = jM+nM
    end do
    iM = iM+nElem
  end do
  write(u6,*) '     -----------------------------------'
  write(u6,*)
  write(u6,*)
  write(u6,*) '     Total Electric Field (cartesian)'
  write(u6,*) '     --------------------------------'
  iM = 1
  do l=0,lMax
    nElem = nTri_Elem1(l)
    jM = iM
    do iElem=1,nElem,7
      nM = min(7,nElem-iElem+1)
      write(u6,'(8X,7ES14.5)') (VTot(i),i=jM,jM+nM-1)
      jM = jM+nM
    end do
    iM = iM+nElem
  end do
  write(u6,*) '     -----------------------------------'
  write(u6,*)
end if

call mma_deallocate(VTot)
call mma_deallocate(QTot)

end subroutine RFmltp
