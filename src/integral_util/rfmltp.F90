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

use rctfld_module, only: Eps, EpsInf, lMax, lRF, lRFCav, MM, NonEq_ref, PCM, rds
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: u6

implicit none
real*8, allocatable :: VTot(:), QTot(:)
integer nComp

if (.not. lRF) return
nComp = (lMax+1)*(lMax+2)*(lMax+3)/6
call mma_allocate(VTot,nComp,Label='VTot')
call mma_allocate(QTot,nComp,Label='QTot')

call RFmltp_Internal(MM,nComp)

call mma_deallocate(VTot)
call mma_deallocate(QTot)

contains

subroutine RFmltp_Internal(Qs,nComp)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!                                                                      *
!     modified by M. P. Fuelscher, 94/04/28                            *
!***********************************************************************

  use Constants, only: Zero, One, Half

  implicit none
  integer nComp
  real*8 Qs(nComp,2)
  integer l, nElem, iM, jM, iElem, iOff, i, nM
  real*8 ESolv, dESolv
  real*8, external :: DDot_

  if (lRF .and. (.not. PCM) .and. lRFCav) then
    call dcopy_(nComp,Qs(1,1),1,QTot,1)
    call DaXpY_(nComp,One,Qs(1,2),1,QTot,1)
#   ifdef _DEBUGPRINT_
    call RecPrt('Total Multipole Moments',' ',QTot,1,nComp)
#   endif
    call dcopy_(nComp,QTot,1,VTot,1)
    ! Compute the electric field due to the total charge distribution.
    call AppFld(VTot,rds,Eps,lMax,EpsInf,NonEq_ref)
#   ifdef _DEBUGPRINT_
    call RecPrt('Total Electric Field',' ',VTot,1,nComp)
#   endif
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    write(u6,*) '     Multipole analysis of the contributions to the dielectric solvation energy'
    write(u6,*)
    write(u6,*) '     --------------------------------------'
    write(u6,*) '        l             dE '
    write(u6,*) '     --------------------------------------'
    Esolv = Zero
    iOff = 1
    do l=0,lMax
      nElem = (l+1)*(l+2)/2
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
      nElem = (l+1)*(l+2)/2
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
      nElem = (l+1)*(l+2)/2
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

end subroutine RFmltp_Internal

end subroutine RFmltp
