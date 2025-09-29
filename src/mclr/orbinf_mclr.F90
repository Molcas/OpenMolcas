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
! Copyright (C) 1991, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ORBINF_MCLR(NIRREP,NRAS1,NRAS2,NRAS3,MXR4tp)
! Obtain information about orbitals from shell information
!
! =====
! input
! =====
! Shell and symmetry information in /LUCINP/
!
! ======
! Output
! ======
! Orbital information in /ORBINP/

use MCLR_Data, only: IBSO, IBTSOB, IREOTS, ISMFTO, ITSOB, NACOB, NOBPT, NOBPTS, NOCOB, NORB1, NORB2, NORB3, NTOOB, NTSOB
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NIRREP, NRAS1(NIRREP), NRAS2(NIRREP), NRAS3(NIRREP), MXR4tp
integer(kind=iwp) :: NNOBPT
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I
#endif
integer(kind=iwp), allocatable :: NRSOBS(:,:), NTOOBS(:)

!***********************************************
!                                              *
! Part 1 : From shell format to orbital format *
!                                              *
!***********************************************

! 2 : Shell information to orbital information for each group of orbital

call mma_allocate(NTOOBS,NIRREP,Label='NTOOBS')
call mma_allocate(NRSOBS,NIRREP,3,Label='NRSOBS')

! RAS1, RAS2, RAS3
NRSOBS(:,1) = NRAS1(1:NIRREP)
NRSOBS(:,2) = NRAS2(1:NIRREP)
NRSOBS(:,3) = NRAS3(1:NIRREP)
NORB1 = sum(NRSOBS(:,1))
NORB2 = sum(NRSOBS(:,2))
NORB3 = sum(NRSOBS(:,3))

! Active, occupied  and total number of orbitals
NACOB = NORB1+NORB2+NORB3
NOCOB = NACOB
NTOOBS(:) = NRSOBS(:,1)+NRSOBS(:,2)+NRSOBS(:,3)
NTOOB = sum(NTOOBS)

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ORBINF_MCLR speaking'
write(u6,*) ' ===================='
write(u6,*) ' Number of orbitals per symmetry + total'
write(u6,'(1X,A,10I4,8X,I3)') '     Ras1             ',(NRSOBS(I,1),I=1,NIRREP),NORB1
write(u6,'(1X,A,10I4,8X,I3)') '     Ras2             ',(NRSOBS(I,2),I=1,NIRREP),NORB2
write(u6,'(1X,A,10I4,8X,I3)') '     Ras3             ',(NRSOBS(I,3),I=1,NIRREP),NORB3
write(u6,'(1X,A,10I4,8X,I3)') '     Total            ',(NTOOBS(I),I=1,NIRREP),NTOOB
#endif

!*******************************************
!                                          *
! Part 2 : Reordering arrays for orbitals  *
!                                          *
!*******************************************
NNOBPT = size(NOBPT)
call ORBORD(NIRREP,MXR4TP,NRSOBS,NTOOBS,IREOTS,ISMFTO,IBSO,NTSOB,IBTSOB,ITSOB,NOBPTS,NNOBPT,NOBPT)

call mma_deallocate(NTOOBS)
call mma_deallocate(NRSOBS)

end subroutine ORBINF_MCLR
