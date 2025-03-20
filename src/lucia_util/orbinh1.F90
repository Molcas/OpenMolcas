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
! Copyright (C) 1995,2000, Jeppe Olsen                                 *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ORBINH1(IORBINH1,IORBINH1_NOCCSYM,NTOOBS,NTOOB,NSMOB)
! Obtain array of 2 orbital indices,
! for symmetry packed matrices
!
! IORBINH1 : Lower half packed
! IORBINH1_NOCCSYM : Complete blocks
!
! resulting indices are with respect to start of given symmetry block
! while input orbital indices are absolute and in symmetry order
!
! Jeppe Olsen, March 1995
!              ORBINH1_NOCCSYM added August 2000

use Index_Functions, only: iTri
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NSMOB, NTOOBS(NSMOB), NTOOB
integer(kind=iwp), intent(out) :: IORBINH1(NTOOB,NTOOB), IORBINH1_NOCCSYM(NTOOB,NTOOB)
integer(kind=iwp) :: I_ABS, INDEX_NOCCSYM, INDX, IOFF, IORB, ISM, JABS, JOFF, JORB, JSM

!write(u6,*) ' ORBINH1 speaking'
!write(u6,*) ' NSMOB NTOOB ',NSMOB,NTOOB
!write(u6,*) ' NTOOBS'
!call IWRTMA(NTOOBS,1,NSMOB,1,NSMOB)
! To eliminate annoying and incorrect compiler warnings
IOFF = 0
JOFF = 0
INDX = 0

! Loop over symmetries of orbitals

do ISM=1,NSMOB
  if (ISM == 1) then
    IOFF = 1
  else
    IOFF = IOFF+NTOOBS(ISM-1)
  end if
  do JSM=1,NSMOB
    if (JSM == 1) then
      JOFF = 1
    else
      JOFF = JOFF+NTOOBS(JSM-1)
    end if
    !write(u6,*) ' ISM JSM IOFF JOFF',ISM,JSM,IOFF,JOFF
    do IORB=1,NTOOBS(ISM)
      I_ABS = IOFF-1+IORB
      do JORB=1,NTOOBS(JSM)
        JABS = JOFF-1+JORB
        !write(u6,*) ' IORB JORB I_ABS JABS ',IORB,JORB,I_ABS,JABS
        if (ISM > JSM) then
          INDX = (IORB-1)*NTOOBS(JSM)+JORB
        else if (ISM == JSM) then
          INDX = iTri(IORB,JORB)
        else if (ISM < JSM) then
          INDX = (JORB-1)*NTOOBS(ISM)+IORB
        end if
        INDEX_NOCCSYM = (IORB-1)*NTOOBS(JSM)+JORB
        IORBINH1(I_ABS,JABS) = INDX
        IORBINH1_NOCCSYM(I_ABS,JABS) = INDEX_NOCCSYM
      end do
    end do
    ! End of loops over orbital indices
  end do
end do
! End of loop over orbital symmetries

#ifdef _DEBUGPRINT_
write(u6,*) ' IORBINH1 matrix delivered from ORBINH1'
call IWRTMA(IORBINH1,NTOOB,NTOOB,NTOOB,NTOOB)
#endif

end subroutine ORBINH1
