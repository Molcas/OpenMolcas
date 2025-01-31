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

implicit real*8(A-H,O-Z)
! Input
dimension NTOOBS(NSMOB)
! output
dimension IORBINH1(NTOOB,NTOOB), IORBINH1_NOCCSYM(NTOOB,NTOOB)

!write(6,*) ' ORBINH1 speaking'
!write(6,*) ' NSMOB NTOOB ',NSMOB,NTOOB
!write(6,*) ' NTOOBS'
!call IWRTMA(NTOOBS,1,NSMOB,1,NSMOB)
! To eliminate annoying and incorrect compiler warnings
IOFF = 0
JOFF = 0
INDEX = 0

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
    !write(6,*) ' ISM JSM IOFF JOFF',ISM,JSM,IOFF,JOFF
    do IORB=1,NTOOBS(ISM)
      IABS = IOFF-1+IORB
      do JORB=1,NTOOBS(JSM)
        JABS = JOFF-1+JORB
        !write(6,*) ' IORB JORB IABS JABS ',IORB,JORB,IABS,JABS
        if (ISM > JSM) then
          INDEX = (IORB-1)*NTOOBS(JSM)+JORB
        else if (ISM == JSM) then
          if (IORB >= JORB) then
            INDEX = IORB*(IORB-1)/2+JORB
          else
            INDEX = JORB*(JORB-1)/2+IORB
          end if
        else if (ISM < JSM) then
          INDEX = (JORB-1)*NTOOBS(ISM)+IORB
        end if
        INDEX_NOCCSYM = (IORB-1)*NTOOBS(JSM)+JORB
        IORBINH1(IABS,JABS) = INDEX
        IORBINH1_NOCCSYM(IABS,JABS) = INDEX_NOCCSYM
      end do
    end do
    ! End of loops over orbital indices
  end do
end do
! End of loop over orbital symmetries

NTEST = 0
if (NTEST >= 100) then
  write(6,*) ' IORBINH1 matrix delivered from ORBINH1'
  call IWRTMA(IORBINH1,NTOOB,NTOOB,NTOOB,NTOOB)
end if

end subroutine ORBINH1
