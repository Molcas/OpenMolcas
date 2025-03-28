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
subroutine ORBINF_MCLR(NIRREP,NSMOB,NRAS1,NRAS2,NRAS3,MXR4tp)
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

use MCLR_Data, only: NORB1, NORB2, NORB3, NORB4, NINOB, NDEOB, NACOB, NOCOB, NTOOB, NRSOBS, ITOOBS, NORB0, IBSO, IBTSOB, IOBPTS, &
                     IOSPIR, IREOST, IREOTS, ISMFSO, ISMFTO, ITPFTO, ITSOB, NACOBS, NDEOBS, NINOBS, NOBPT, NOBPTS, NOCOBS, NOSPIR, &
                     NR0OBS, NR4OBS, NTOOBS, NTSOB
use DetDim, only: MXPIRR, MXPOBS, MXPR4T
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer NIRREP, NSMOB, MXR4tp
integer NRAS1(NIRREP), NRAS2(NIRREP), NRAS3(NIRREP)
! Local variables
integer IPIRR, IRREP, ISM, IISM, IPR4T, ISMOB
#ifdef _DEBUGPRINT_
integer I
#endif

!***********************************************
!                                              *
! Part 1 : From shell format to orbital format *
!                                              *
!***********************************************
do IPIRR=1,MXPIRR
  NOSPIR(IPIRR) = 1
  IOSPIR(1,IPIRR) = IPIRR
end do

! 2 : Shell information to orbital information for each group of orbital

! RAS1
call iCopy(3*MXPOBS,[0],0,NRSOBS,1)
NORB1 = 0
do IRREP=1,NIRREP
  do ISM=1,NOSPIR(IRREP)
    IISM = IOSPIR(ISM,IRREP)
    NRSOBS(IISM,1) = NRSOBS(IISM,1)+NRAS1(IRREP)
    NORB1 = NORB1+NRAS1(IRREP)
  end do
end do
! RAS2

NORB2 = 0
do IRREP=1,NIRREP
  do ISM=1,NOSPIR(IRREP)
    IISM = IOSPIR(ISM,IRREP)
    NRSOBS(IISM,2) = NRSOBS(IISM,2)+NRAS2(IRREP)
    NORB2 = NORB2+NRAS2(IRREP)
  end do
end do

! RAS3
NORB3 = 0
do IRREP=1,NIRREP
  do ISM=1,NOSPIR(IRREP)
    IISM = IOSPIR(ISM,IRREP)
    NRSOBS(IISM,3) = NRSOBS(IISM,3)+NRAS3(IRREP)
    NORB3 = NORB3+NRAS3(IRREP)
  end do
end do
! Inactive, RAS0, RAS4, deleted
NORB4 = 0
NORB0 = 0
NINOB = 0
NDEOB = 0
do IRREP=1,NIRREP
  do ISM=1,NOSPIR(IRREP)
    IISM = IOSPIR(ISM,IRREP)
    NINOBS(IISM) = 0
    NDEOBS(IISM) = 0
    NR0OBS(1,IISM) = 0
    do IPR4T=1,MXPR4T
      NR4OBS(IISM,IPR4T) = 0
    end do
  end do
end do
! Active, occupied  and total number of orbitals
NACOB = NORB1+NORB2+NORB3
NOCOB = NACOB+NINOB+NORB0+NORB4
NTOOB = NOCOB+NDEOB
do ISMOB=1,NSMOB
  NACOBS(ISMOB) = NRSOBS(ISMOB,1)+NRSOBS(ISMOB,2)+NRSOBS(ISMOB,3)
  NOCOBS(ISMOB) = NACOBS(ISMOB)
  NTOOBS(ISMOB) = NACOBS(ISMOB)
end do

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ORBINF_MCLR speaking'
write(u6,*) ' ===================='
write(u6,*) ' Number of orbitals per symmetry + total'
write(u6,'(1X,A,10I4,8X,I3)') '     Ras1             ',(NRSOBS(I,1),I=1,NSMOB),NORB1
write(u6,'(1X,A,10I4,8X,I3)') '     Ras2             ',(NRSOBS(I,2),I=1,NSMOB),NORB2
write(u6,'(1X,A,10I4,8X,I3)') '     Ras3             ',(NRSOBS(I,3),I=1,NSMOB),NORB3
write(u6,'(1X,A,10I4,8X,I3)') '     Active           ',(NACOBS(I),I=1,NSMOB),NACOB
write(u6,'(1X,A,10I4,8X,I3)') '     Total            ',(NTOOBS(I),I=1,NSMOB),NTOOB
#endif
! Offsets for orbitals of given symmetry
ITOOBS(1) = 1
do ISMOB=2,NSMOB
  ITOOBS(ISMOB) = ITOOBS(ISMOB-1)+NTOOBS(ISMOB-1)
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Offsets for orbital of given symmetry'
call IWRTMA(ITOOBS,1,NSMOB,1,NSMOB)
#endif
!*******************************************
!                                          *
! Part 2 : Reordering arrays for orbitals  *
!                                          *
!*******************************************
call ORBORD(NSMOB,MXPOBS,MXR4TP,NDEOBS,NINOBS,NR0OBS,NACOBS,NRSOBS,NR4OBS,NTOOBS,IREOST,IREOTS,ISMFTO,IBSO,NTSOB,IBTSOB,ITSOB, &
            NOBPTS,IOBPTS,MXPR4T,ISMFSO,ITPFTO,NOBPT)

end subroutine ORBINF_MCLR
