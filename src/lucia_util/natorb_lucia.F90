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
! Copyright (C) 1994,1998, Jeppe Olsen                                 *
!               2017, Roland Lindh                                     *
!***********************************************************************

subroutine NATORB_LUCIA(RHO1,NSMOB,NTOPSM,NACPSM,NINPSM,ISTOB,XNAT,RHO1SM,OCCNUM,NACOB,SCR)
! Obtain natural orbitals in symmetry blocks
!
! Jeppe Olsen, June 1994
!              Modification, Oct 94
!              Last modification, Feb. 1998 (reorder deg eigenvalues)

use Index_Functions, only: nTri_Elem
use lucia_data, only: IPRDEN
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NSMOB, NTOPSM(NSMOB), NACPSM(NSMOB), NINPSM(NSMOB), ISTOB(*), NACOB
real(kind=wp), intent(in) :: RHO1(NACOB,NACOB)
real(kind=wp), intent(_OUT_) :: XNAT(*), RHO1SM(*), OCCNUM(*), SCR(*)
integer(kind=iwp) :: I, IMTOFF, IOB, IOBOFF, IOBP, ISMOB, JOB, JOBP, LOB
real(kind=wp) :: SS, TESTY, XI1I, XI1I1, XII, XII1

! To get rid of annoying and incorrect compiler warnings
IOBOFF = 0
IMTOFF = 0
! IOBOFF : Offset for active orbitals in symmetry order
do ISMOB=1,NSMOB
  LOB = NACPSM(ISMOB)
  if (ISMOB == 1) then
    IOBOFF = NINPSM(1)+1
    IMTOFF = 1
  else
    IOBOFF = IOBOFF+NTOPSM(ISMOB-1)-NINPSM(ISMOB-1)+NINPSM(ISMOB)
    IMTOFF = IMTOFF+NACPSM(ISMOB-1)**2
  end if
  LOB = NACPSM(ISMOB)

  ! Extract symmetry block of density matrix

  do IOB=IOBOFF,IOBOFF+LOB-1
    do JOB=IOBOFF,IOBOFF+LOB-1
      ! Corresponding type indices
      IOBP = ISTOB(IOB)
      JOBP = ISTOB(JOB)
      RHO1SM(IMTOFF-1+(JOB-IOBOFF)*LOB+IOB-IOBOFF+1) = RHO1(IOBP,JOBP)
    end do
  end do

  if (IPRDEN >= 2) then
    write(u6,*)
    write(u6,*) ' Density matrix for symmetry  = ',ISMOB
    write(u6,*) ' ======================================='
    write(u6,*)
    call WRTMAT(RHO1SM(IMTOFF),LOB,LOB,LOB,LOB)
  end if
  ! Pack and diagonalize
  !    TRIPAK(AUTPAK,APAK,MATDIM,NDIM)
  call TRIPAK(RHO1SM(IMTOFF),SCR,LOB,LOB)
  ! scale with -1 to get highest occupation numbers as first eigenvectors
  SCR(1:nTri_Elem(LOB)) = -SCR(1:nTri_Elem(LOB))
  call unitmat(XNAT(IMTOFF),LOB)
  call NIDiag(SCR,XNAT(IMTOFF),LOB,LOB)
  call JACORD(SCR,XNAT(IMTOFF),LOB,LOB)

  do I=1,LOB
    OCCNUM(IOBOFF-1+I) = -SCR(nTri_Elem(I))
  end do
  ! Order the degenerate eigenvalues so diagonal terms are maximized
  TESTY = 1.0e-11_wp
  do IOB=2,LOB
    if (abs(OCCNUM(IOBOFF-1+IOB)-OCCNUM(IOBOFF-2+IOB)) <= TESTY) then
      XII = abs(XNAT(IMTOFF-1+(IOB-1)*LOB+IOB))
      XI1I1 = abs(XNAT(IMTOFF-1+(IOB-1-1)*LOB+IOB-1))
      XII1 = abs(XNAT(IMTOFF-1+(IOB-1-1)*LOB+IOB))
      XI1I = abs(XNAT(IMTOFF-1+(IOB-1)*LOB+IOB-1))

      if ((XI1I > XII) .and. (XII1 > XI1I1)) then
        ! interchange orbital IOB and IOB -1
        call dSwap_(LOB,XNAT(IMTOFF+(IOB-1)*LOB),1,XNAT(IMTOFF+(IOB-1-1)*LOB),1)
        SS = OCCNUM(IOBOFF-1+IOB-1)
        OCCNUM(IOBOFF-1+IOB-1) = OCCNUM(IOBOFF-1+IOB)
        OCCNUM(IOBOFF-1+IOB) = SS
        if (IPRDEN >= 1) write(u6,*) ' Orbitals interchanged ',IOBOFF-1+IOB,IOBOFF-2+IOB
      end if
    end if
  end do

  if (IPRDEN >= 1) then
    write(u6,*)
    write(u6,*) ' Natural occupation numbers for symmetry = ',ISMOB
    write(u6,*) ' ==================================================='
    write(u6,*)
    call WRTMAT(OCCNUM(IOBOFF),1,LOB,1,LOB)
    if (IPRDEN >= 2) then
      write(u6,*)
      write(u6,*) ' Corresponding Eigenvectors'
      write(u6,*)
      call WRTMAT(XNAT(IMTOFF),LOB,LOB,LOB,LOB)
    end if
  end if
end do
! (End of loop over orbital symmetries)

end subroutine NATORB_LUCIA
