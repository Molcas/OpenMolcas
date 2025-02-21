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
! Copyright (C) 1994, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ORBORD_GAS(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB,NGSOBT,NTOOBS,NTOOB,IREOST,IREOTS,ISFTO,IBSO,NOBPTS,IOBPTS,ISFSO,NOBPT)
! Obtain Reordering arrays for orbitals
! (See note below for assumed ordering)
!
! GAS version
!
! =====
! Input
! =====
!  NSMOB  : Number of orbital symmetries
!  MXPOBS : Max number of orbital symmetries allowed by program
!  MXPNGAS: Max number of GAS spaces allowed by program
!  NGAS   : Number of GAS spaces
!  NGSOB  : Number of GAS orbitals per symmetry and space
!  NGSOBT : Number of GAS orbitals per space
!  NTOOBS : Number of orbitals per symmetry,all types
!
! ======
! Output
! ======
!  IREOST : Reordering array symmetry => type
!  IREOTS : Reordering array type     => symmetry
!  ISFTO  : Symmetry array for type ordered orbitals
!  IBSO   : First orbital of given symmetry (symmetry ordered)
!  NOBPTS : Number of orbitals per subtype and symmetry
!  IOBPTS : Off sets for orbitals of given subtype and symmetry
!           ordered according to input integrals
!
! ISFSO  : Symmetry of orbitals, symmetry ordereing
!
! Jeppe Olsen, Winter 1994

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NSMOB, MXPOBS, MXPNGAS, NGAS, NGSOB(MXPOBS,MXPNGAS), NGSOBT(MXPNGAS), NTOOBS(NSMOB), NTOOB
integer(kind=iwp), intent(out) :: IREOST(NTOOB), IREOTS(NTOOB), ISFTO(NTOOB), IBSO(NSMOB), NOBPTS(MXPNGAS,MXPOBS), &
                                  IOBPTS(MXPNGAS,MXPOBS), ISFSO(NTOOB), NOBPT(MXPNGAS)
integer(kind=iwp) :: IADD, IBSSM, IGAS, IOFF, IORB, ISM, ISTOFF, ISYM, ITSOFF, NPREV

! ==========================
! Note on order of orbitals
! ==========================
!
! The orbitals are supposed to be imported ordered symmetry-type
! ordered as
!
! Loop over symmetries of orbitals
!  Loop over GAS spaces
!   Loop over orbitals of this sym and GAS
!   End of Loop over orbitals
!  End of Loop over Gas spaces
! End of loop over symmetries
!
! Internally the orbitals are reordered to type symmetry order
! where the outer loop is over types and the inner loop is
! over symmetries, i.e.
!
! Loop over GAS spaces
!  Loop over symmetries of orbitals
!   Loop over orbitals of this sym and GAS
!   End of Loop over orbitals
!  End of loop over symmetries
! End of Loop over Gas spaces
!
! 1:  Construct ISFTO, IREOST,IREOTS,NOBPTS,IOBPTS
!
! To get rid of annoying and incorrect compiler warnings
IBSSM = 0

ITSOFF = 1
do IGAS=1,NGAS
  do ISYM=1,NSMOB
    if (ISYM == 1) then
      IBSSM = 1
    else
      IBSSM = IBSSM+NTOOBS(ISYM-1)
    end if
    NPREV = sum(NGSOB(ISYM,1:IGAS-1))
    IADD = 0
    NOBPTS(IGAS,ISYM) = NGSOB(ISYM,IGAS)
    IOBPTS(IGAS,ISYM) = ITSOFF
    !NOBPTS(ISYM,IGAS) = NGSOB(ISYM,IGAS)
    !IOBPTS(ISYM,IGAS) = ITSOFF
    do IORB=ITSOFF,ITSOFF+NGSOB(ISYM,IGAS)-1
      IADD = IADD+1
      IREOTS(IORB) = IBSSM-1+NPREV+IADD
      IREOST(IBSSM-1+NPREV+IADD) = IORB
      !ITFTO(IORB) = IGAS
      ISFTO(IORB) = ISYM
    end do
    ITSOFF = ITSOFF+NGSOB(ISYM,IGAS)
  end do
end do

! 2 : ISFSO,ITFSO

ISTOFF = 1
do ISYM=1,NSMOB
  do IGAS=1,NGAS
    ISFSO(ISTOFF:ISTOFF+NGSOB(ISYM,IGAS)-1) = ISYM
    !ITFSO(ISTOFF:ISTOFF+NGSOB(ISYM,IGAS)-1) = IGAS
    ISTOFF = ISTOFF+NGSOB(ISYM,IGAS)
  end do
end do

! 3 IBSO, NOBPT

IOFF = 1
do ISM=1,NSMOB
  IBSO(ISM) = IOFF
  IOFF = IOFF+NTOOBS(ISM)
end do
NOBPT(1:NGAS) = NGSOBT(1:NGAS)

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' =================='
write(u6,*) ' Output from ORBORD'
write(u6,*) ' =================='
write(u6,*)
write(u6,*) ' Symmetry of orbitals, type ordered'
call IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
write(u6,*) ' Symmetry => type reordering array'
call IWRTMA(IREOST,1,NTOOB,1,NTOOB)
write(u6,*) ' Type => symmetry reordering array'
call IWRTMA(IREOTS,1,NTOOB,1,NTOOB)
write(u6,*) ' IBSO array'
call IWRTMA(IBSO,1,NSMOB,1,NSMOB)

write(u6,*) ' NOBPTS'
call IWRTMA(NOBPTS,NGAS,NSMOB,MXPNGAS,MXPOBS)
write(u6,*) ' NOBPT'
call IWRTMA(NOBPT,NGAS,1,MXPNGAS,1)
write(u6,*) ' IOBPTS'
call IWRTMA(IOBPTS,NGAS,NSMOB,MXPNGAS,MXPOBS)

write(u6,*) ' ISFTO array :'
call IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
!write(u6,*) ' ITFSO array :'
!call IWRTMA(ITFSO,1,NTOOB,1,NTOOB)

write(u6,*) ' ISFSO array :'
call IWRTMA(ISFSO,1,NTOOB,1,NTOOB)
!write(u6,*) ' ITFTO array :'
!call IWRTMA(ITFTO,1,NTOOB,1,NTOOB)
#endif

end subroutine ORBORD_GAS
