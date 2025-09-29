!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ADS1(NK,I1,XI1S,LI1,IORB,LORB,ICLS,ISM,IMAPO,IMAPS,IMPL,IMPO,IMPF,LMAP,IEL1,IEL3,I1EL1,I1EL3,ISSO,I1SSO,N1SSO,NOCTP, &
                N1OCTP,NORB1,NORB2,ORBSM,KMAX,KMIN,IEND)
! Obtain I1(KSTR) = +/- A+ IORB !KSTR>
!
! KSTR is restricted to strings with relative numbers in the
! range KMAX to KMIN

! =====
! Input
! =====
! IORB        : Firat orbital to be added
! LORB        : Number of orbitals to be added : IORB to IORB-1+LORB
!               are used. They must all be in the same TS group
! ICLS,ISM    : Class and symmetry of string with added electron
! IMAPO,IMAPS : map from Kstrings to Istrings
! IEL1(3)     : Number of electrons in RAS1(3) for I strings
! I1EL1(3)    : Number of electrons in RAS1(3) for K strings
! ISSO        : TS symmetry offset for I strings
! I1SSO       : TS symmetry offset for K strings
! N1SSO       : Number of TS strings for K strings
! NOCTP       : Number of occupation types for I strings
! N1OCTP      : Number of occupation types for K strings
! NORB1(2,3)  : Number of RAS1(2,3) orbitals
! IORBSM      : Orbital symmety array
! KMAX        : Largest allowed relative number for K strings
!               If Kmax is set to -1 all strings are searched
! KMIN        : Smallest allowed relative number for K strings
!
! ======
! Output
! ======
! NK              : Number of K strings
! I1(KSTR,JORB)   : /= 0 => a + JORB !KSTR> = +/-!ISTR>
! XI1S(KSTR,JORB) : above +/-
!                 : == 0    a + JORB !KSTR> = 0
! Offset is KMIN

use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: NK, IEND
integer(kind=iwp), intent(_OUT_) :: I1(*)
real(kind=wp), intent(_OUT_) :: XI1S(*)
integer(kind=iwp), intent(in) :: LI1, IORB, LORB, ICLS, ISM, IMAPO(*), IMAPS(*), IMPL(*), IMPO(*), IMPF, LMAP, IEL1(*), IEL3(*), &
                                 I1EL1(*), I1EL3(*), NOCTP, ISSO(NOCTP,*), N1OCTP, I1SSO(N1OCTP,*), N1SSO(N1OCTP,*), NORB1, NORB2, &
                                 ORBSM(*), KMAX, KMIN
integer(kind=iwp) :: IIIORB, IIORB, IOFF, IORBR, ISTR, KEL1, KEL3, KEND, KKTYPE, KOFF, KREL, KSM, KSTR, KSUB, KTYPE, LDIM
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IIORBR
#endif
logical(kind=iwp) :: Skip

LDIM = 0 ! dummy initialize
#ifdef _DEBUGPRINT_
write(u6,*) ' =============='
write(u6,*) ' ADSTS speaking'
write(u6,*) ' =============='
write(u6,*) ' IORB,ISM,ICLS',IORB,ISM,ICLS
write(u6,*) ' IMPF, LMAP ',IMPF,LMAP
write(u6,*) ' N1SSO :'
call IWRTMA(N1SSO,N1OCTP,8,N1OCTP,8)
#endif
NK = KMAX-KMIN+1
! Type of kstrings
if (IORB <= NORB1) then
  KEL1 = IEL1(ICLS)-1
  KEL3 = IEL3(ICLS)
else if (IORB <= NORB1+NORB2) then
  KEL1 = IEL1(ICLS)
  KEL3 = IEL3(ICLS)
else
  KEL1 = IEL1(ICLS)
  KEL3 = IEL3(ICLS)-1
end if
KTYPE = 0
!write(u6,*) ' N1OCTP ',N1OCTP
do KKTYPE=1,N1OCTP
  if ((I1EL1(KKTYPE) == KEL1) .and. (I1EL3(KKTYPE) == KEL3)) KTYPE = KKTYPE
end do
!write(u6,*) ' kel1 kel3 ktype ',KEL1,KEL3,KTYPE
Skip = .false.
if (KTYPE == 0) then
  NK = 0
  IEND = 1
  Skip = .true.
else
  ! Symmetry of K strings
  KSM = Mul(ORBSM(IORB),ISM)
  if (KSM == 0) then
    NK = 0
    IEND = 1
    Skip = .true.
  end if
end if
if (.not. Skip) then
  KOFF = I1SSO(KTYPE,KSM)
  !write(u6,*) ' KTYPE KSM ',KTYPE,KSM
  if (KMAX == -1) then
    KEND = N1SSO(KTYPE,KSM)
  else
    KEND = min(N1SSO(KTYPE,KSM),KMAX)
  end if
  if (KEND < N1SSO(KTYPE,KSM)) then
    IEND = 0
  else
    IEND = 1
  end if
  NK = KEND-KMIN+1
  if (KMAX == -1) then
    LDIM = NK
  else
    LDIM = LI1
  end if
  !if (KMAX == -1) write(u6,*) ' KMAX = -1, LDIM=',LDIM
  IOFF = ISSO(ICLS,ISM)
  KSUB = KOFF+KMIN-2
  do IIORB=IORB,IORB+LORB-1
    IORBR = IIORB-IORB+1
    do KSTR=KOFF+KMIN-1,KOFF+KEND-1
      !write(u6,*) ' KSTR = ',KSTR
      KREL = KSTR-KSUB

      ISTR = 0
      if (IMPF == 1) then
        if (IMAPO((KSTR-1)*LMAP+IIORB) == IIORB) ISTR = IMAPS((KSTR-1)*LMAP+IIORB)
      else
        !write(u6,*) ' IMPL = ',IMPL(KSTR)
        !write(u6,*) ' IMPO = ',IMPO(KSTR)
        do IIIORB=1,IMPL(KSTR)
          if (IMAPO(IMPO(KSTR)-1+IIIORB) == IIORB) ISTR = IMAPS(IMPO(KSTR)-1+IIIORB)
        end do
      end if
      if (ISTR > 0) then
        I1(KREL+(IORBR-1)*LDIM) = ISTR-IOFF+1
        XI1S(KREL+(IORBR-1)*LDIM) = One
      else if (ISTR < 0) then
        I1(KREL+(IORBR-1)*LDIM) = -ISTR-IOFF+1
        XI1S(KREL+(IORBR-1)*LDIM) = -One
      else if (ISTR == 0) then
        I1(KREL+(IORBR-1)*LDIM) = 0
        XI1S(KREL+(IORBR-1)*LDIM) = Zero
      end if
    end do
  end do
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Output from ASTR'
write(u6,*) ' ================'
write(u6,*) ' Number of K strings accessed ',NK
if (NK /= 0) then
  do IIORB=IORB,IORB+LORB-1
    IIORBR = IIORB-IORB
    write(u6,*) ' Info for orbital ',IIORB
    write(u6,*) ' Excited strings and sign'
    call IWRTMA(I1(IIORBR*LDIM+1),1,NK,1,NK)
    call WRTMAT(XI1S(IIORBR*LDIM+1),1,NK,1,NK)
  end do
end if
#endif

end subroutine ADS1
