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

subroutine ADS1(NK,I1,XI1S,LI1,IORB,LORB,ICLS,ISM,IMAPO,IMAPS,IMPL,IMPO,IMPF,LMAP,IEL1,IEL3,I1EL1,I1EL3,ISSO,NSSO,I1SSO,N1SSO, &
                NOCTP,N1OCTP,NORB1,NORB2,NORB3,ORBSM,NORB,KMAX,KMIN,IEND)
! Obtain I1(KSTR) = +/- A+ IORB !KSTR>
!
! KSTR is restricted to strings with relative numbers in the
! range KMAX to KMIN
! =====
! Input
! =====
! IORB : Firat orbital to be added
! LORB : Number of orbitals to be added : IORB to IORB-1+LORB
!        are used. They must all be in the same TS group
! ICLS,ISM : Class and symmetry of string with added electron
! IMAPO,IMAPS : map from Kstrings to Istrings
! IEL1(3) : Number of electrons in RAS1(3) for I strings
! I1EL1(3) : Number of electrons in RAS1(3) for K strings
! ISSO : TS symmetry offset for I strings
! NSSO : Number of TS strings for I strings
! I1SSO : TS symmetry offset for K strings
! N1SSO : Number of TS strings for K strings
! NOCTP : Number of occupation types for I strings
! N1OCTP : Number of occupation types for K strings
! NORB1(2,3) : Number of RAS1(2,3) orbitals
! IORBSM : Orbital symmety array
! NORB : Number of active  orbitals
! KMAX : Largest allowed relative number for K strings
!        If Kmax is set to -1 all strings are searched
! KMIN : Smallest allowed relative number for K strings
!
! ======
! Output
! ======
!
! NK      : Number of K strings
! I1(KSTR,JORB) : ne. 0 => a + JORB !KSTR> = +/-!ISTR>
! XI1S(KSTR,JORB) : above +/-
!          : eq. 0    a + JORB !KSTR> = 0
! Offset is KMIN

use Symmetry_Info, only: Mul
implicit real*8(A-H,O-Z)
! Input
integer IEL1(*), IEL3(*), I1EL1(*), I1EL3(*)
integer ISSO(NOCTP,*), NSSO(NOCTP,*)
integer I1SSO(N1OCTP,*), N1SSO(N1OCTP,*)
integer ORBSM(*)
!integer IMAPO(NORB,*), IMAPS(NORB,*)
integer IMAPO(*), IMAPS(*)
integer IMPL(*), IMPO(*)
! Output
integer I1(*)
dimension XI1S(*)

LDIM = 0 ! dummy initialize
NTEST = 000
if (NTEST /= 0) then
  write(6,*) ' =============='
  write(6,*) ' ADSTS speaking'
  write(6,*) ' =============='
  write(6,*) ' IORB,ISM,ICLS',IORB,ISM,ICLS
  write(6,*) ' IMPF, LMAP ',IMPF,LMAP
  write(6,*) ' N1SSO :'
  call IWRTMA(N1SSO,N1OCTP,8,N1OCTP,8)
end if
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
!write(6,*) ' N1OCTP ',N1OCTP
do KKTYPE=1,N1OCTP
  if ((I1EL1(KKTYPE) == KEL1) .and. (I1EL3(KKTYPE) == KEL3)) KTYPE = KKTYPE
end do
!write(6,*) ' kel1 kel3 ktype ',KEL1,KEL3,KTYPE
if (KTYPE == 0) then
  NK = 0
  IEND = 1
  goto 101
end if
! Symmetry of K strings
KSM = Mul(ORBSM(IORB),ISM)
if (KSM == 0) then
  NK = 0
  IEND = 1
  goto 101
end if
KOFF = I1SSO(KTYPE,KSM)
!? write(6,*) ' KTYPE KSM ',KTYPE,KSM
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
!? if (KMAX == -1) write(6,*) ' KMAX = -1, LDIM=',LDIM
IOFF = ISSO(ICLS,ISM)
KSUB = KOFF+KMIN-2
do IIORB=IORB,IORB+LORB-1
  IORBR = IIORB-IORB+1
  do KSTR=KOFF+KMIN-1,KOFF+KEND-1
    !write(6,*) ' KSTR = ',KSTR
    KREL = KSTR-KSUB

    ISTR = 0
    if (IMPF == 1) then
      if (IMAPO((KSTR-1)*LMAP+IIORB) == IIORB) ISTR = IMAPS((KSTR-1)*LMAP+IIORB)
    else
      !write(6,*) ' IMPL = ',IMPL(KSTR)
      !write(6,*) ' IMPO = ',IMPO(KSTR)
      do IIIORB=1,IMPL(KSTR)
        if (IMAPO(IMPO(KSTR)-1+IIIORB) == IIORB) ISTR = IMAPS(IMPO(KSTR)-1+IIIORB)
      end do
    end if
    if (ISTR > 0) then
      I1(KREL+(IORBR-1)*LDIM) = ISTR-IOFF+1
      XI1S(KREL+(IORBR-1)*LDIM) = 1.0d0
    else if (ISTR < 0) then
      I1(KREL+(IORBR-1)*LDIM) = -ISTR-IOFF+1
      XI1S(KREL+(IORBR-1)*LDIM) = -1.0d0
    else if (ISTR == 0) then
      I1(KREL+(IORBR-1)*LDIM) = 0
      XI1S(KREL+(IORBR-1)*LDIM) = 0.0d0
    end if
  end do
end do
101 continue

if (NTEST > 0) then
  write(6,*) ' Output from ASTR'
  write(6,*) ' ================'
  write(6,*) ' Number of K strings accessed ',NK
  if (NK /= 0) then
    do IIORB=IORB,IORB+LORB-1
      IIORBR = IIORB-IORB+1
      write(6,*) ' Info for orbital ',IIORB
      write(6,*) ' Excited strings and sign'
      call IWRTMA(I1(1+(IIORBR-1)*LDIM),1,NK,1,NK)
      call WRTMAT(XI1S(1+(IIORBR-1)*LDIM),1,NK,1,NK)
    end do
  end if
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(NSSO)
  call Unused_integer(NORB3)
  call Unused_integer(NORB)
end if

end subroutine ADS1
