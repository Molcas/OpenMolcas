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
subroutine ADADS1(NK,I1,XI1S,IOBSM,IOBTP,IOBOFF,NIOB,JOBSM,JOBTP,JOBOFF,NJOB,IJORD,NKDIM,ICLS,ISM,I2MAPO,I2MAPS,I2MPF,L2MP,I2MPO, &
                  I2MPL,I1MAPO,I1MAPS,I1MPF,L1MP,I1MPO,I1MPL,IEL1,IEL3,I2EL1,I2EL3,ISSO,I2SSO,N2SSO,NOCTP,N2OCTP,NORB,KMAX,KMIN, &
                  IEND)
! Obtain I1(KSTR) = +/- A+ IORB !KSTR>
!
! KSTR is restricted to strings with relative numbers in the
! range KMAX to KMIN
! =====
! Input
! =====
! ICLS,ISM : Class and symmetry of string with added electron
! I2MAPO,I2MAPS : map N-2 strings to N-1 strings
! I1MAPO,I1MAPS : map N-1 strings to  strings
! IEL1(3) : Number of electrons in RAS1(3) for I strings
! I2EL1(3) : Number of electrons in RAS1(3) for K strings
! ISSO : TS symmetry offset for I strings
! I2SSO : TS symmetry offset for K strings
! N2SSO : Number of TS strings for K strings
! NOCTP : Number of occupation types for I strings
! N2OCTP : Number of occupation types for K strings
! IORBSM : Orbital symmety array
! NORB : Number of active  orbitals
! KMAX : Largest allowed relative number for K strings
! KMIN : Smallest allowed relative number for K strings
!
! ======
! Output
! ======
!
! NK         : Number of K strings
! I1(KSTR)   : /= 0 => a + IORB a+ JORB !KSTR> = +/-!ISTR>
! XI1S(KSTR) : above +/-
!            : == 0    a + IORB !KSTR> = 0
! Offset is KMIN
! IEND : = 0 => end of N-2 strings has not been encountered
! IEND : = 1 => end of N-2 strings has     been encountered

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: NK, IEND
integer(kind=iwp), intent(in) :: IOBSM, IOBTP, IOBOFF, NIOB, JOBSM, JOBTP, JOBOFF, NJOB, IJORD, NKDIM, ICLS, ISM, I2MAPO(*), &
                                 I2MAPS(*), I2MPF, L2MP, I2MPO(*), I2MPL(*), I1MAPO(*), I1MAPS(*), I1MPF, L1MP, I1MPO(*), &
                                 I1MPL(*), IEL1(*), IEL3(*), I2EL1(*), I2EL3(*), NOCTP, ISSO(NOCTP,*), N2OCTP, I2SSO(N2OCTP,*), &
                                 N2SSO(N2OCTP,*), NORB, KMAX, KMIN
integer(kind=iwp), intent(_OUT_) :: I1(NKDIM,*)
real(kind=wp), intent(_OUT_) :: XI1S(NKDIM,*)
integer(kind=iwp) :: i, IFST, IIORB, ij, IJOFF, IKSM, IMIN, IOFF, iorb, ISTR, J, JJORB, JKSM, JKSTR, JORB, KEL1, KEL3, KEND, &
                     KKTYPE, KOFF, KSM, KSTR, KTYPE, NIJ
real(kind=wp) :: SGN
logical(kind=iwp) :: Skip

NIJ = 0 ! dummy initialize
#ifdef _DEBUGPRINT_
write(u6,*) ' ================'
write(u6,*) ' Info from ADADS1'
write(u6,*) ' ================'
write(u6,*) ' Iobsm Iobtp Ioboff, Niob',IOBSM,IOBTP,IOBOFF,NIOB
write(u6,*) ' Jobsm Jobtp Joboff, NJob',JOBSM,JOBTP,JOBOFF,NJOB
write(u6,*) ' icls ism',ICLS,ISM
write(u6,*) ' I1MPF, I2MPF ',I1MPF,I2MPF
#endif

NK = KMAX-KMIN+1
! Type of kstrings
if (IOBTP == 1) then
  KEL1 = IEL1(ICLS)-1
  KEL3 = IEL3(ICLS)
else if (IOBTP == 2) then
  KEL1 = IEL1(ICLS)
  KEL3 = IEL3(ICLS)
else
  KEL1 = IEL1(ICLS)
  KEL3 = IEL3(ICLS)-1
end if
if (JOBTP == 1) then
  KEL1 = KEL1-1
  KEL3 = KEL3
else if (JOBTP == 2) then
  KEL1 = KEL1
  KEL3 = KEL3
else
  KEL1 = KEL1
  KEL3 = KEL3-1
end if
!write(u6,*) ' icls ',icls
!write(u6,*) ' iel1 iel3 ',iel1(icls),iel3(icls)
!write(u6,*) ' kel1 kel3 ',kel1,kel3
KTYPE = 0
do KKTYPE=1,N2OCTP
  if ((I2EL1(KKTYPE) == KEL1) .and. (I2EL3(KKTYPE) == KEL3)) KTYPE = KKTYPE
end do
!write(u6,*) ' ktype ',ktype
Skip = .false.
if (KTYPE == 0) then
  NK = 0
  IEND = 1
  Skip = .true.
else
  ! Symmetry of K strings
  KSM = 0
  JKSM = Mul(IOBSM,ISM)
  if (JKSM == 0) then
    IKSM = Mul(JOBSM,ISM)
    if (IKSM == 0) then
      NK = 0
      IEND = 1
      Skip = .true.
    else
      KSM = Mul(IOBSM,IKSM)
    end if
  else
    KSM = Mul(JOBSM,JKSM)
  end if
  !write(u6,*) ' JOBSM,KSM,JKSM ',JOBSM,KSM,JKSM
  if (KSM == 0) then
    NK = 0
    IEND = 1
    Skip = .true.
  end if
end if
if (.not. Skip) then
  KOFF = I2SSO(KTYPE,KSM)
  KEND = min(KMAX,N2SSO(KTYPE,KSM))
  if (KEND == N2SSO(KTYPE,KSM)) then
    IEND = 1
  else
    IEND = 0
  end if
  NK = KEND-KMIN+1
  IOFF = ISSO(ICLS,ISM)
  ! Loop over iorb,jorb
  if (IJORD == 0) then
    NIJ = NIOB*NJOB
  else
    NIJ = nTri_Elem(NIOB)
  end if
  if (NKDim > 0) then
    I1(:,1:NIJ) = 0
    XI1S(:,1:NIJ) = Zero
  end if

  do J=1,NJOB
    if (IJORD == 1) then
      IMIN = J
    else
      IMIN = 1
    end if
    !do I=IMIN,NIOB
    !  IJ = IJ+1
    !  do IJ=1,NIJ
    !    call NXTIJ(I,J,NIOB,NJOB,IJORD,NONEW)
    !    IORB = IOBOFF-1+I
    JORB = JOBOFF-1+J
    if (IJORD == 1) then
      IJOFF = (J-1)*NIOB-nTri_Elem(J-2)+1-imin+1
    else
      IJOFF = (J-1)*NIOB+1
    end if
    do KSTR=KOFF+KMIN-1,KOFF+KEND-1
      ! N-2 => N-1
      JKSTR = 0
      if (I2MPF == 1) then
        if (I2MAPO((KSTR-1)*L2MP+JORB) == JORB) JKSTR = I2MAPS((KSTR-1)*L2MP+JORB)
      else if (I2MPF == 0) then
        IFST = max(1,JORB+I2MPL(KSTR)-NORB)
        do JJORB=IFST,min(JORB,I2MPL(KSTR))
          if (I2MAPO(I2MPO(KSTR)-1+JJORB) == JORB) JKSTR = I2MAPS(I2MPO(KSTR)-1+JJORB)
        end do
      end if
      if (JKSTR == 0) cycle
      if (JKSTR > 0) then
        SGN = One
      else
        JKSTR = -JKSTR
        SGN = -One
      end if
      ! N-1 => N
      do i=imin,niob
        ij = ijoff+i-1
        iorb = ioboff-1+i
        ISTR = 0
        if (I1MPF == 1) then
          if (I1MAPO((JKSTR-1)*L1MP+IORB) == IORB) ISTR = I1MAPS((JKSTR-1)*L1MP+IORB)
        else if (I1MPF == 0) then
          IFST = max(1,IORB+NORB-I1MPL(JKSTR))
          do IIORB=IFST,min(IORB,I1MPL(JKSTR))
            if (I1MAPO(I1MPO(JKSTR)-1+IIORB) == IORB) ISTR = I1MAPS(I1MPO(JKSTR)-1+IIORB)
          end do
        end if
        if (ISTR == 0) cycle
        ! Synthesis
        if (ISTR > 0) then
          I1(KSTR-KOFF-KMIN+2,IJ) = ISTR-IOFF+1
          XI1S(KSTR-KOFF-KMIN+2,IJ) = SGN
        else
          I1(KSTR-KOFF-KMIN+2,IJ) = -ISTR-IOFF+1
          XI1S(KSTR-KOFF-KMIN+2,IJ) = -SGN
        end if
      end do

    end do

    !end do
  end do
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Output from ADADS1'
write(u6,*) ' =================='
write(u6,*) ' Number of K strings accessed ',NK
if (NK /= 0) then
  do IJ=1,NIJ
    write(u6,*) ' Excited strings and sign for ij = ',IJ
    call IWRTMA(I1(1,IJ),1,NK,1,NK)
    call WRTMAT(XI1S(1,IJ),1,NK,1,NK)
  end do
end if
#endif

end subroutine ADADS1
