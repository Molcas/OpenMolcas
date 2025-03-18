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

subroutine ADADS1(NK,I1,XI1S,IOBSM,IOBTP,IOBOFF,NIOB,JOBSM,JOBTP,JOBOFF,NJOB,IJORD,NKDIM,ICLS,ISM,I2MAPO,I2MAPS,I2MPF,L2MP,I2MPO, &
                  I2MPL,I1MAPO,I1MAPS,I1MPF,L1MP,I1MPO,I1MPL,IEL1,IEL3,I2EL1,I2EL3,ISSO,NSSO,I2SSO,N2SSO,NOCTP,N2OCTP,NORB1,NORB2, &
                  NORB3,NORB,KMAX,KMIN,IEND)
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
! NSSO : Number of TS strings for I strings
! I2SSO : TS symmetry offset for K strings
! N2SSO : Number of TS strings for K strings
! NOCTP : Number of occupation types for I strings
! N2OCTP : Number of occupation types for K strings
! NORB1(2,3) : Number of RAS1(2,3) orbitals
! IORBSM : Orbital symmety array
! NORB : Number of active  orbitals
! KMAX : Largest allowed relative number for K strings
! KMIN : Smallest allowed relative number for K strings
!
! ======
! Output
! ======
!
! NK      : Number of K strings
! I1(KSTR) : ne. 0 => a + IORB a+ JORB !KSTR> = +/-!ISTR>
! XI1S(KSTR) : above +/-
!          : eq. 0    a + IORB !KSTR> = 0
! Offset is KMIN
! IEND : = 0 => end of N-2 strings has not been encountered
! IEND : = 1 => end of N-2 strings has     been encountered

use Symmetry_Info, only: Mul

implicit real*8(A-H,O-Z)
!.Input
integer IEL1(*), IEL3(*), I2EL1(*), I2EL3(*)
integer ISSO(NOCTP,*), NSSO(NOCTP,*)
integer I2SSO(N2OCTP,*), N2SSO(N2OCTP,*)
integer I2MAPO(*), I2MAPS(*)
integer I1MAPO(*), I1MAPS(*)
integer I2MPL(*), I2MPO(*)
!eaw
integer I1mpo(*), i1mpl(*)
!eaw
!.Output
integer I1(NKDIM,*)
dimension XI1S(NKDIM,*)

NIJ = 0 ! dummy initialize
iprstr = 0
NTEST = IPRSTR !100
if (NTEST > 0) then
  write(6,*) ' ================'
  write(6,*) ' Info from ADADS1'
  write(6,*) ' ================'
  write(6,*) ' Iobsm Iobtp Ioboff, Niob',IOBSM,IOBTP,IOBOFF,NIOB
  write(6,*) ' Jobsm Jobtp Joboff, NJob',JOBSM,JOBTP,JOBOFF,NJOB
  write(6,*) ' icls ism',ICLS,ISM
  write(6,*) ' I1MPF, I2MPF ',I1MPF,I2MPF
end if

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
!? write(6,*) ' icls ',icls
!? write(6,*) ' iel1 iel3 ',iel1(icls),iel3(icls)
!? write(6,*) ' kel1 kel3 ',kel1,kel3
KTYPE = 0
do KKTYPE=1,N2OCTP
  if ((I2EL1(KKTYPE) == KEL1) .and. (I2EL3(KKTYPE) == KEL3)) KTYPE = KKTYPE
end do
!? write(6,*) ' ktype ',ktype
if (KTYPE == 0) then
  NK = 0
  IEND = 1
  goto 101
end if
! Symmetry of K strings
JKSM = Mul(IOBSM,ISM)
if (JKSM == 0) then
  IKSM = Mul(JOBSM,ISM)
  if (IKSM == 0) then
    NK = 0
    IEND = 1
    goto 101
  else
    KSM = Mul(IOBSM,IKSM)
  end if
else
  KSM = Mul(JOBSM,JKSM)
end if
!? write(6,*) ' JOBSM,KSM,JKSM ',JOBSM,KSM,JKSM
if (KSM == 0) then
  NK = 0
  IEND = 1
  goto 101
end if
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
  NIJ = NIOB*(NIOB+1)/2
end if
if (NKDim > 0) then
  do IJ=1,NIJ
    call ICopy(NKDIM,[0],0,I1(1,IJ),1)
    call dcopy_(NKDIM,[0.0d0],0,XI1S(1,IJ),1)
  end do
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
    IJOFF = (J-1)*NIOB-(J-1)*(J-2)/2+1-imin+1
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
      do JJORB=IFST,min(JORb,I2MPL(KSTR))
        if (I2MAPO(I2MPO(KSTR)-1+JJORB) == JORB) JKSTR = I2MAPS(I2MPO(KSTR)-1+JJORB)
      end do
    end if
    if (JKSTR == 0) goto 100
    if (JKSTR > 0) then
      SIGN = 1.0d0
    else
      JKSTR = -JKSTR
      SIGN = -1.0d0
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
      if (ISTR == 0) goto 99
      ! Synthesis
      if (ISTR > 0) then
        I1(KSTR-KOFF-KMIN+2,IJ) = ISTR-IOFF+1
        XI1S(KSTR-KOFF-KMIN+2,IJ) = SIGN
      else
        I1(KSTR-KOFF-KMIN+2,IJ) = -ISTR-IOFF+1
        XI1S(KSTR-KOFF-KMIN+2,IJ) = -SIGN
      end if
99    continue
    end do

100 continue
  end do

  !end do
end do
101 continue

if (NTEST > 0) then
  write(6,*) ' Output from ADADS1'
  write(6,*) ' =================='
  write(6,*) ' Number of K strings accessed ',NK
  if (NK /= 0) then
    do IJ=1,NIJ
      write(6,*) ' Excited strings and sign for ij = ',IJ
      call IWRTMA(I1(1,IJ),1,NK,1,NK)
      call WRTMAT(XI1S(1,IJ),1,NK,1,NK)
    end do
  end if
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(NSSO)
  call Unused_integer(NORB1)
  call Unused_integer(NORB2)
  call Unused_integer(NORB3)
end if

end subroutine ADADS1
