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
! Copyright (C) 1991-1994, Jeppe Olsen                                 *
!***********************************************************************

subroutine RSBB2BN_MCLR(IASM,IATP,IBSM,IBTP,NIA,NIB,JASM,JATP,JBSM,JBTP,NJA,NJB,IAGRP,IBGRP,IAEL1,IAEL3,JAEL1,JAEL3,IBEL1,IBEL3, &
                        JBEL1,JBEL3,SB,CB,NTSOB,IBTSOB,ITSOB,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,NSMOB,NSMST, &
                        NSMSX,NSMDX,MXPOBS,IUSEAB,ICJKAIB,CJRES,SIRES,S2,ISIGN,ieaw,TimeDep)
! Combined alpha-beta double excitation
! contribution from given C block to given S block
! If IUSAB only half the terms are constructed
! =====
! Input
! =====
!
! IASM,IATP : Symmetry and type of alpha  strings in sigma
! IBSM,IBTP : Symmetry and type of beta   strings in sigma
! JASM,JATP : Symmetry and type of alpha  strings in C
! JBSM,JBTP : Symmetry and type of beta   strings in C
! NIA,NIB : Number of alpha-(beta-) strings in sigma
! NJA,NJB : Number of alpha-(beta-) strings in C
! IAGRP : String group of alpha strings
! IBGRP : String group of beta strings
! IAEL1(3) : Number of electrons in RAS1(3) for alpha strings in sigma
! IBEL1(3) : Number of electrons in RAS1(3) for beta  strings in sigma
! JAEL1(3) : Number of electrons in RAS1(3) for alpha strings in C
! JBEL1(3) : Number of electrons in RAS1(3) for beta  strings in C
! CB   : Input C block
! NTSOB  : Number of orbitals per type and symmetry
! IBTSOB : base for orbitals of given type and symmetry
! IBORB  : Orbitals of given type and symmetry
! NSMOB,NSMST,NSMSX : Number of symmetries of orbitals,strings,
!       single excitations
! MAXK   : Largest number of inner resolution strings treated at simult.
!
! ICJKAIB =1 =>  construct C(Ka,Jb,j) and S(Ka,Ib,i) as intermediate
!                 matrices in order to reduce overhead
!
! ======
! Output
! ======
! SB : updated sigma block
!
! =======
! Scratch
! =======
!
! SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
!              largest number of orbital pairs of given symmetries and
!              types.
! I1, XI1S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! I2, XI2S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! H : Space for two electron integrals
!
! Jeppe Olsen, Winter of 1991
!
! Feb 92 : Loops restructured ; Generation of I2,XI2S moved outside
! October 1993 : IUSEAB added
! January 1994 : Loop restructured + ICJKAIB introduced
! February 1994 : Fetching and adding to transposed blocks

use Symmetry_Info, only: Mul

implicit real*8(A-H,O-Z)
! General input
logical TimeDep
integer NTSOB(3,*), IBTSOB(3,*), ITSOB(*)
! Input
dimension CB(*)
! Output
dimension SB(*)
! Scratch
dimension SSCR(*), CSCR(*), I1(MAXK,*), XI1S(MAXK,*)
dimension I2(MAXK,*), XI2S(MAXK,*)
dimension I3(MAXK,*), XI3S(MAXK,*)
dimension I4(MAXK,*), XI4S(MAXK,*)
dimension XINT(*)
dimension CJRES(*), SIRES(*)
dimension S2(*)
! Local arrays
dimension ITP(3), JTP(3), KTP(3), LTP(3)

ZERO = 0.0d0
ONEM = -1.0d0
!IUSEAB = 0
!ICJKAIB = 1
IROUTE = 1
! Symmetry of allowed excitations
IJSM = Mul(IASM,JASM)
KLSM = Mul(IBSM,JBSM)
if ((IJSM == 0) .or. (KLSM == 0)) goto 9999
! Types of SX that connects the two strings
call SXTYP(NKLTYP,KTP,LTP,IBEL1,IBEL3,JBEL1,JBEL3)
call SXTYP(NIJTYP,ITP,JTP,IAEL1,IAEL3,JAEL1,JAEL3)
if ((NIJTYP == 0) .or. (NKLTYP == 0)) goto 9999
do IJTYP=1,NIJTYP
  ITYP = ITP(IJTYP)
  JTYP = JTP(IJTYP)
  ! TESTTING
  N1IND = 0
  N2IND = 0
  N3IND = 0
  if (ITYP == 1) N1IND = N1IND+1
  if (ITYP == 2) N2IND = N2IND+1
  if (ITYP == 3) N3IND = N3IND+1

  if (JTYP == 1) N1IND = N1IND+1
  if (JTYP == 2) N2IND = N2IND+1
  if (JTYP == 3) N3IND = N3IND+1

  do ISM=1,NSMOB
    JSM = Mul(ISM,IJSM)
    if (JSM == 0) goto 1940
    IOFF = IBTSOB(ITYP,ISM)
    JOFF = IBTSOB(JTYP,JSM)
    NI = NTSOB(ITYP,ISM)
    NJ = NTSOB(JTYP,JSM)
    if ((NI == 0) .or. (NJ == 0)) goto 1940
    ! Loop over batches of KA strings
    KABOT = 1-MAXK
    KATOP = 0
1801 continue
    KABOT = KABOT+MAXK
    KATOP = KATOP+MAXK
    ! Find Ka strings that connect with Ja strings for given group of Jorbs
    call ADST(JOFF,NJ,JATP,JASM,IAGRP,KABOT,KATOP,I1(1,1),XI1S(1,1),MAXK,NKABTC,KAEND)
    call ADST(IOFF,NI,IATP,IASM,IAGRP,KABOT,KATOP,I3(1,1),XI3S(1,1),MAXK,NKABTC,KAEND)
    if (NKABTC == 0) goto 1940
    ! Generate - if required C(Ka,Jb,j)
    if (ICJKAIB /= 0) then
      IIOFF = 1
      LCJ = NJB*NKABTC
      do JJ=1,NJ
        if (JJ == 1) then
          IIOFF = 1
        else
          IIOFF = IIOFF+LCJ
        end if
        call GATRMT(CB,NJA,NJB,CJRES(IIOFF),NKABTC,NJB,I1(1,JJ),XI1S(1,JJ))
      end do

      ! We have now C gathered in the form C(Ka,Jb,j).
      ! If Ka is small, say 1, it can be advantageous to switch
      ! around to C(j,Ka,Jb). This is mediated by the switch IROUTE
      ! IROUTE = 1 : Normal (i.e. old) route,
      ! IROUTE = 2 : New route with j first
      ! IROUTE = 3 : C(Ka,j,Jb)

      if ((NJ >= NKABTC) .and. (NI >= NKABTC)) then
        IROUTE = 2
      else
        IROUTE = 3
      end if
      ! DOES THIS WORK?
      !9805EAW     IROUTE = 1
      if (TimeDep) IROUTE = 1

      if (IROUTE == 2) then
        ! C(Ka,Jb,j) => C(j,Ka,Jb)
        call TRNSPS(LCJ,NJ,CJRES,SIRES)
        call DCOPY_(NJ*LCJ,SIRES,1,CJRES,1)
      end if
      if (IROUTE == 3) then
        ! C(Ka,Jb,j) => C(Ka,j,JB)
        do JB=1,NJB
          do J=1,NJ
            IOFFIN = (J-1)*NJB*NKABTC+(JB-1)*NKABTC+1
            IOFFOUT = (JB-1)*NKABTC*NJ+(J-1)*NKABTC+1
            do KA=1,NKABTC
              SIRES(IOFFOUT-1+KA) = CJRES(IOFFIN-1+KA)
            end do
          end do
        end do
        call DCOPY_(NJ*LCJ,SIRES,1,CJRES,1)
      end if

      call DCOPY_(NIB*NKABTC*NI,[ZERO],0,SIRES,1)

    end if

    do KLTYP=1,NKLTYP
      KTYP = KTP(KLTYP)
      LTYP = LTP(KLTYP)
      ! Testing
      N1IND = 0
      N2IND = 0
      N3IND = 0
      if (ITYP == 1) N1IND = N1IND+1
      if (ITYP == 2) N2IND = N2IND+1
      if (ITYP == 3) N3IND = N3IND+1

      if (JTYP == 1) N1IND = N1IND+1
      if (JTYP == 2) N2IND = N2IND+1
      if (JTYP == 3) N3IND = N3IND+1

      if (KTYP == 1) N1IND = N1IND+1
      if (KTYP == 2) N2IND = N2IND+1
      if (KTYP == 3) N3IND = N3IND+1

      if (LTYP == 1) N1IND = N1IND+1
      if (LTYP == 2) N2IND = N2IND+1
      if (LTYP == 3) N3IND = N3IND+1
      do KSM=1,NSMOB
        IFIRST = 1
        LSM = Mul(KSM,KLSM)
        if (LSM == 0) goto 1930
        KOFF = IBTSOB(KTYP,KSM)
        LOFF = IBTSOB(LTYP,LSM)
        NK = NTSOB(KTYP,KSM)
        NL = NTSOB(LTYP,LSM)
        ! If IUSEAB  is used, only terms with i >= k will be generated so
        IKORD = 0
        if ((IUSEAB == 1) .and. (ISM > KSM)) goto 1930
        if ((IUSEAB == 1) .and. (ISM == KSM) .and. (ITYP > KTYP)) goto 1930
        if ((IUSEAB == 1) .and. (ISM == KSM) .and. (ITYP == KTYP)) IKORD = 1

        if ((NK == 0) .or. (NL == 0)) goto 1930
        ! Loop over batches of KB strings
        KBBOT = 1-MAXK
        KBTOP = 0
1800    continue
        KBBOT = KBBOT+MAXK
        KBTOP = KBTOP+MAXK
        ! obtain cb(KA,KB,jl) =  sum(JA,JB)<KA!a la!JA><KB!a jb !JB>C(JA,JB)

        call ADST(LOFF,NL,JBTP,JBSM,IBGRP,KBBOT,KBTOP,I2(1,1),XI2S(1,1),MAXK,NKBBTC,KBEND)
        call ADST(KOFF,NK,IBTP,IBSM,IBGRP,KBBOT,KBTOP,I4(1,1),XI4S(1,1),MAXK,NKBBTC,KBEND)
        if (NKBBTC == 0) goto 1930

        ! Modern low copy version
        if (IFIRST == 1) then
          IXCHNG = 0
          if (IROUTE == 1) then
            ! Integrals stored as (j l i k)
            !write(6,*) 'Timedep in rsbb2bn_mclr;',TimeDep
            if (TimeDep) then
              call GETINT_td(XINT,JTYP,JSM,ITYP,ISM,LTYP,LSM,KTYP,KSM,0,0,iroute,ieaw)
            else
              call GETINT_MCLR(XINT,JTYP,JSM,ITYP,ISM,LTYP,LSM,KTYP,KSM,IXCHNG,0,0,0,ieaw)
            end if
          else if (IROUTE == 2) then
            ! Integrals stored as (i j k l)
            if (TimeDep) then
              call GETINT_td(XINT,ITYP,ISM,JTYP,JSM,KTYP,KSM,LTYP,LSM,0,0,iroute,ieaw)
            else
              call GETINT_MCLR(XINT,ITYP,ISM,JTYP,JSM,KTYP,KSM,LTYP,LSM,IXCHNG,0,0,1,ieaw)
            end if
          else if (IROUTE == 3) then
            ! Integrals stored as (j i k l)
            if (TimeDep) then
              call GETINT_td(XINT,JTYP,JSM,ITYP,ISM,KTYP,KSM,LTYP,LSM,0,0,iroute,ieaw)
            else
              call GETINT_MCLR(XINT,JTYP,JSM,ITYP,ISM,KTYP,KSM,LTYP,LSM,IXCHNG,0,0,1,ieaw)
            end if
          end if
          if (ISIGN == -1) call DSCAL_(NI*NJ*NK*NL,ONEM,XINT,1)
          IFIRST = 0
        end if
        call SKICKJ_MCLR(SIRES,CJRES,NKABTC,NIB,NJB,NKBBTC,XINT,NI,NJ,NK,NL,MAXK,I4,XI4S,I2,XI2S,IKORD,IDUM,iXDUM,XDUM,IROUTE)

        if (KBEND == 0) goto 1800
        ! End of loop over partitioning of beta strings
1930    continue
      end do
    end do
    ! Scatter out from s(Ka,Ib,i)
    ! Restore order !!
    if (IROUTE == 2) then
      call TRNSPS(NI,NIB*NKABTC,SIRES,CJRES)
      call DCOPY_(NI*NIB*NKABTC,CJRES,1,SIRES,1)
    end if
    if (IROUTE == 3) then
      do JB=1,NIB
        do J=1,NI
          IOFFIN = (J-1)*NIB*NKABTC+(JB-1)*NKABTC+1
          IOFFOUT = (JB-1)*NKABTC*NI+(J-1)*NKABTC+1
          do KA=1,NKABTC
            CJRES(IOFFIN-1+KA) = SIRES(IOFFOUT-1+KA)
          end do
        end do
      end do
      call DCOPY_(NI*NIB*NKABTC,CJRES,1,SIRES,1)
    end if
    if (ICJKAIB == 1) then
      do II=1,NI
        call SCARMT(SIRES((II-1)*NKABTC*NIB+1),NKABTC,NIB,SB,NIA,NIB,I3(1,II),XI3S(1,II))

      end do
    end if
    if (KAEND == 0) goto 1801
    ! End of loop over partitioning of alpha strings
1940 continue
  end do
end do

9999 continue

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(ITSOB)
  call Unused_real_array(SSCR)
  call Unused_real_array(CSCR)
  call Unused_integer(NSMST)
  call Unused_integer(NSMSX)
  call Unused_integer(NSMDX)
  call Unused_integer(MXPOBS)
  call Unused_real_array(S2)
end if

end subroutine RSBB2BN_MCLR
