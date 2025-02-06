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
! Copyright (C) 1996, Jeppe Olsen                                      *
!***********************************************************************

subroutine GSBBD2A_LUCIA(RHO2,RHO2S,RHO2A,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NGAS,ISEL,ICEL,SB,CB,ADSXA,STSTSX,SXDXSX, &
                         MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,X,NSMOB,NSMST,MXPOBS,SCLFAC,IPACK)
! SUBROUTINE GSBBD2A_LUCIA --> 37
!
! Contributions to two-electron density matrix from column excitations
!
! GAS version, '96, Jeppe Olsen
!
! =====
! Input
! =====
! RHO2  : two body density matrix to be updated
! NACOB : Number of active orbitals
! ISCSM,ISCTP : Symmetry and type of sigma columns
! ICCSM,ICCTP : Symmetry and type of C     columns
! IGRP : String group of columns
! NROW : Number of rows in S and C block
! NGAS : Number of active spaces
! ISEL : Number of electrons per AS for S block
! ICEL : Number of electrons per AS for C block
! CB   : Input C block
! ADASX : sym of a+, a => sym of a+a
! ADSXA : sym of a+, a+a => sym of a
! STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
! MXPNGAS : Max number of AS spaces ( program parameter )
! NOBPTS  : Number of orbitals per type and symmetry
! IOBPTS : base for orbitals of given type and symmetry
! IBORB  : Orbitals of given type and symmetry
! NSMOB,NSMST,NSMDX : Number of symmetries of orbitals, strings, double excitations
! MAXI   : Largest Number of ' spectator strings 'treated simultaneously
! MAXK   : Largest number of inner resolution strings treated at simult.
! IPACK  : Should we pack the density?
!
! ======
! Output
! ======
! RHO2, RHO2S, RHO2A : Updated density block
!
! =======
! Scratch
! =======
!
! SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
!              largest number of orbital pairs of given symmetries and
!              types.
! I1, XI1S,  : For holding creations/annihilations type and symmetry
!
! Jeppe Olsen, Fall of 96

use Para_Info, only: MyRank, nProcs
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NACOB, ISCSM, ISCTP, ICCSM, ICCTP, IGRP, NROW, NGAS, ISEL(NGAS), ICEL(NGAS), MXPOBS, ADSXA(MXPOBS,2*MXPOBS), &
                     NSMST, STSTSX(NSMST,NSMST), SXDXSX(2*MXPOBS,4*MXPOBS), MXPNGAS, NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*), MAXI, &
                     MAXK, I1(MAXK,*), NSMOB
real(kind=wp) :: RHO2(*), RHO2S(*), RHO2A(*), SB(*), CB(*), SSCR(*), CSCR(*), XI1S(MAXK,*), X(*), SCLFAC
logical(kind=iwp) :: IPACK
integer(kind=iwp) :: I, I1IK, I1JL, IBOT, IDXSM, IDXTP, IFIRST, IFRST, II12, IIK, IIKE, IIPART, IJL, IJLE, IKBOFF, IKOBSM, IKOFF, &
                     IKSM, IOFF, ISM, ITOP, ITP(256), ITYP, J, JFRST, JLBOFF, JLOBSM, JLOFF, JLSM, JOFF, JSM, JTP(256), JTYP, K, &
                     K12, KBOT, KEND, KFRST, KOFF, KSM, KTOP, KTP(256), KTYP, L, LDUMMY, LOFF, LSM, LTP(256), LTYP, NDXTP, NI, &
                     NIBTC, NIK, NJ, NJL, NK, NKBTC, NL, NONEW, NPART, NPARTSZ, NTEST
real(kind=wp) :: FACTOR

NTEST = 0
if (NTEST >= 1000) then
  write(u6,*)
  write(u6,*) ' ================='
  write(u6,*) ' GSBBD2A in action'
  write(u6,*) ' ================='
  write(u6,*)
  write(u6,*) ' Occupation of active left strings'
  call IWRTMA(ISEL,1,NGAS,1,NGAS)
  write(u6,*) ' Occupation of active Right strings'
  call IWRTMA(ICEL,1,NGAS,1,NGAS)
end if

IFRST = 1
JFRST = 1

!SVC: determine optimum number of partitions as the lowest multiple of
!     NPROCS that satisfies a block size smaller than MAXI:
NPART = 0
do
  NPART = NPART+NPROCS
  NPARTSZ = max(NROW-1,0)/NPART+1
  if (NPARTSZ <= MAXI) exit
end do

! Type of single excitations that connects the two column strings
call DXTYP_GAS(NDXTP,ITP,JTP,KTP,LTP,NGAS,ISEL,ICEL)
! Symmetry of Double excitation that connects IBSM and JBSM
! For general use : STSTSX => STSTDX
IDXSM = STSTSX(ISCSM,ICCSM)
if (IDXSM == 0) goto 2001
if (NTEST >= 1000) write(u6,*) ' ISCSM,ICCSM ',ISCSM,ICCSM
do IDXTP=1,NDXTP
  ITYP = ITP(IDXTP)
  JTYP = JTP(IDXTP)
  KTYP = KTP(IDXTP)
  LTYP = LTP(IDXTP)
  if (NTEST >= 1000) write(u6,*) ' ITYP JTYP KTYP LTYP ',ITYP,JTYP,KTYP,LTYP
  do IKOBSM=1,NSMOB
    JLOBSM = SXDXSX(IKOBSM,IDXSM)
    if (JLOBSM == 0) goto 1950
    ! types + symmetries defined => K strings are defined
    KFRST = 1
    ! Loop over of symmetry of i orbitals
    do ISM=1,NSMOB
      KSM = ADSXA(ISM,IKOBSM)
      NI = NOBPTS(ITYP,ISM)
      NK = NOBPTS(KTYP,KSM)
      if ((NI == 0) .or. (NK == 0)) goto 1940
      ! Loop over batches of j orbitals
      do JSM=1,NSMOB
        LSM = ADSXA(JSM,JLOBSM)
        NJ = NOBPTS(JTYP,JSM)
        NL = NOBPTS(LTYP,LSM)
        if ((NJ == 0) .or. (NL == 0)) goto 1930

        IOFF = IOBPTS(ITYP,ISM)
        JOFF = IOBPTS(JTYP,JSM)
        KOFF = IOBPTS(KTYP,KSM)
        LOFF = IOBPTS(LTYP,LSM)

        if (IOFF < KOFF) goto 1930
        if (JOFF < LOFF) goto 1930

        ! ==============================================================
        !                    Use N-2 projection method
        ! ==============================================================

        IFIRST = 1
        ! Loop over batches of I strings
        call SETVEC(X,Zero,NI*NJ*NK*NL)
        do IIPART=1+MYRANK,NPART,NPROCS
          IBOT = 1+(IIPART-1)*NPARTSZ
          ITOP = min(IBOT+NPARTSZ-1,NROW)
          NIBTC = ITOP-IBOT+1
          if (NIBTC <= 0) exit
          ! Loop over batches of intermediate strings
          KBOT = 1-MAXK
          KTOP = 0
1800      continue
          KBOT = KBOT+MAXK
          KTOP = KTOP+MAXK

          ! =======================================================
          !
          ! obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
          !
          ! =======================================================

          JLBOFF = 1
          if ((JSM == LSM) .and. (JTYP == LTYP)) then
            NJL = NJ*(NJ+1)/2
            JLSM = 1
          else
            NJL = NJ*NL
            JLSM = 0
          end if
          ! Obtain all double excitations from this group of K strings
          II12 = 1
          K12 = 1
          call ADADST_GAS(1,JSM,JTYP,NJ,1,LSM,LTYP,NL,ICCTP,ICCSM,IGRP,KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,II12,K12, &
                          SCLFAC)

          JFRST = 0
          KFRST = 0

          if (NKBTC == 0) goto 1930
          ! Loop over jl in TS classes
          J = 0
          L = 1

          do IJL=1,NJL
            call NXTIJ(J,L,NJ,NL,JLSM,NONEW)
            I1JL = (L-1)*NJ+J
            ! JAN28
            if (JLSM /= 0) then
              IJLE = J*(J-1)/2+L
            else
              IJLE = IJL
            end if
            ! JAN28
            ! CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
            !JLOFF = (JLBOFF-1+IJL-1)*NKBTC*NIBTC+1
            JLOFF = (JLBOFF-1+IJLE-1)*NKBTC*NIBTC+1
            if ((JLSM == 1) .and. (J == L)) then
              ! a+j a+j gives trivially zero
              call SETVEC(CSCR(JLOFF),Zero,NKBTC*NIBTC)
            else
              call MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,I1(1,I1JL),XI1S(1,I1JL))
            end if
          end do

          ! =======================================================
          !
          ! obtain sb(KB,IA,ik) = sum(IB)<KB!a kb a ib !IB>S(IA,IB)
          !
          ! =======================================================

          IKBOFF = 1
          if ((ISM == KSM) .and. (ITYP == KTYP)) then
            NIK = NI*(NI+1)/2
            IKSM = 1
          else
            NIK = NI*NK
            IKSM = 0
          end if
          ! Obtain all double excitations from this group of K strings
          II12 = 2
          K12 = 1
          if (IFRST == 1) KFRST = 1
          call ADADST_GAS(1,ISM,ITYP,NI,1,KSM,KTYP,NK,ISCTP,ISCSM,IGRP,KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,IFRST,KFRST,II12,K12,One)

          IFRST = 0
          KFRST = 0

          if (NKBTC == 0) goto 1930
          ! Loop over jl in TS classes
          I = 0
          K = 1

          do IIK=1,NIK
            call NXTIJ(I,K,NI,NK,IKSM,NONEW)
            I1IK = (K-1)*NI+I
            ! JAN28
            if (IKSM /= 0) then
              IIKE = I*(I-1)/2+K
            else
              IIKE = IIK
            end if
            ! JAN28
            ! SB(IA,KB,ik) = +/-S(IA,a+ka+iIA)
            !IKOFF = (IKBOFF-1+IIK-1)*NKBTC*NIBTC+1
            IKOFF = (IKBOFF-1+IIKE-1)*NKBTC*NIBTC+1
            if ((IKSM == 1) .and. (I == K)) then
              ! a+j a+j gives trivially zero
              call SETVEC(SSCR(IKOFF),Zero,NKBTC*NIBTC)
            else
              call MATCG(SB,SSCR(IKOFF),NROW,NIBTC,IBOT,NKBTC,I1(1,I1IK),XI1S(1,I1IK))
            end if
          end do

          ! ==================================================================
          !
          ! RHO2C(ik,jl)  = RHO2C(ik,jl) - sum(Ia,Kb)SB(Ia,Kb,ik)*CB(Ia,Kb,jl)
          !
          ! ==================================================================

          ! The minus ??
          !
          ! Well, the density matrices are constructed as
          !
          ! <I!a+i a+k aj al!> = -sum(K) <I!a+ia+k!K><J!aj al!K>, and
          ! the latter matrices are the ones we are constructing

          IOFF = IOBPTS(ITYP,ISM)
          JOFF = IOBPTS(JTYP,JSM)
          KOFF = IOBPTS(KTYP,KSM)
          LOFF = IOBPTS(LTYP,LSM)
          !if ((IOFF == 3) .and. (JOFF == 3) .and. (KOFF == 4) .and. (LOFF == 4)) NTEST = 5000
          LDUMMY = NKBTC*NIBTC
          if (NTEST >= 2000) then
            write(u6,*) ' CSCR matrix'
            call WRTMAT(CSCR,LDUMMY,NJL,LDUMMY,NJL)
            write(u6,*) ' SSCR matrix'
            call WRTMAT(SSCR,LDUMMY,NIK,LDUMMY,NIK)
          end if

          if (IFIRST == 1) then
            FACTOR = Zero
          else
            FACTOR = One
          end if
          LDUMMY = NKBTC*NIBTC
          !    MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP)
          call MATML7(X,SSCR,CSCR,NIK,NJL,LDUMMY,NIK,LDUMMY,NJL,FACTOR,-One,1)
          IFIRST = 0
          if (NTEST >= 2000) then
            write(u6,*) ' Updated X matrix IK,JL,IK,JL',NIK,NJL,NIK,NJL
            call WRTMAT(X,NIK,NJL,NIK,NJL)
          end if

          if (KEND == 0) goto 1800
          ! End of loop over partitionings of resolution strings
        end do
        ! Rho2(ik,jl) has been constructed for ik,jl belonging to
        ! Scatter out to density matrix
        IOFF = IOBPTS(ITYP,ISM)
        JOFF = IOBPTS(JTYP,JSM)
        KOFF = IOBPTS(KTYP,KSM)
        LOFF = IOBPTS(LTYP,LSM)
        !write(u6,*) 'I, J, K, L offsets & IPACK :',IOFF,JOFF,KOFF,LOFF,Ipack
        call ADTOR2(RHO2,RHO2S,RHO2A,X,1,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NACOB,IPACK)
        !    ADTOR2(RHO2,RHO2T,ITYPE,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)

        !write(u6,*) ' updated density matrix A',' norb = 4 ',norb
        !write(u6,*) ' offset ','IOFF,JOFF,KOFF,LOFF',IOFF,JOFF,KOFF,LOFF
        !call prsym(rho2s,NORB*(NORB+1)/2)

1930    continue
      end do
1940  continue
    end do
1950 continue
  end do
end do
2001 continue

end subroutine GSBBD2A_LUCIA
