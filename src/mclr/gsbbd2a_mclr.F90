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

!#define _DEBUGPRINT_
subroutine GSBBD2A_MCLR(RHO2,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NGAS,ISEL,ICEL,SB,CB,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR, &
                        CSCR,I1,XI1S,X,NSM)
! Contributions to two-electron density matrix from column excitations
!
! GAS version, '96, Jeppe Olsen
!
! =====
! Input
! =====
! RHO2        : two body density matrix to be updated
! NACOB       : Number of active orbitals
! ISCSM,ISCTP : Symmetry and type of sigma columns
! ICCSM,ICCTP : Symmetry and type of C     columns
! IGRP        : String group of columns
! NROW        : Number of rows in S and C block
! NGAS        : Number of active spaces
! ISEL        : Number of electrons per AS for S block
! ICEL        : Number of electrons per AS for C block
! CB          : Input C block
! MXPNGAS     : Max number of AS spaces (program parameter)
! NOBPTS      : Number of orbitals per type and symmetry
! IOBPTS      : base for orbitals of given type and symmetry
! IBORB       : Orbitals of given type and symmetry
! NSM         : Number of symmetries of orbitals
! MAXI        : Largest Number of "spectator strings" treated simultaneously
! MAXK        : Largest number of inner resolution strings treated at simult.
!
! ======
! Output
! ======
! RHO2 : Updated density block
!
! =======
! Scratch
! =======
! SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the largest
!              number of orbital pairs of given symmetries and types.
! I1, XI1S   : For holding creations/annihilations type and symmetry
!
! Jeppe Olsen, Fall of 96

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: RHO2(*)
integer(kind=iwp), intent(in) :: NACOB, ISCSM, ISCTP, ICCSM, ICCTP, IGRP, NROW, NGAS, ISEL(NGAS), ICEL(NGAS), MXPNGAS, &
                                 NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*), MAXI, MAXK, NSM
real(kind=wp), intent(in) :: SB(*), CB(*)
real(kind=wp), intent(_OUT_) :: SSCR(*), CSCR(*), XI1S(MAXK,*), X(*)
integer(kind=iwp), intent(_OUT_) :: I1(MAXK,*)
integer(kind=iwp) :: I, IBOT, IDXSM, IDXTP, IFIRST, IIK, IIKE, IJL, IJLE, IKBOFF, IKOBSM, IKOFF, IKSM, IOFF, IPART, ISM, ITOP, &
                     ITP(3*3), ITYP, J, JLBOFF, JLOBSM, JLOFF, JLSM, JOFF, JSM, JTP(3*3), JTYP, K, KBOT, KEND, KOFF, KSM, KTOP, &
                     KTP(3*3), KTYP, L, LDUMMY, LOFF, LSM, LTP(3*3), LTYP, NDXTP, NI, NIBTC, NIK, NJ, NJL, NK, NKBTC, nkStref, NL, &
                     NONEW, NPART
real(kind=wp) :: FACTOR

! Type of single excitations that connects the two column strings
call DXTYP_GAS(NDXTP,ITP,JTP,KTP,LTP,3,ISEL,ICEL)
! Symmetry of Double excitation that connects IBSM and JBSM
IDXSM = Mul(ISCSM,ICCSM)
if (IDXSM == 0) return
do IDXTP=1,NDXTP
  ITYP = ITP(IDXTP)
  JTYP = JTP(IDXTP)
  KTYP = KTP(IDXTP)
  LTYP = LTP(IDXTP)
  do IKOBSM=1,NSM
    JLOBSM = Mul(IKOBSM,IDXSM)
    if (JLOBSM == 0) cycle
    ! types + symmetries defined => K strings are defined
    !        KFRST = 1
    ! Loop over of symmetry of i orbitals
    do ISM=1,NSM
      KSM = Mul(ISM,IKOBSM)
      NI = NOBPTS(ITYP,ISM)
      NK = NOBPTS(KTYP,KSM)
      IOFF = IOBPTS(ITYP,ISM)
      KOFF = IOBPTS(KTYP,KSM)
      if ((NI == 0) .or. (NK == 0)) cycle
      ! Loop over batches of j orbitals
      outer: do JSM=1,NSM
        IFIRST = 1
        LSM = Mul(JSM,JLOBSM)
        NJ = NOBPTS(JTYP,JSM)
        NL = NOBPTS(LTYP,LSM)
        JOFF = IOBPTS(JTYP,JSM)
        LOFF = IOBPTS(LTYP,LSM)
        if (IOFF < KOFF) cycle
        if (JOFF < LOFF) cycle
        if ((NJ == 0) .or. (NL == 0)) cycle outer

        ! ==============================================================
        !         Use N-2 projection method
        ! ==============================================================

        ! Loop over batches of I strings
        NPART = NROW/MAXI
        if (NPART*MAXI /= NROW) NPART = NPART+1
#       ifdef _DEBUGPRINT_
        write(u6,*) ' NROW, MAXI NPART ',NROW,MAXI,NPART
#       endif
        do IPART=1,NPART
          IBOT = 1+(IPART-1)*MAXI
          ITOP = min(IBOT+MAXI-1,NROW)
          NIBTC = ITOP-IBOT+1
          ! Loop over batches of intermediate strings
          KBOT = 1-MAXK
          KTOP = 0
          do
            KBOT = KBOT+MAXK
            KTOP = KTOP+MAXK

            ! =======================================================
            !
            ! obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
            !
            ! =======================================================

            nkStref = maxk  ! ????????
            JLBOFF = 1
            if ((JSM == LSM) .and. (JTYP == LTYP)) then
              NJL = nTri_Elem(NJ)
              JLSM = 1
            else
              NJL = NJ*NL
              JLSM = 0
            end if
            ! Obtain all double excitations from this group of K strings
            lOFF = IOBPTS(lTYP,lSM)
            jOFF = IOBPTS(jTYP,jSM)
            call ADADST(JTYP,JSM,JOFF,NJ,LTYP,LSM,LOFF,NL,jlsm,ICCTP,ICCSM,IGRP,KBOT,KTOP,I1,XI1S,NKBTC,nkstref,KEND)

            if (NKBTC == 0) cycle outer
            ! Loop over jl in TS classes
            J = 0
            L = 1

            do IJL=1,NJL
              call NXTIJ(J,L,NJ,NL,JLSM,NONEW)
              !I1JL = (J-1)*NJ+L
              ! CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
              ! JAN28
              if (JLSM /= 0) then
                IJLE = nTri_Elem(J-1)+L
              else
                IJLE = IJL
              end if

              JLOFF = (JLBOFF-1+IJLE-1)*NKBTC*NIBTC+1
              if ((JLSM == 1) .and. (J == L)) then
                ! a+j a+j gives trivially zero
                CSCR(JLOFF:JLOFF+NKBTC*NIBTC-1) = Zero
              else
                !EAW BEGIN 970407
                !call MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,I1(1,I1JL),XI1S(1,I1JL))
                call MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,I1(:,IJL),XI1S(:,IJL))
                !EAW END
              end if
            end do

            ! =======================================================
            !
            ! obtain sb(KB,IA,ik) = sum(IB)<KB!a kb a ib !IB>S(IA,IB)
            !
            ! =======================================================

            IKBOFF = 1
            if ((ISM == KSM) .and. (ITYP == KTYP)) then
              NIK = nTri_Elem(NI)
              IKSM = 1
            else
              NIK = NI*NK
              IKSM = 0
            end if
            ! Obtain all double excitations from this group of K strings
            KOFF = IOBPTS(KTYP,KSM)
            iOFF = IOBPTS(iTYP,iSM)
            !if (IFRST == 1) KFRST = 1
            call ADADST(ITYP,ISM,IOFF,NI,KTYP,KSM,KOFF,NK,IKSM,ISCTP,ISCSM,IGRP,KBOT,KTOP,I1,XI1S,NKBTC,NKSTREF,KEND)

            if (NKBTC == 0) cycle outer
            ! Loop over jl in TS classes
            I = 0
            K = 1

            do IIK=1,NIK
              call NXTIJ(I,K,NI,NK,IKSM,NONEW)
              !I1IK = (K-1)*NI+I
              ! JAN28
              if (IKSM /= 0) then
                IIKE = nTri_Elem(I-1)+K
              else
                IIKE = IIK
              end if
              ! JAN28
              ! SB(IA,KB,ik) = +/-S(IA,a+ka+iIA)
              IKOFF = (IKBOFF-1+IIKE-1)*NKBTC*NIBTC+1
              if ((IKSM == 1) .and. (I == K)) then
                ! a+j a+j gives trivially zero
                SSCR(IKOFF:IKOFF+NKBTC*NIBTC-1) = Zero
              else
                !EAW-BEGIN 970407
                !call MATCG(SB,SSCR(IKOFF),NROW,NIBTC,IBOT,NKBTC,I1(1,I1IK),XI1S(1,I1IK))
                call MATCG(SB,SSCR(IKOFF),NROW,NIBTC,IBOT,NKBTC,I1(:,IIK),XI1S(:,IIK))
                !EAW-END
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

            LDUMMY = NKBTC*NIBTC

            if (IFIRST == 1) then
              FACTOR = Zero
            else
              FACTOR = One
            end if
            LDUMMY = NKBTC*NIBTC
            call DGEMM_('T','N',NIK,NJL,LDUMMY,-One,SSCR,max(1,LDUMMY),CSCR,max(1,LDUMMY),FACTOR,X,max(1,NIK))
            IFIRST = 0

            if (KEND /= 0) exit
          end do
          ! End of loop over partitionings of resolution strings
        end do
        ! Rho2(ik,jl) has been constructed for ik,jl belonging to
        ! Scatter out to density matrix
        IOFF = IOBPTS(ITYP,ISM)
        JOFF = IOBPTS(JTYP,JSM)
        KOFF = IOBPTS(KTYP,KSM)
        LOFF = IOBPTS(LTYP,LSM)
        call ADTOR2_MCLR(RHO2,X,1,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NACOB)
      end do outer
    end do
  end do
end do

end subroutine GSBBD2A_MCLR
