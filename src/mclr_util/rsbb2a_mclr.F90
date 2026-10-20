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

subroutine RSBB2A_MCLR(ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,ISEL1,ISEL3,ICEL1,ICEL3,SB,CB,NTSOB,IBTSOB,MAXI,MAXK,SSCR,CSCR,I1, &
                       XI1S,XINT,NSM,SGN,NOPART,TimeDep,ieaw)
! two electron excitations on column strings
!
! =====
! Input
! =====
!
! ISCSM,ISCTP : Symmetry and type of sigma columns
! ICCSM,ICCTP : Symmetry and type of C     columns
! IGRP        : String group of columns
! NROW        : Number of rows in S and C block
! ISEL1(3)    : Number of electrons in RAS1(3) for S block
! ICEL1(3)    : Number of electrons in RAS1(3) for C block
! CB          : Input C block
! NTSOB       : Number of orbitals per type and symmetry
! IBTSOB      : base for orbitals of given type and symmetry
! IBORB       : Orbitals of given type and symmetry
! NSM         : Number of symmetries of orbitals, single excitations
! MAXI        : Largest Number of "spectator strings" treated simultaneously
! MAXK        : Largest number of inner resolution strings treated at simult.
!
! ======
! Output
! ======
! SB : updated sigma block
!
! =======
! Scratch
! =======
! SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the largest
!              number of orbital pairs of given symmetries and types.
! I1, XI1S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! XINT       : Space for two electron integrals
!
! Jeppe Olsen, Winter of 1991

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Str_Info, only: NOCTYP, STR
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ISCSM, ISCTP, ICCSM, ICCTP, IGRP, NROW, ISEL1, ISEL3, ICEL1, ICEL3, NTSOB(3,*), IBTSOB(3,*), &
                                 MAXI, MAXK, NSM, NOPART, ieaw
real(kind=wp), intent(inout) :: SB(*)
real(kind=wp), intent(in) :: CB(*), SGN
real(kind=wp), intent(_OUT_) :: SSCR(*), CSCR(*), XI1S(*), XINT(*)
integer(kind=iwp), intent(_OUT_) :: I1(*)
logical(kind=iwp), intent(in) :: TimeDep
integer(kind=iwp) :: I, IBOT, IDXSM, IDXTYP, IFIRST, IJL, IK, IKOFF, IKPSM, IKSM, IOFF, IPART, ISBOFF, ISM, ITOP, ITP(36), ITYP, &
                     IXCHNG, J, JLOFF, JLOFF2, JLPSM, JLSM, JOFF, JSM, JTP(36), JTYP, K, K1GRP, K1SM, K1TP, K2GRP, K2SM, K2TP, &
                     KBOT, KEND, KOFF, KSM, KTOP, KTP(36), KTYP, L, LIKB, LOFF, LSM, LTP(36), LTYP, NDXTYP, NI, NIBTC, NIK, NJ, &
                     NJL, NK, NKBTC, NKSTR, NKSTREF, NL, NONEW, NPART
real(kind=wp) :: FACTORAB, FACTORC

! Types of DX that connects the two strings

!write(u6,*) 'ieaw in rsbb2a_mclr',ieaw
IDXSM = Mul(ISCSM,ICCSM)
if (IDXSM == 0) return
call DXTYP(NDXTYP,ITP,JTP,KTP,LTP,ISEL1,ISEL3,ICEL1,ICEL3)
do IDXTYP=1,NDXTYP
  ITYP = ITP(IDXTYP)
  JTYP = JTP(IDXTYP)
  KTYP = KTP(IDXTYP)
  LTYP = LTP(IDXTYP)
  ! Type of intermediate strings
  call NEWTYP_MCLR(IGRP,ICCTP,[1],[JTYP],1,K1GRP,K1TP)
  call NEWTYP_MCLR(K1GRP,K1TP,[1],[LTYP],1,K2GRP,K2TP)
  if (K2TP <= 0) cycle
  ! Symmetry of allowed Double excitation,loop over excitations
  do IKSM=1,NSM
    JLSM = Mul(IKSM,IDXSM)
    if (JLSM == 0) cycle
    do ISM=1,NSM
      ! Works only for D2h
      KSM = Mul(ISM,IKSM)
      if (KSM == 0) cycle
      ! sym of intermediate strings

      K1SM = Mul(ISM,ISCSM)
      K2SM = Mul(KSM,K1SM)
      ! Intermediate K strings are of type K2TP and Sym K2Sm

      NKSTR = Str(K2GRP)%NSTSO((K2SM-1)*NOCTYP(K2GRP)+K2TP)
      if (NOPART == 0) then
        NKSTREF = min(NKSTR,MAXK)
      else
        NKSTREF = NKSTR
      end if
      IOFF = IBTSOB(ITYP,ISM)
      KOFF = IBTSOB(KTYP,KSM)
      NI = NTSOB(ITYP,ISM)
      NK = NTSOB(KTYP,KSM)
      if (KOFF > IOFF) cycle
      if ((ISM == KSM) .and. (ITYP == KTYP)) then
        IKPSM = 1
        NIK = nTri_Elem(NI)
      else
        IKPSM = 0
        NIK = NI*NK
      end if
      if (NOPART == 1) SSCR(1:NKSTREF*NROW*NIK) = Zero
      do JSM=1,NSM
        LSM = Mul(JSM,JLSM)
        if (LSM == 0) cycle
        JOFF = IBTSOB(JTYP,JSM)
        LOFF = IBTSOB(LTYP,LSM)
        if (LOFF > JOFF) cycle
        NJ = NTSOB(JTYP,JSM)
        NL = NTSOB(LTYP,LSM)
        if ((JSM == LSM) .and. (JTYP == LTYP)) then
          JLPSM = 1
          NJL = nTri_Elem(NJ)
        else
          JLPSM = 0
          NJL = NJ*NL
        end if
        if ((NI == 0) .or. (NJ == 0) .or. (NK == 0) .or. (NL == 0)) cycle
        IFIRST = 1
        ! Loop over batches of I strings
        if (NOPART == 0) then
          NPART = NROW/MAXI
          if (NPART*MAXI /= NROW) NPART = NPART+1
        else
          NPART = 1
        end if
        outer: do IPART=1,NPART
          IBOT = 1+(IPART-1)*MAXI
          if (NOPART == 0) then
            ITOP = min(IBOT+MAXI-1,NROW)
          else
            ITOP = NROW
          end if
          NIBTC = ITOP-IBOT+1

          ! Loop over batches of intermediate strings

          KBOT = 1-MAXK
          KTOP = 0
          do
            if (NOPART == 0) then
              KBOT = KBOT+MAXK
              KTOP = KTOP+MAXK
            else
              KBOT = 1
              KTOP = NKSTREF
            end if

            ! obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)

            ! Generate arrays a+ja+l |kstr> for kstr in current interval
            call ADADST(JTYP,JSM,JOFF,NJ,LTYP,LSM,LOFF,NL,JLPSM,ICCTP,ICCSM,IGRP,KBOT,KTOP,I1,XI1S,NKBTC,NKSTREF,KEND)
            if (NKBTC == 0) exit outer
            J = 0
            L = 1
            do IJL=1,NJL
              call NXTIJ(J,L,NJ,NL,JLPSM,NONEW)
              ! CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
              JLOFF = (IJL-1)*NKBTC*NIBTC+1
              JLOFF2 = (IJL-1)*NKSTREF+1
              call MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,I1(JLOFF2),XI1S(JLOFF2))
            end do
            !=============================================
            ! SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
            !==============================================
            ! Obtain two electron integrals (ij!kl)-(il!kj)
            if (IFIRST == 1) then
              IXCHNG = 1
              !write(u6,*) 'TimeDep in rsbb2a_mclr is:',TimeDep
              if (TimeDep) then
                call GETINT_td(XINT,ITYP,ISM,JTYP,JSM,KTYP,KSM,LTYP,LSM,IKPSM,JLPSM,4,ieaw)
              else
                !write(u6,*) 'I call getint not getint_td'
                call GETINT_MCLR(XINT,ITYP,ISM,JTYP,JSM,KTYP,KSM,LTYP,LSM,IXCHNG,IKPSM,JLPSM,0,0)
              end if
            end if
            IFIRST = 0
            ! and now, to the work
            LIKB = NIBTC*NKBTC
            if (NOPART == 1) then
              FACTORC = One
            else
              FACTORC = Zero
            end if
            FACTORAB = One
            call DGEMM_('N','T',LIKB,NIK,NJL,FACTORAB,CSCR,LIKB,XINT,NIK,FACTORC,SSCR,LIKB)
            ! ============================
            ! Loop over ik and scatter out
            ! ============================

            ! Generate arrays a+i a+k !kstr>
            if (NOPART == 0) then
              call ADADST(ITYP,ISM,IOFF,NI,KTYP,KSM,KOFF,NK,IKPSM,ISCTP,ISCSM,IGRP,KBOT,KTOP,I1,XI1S,NKBTC,NKSTREF,KEND)
              I = 0
              K = 1
              do IK=1,NIK
                call NXTIJ(I,K,NI,NK,IKPSM,NONEW)
                ISBOFF = 1+(IK-1)*NIBTC*NKBTC
                IKOFF = (IK-1)*NKSTREF+1
                if (SGN == -One) XI1S(IKOFF:IKOFF+NKSTREF-1) = -XI1S(IKOFF:IKOFF+NKSTREF-1)
                call MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,NKBTC,I1(IKOFF),XI1S(IKOFF))
              end do
            end if

            if ((KEND /= 0) .or. (NOPART /= 0)) exit
          end do
          ! End of loop over partitionings of resolution strings
        end do outer
      end do
      ! End of loop over JSM
      if (NOPART == 1) then
        ! processing of ik dependent terms is done indepDently of jl dependent terms:
        call ADADST(ITYP,ISM,IOFF,NI,KTYP,KSM,KOFF,NK,IKPSM,ISCTP,ISCSM,IGRP,KBOT,KTOP,I1,XI1S,NKBTC,NKSTREF,KEND)
        I = 0
        K = 1
        do IK=1,NIK
          call NXTIJ(I,K,NI,NK,IKPSM,NONEW)
          ISBOFF = 1+(IK-1)*NIBTC*NKBTC
          IKOFF = (IK-1)*NKSTREF+1
          ! Well, someplace the minus must come in
          if (SGN == -One) XI1S(IKOFF:IKOFF+NKSTREF-1) = -XI1S(IKOFF:IKOFF+NKSTREF-1)
          call MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,NKBTC,I1(IKOFF),XI1S(IKOFF))
        end do
      end if

    end do
  end do
end do

end subroutine RSBB2A_MCLR
