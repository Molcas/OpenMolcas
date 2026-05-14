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

subroutine RSBB2A(ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NGAS,ISOC,ICOC,SB,CB,NOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,XINT,NSMOB,NSMST, &
                  SCLFAC,IPHGAS)
! SUBROUTINE RSBB2A --> 46
!
! two electron excitations on column strings
!
! =====
! Input
! =====
!
! ISCSM,ISCTP : Symmetry and type of sigma columns
! ICCSM,ICCTP : Symmetry and type of C     columns
! IGRP : String group of columns
! NROW : Number of rows in S and C block
! ISEL1(3) : Number of electrons in RAS1(3) for S block
! ICEL1(3) : Number of electrons in RAS1(3) for C block
! CB   : Input C block
! NTSOB  : Number of orbitals per type and symmetry
! IBTSOB : base for orbitals of given type and symmetry
! IBORB  : Orbitals of given type and symmetry
! NSMOB,NSMST : Number of symmetries of orbitals, strings
! MAXI   : Largest number of "spectator strings" treated simultaneously
! MAXK   : Largest number of inner resolution strings treated at simult.
!
! ======
! Output
! ======
!
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
! XINT : Space for two electron integrals
!
! Jeppe Olsen, Winter of 1991

use Symmetry_Info, only: Mul
use Index_Functions, only: nTri_Elem
use Para_Info, only: MyRank, nProcs
use lucia_data, only: MXPNGAS, MXPTSOB
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ISCSM, ISCTP, ICCSM, ICCTP, IGRP, NROW, NGAS, ISOC(NGAS), ICOC(NGAS), NSMST, NOBPTS(MXPNGAS,*), &
                                 MAXI, MAXK, NSMOB, IPHGAS(NGAS)
real(kind=wp), intent(inout) :: SB(*)
real(kind=wp), intent(in) :: CB(*), SCLFAC
real(kind=wp), intent(_OUT_) :: SSCR(*), CSCR(*), XI1S(MAXK,*), XINT(*)
integer(kind=iwp), intent(_OUT_) :: I1(MAXK,*)
integer(kind=iwp) :: I, I1JL, I4_AC(4), I4_REO(4), I4_TP(4), IAC, IBOT, ICOUL, IDXSM, IDXTYP, IFIRST, IFRST, II12, IIPART, IJKL, &
                     IJL, IK, IKBOFF, IKBT(3,8), IKBTC, IKOBSM, IKOFF, IKPAIR, IKSM, IKSMBT(2,8), IOBSM, ISBOFF, ISCR(4), ISM, &
                     ISM_ORIG, ITOP, ITP(256), ITPSM_ORIG, ITYP, ITYP_ORIG, IXCHNG, J, JAC, JFRST, JL, JLBOFF, JLBT(3,8), JLBTC, &
                     JLOBSM, JLOFF, JLPAIR, JLSM, JLSMBT(2,8), JSM, JSM_ORIG, JTP(256), JTPSM_ORIG, JTYP, JTYP_ORIG, K, K12, KAC, &
                     KBOT, KEND, KFRST, KSM, KSM_ORIG, KTOP, KTP(256), KTPSM_ORIG, KTYP, KTYP_ORIG, L, LAC, LENGTH, LIKB, LSM, &
                     LSM_ORIG, LTP(256), LTPSM_ORIG, LTYP, LTYP_ORIG, MI, MJ, MK, ML, MXPAIR, NBLK, NBLKT, NDXTYP, NI, NIBTC, &
                     NIJKL1, NIK, NIKBT, NIKT, NJ, NJL, NJLBT, NJLT, NK, NKBTC, NL, NONEW, NPART, NPARTSZ
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: II, JIKBT, JJLBT
#endif
real(kind=wp) :: FACTORAB, FACTORC, FACX
real(kind=wp), allocatable :: SCR(:)

call mma_allocate(SCR,MXPTSOB**4,Label='SCR')
#ifdef _DEBUGPRINT_
write(u6,*) ' ==============='
write(u6,*) ' RSBB2A speaking'
write(u6,*) ' ==============='
write(u6,*) ' ISOC and ICOC :'
call IWRTMA(ISOC,1,NGAS,1,NGAS)
call IWRTMA(ICOC,1,NGAS,1,NGAS)
#endif
IFRST = 1
JFRST = 1
FACX = Zero

! Types of DX that connects the two strings

IDXSM = Mul(ISCSM,ICCSM)
if (IDXSM /= 0) then
  ! Connecting double excitations
  call DXTYP2_GAS(NDXTYP,ITP,JTP,KTP,LTP,NGAS,ISOC,ICOC,IPHGAS)
  do IDXTYP=1,NDXTYP
    ITYP = ITP(IDXTYP)
    JTYP = JTP(IDXTYP)
    KTYP = KTP(IDXTYP)
    LTYP = LTP(IDXTYP)
    ! Is this combination of types allowed
    !IJKL_ACT = 1
    !if (IJKL_ACT == 0) cycle

    !write(u6,*) ' test inserted in RSBB2A'
    !NPTOT = 0
    !if (ITYP == 3) NPTOT = NPTOT+1
    !if (JTYP == 3) NPTOT = NPTOT+1
    !if (KTYP == 3) NPTOT = NPTOT+1
    !if (LTYP == 3) NPTOT = NPTOT+1
    !if (NPTOT == 3) cycle
    ITYP_ORIG = ITYP
    JTYP_ORIG = JTYP
    KTYP_ORIG = KTYP
    LTYP_ORIG = LTYP

#   ifdef _DEBUGPRINT_
    write(u6,*) ' ITYP_ORIG, JTYP_ORIG, KTYP_ORIG, LTYP_ORIG',ITYP_ORIG,JTYP_ORIG,KTYP_ORIG,LTYP_ORIG
#   endif
    NIJKL1 = 0
    if (ITYP == 1) NIJKL1 = NIJKL1+1
    if (JTYP == 1) NIJKL1 = NIJKL1+1
    if (KTYP == 1) NIJKL1 = NIJKL1+1
    if (LTYP == 1) NIJKL1 = NIJKL1+1
    ! Optimal ordering of operators
    I4_AC(1) = 2
    I4_AC(2) = 2
    I4_AC(3) = 1
    I4_AC(4) = 1
    I4_TP(1) = ITYP
    I4_TP(2) = KTYP
    I4_TP(3) = LTYP
    I4_TP(4) = JTYP
    !if (IUSE_PH == 1) then
    !  NOP = 4
    !  call ALG_ROUTERX(ISOC,JSOC,NOP,I4_TP,I4_AC,I4_REO,SIGN4)
    !else
    I4_REO(:) = [1,2,3,4]
    !end if
    ! Type of operators : TP and AC
    do IJKL=1,4
      !ISCR(I4_REO(IJKL)) = I4_TP(IJKL)
      ISCR(IJKL) = I4_TP(I4_REO(IJKL))
    end do
    do IJKL=1,4
      I4_TP(IJKL) = ISCR(IJKL)
    end do
    ITYP = I4_TP(1)
    KTYP = I4_TP(2)
    LTYP = I4_TP(3)
    JTYP = I4_TP(4)
    do IJKL=1,4
      ISCR(IJKL) = I4_AC(I4_REO(IJKL))
    end do
    do IJKL=1,4
      I4_AC(IJKL) = ISCR(IJKL)
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) ' I4_AC, IT_TP  defined'
    write(u6,*) ' I4_AC, I4_TP'
    call IWRTMA(I4_AC,1,4,1,4)
    call IWRTMA(I4_TP,1,4,1,4)
#   endif

    !SVC: determine optimum number of partitions as the lowest multiple of
    !     NPROCS that satisfies a block size smaller than MAXI:
    NPART = 0
    do
      NPART = NPART+NPROCS
      NPARTSZ = max(NROW-1,0)/NPART+1
      if (NPARTSZ <= MAXI) exit
    end do

    if (I4_AC(1) == I4_AC(2)) then

      ! a+ a+ a a or a a a+ a+
      ! Largest possible number of orbital pairs
      MI = 0
      MJ = 0
      MK = 0
      ML = 0
      do IOBSM=1,NSMST
        MI = max(MI,NOBPTS(ITYP,IOBSM))
        MJ = max(MJ,NOBPTS(JTYP,IOBSM))
        MK = max(MK,NOBPTS(KTYP,IOBSM))
        ML = max(ML,NOBPTS(LTYP,IOBSM))
      end do
      MXPAIR = max(MI*MK,MJ*ML)
      ! Largest posssible
      ! Symmetry of allowed Double excitation,loop over excitations
      do IKOBSM=1,NSMOB
        JLOBSM = Mul(IKOBSM,IDXSM)
#       ifdef _DEBUGPRINT_
        write(u6,*) ' IKOBSM,JLOBSM',IKOBSM,JLOBSM
#       endif
        if (JLOBSM == 0) cycle
        ! types + symmetries defined => K strings are defined
        KFRST = 1

        ! Number of batchs of symmetry pairs IK

        LENGTH = 0
        NIKBT = 0
        NBLK = 0
        NBLKT = 0
        do ISM=1,NSMOB
          KSM = Mul(ISM,IKOBSM)
          NI = NOBPTS(ITYP,ISM)
          NK = NOBPTS(KTYP,KSM)
#         ifdef _DEBUGPRINT_
          write(u6,*) ' NI, NK',NI,NK
#         endif

          if ((ISM == KSM) .and. (ITYP == KTYP)) then
            NIK = nTri_Elem(NI)
          else if ((ITYP > KTYP) .or. ((ITYP == KTYP) .and. (ISM > KSM))) then
            NIK = NI*NK
          else
            NIK = 0
          end if
          if (NIK /= 0) then
            NBLKT = NBLKT+1
            if (LENGTH+NIK > MXPAIR) then
              ! The present batch is complete
              NIKBT = NIKBT+1
              IKBT(1,NIKBT) = NBLKT-NBLK
              IKBT(2,NIKBT) = NBLK
              IKBT(3,NIKBT) = LENGTH
              LENGTH = 0
              NBLK = 0
            end if
            NBLK = NBLK+1
            LENGTH = LENGTH+NIK
            IKSMBT(1,NBLKT) = ISM
            IKSMBT(2,NBLKT) = KSM
          end if
        end do
        ! The last batch
        if (NBLK /= 0) then
          NIKBT = NIKBT+1
          IKBT(1,NIKBT) = NBLKT-NBLK+1
          IKBT(2,NIKBT) = NBLK
          IKBT(3,NIKBT) = LENGTH
        end if

#       ifdef _DEBUGPRINT_
        write(u6,*) ' ITYP, KTYP, IKOBSM,  NIKBT = ',ITYP,KTYP,IKOBSM,NIKBT
        write(u6,*) ' IKBT : Offset, number, length'
        do JIKBT=1,NIKBT
          write(u6,'(3i3)') (IKBT(II,JIKBT),II=1,3)
        end do
        write(u6,*) ' IKSMBT'
        call IWRTMA(IKSMBT,2,NBLKT,2,8)
#       endif

        ! Number of batches of symmetry pairs JL

        LENGTH = 0
        NJLBT = 0
        NBLK = 0
        NBLKT = 0
        do JSM=1,NSMOB
          LSM = Mul(JSM,JLOBSM)
          NJ = NOBPTS(JTYP,JSM)
          NL = NOBPTS(LTYP,LSM)

          if ((JSM == LSM) .and. (JTYP == LTYP)) then
            NJL = nTri_Elem(NJ)
          else if ((JTYP > LTYP) .or. ((JTYP == LTYP) .and. (JSM > LSM))) then
            NJL = NJ*NL
          else
            NJL = 0
          end if
          if (NJL /= 0) then
            NBLKT = NBLKT+1
            if (LENGTH+NJL > MXPAIR) then
              ! The present batch is complete
              NJLBT = NJLBT+1
              JLBT(1,NJLBT) = NBLKT-NBLK
              JLBT(2,NJLBT) = NBLK
              JLBT(3,NJLBT) = LENGTH
              LENGTH = 0
              NBLK = 0
            end if
            NBLK = NBLK+1
            LENGTH = LENGTH+NJL
            JLSMBT(1,NBLKT) = JSM
            JLSMBT(2,NBLKT) = LSM
          end if
        end do
        ! The last batch
        if (NBLK /= 0) then
          NJLBT = NJLBT+1
          JLBT(1,NJLBT) = NBLKT-NBLK+1
          JLBT(2,NJLBT) = NBLK
          JLBT(3,NJLBT) = LENGTH
        end if

#       ifdef _DEBUGPRINT_
        write(u6,*) ' JTYP, LTYP, JLOBSM,  NJLBT = ',JTYP,LTYP,JLOBSM,NJLBT
        write(u6,*) ' JLBT : Offset, number, length'
        do JJLBT=1,NJLBT
          write(u6,'(3i3)') (JLBT(II,JJLBT),II=1,3)
        end do
        write(u6,*) ' JLSMBT'
        call IWRTMA(JLSMBT,2,NBLKT,2,8)
#       endif

        ! Loop over batches of IK strings
        do IKBTC=1,NIKBT
#         ifdef _DEBUGPRINT_
          write(u6,*) ' IKBTC = ',IKBTC
#         endif
          ! Loop over batches of JL strings
          do JLBTC=1,NJLBT
            IFIRST = 1
            ! Loop over batches of I strings
            do IIPART=1+MYRANK,NPART,NPROCS
              IBOT = 1+(IIPART-1)*NPARTSZ
              ITOP = min(IBOT+NPARTSZ-1,NROW)
              NIBTC = ITOP-IBOT+1
              if (NIBTC <= 0) exit
              ! Loop over batches of intermediate strings
              KBOT = 1-MAXK
              KTOP = 0
              outer: do
                KBOT = KBOT+MAXK
                KTOP = KTOP+MAXK

                JLBOFF = 1
                NJLT = JLBT(3,JLBTC)
                do JLPAIR=1,JLBT(2,JLBTC)
                  JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                  LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)
                  NJ = NOBPTS(JTYP,JSM)
                  NL = NOBPTS(LTYP,LSM)
                  if ((JSM == LSM) .and. (JTYP == LTYP)) then
                    NJL = nTri_Elem(NJ)
                    JLSM = 1
                  else
                    NJL = NJ*NL
                    JLSM = 0
                  end if

                  ! obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)

                  ! Obtain all double excitations from this group of K strings
                  II12 = 1
                  K12 = 1
                  !write(u6,*) ' Before ADAADAST'
                  ! Creation / annihilation maps, conjugated of above
                  if (I4_AC(4) == 1) then
                    JAC = 2
                  else
                    JAC = 1
                  end if
                  if (I4_AC(3) == 1) then
                    LAC = 2
                  else
                    LAC = 1
                  end if
                  call ADAADAST_GAS(1,JSM,JTYP,NJ,JAC,1,LSM,LTYP,NL,LAC,ICCTP,ICCSM,IGRP,KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,JFRST, &
                                    KFRST,II12,K12,SCLFAC)

                  JFRST = 0
                  KFRST = 0

                  if  (NKBTC == 0) exit outer
                  ! Loop over jl in TS classes
                  J = 0
                  L = 1

                  do IJL=1,NJL
                    call NXTIJ(J,L,NJ,NL,JLSM,NONEW)
                    I1JL = (L-1)*NJ+J
                    ! CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
                    JLOFF = (JLBOFF-1+IJL-1)*NKBTC*NIBTC+1
                    if ((JLSM == 1) .and. (J == L)) then
                      ! a+j a+j gives trivially zero
                      CSCR(JLOFF:JLOFF+NKBTC*NIBTC-1) = Zero
                    else
                      call MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,I1(:,I1JL),XI1S(:,I1JL))
                    end if
                  end do

                  JLBOFF = JLBOFF+NJL
                end do

                ! (End of loop over jlpair in batch)
                !=============================================
                ! SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
                !=============================================
                ! Obtain two electron integrals xint(ik,jl) = (ij!kl)-(il!kj)
                if (IFIRST == 1) then
                  IXCHNG = 1
                  ! Obtain integrals in ik batch
                  NIKT = IKBT(3,IKBTC)
                  NJLT = JLBT(3,JLBTC)
                  JLOFF = 1
                  do JLPAIR=1,JLBT(2,JLBTC)
                    IKOFF = 1
                    do IKPAIR=1,IKBT(2,IKBTC)

                      ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                      KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                      JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                      LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)

                      if ((ISM == KSM) .and. (ITYP == KTYP)) then
                        IKSM = 1
                        NIK = nTri_Elem(NOBPTS(ITYP,ISM))
                      else
                        IKSM = 0
                        NIK = NOBPTS(ITYP,ISM)*NOBPTS(KTYP,KSM)
                      end if

                      if ((JSM == LSM) .and. (JTYP == LTYP)) then
                        JLSM = 1
                        NJL = nTri_Elem(NOBPTS(JTYP,JSM))
                      else
                        JLSM = 0
                        NJL = NOBPTS(JTYP,JSM)*NOBPTS(LTYP,LSM)
                      end if
                      ! ===============================================================
                      ! Required form of integrals : Coulomb - Exchange of just Coulomb
                      ! ===============================================================
                      ICOUL = 0
                      ! Use coulomb - exchange
                      IXCHNG = 1
                      ! fetch integrals
                      ! Full conjugation symmetry, do do not worry
                      call GETINT(SCR,ITYP,ISM,JTYP,JSM,KTYP,KSM,LTYP,LSM,IXCHNG,IKSM,JLSM,ICOUL)
                      ! End if similarity transformed Hamiltonian is used
                      do JL=1,NJL
                        XINT(IKOFF+(JLOFF-1+JL-1)*NIKT:IKOFF+(JLOFF-1+JL-1)*NIKT+NIK-1) = SCR((JL-1)*NIK+1:JL*NIK)
                      end do
                      IKOFF = IKOFF+NIK
                    end do
                    JLOFF = JLOFF+NJL
                  end do
                end if
                ! End if integrals should be fetched
                IFIRST = 0
                ! and now, to the work
                LIKB = NIBTC*NKBTC
#               ifdef _DEBUGPRINT_
                write(u6,*) ' Integral block'
                call WRTMAT(XINT,NIKT,NJLT,NIKT,NJLT)
                write(u6,*) ' CSCR matrix'
                call WRTMAT(CSCR,LIKB,NJLT,LIKB,NJLT)
#               endif

                !!MXACIJO = MXACIJ
                !MXACIJ = MAX(MXACIJ,LIKB*NJLT,LIKB*NIKT)
                !!if (MXACIJ > MXACIJO) then
                !!  write(u6,*) ' New max MXACIJ = ', MXACIJ
                !!  write(u6,*) ' ISCTP,ICCTP', ISCTP,ICCTP
                !!  write(u6,*) ' ITYP,JTYP,KTYP,LTYP',ITYP,JTYP,KTYP,LTYP
                !!  write(u6,*) 'NIJT, NJLT, NIBTC NKBTC',NIJT,NJLT,NIBTC,NKBTC
                !!end if

                FACTORC = Zero
                FACTORAB = One
                call MATML7(SSCR,CSCR,XINT,LIKB,NIKT,LIKB,NJLT,NIKT,NJLT,FACTORC,FACTORAB,2)
#               ifdef _DEBUGPRINT_
                write(u6,*) ' SSCR matrix'
                call WRTMAT(SSCR,LIKB,NIKT,LIKB,NIKT)
#               endif
                ! ============================
                ! Loop over ik and scatter out
                ! ============================
                ! Generate double excitations from K strings
                ! I strings connected with K strings in batch <I!a+i a+k!K)
                II12 = 2

                IKBOFF = 1
                do IKPAIR=1,IKBT(2,IKBTC)
                  ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                  KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                  NI = NOBPTS(ITYP,ISM)
                  NK = NOBPTS(KTYP,KSM)
                  if ((ISM == KSM) .and. (ITYP == KTYP)) then
                    NIK = nTri_Elem(NI)
                    IKSM = 1
                  else
                    NIK = NI*NK
                    IKSM = 0
                  end if
                  if (IFRST == 1) KFRST = 1

                  IAC = I4_AC(1)
                  KAC = I4_AC(2)

                  call ADAADAST_GAS(1,ISM,ITYP,NI,IAC,1,KSM,KTYP,NK,KAC,ISCTP,ISCSM,IGRP,KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,IFRST, &
                                    KFRST,II12,K12,One)

                  IFRST = 0
                  KFRST = 0

                  I = 0
                  K = 1
                  do IK=1,NIK
                    call NXTIJ(I,K,NI,NK,IKSM,NONEW)
                    IKOFF = (K-1)*NI+I
                    ISBOFF = 1+(IKBOFF-1+IK-1)*NIBTC*NKBTC
                    if ((IKSM == 1) .and. (I == k)) then
                      ! a+ i a+i gives trivially zero
                    else
                      call MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,NKBTC,I1(:,IKOFF),XI1S(:,IKOFF))
                    end if
                  end do
                  IKBOFF = IKBOFF+NIK

                end do
                ! End of loop over IKPAIRS in batch

                if (KEND /= 0) exit outer
              end do outer
              ! End of loop over partitionings of resolution strings
            end do
            ! End of loop over partitionings of I strings
          end do
          ! End of loop over batches of JL
        end do
        ! End of loop over batches of IK
      end do
      ! End of loop over IKOBSM

    else if (.not. (I4_AC(1) == I4_AC(2))) then

      ! Three types of operators :
      ! a+ a  a+ a
      ! a+ a  a  a+
      ! a  a+ a+ a
      !
      ! The first end up with
      ! -a+ i ak a+l aj X2(ik,jl)
      !
      ! Number two and three end up with
      ! -a i a k a l aj XC(ik,jl)  (In coulomb form)

      JLSM = 0
      IKSM = 0
      ! Symmetry of allowed Double excitation,loop over excitations
      do IKOBSM=1,NSMOB
        JLOBSM = Mul(IKOBSM,IDXSM)
        if (JLOBSM == 0) cycle
        ! types + symmetries defined => K strings are defined
        KFRST = 1
        do ISM=1,NSMOB
          KSM = Mul(ISM,IKOBSM)
          do JSM=1,NSMOB
            LSM = Mul(JSM,JLOBSM)
#           ifdef _DEBUGPRINT_
            write(u6,*) ' ISM KSM LSM JSM',ISM,KSM,LSM,JSM
#           endif
            ISCR(I4_REO(1)) = ISM
            ISCR(I4_REO(2)) = KSM
            ISCR(I4_REO(3)) = LSM
            ISCR(I4_REO(4)) = JSM

            ISM_ORIG = ISCR(1)
            KSM_ORIG = ISCR(2)
            LSM_ORIG = ISCR(3)
            JSM_ORIG = ISCR(4)

            !do ISM_ORIG=1,NSMOB
            !  KSM_ORIG = Mul(ISM_ORIG,IKOBSM)
            !  do JSM_ORIG=1,NSMOB
            !    LSM_ORIG = Mul(JSM_ORIG,JLOBSM)
            !
            !    ISCR(1) = ISM_ORIG
            !    ISCR(2) = KSM_ORIG
            !    ISCR(3) = LSM_ORIG
            !    ISCR(4) = JSM_ORIG
            !
            !    ISM = ISCR(I4_REO(1))
            !    KSM = ISCR(I4_REO(2))
            !    LSM = ISCR(I4_REO(3))
            !    JSM = ISCR(I4_REO(4))

            NI = NOBPTS(ITYP,ISM)
            NJ = NOBPTS(JTYP,JSM)
            NK = NOBPTS(KTYP,KSM)
            NL = NOBPTS(LTYP,LSM)
            NIK = NI*NK
            NJL = NJ*NL
            if ((NIK == 0) .or. (NJL == 0)) cycle

            ITPSM_ORIG = (ITYP_ORIG-1)*NSMOB+ISM_ORIG
            JTPSM_ORIG = (JTYP_ORIG-1)*NSMOB+JSM_ORIG
            KTPSM_ORIG = (KTYP_ORIG-1)*NSMOB+KSM_ORIG
            LTPSM_ORIG = (LTYP_ORIG-1)*NSMOB+LSM_ORIG

            if ((ITPSM_ORIG >= KTPSM_ORIG) .and. (JTPSM_ORIG >= LTPSM_ORIG)) then

              IFIRST = 1
              ! Loop over batches of I strings
              do IIPART=1+MYRANK,NPART,NPROCS
                IBOT = 1+(IIPART-1)*NPARTSZ
                ITOP = min(IBOT+NPARTSZ-1,NROW)
                NIBTC = ITOP-IBOT+1
                if (NIBTC <= 0) exit
                ! Loop over batches of intermediate strings
                KBOT = 1-MAXK
                KTOP = 0
                do
                  KBOT = KBOT+MAXK
                  KTOP = KTOP+MAXK

                  ! obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)

                  ! Obtain all double excitations from this group of K strings
                  II12 = 1
                  K12 = 1
                  ! Creation / annihilation maps, conjugated of above
                  if (I4_AC(4) == 1) then
                    JAC = 2
                  else
                    JAC = 1
                  end if
                  if (I4_AC(3) == 1) then
                    LAC = 2
                  else
                    LAC = 1
                  end if
                  !KFRST = 1
                  call ADAADAST_GAS(1,JSM,JTYP,NJ,JAC,1,LSM,LTYP,NL,LAC,ICCTP,ICCSM,IGRP,KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,JFRST, &
                                    KFRST,II12,K12,SCLFAC)

                  JFRST = 0
                  KFRST = 0

                  if (NKBTC == 0) exit
                  ! Loop over jl in TS classes and gather
                  J = 0
                  L = 1
                  do IJL=1,NJL
                    call NXTIJ(J,L,NJ,NL,JLSM,NONEW)
                    I1JL = (L-1)*NJ+J
                    ! CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
                    JLOFF = (IJL-1)*NKBTC*NIBTC+1
                    call MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,I1(:,I1JL),XI1S(:,I1JL))
                  end do

                  !=============================================
                  ! SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
                  !=============================================
                  ! Obtain two electron integrals as xint(ik,jl) = (ij!kl)-(il!kj)
                  IKSM = 0
                  JLSM = 0
                  if (IFIRST == 1) then
                    if (I4_AC(1) == I4_AC(3)) then
                      ! a+ a a+ a
                      ICOUL = 2
                    else
                      ! a+ a a a+ or a+ a a a+
                      ICOUL = 1
                    end if
                    ! Use coulomb - exchange or just coulomb integrals ?
                    if ((ITPSM_ORIG == KTPSM_ORIG) .and. (JTPSM_ORIG == LTPSM_ORIG)) then
                      ! No use of exchange
                      IXCHNG = 0
                    else
                      ! Exchange used, combines two terms
                      IXCHNG = 1
                    end if
                    if ((ITPSM_ORIG /= KTPSM_ORIG) .and. (JTPSM_ORIG /= LTPSM_ORIG)) then
                      ! Exchange used, combines four terms
                      FACX = -One
                    else
                      FACX = -Half
                    end if
#                   ifdef _DEBUGPRINT_
                    write(u6,*) ' ITPSM_ORIG,KTPSM_ORIG,JTPSM_ORIG,LTPSM_ORIG,FACX',ITPSM_ORIG,KTPSM_ORIG,JTPSM_ORIG,LTPSM_ORIG,FACX
#                   endif
                    ! fetch integrals
                    ! we want the operator in the form a+i ak a+l aj ((ij!lk)-(ik!lj))
                    if (ICOUL == 2) then
                      ! Obtain X2(ik,lj) = (ij!lk)
                      call GETINT(XINT,ITYP,ISM,JTYP,JSM,LTYP,LSM,KTYP,KSM,IXCHNG,IKSM,JLSM,ICOUL)
                    else if (ICOUL == 1) then
                      call GETINT(XINT,ITYP,ISM,KTYP,KSM,JTYP,JSM,LTYP,LSM,IXCHNG,IKSM,JLSM,ICOUL)
                    end if

                  end if
                  ! End if integrals should be fetched
                  IFIRST = 0
                  ! and now, to the work
                  LIKB = NIBTC*NKBTC
#                 ifdef _DEBUGPRINT_
                  write(u6,*) ' Integral block'
                  call WRTMAT(XINT,NIK,NJL,NIK,NJL)
                  write(u6,*) ' CSCR matrix'
                  call WRTMAT(CSCR,LIKB,NJL,LIKB,NJL)
#                 endif

                  !!MXACIJO = MXACIJ
                  !MXACIJ = MAX(MXACIJ,LIKB*NJL,LIKB*NIK)
                  !!if (MXACIJ > MXACIJO) then
                  !!  write(u6,*) ' New max MXACIJ = ', MXACIJ
                  !!  write(u6,*) ' ISCTP,ICCTP', ISCTP,ICCTP
                  !!  write(u6,*) ' ITYP,JTYP,KTYP,LTYP',ITYP,JTYP,KTYP,LTYP
                  !!  write(u6,*) 'NIJ NJL NIBTC NKBTC',NIJ,NJL,NIBTC,NKBTC
                  !!end if

                  FACTORC = Zero
                  FACTORAB = FACX
                  call MATML7(SSCR,CSCR,XINT,LIKB,NIK,LIKB,NJL,NIK,NJL,FACTORC,FACTORAB,2)
#                 ifdef _DEBUGPRINT_
                  write(u6,*) ' SSCR matrix'
                  call WRTMAT(SSCR,LIKB,NIK,LIKB,NIK)
#                 endif
                  ! ============================
                  ! Loop over ik and scatter out
                  ! ============================
                  ! Generate double excitations from K strings
                  ! I strings connected with K strings in batch <I!a+i a+k!K)
                  II12 = 2

                  if (IFRST == 1) KFRST = 1

                  IAC = I4_AC(1)
                  KAC = I4_AC(2)

                  !KFRST = 1
                  call ADAADAST_GAS(1,ISM,ITYP,NI,IAC,1,KSM,KTYP,NK,KAC,ISCTP,ISCSM,IGRP,KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,IFRST, &
                                  KFRST,II12,K12,One)

                  IFRST = 0
                  KFRST = 0

                  I = 0
                  K = 1
                  do IK=1,NIK
                    call NXTIJ(I,K,NI,NK,IKSM,NONEW)
                    IKOFF = (K-1)*NI+I
                    ISBOFF = 1+(IK-1)*NIBTC*NKBTC
                    call MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,NKBTC,I1(:,IKOFF),XI1S(:,IKOFF))
                  end do
                  !write(u6,*) ' first element of updated SB',SB(1)

                  if (KEND /= 0) exit
                end do
                ! End of loop over partitionings of resolution strings
              end do
              ! End of loop over batches of I strings
            end if
            ! End of if I >= K, J >= L
          end do
          ! End of loop over KSM
        end do
        ! End of loop over ISM
      end do
    end if
    ! End of a+ a+ a a/a a a+ a+ versus a+ a a+ a switch

  end do

end if

call mma_deallocate(SCR)

end subroutine RSBB2A
