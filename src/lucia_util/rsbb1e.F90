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
! Copyright (C) 1991,1997, Jeppe Olsen                                 *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine RSBB1E(ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NGAS,ISEL,ICEL,SB,CB,NOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,H,NSMOB,MOC, &
                  SCLFAC,IPHGAS)
! SUBROUTINE RSBB1E --> 33
!
! One electron excitations on column strings
! If IH2TRM /= 0 then the diagonal and one-electron
! excitations arising from the two body operator is also included
!
! =====
! Input
! =====
!
! ISCSM,ISCTP : Symmetry and type of sigma columns
! ISBCTP : Base for sigma column types
! ICCSM,ICCTP : Symmetry and type of C     columns
! ICBCTP : Base for C     column types
! IGRP : String group of columns
! NROW : Number of rows in S and C block
! NGAS : Number of active sets
! ISEL : Occupation in each active set for sigma block
! ICEL : Occupation in each active set for C     block
! CB   : Input C block
! NTSOB  : Number of orbitals per type and symmetry
! IBTSOB : base for orbitals of given type and symmetry
! IBORB  : Orbitals of given type and symmetry
! NSMOB  : Number of symmetries of orbitals
! MAXI   : Largest Number of "spectator strings" treated simultaneously
! MAXK   : Largest number of inner resolution strings treated at simult.
!
! MOC  : Use MOC method (instead of N-1 resolution method)
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
! H : Space for one electron integrals
!
! Jeppe Olsen, Winter of 1991
!              IUSE_PH added winter of 97

use Symmetry_Info, only: Mul
use Para_Info, only: MyRank, nProcs
use lucia_data, only: MXPNGAS, MXPTSOB
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ISCSM, ISCTP, ICCSM, ICCTP, IGRP, NROW, NGAS, ISEL(NGAS), ICEL(NGAS), NOBPTS(MXPNGAS,*), MAXI, &
                                 MAXK, NSMOB, MOC, IPHGAS(NGAS)
real(kind=wp), intent(_OUT_) :: SB(*), SSCR(*), CSCR(*), H(*)
real(kind=wp), intent(in) :: CB(*), SCLFAC
integer(kind=iwp), intent(inout) :: I1(*), I2(*)
real(kind=wp), intent(inout) :: XI1S(*), XI2S(*)
integer(kind=iwp) :: IBOT, ICGOFF, ICGRP(16), IDOCOMP, IIORB, IIPART, IJ_AC(2), IJ_DIM(2), IJ_REO(2), IJ_SM(2), IJ_TP(2), IJSM, &
                     IJTP, ISBOFF, ISGRP(16), ISM, ITOP, ITP(16), ITYP, IXXX, JJORB, JSM, JTP(16), JTYP, KACT, KBOT, KEND, KTOP, &
                     LKABTC, NIBTC, NIK, NIORB, NIPART, NIPARTSZ, NJORB, NKAEFF, NKASTR, NSXTP
real(kind=wp) :: FACTORAB, FACTORC, HSCR(MXPTSOB*MXPTSOB), SCLFACS, SIGNIJ

!MOC = 1
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ======================='
write(u6,*) ' Information from RSBB1E'
write(u6,*) ' ======================='
write(u6,*)
write(u6,*) ' RSBB1E : MOC ',MOC
write(u6,*) ' ISEL :'
call IWRTMA(ISEL,1,NGAS,1,NGAS)
write(u6,*) ' ICEL :'
call IWRTMA(ICEL,1,NGAS,1,NGAS)
#endif

! Number of partitionings over column strings
!SVC: determine optimum number of partitions as the lowest multiple of
!     NPROCS that satisfies a block size smaller than MAXI:
NIPART = 0
do
  NIPART = NIPART+NPROCS
  NIPARTSZ = max(NROW-1,0)/NIPART+1
  if (NIPARTSZ <= MAXI) exit
end do

! Obtain groups
!    GET_SPGP_INF(ISPGP,ITP,IGRP)
call GET_SPGP_INF(ICCTP,IGRP,ICGRP)
call GET_SPGP_INF(ISCTP,IGRP,ISGRP)

! Types of single excitations that connect ISEL and ICEL
call SXTYP2_GAS(NSXTP,ITP,JTP,NGAS,ISEL,ICEL,IPHGAS)
! Symmetry of single excitation that connects IBSM and JBSM
IJSM = Mul(ISCSM,ICCSM)
if (IJSM /= 0) then
  do IJTP=1,NSXTP
    ITYP = ITP(IJTP)
    JTYP = JTP(IJTP)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' ITYP JTYP ',ITYP,JTYP
#   endif
    ! Is this combination of types allowed
    !IJ_ACT = 1
    !if (IJ_ACT == 0) cycle
    ! Hvilken vej skal vi valge,
    !NOP = 2
    IJ_AC(1) = 2
    IJ_AC(2) = 1
    IJ_TP(1) = ITYP
    IJ_TP(2) = JTYP
    !if (IUSE_PH == 1) then
    !  call ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TP,IJ_AC,IJ_REO,SIGNIJ)
    !else
    IJ_REO(1) = 1
    IJ_REO(2) = 2
    SIGNIJ = One
    !end if

    if (IJ_REO(1) == 1) then

      IJ_TP(1) = ITYP
      IJ_TP(2) = JTYP
    else
      IXXX = IJ_AC(1)
      IJ_AC(1) = IJ_AC(2)
      IJ_AC(2) = IXXX

      !ISCR(1) = ITYP
      !ISCR(2) = JTYP
      IJ_TP(1) = JTYP
      IJ_TP(2) = ITYP
    end if

    !ISCR(1) = IJ_AC(1)
    !ISCR(2) = IJ_AC(2)
    !IJ_AC(1) = ISCR(IJ_REO(1))
    !IJ_AC(2) = ISCR(IJ_REO(2))

    ! nasty code to avoid optimization
    !if (iscr(1) == -1000) write(u6,*) IJ_TP,IJ_REO
    !
    !ISCR(1) = ITYP
    !ISCR(2) = JTYP
    !IJ_TP(1) = ISCR(IJ_REO(1))
    !IJ_TP(2) = ISCR(IJ_REO(2))

    do ISM=1,NSMOB
      JSM = Mul(ISM,IJSM)
      ! New intermediate strings will be accessed so
      if (JSM == 0) cycle
#     ifdef _DEBUGPRINT_
      write(u6,*) ' ISM JSM ',ISM,JSM
#     endif
      NIORB = NOBPTS(ITYP,ISM)
      NJORB = NOBPTS(JTYP,JSM)
      ! Reorder

      !ISCR(1) = ISM
      !ISCR(2) = JSM
      !IJ_SM(1) = ISCR(IJ_REO(1))
      !IJ_SM(2) = ISCR(IJ_REO(2))

      !ISCR(1) = NIORB
      !ISCR(2) = NJORB
      !IJ_DIM(1) = ISCR(IJ_REO(1))
      !IJ_DIM(2) = ISCR(IJ_REO(2))

      if (IJ_REO(1) == 1) then
        IJ_SM(1) = ISM
        IJ_SM(2) = JSM
        IJ_DIM(1) = NIORB
        IJ_DIM(2) = NJORB
      else
        IJ_SM(1) = JSM
        IJ_SM(2) = ISM
        IJ_DIM(1) = NJORB
        IJ_DIM(2) = NIORB
      end if

      if ((NIORB == 0) .or. (NJORB == 0)) cycle
      ! Fetch integrals : For CI-transformations using RSBB1E
      ! most of the blocks vanishes
      ! Obtain one electron integrals (ISM,ITP,JSM,JTP) transposed
      if (IJ_REO(1) == 1) then
        ! obtain integrals h(j,i)
        call GETH1(HSCR,IJ_SM(1),IJ_TP(1),IJ_SM(2),IJ_TP(2))
        call TRNSPS(IJ_DIM(1),IJ_DIM(2),HSCR,H)
      else
        ! Obtain integrals h(i,j)
        call GETH1(H,IJ_SM(2),IJ_TP(2),IJ_SM(1),IJ_TP(1))
      end if
      !OLD XNORM = dDot_(IJ_DIM(1)*IJ_DIM(2),H,1,H,1)
      !OLD if (XNORM == 0) cycle
      if (MOC == 0) then

        ! ==============================================================
        !                    Use N-1 resolution method
        ! ==============================================================

        ! Obtain annihilation/creation maps for all K strings

        ! For operator connecting to |Ka> and |Ja> i.e. operator 2
        SCLFACS = SIGNIJ*SCLFAC
#       ifdef _DEBUGPRINT_
        write(u6,*) ' IJ_SM,IJ_TP,IJ_AC',IJ_SM(2),IJ_TP(2),IJ_AC(2)
#       endif
        call ADAST_GAS(IJ_SM(2),IJ_TP(2),NGAS,ICGRP,ICCSM,I1,XI1S,NKASTR,KACT,SCLFACS,IJ_AC(1))
        ! For operator connecting |Ka> and |Ia>, i.e. operator 1
        call ADAST_GAS(IJ_SM(1),IJ_TP(1),NGAS,ISGRP,ISCSM,I2,XI2S,NKASTR,KACT,One,IJ_AC(1))
        ! Compress list to common nonvanishing elements
        IDOCOMP = 1
        if (IDOCOMP == 1) then
          call COMPRS2LST(I1,XI1S,IJ_DIM(2),I2,XI2S,IJ_DIM(1),NKASTR,NKAEFF)
        else
          NKAEFF = NKASTR
        end if
        ! Loop over partitionings of the row strings
        ! Loop over partitionings of N-1 strings
        KBOT = 1-MAXK
        KTOP = 0
        do
          KBOT = KBOT+MAXK
          KTOP = min(KTOP+MAXK,NKAEFF)
          if (KTOP == NKAEFF) then
            KEND = 1
          else
            KEND = 0
          end if
          LKABTC = KTOP-KBOT+1
          ! This is the place to start over partitioning of I strings
          do IIPART=1+MYRANK,NIPART,NPROCS
            IBOT = (IIPART-1)*MAXI+1
            ITOP = min(IBOT+MAXI-1,NROW)
            NIBTC = ITOP-IBOT+1
            if (NIBTC <= 0) exit
            ! Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
            do JJORB=1,IJ_DIM(2)
              ICGOFF = 1+(JJORB-1)*LKABTC*NIBTC
              call MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,LKABTC,I1(KBOT+(JJORB-1)*NKASTR),XI1S(KBOT+(JJORB-1)*NKASTR))
            end do
            ! Obtain one electron integrals (ISM,ITP,JSM,JTP) transposed
            !call GETH1(HSCR,IJ_SM(1),IJ_TP(1),IJ_SM(2),IJ_TP(2))
            !call TRNSPS(IJ_DIM(1),IJ_DIM(2),HSCR,H)
            ! Problems when HOLE switches blocks around ?
            !call GETH1(H,IJ_SM(2),IJ_TP(2),IJ_SM(1),IJ_TP(1))
#           ifdef _DEBUGPRINT_
            write(u6,*) ' RSBB1E H BLOCK'
            call WRTMAT(H,IJ_DIM(2),IJ_DIM(1),IJ_DIM(2),IJ_DIM(1))
#           endif
            ! Sscr(I,K,i) = CSCR(I,K,j)*h(j,i)
            NIK = NIBTC*LKABTC
            FACTORC = Zero
            FACTORAB = One
#           ifdef _DEBUGPRINT_
            write(u6,*) ' CSCR array,NIK X NJORB array'
            call WRTMAT(CSCR,NIK,IJ_DIM(2),NIK,IJ_DIM(2))
#           endif
            call MATML7(SSCR,CSCR,H,NIK,IJ_DIM(1),NIK,IJ_DIM(2),IJ_DIM(2),IJ_DIM(1),FACTORC,FACTORAB,0)
#           ifdef _DEBUGPRINT_
            write(u6,*) ' SSCR array,NIK X NIORB array'
            call WRTMAT(SSCR,NIK,IJ_DIM(1),NIK,IJ_DIM(1))
#           endif
            ! S(I,a+ K) =  S(I, a+ K) + sgn*Sscr(I,K,i)
            do IIORB=1,IJ_DIM(1)
              ISBOFF = 1+(IIORB-1)*LKABTC*NIBTC
              call MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,LKABTC,I2(KBOT+(IIORB-1)*NKASTR),XI2S(KBOT+(IIORB-1)*NKASTR))
            end do

          end do
          ! end of this K partitioning
          if (KEND /= 0) exit
        end do
        ! End of loop over I partitioninigs
      end if
      ! (End of algorithm switch)
    end do
    ! (end of loop over symmetries)
  end do
end if

end subroutine RSBB1E
