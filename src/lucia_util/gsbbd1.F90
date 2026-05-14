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
! Copyright (C) 1991,1995,1998, Jeppe Olsen                            *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine GSBBD1(RHO1,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NGAS,ISEL,ICEL,SB,CB,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1, &
                  XI1S,I2,XI2S,NSMOB,RHO1S,SCLFAC,IPHGAS,IDOSRHO1,SRHO1,IAB)
! SUBROUTINE GSBBD1 --> 40
!
! Contributions to one electron density matrix from column excitations
!
! GAS version, August 95, Jeppe Olsen
! Particle-Hole version of Jan. 98
!
! =====
! Input
! =====
! RHO1  : One body density matrix to be updated
! NACOB : Number of active orbitals
! ISCSM,ISCTP : Symmetry and type of sigma columns
! ICCSM,ICCTP : Symmetry and type of C     columns
! IGRP : String group of columns
! NROW : Number of rows in S and C block
! NGAS : Number of active spaces
! ISEL : Number of electrons per AS for S block
! ICEL : Number of electrons per AS for C block
! CB   : Input C block
! MXPNGAS : Max number of AS spaces (program parameter)
! NOBPTS  : Number of orbitals per type and symmetry
! IOBPTS : base for orbitals of given type and symmetry
! IBORB  : Orbitals of given type and symmetry
! NSMOB  : Number of symmetries of orbitals
! MAXI   : Largest Number of "spectator strings" treated simultaneously
! MAXK   : Largest number of inner resolution strings treated at simult.
!
! ======
! Output
! ======
! RHO1 : Updated density block
!
! =======
! Scratch
! =======
!
! SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
!              largest number of orbital pairs of given symmetries and
!              types.
! I1, XI1S   : MAXK*Max number of orbitals of given type and symmetry
! I2, XI2S   : MAXK*Max number of orbitals of given type and symmetry
!              type and symmetry
! RHO1S : Space for one electron density
!
! Jeppe Olsen, Winter of 1991
! Updated for GAS, August '95

use Symmetry_Info, only: Mul
use Para_Info, only: MyRank, nProcs
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: RHO1(*), XI1S(*), XI2S(*), RHO1S(*), SRHO1(*)
integer(kind=iwp), intent(in) :: NACOB, ISCSM, ISCTP, ICCSM, ICCTP, IGRP, NROW, NGAS, ISEL(NGAS), ICEL(NGAS), MXPNGAS, &
                                 NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*), MAXI, MAXK, NSMOB, IPHGAS(*), IDOSRHO1, IAB
real(kind=wp), intent(in) :: SB(*), CB(*), SCLFAC
real(kind=wp), intent(_OUT_) :: SSCR(*), CSCR(*)
integer(kind=iwp), intent(inout) :: I1(*), I2(*)
integer(kind=iwp) :: IBIORB, IBJORB, IBOT, ICGOFF, ICGRP(16), IDOCOMP, IIORB, IIPART, IJ_AC(2), IJ_DIM(2), IJ_OFF(2), IJ_REO(2), &
                     IJ_SM(2), IJ_TP(2), IJSM, IJTP, IORB, ISGOFF, ISGRP(16), ISM, ITOP, ITP(256), ITYP, IXXX, JJORB, JORB, JSM, &
                     JTP(256), JTYP, KACT, KBOT, KEND, KTOP, LKABTC, NIBTC, NIORB, NIPART, NIPARTSZ, NJORB, NKAEFF, NKASTR, NKI, &
                     NSXTP
real(kind=wp) :: FACTORAB, FACTORC, SCLFACS, SIGNIJ, XAB

! Add or subtract for spindensity
if (IAB == 1) then
  XAB = One
else
  XAB = -One
end if
! Local arrays
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ================'
write(u6,*) ' GSBBu61 in action'
write(u6,*) ' ================'
write(u6,*)
write(u6,*) ' Occupation of active left strings'
call IWRTMA(ISEL,1,NGAS,1,NGAS)
write(u6,*) ' Occupation of active Right strings'
call IWRTMA(ICEL,1,NGAS,1,NGAS)
write(u6,*) ' ISCSM, ICCSM = ',ISCSM,ICCSM

write(u6,*) ' GSBBD1, sclfac ',SCLFAC
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

! Groups defining supergroups
!    GET_SPGP_INF(ISPGP,ITP,IGRP)
call GET_SPGP_INF(ICCTP,IGRP,ICGRP)
call GET_SPGP_INF(ISCTP,IGRP,ISGRP)

! Type of single excitations that connects the two column strings
call SXTYP2_GAS(NSXTP,ITP,JTP,NGAS,ISEL,ICEL,IPHGAS)
! Symmetry of single excitation that connects IBSM and JBSM
IJSM = Mul(ISCSM,ICCSM)
#ifdef _DEBUGPRINT_
write(u6,*) ' ISCSM,ICCSM IJSM ',ISCSM,ICCSM,IJSM
#endif
if (IJSM /= 0) then
  do IJTP=1,NSXTP
    ITYP = ITP(IJTP)
    JTYP = JTP(IJTP)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' ITYP JTYP ',ITYP,JTYP
#   endif
    ! Hvilken vej skal vi valge,
    ! Mi pojdem drugim putem (C)
    ! VV: the code below confuses Absoft compiler and was rewritten.
    !NOP = 2
    IJ_AC(1) = 2
    IJ_AC(2) = 1
    IJ_TP(1) = ITYP
    IJ_TP(2) = JTYP
    !if (IUSE_PH == 1) call ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TP,IJ_AC,IJ_REO,SIGNIJ)
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
    !ISCR(1) = ITYP
    !ISCR(2) = JTYP
    !IJ_TP(1) = ISCR(IJ_REO(1))
    !IJ_TP(2) = ISCR(IJ_REO(2))

    do ISM=1,NSMOB
      ! new i and j so new intermediate strings

      JSM = Mul(ISM,IJSM)
      if (JSM == 0) cycle
#     ifdef _DEBUGPRINT_
      write(u6,*) ' ISM JSM ',ISM,JSM
#     endif
      NIORB = NOBPTS(ITYP,ISM)
      NJORB = NOBPTS(JTYP,JSM)
      IBIORB = IOBPTS(ITYP,ISM)
      IBJORB = IOBPTS(JTYP,JSM)
      ! Reorder
      !ISCR(1) = ISM
      !ISCR(2) = JSM
      !IJ_SM(1) = ISCR(IJ_REO(1))
      !IJ_SM(2) = ISCR(IJ_REO(2))

      !ISCR(1) = NIORB
      !ISCR(2) = NJORB
      !IJ_DIM(1) = ISCR(IJ_REO(1))
      !IJ_DIM(2) = ISCR(IJ_REO(2))

      !ISCR(1) = IBIORB
      !ISCR(2) = IBJORB
      !IJ_OFF(1) = ISCR(IJ_REO(1))
      !IJ_OFF(2) = ISCR(IJ_REO(2))

      if (IJ_REO(1) == 1) then
        IJ_SM(1) = ISM
        IJ_SM(2) = JSM
        IJ_DIM(1) = NIORB
        IJ_DIM(2) = NJORB
        IJ_OFF(1) = IBIORB
        IJ_OFF(2) = IBJORB
      else
        IJ_SM(1) = JSM
        IJ_SM(2) = ISM
        IJ_DIM(1) = NJORB
        IJ_DIM(2) = NIORB
        IJ_OFF(1) = IBJORB
        IJ_OFF(2) = IBIORB
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) ' NIORB NJORB ',NIORB,NJORB
#     endif
      if ((NIORB == 0) .or. (NJORB == 0)) cycle

      ! For operator connecting to |Ka> and |Ja> i.e. operator 2
      SCLFACS = SCLFAC*SIGNIJ
#     ifdef _DEBUGPRINT_
      write(u6,*) ' IJ_SM,IJ_TP,IJ_AC',IJ_SM(2),IJ_TP(2),IJ_AC(2)
#     endif
      call ADAST_GAS(IJ_SM(2),IJ_TP(2),NGAS,ICGRP,ICCSM,I1,XI1S,NKASTR,KACT,SCLFACS,IJ_AC(1))
      ! For operator connecting |Ka> and |Ia>, i.e. operator 1
      if (NKASTR == 0) cycle
      call ADAST_GAS(IJ_SM(1),IJ_TP(1),NGAS,ISGRP,ISCSM,I2,XI2S,NKASTR,KACT,One,IJ_AC(1))
      if (NKASTR == 0) cycle
      ! Compress list to common nonvanishing elements
      IDOCOMP = 1
      if (IDOCOMP == 1) then
        call COMPRS2LST(I1,XI1S,IJ_DIM(2),I2,XI2S,IJ_DIM(1),NKASTR,NKAEFF)
      else
        NKAEFF = NKASTR
      end if
      !write(u6,*) ' NKAEFF NKASTR',NKAEFF,NKASTR

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
          IBOT = (IIPART-1)*NIPARTSZ+1
          ITOP = min(IBOT+NIPARTSZ-1,NROW)
          NIBTC = ITOP-IBOT+1
          if (NIBTC <= 0) exit
          ! Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
          do JJORB=1,IJ_DIM(2)
            ICGOFF = 1+(JJORB-1)*LKABTC*NIBTC
            call MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,LKABTC,I1(KBOT+(JJORB-1)*NKASTR),XI1S(KBOT+(JJORB-1)*NKASTR))
          end do
          ! Obtain SSCR(I,K,IORB) = SUM(I)<K!A IORB!J>S(I,J)
          do IIORB=1,IJ_DIM(1)
            ! Gather S Block
            ISGOFF = 1+(IIORB-1)*LKABTC*NIBTC
            call MATCG(SB,SSCR(ISGOFF),NROW,NIBTC,IBOT,LKABTC,I2(KBOT+(IIORB-1)*NKASTR),XI2S(KBOT+(IIORB-1)*NKASTR))
          end do

          ! And then the hard work
          NKI = LKABTC*NIBTC
#         ifdef _DEBUGPRINT_
          write(u6,*) ' CSCR and SSCR'
          call WRTMAT(CSCR,IJ_DIM(2),NKI,IJ_DIM(2),NKI)
          call WRTMAT(SSCR,IJ_DIM(1),NKI,IJ_DIM(1),NKI)
#         endif
          FACTORC = Zero
          FACTORAB = One
          call MATML7(RHO1S,SSCR,CSCR,IJ_DIM(1),IJ_DIM(2),NKI,IJ_DIM(1),NKI,IJ_DIM(2),FACTORC,FACTORAB,1)

#         ifdef _DEBUGPRINT_
          write(u6,*) ' Block to one-body density'
          call WRTMAT(RHO1S,IJ_DIM(1),IJ_DIM(2),IJ_DIM(1),IJ_DIM(2))
#         endif
          ! Scatter out to complete matrix
          do JJORB=1,IJ_DIM(2)
            JORB = IJ_OFF(2)-1+JJORB
            do IIORB=1,IJ_DIM(1)
              IORB = IJ_OFF(1)-1+IIORB
              RHO1((JORB-1)*NACOB+IORB) = RHO1((JORB-1)*NACOB+IORB)+RHO1S((JJORB-1)*IJ_DIM(1)+IIORB)
              if (IDOSRHO1 == 1) SRHO1((JORB-1)*NACOB+IORB) = SRHO1((JORB-1)*NACOB+IORB)+XAB*RHO1S((JJORB-1)*IJ_DIM(1)+IIORB)
            end do
          end do
          ! End of hard work

        end do
        ! end of this I partitioning
        ! end of this K partitioning
        if (KEND /= 0) exit
      end do
      ! End of loop over I partitionings
    end do
    ! (end of loop over symmetries)
  end do
end if

!stop ' enforced stop in RSBBD1'

end subroutine GSBBD1
