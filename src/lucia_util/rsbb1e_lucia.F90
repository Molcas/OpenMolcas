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

subroutine RSBB1E_LUCIA(ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NGAS,ISEL,ICEL,SB,CB,ADSXA,STSTSX,NOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2, &
                        XI2S,H,NSMOB,NSMST,NSMSX,MOC,MXSXST,IH2TRM,SCLFAC,IUSE_PH,IPHGAS,NTESTG)
! SUBROUTINE RSBB1E_LUCIA --> 33
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
! ADASX : sym of a+, a => sym of a+a
! ADSXA : sym of a+, a+a => sym of a
! SXSTST : Sym of sx,!st> => sym of sx !st>
! STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
! NTSOB  : Number of orbitals per type and symmetry
! IBTSOB : base for orbitals of given type and symmetry
! IBORB  : Orbitals of given type and symmetry
! NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
!       single excitations, double excitations
! MAXI   : Largest Number of ' spectator strings 'treated simultaneously
! MAXK   : Largest number of inner resolution strings treated at simult.
!
! MOC  : Use MOC method ( instead of N-1 resolution method )
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

use Constants, only: One
use Para_Info, only: MyRank, nProcs
use lucia_data, only: MXPOBS, MXPNGAS, MXPTSOB

implicit none
integer ISCSM, ISCTP, ICCSM, ICCTP, IGRP, NROW, NGAS, MAXI, MAXK, NSMOB, NSMST, NSMSX, MOC, MXSXST, IH2TRM, IUSE_PH, NTESTG
real*8 SCLFAC
! General input
integer ADSXA(MXPOBS,2*MXPOBS), STSTSX(NSMST,NSMST)
integer NOBPTS(MXPNGAS,*)
integer IPHGAS(NGAS)
! Specific Input
integer ISEL(NGAS), ICEL(NGAS)
real*8 CB(*)
! Output
real*8 SB(*)
! Scatch
real*8 SSCR(*), CSCR(*), XI1S(*), H(*), XI2S(*)
integer I1(*), I2(*)
! Local arrays ( assume MPNGAS = 16 ) !!!
integer ITP(16), JTP(16)
integer ISGRP(16), ICGRP(16)
! For transposing integral block
real*8 HSCR(MXPTSOB*MXPTSOB)
integer IJ_REO(2), IJ_DIM(2), IJ_SM(2), IJ_TP(2), IJ_AC(2)
! Type of single excitations that connects the two column strings
integer NTESTL, NTEST, NIPART, NIPARTSZ, IFRST, IJSM, IJTP, NSXTP, ITYP, JTYP, IXXX, ISM, JSM, KFRST, NIORB, NJORB, IDOCOMP, &
        NKAEFF, NKASTR, KBOT, KTOP, KEND, LKABTC, IIPART, IBOT, ITOP, NIBTC, JJORB, ICGOFF, NIK, IIORB, ISBOFF, IEND, KACT
real*8 SIGNIJ, SCLFACS, FACTORC, FACTORAB

!MOC = 1
NTESTL = 0
NTEST = max(NTESTL,NTESTG)
if (NTEST >= 500) then
  write(6,*)
  write(6,*) ' ======================='
  write(6,*) ' Information from RSBB1E'
  write(6,*) ' ======================='
  write(6,*)
  write(6,*) ' RSBB1E : MOC,IH2TRM,IUSE_PH ',MOC,IH2TRM,IUSE_PH
  write(6,*) ' ISEL :'
  call IWRTMA(ISEL,1,NGAS,1,NGAS)
  write(6,*) ' ICEL :'
  call IWRTMA(ICEL,1,NGAS,1,NGAS)
end if

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

IFRST = 1
! Types of single excitations that connect ISEL and ICEL
call SXTYP2_GAS(NSXTP,ITP,JTP,NGAS,ISEL,ICEL,IPHGAS)
! Symmetry of single excitation that connects IBSM and JBSM
IJSM = STSTSX(ISCSM,ICCSM)
if (IJSM == 0) goto 1001
do IJTP=1,NSXTP
  ITYP = ITP(IJTP)
  JTYP = JTP(IJTP)
  if (NTEST >= 2000) write(6,*) ' ITYP JTYP ',ITYP,JTYP
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
  SIGNIJ = 1.0d0
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
  !if (iscr(1) == -1000) write(6,*) IJ_TP,IJ_REO
  !
  !ISCR(1) = ITYP
  !ISCR(2) = JTYP
  !IJ_TP(1) = ISCR(IJ_REO(1))
  !IJ_TP(2) = ISCR(IJ_REO(2))

  do ISM=1,NSMOB
    JSM = ADSXA(ISM,IJSM)
    ! New intermediate strings will be accessed so
    KFRST = 1
    if (JSM == 0) goto 800
    if (NTEST >= 2000) write(6,*) ' ISM JSM ',ISM,JSM
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

    if ((NIORB == 0) .or. (NJORB == 0)) goto 800
      ! Fetch integrals : For CI-transformations using RSBB1E
      ! most of the blocks vanishes
      ! Obtain one electron integrals (ISM,ITP,JSM,JTP) transposed
    if (IJ_REO(1) == 1) then
      ! obtain integrals h(j,i)
      call GETH1(HSCR,IJ_SM(1),IJ_TP(1),IJ_SM(2),IJ_TP(2))
      call TRPMAT(HSCR,IJ_DIM(1),IJ_DIM(2),H)
    else
      ! Obtain integrals h(i,j)
      call GETH1(H,IJ_SM(2),IJ_TP(2),IJ_SM(1),IJ_TP(1))
    end if
    !OLD XNORM = INPROD(H,H,IJ_DIM(1)*IJ_DIM(2))
    !OLD if (XNORM == 0) goto 800
    if (MOC == 0) then

      ! ================================================================
      !                    Use N-1 resolution method
      ! ================================================================

      ! Obtain annihilation/creation maps for all K strings

      ! For operator connecting to |Ka> and |Ja> i.e. operator 2
      SCLFACS = SIGNIJ*SCLFAC
      if (NTEST >= 1000) write(6,*) ' IJ_SM,IJ_TP,IJ_AC',IJ_SM(2),IJ_TP(2),IJ_AC(2)
      call ADAST_GAS(IJ_SM(2),IJ_TP(2),NGAS,ICGRP,ICCSM,I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFACS,IJ_AC(1))
      ! For operator connecting |Ka> and |Ia>, i.e. operator 1
      call ADAST_GAS(IJ_SM(1),IJ_TP(1),NGAS,ISGRP,ISCSM,I2,XI2S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJ_AC(1))
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
700   continue
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
        !call TRPMAT(HSCR,IJ_DIM(1),IJ_DIM(2),H)
        ! Problems when HOLE switches blocks around ?
        !call GETH1(H,IJ_SM(2),IJ_TP(2),IJ_SM(1),IJ_TP(1))
        if (NTEST >= 1000) then
          write(6,*) ' RSBB1E H BLOCK'
          call WRTMAT(H,IJ_DIM(2),IJ_DIM(1),IJ_DIM(2),IJ_DIM(1))
        end if
        ! Sscr(I,K,i) = CSCR(I,K,j)*h(j,i)
        NIK = NIBTC*LKABTC
        FACTORC = 0.0d0
        FACTORAB = 1.0d0
        if (NTEST >= 2000) then
          write(6,*) ' CSCR array,NIK X NJORB array'
          call WRTMAT(CSCR,NIK,IJ_DIM(2),NIK,IJ_DIM(2))
        end if
        call MATML7(SSCR,CSCR,H,NIK,IJ_DIM(1),NIK,IJ_DIM(2),IJ_DIM(2),IJ_DIM(1),FACTORC,FACTORAB,0)
        if (NTEST >= 2000) then
          write(6,*) ' SSCR array,NIK X NIORB array'
          call WRTMAT(SSCR,NIK,IJ_DIM(1),NIK,IJ_DIM(1))
        end if
        ! S(I,a+ K) =  S(I, a+ K) + sgn*Sscr(I,K,i)
        do IIORB=1,IJ_DIM(1)
          ISBOFF = 1+(IIORB-1)*LKABTC*NIBTC
          call MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,LKABTC,I2(KBOT+(IIORB-1)*NKASTR),XI2S(KBOT+(IIORB-1)*NKASTR))
        end do

      end do
      ! end of this K partitioning
      if (KEND == 0) goto 700
      ! End of loop over I partitioninigs
    end if
    ! (End of algorithm switch)
800 continue
  end do
  ! (end of loop over symmetries)
end do
1001 continue

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(NSMSX)
  call Unused_integer(MXSXST)
end if

end subroutine RSBB1E_LUCIA
