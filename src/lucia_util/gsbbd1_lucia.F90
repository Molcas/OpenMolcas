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
subroutine GSBBD1_LUCIA(RHO1,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NGAS,ISEL,ICEL,SB,CB,ADSXA,SXSTST,STSTSX,MXPNGAS,NOBPTS, &
                        IOBPTS,ITSOB,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,H,NSMOB,NSMST,NSMSX,MXPOBS,RHO1S,SCLFAC,IUSE_PH,IPHGAS, &
                        IDOSRHO1,SRHO1,IAB)
! SUBROUTINE GSBBD1_LUCIA --> 40
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
! ADASX : sym of a+, a => sym of a+a
! ADSXA : sym of a+, a+a => sym of a
! SXSTST : Sym of sx,!st> => sym of sx !st>
! STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
! MXPNGAS : Max number of AS spaces ( program parameter )
! NOBPTS  : Number of orbitals per type and symmetry
! IOBPTS : base for orbitals of given type and symmetry
! IBORB  : Orbitals of given type and symmetry
! NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
!       single excitations, double excitations
! MAXI   : Largest Number of ' spectator strings 'treated simultaneously
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

use Para_Info, only: MyRank, nProcs

implicit real*8(A-H,O-Z)
! General input
integer ADSXA(MXPOBS,2*MXPOBS), SXSTST(NSMSX,NSMST), STSTSX(NSMST,NSMST)
integer NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*), ITSOB(*)
integer IPHGAS(*)
! Input
integer ISEL(NGAS), ICEL(NGAS)
dimension CB(*), SB(*)
! Output
dimension RHO1(*), SRHO1(*)
! Scratch
dimension SSCR(*), CSCR(*), RHO1S(*)
dimension I1(*), XI1S(*)
dimension I2(*), XI2S(*)
! Local arrays ( assume MPNGAS = 16 ) !!!
dimension ITP(16*16), JTP(16*16)
dimension IJ_REO(2), IJ_DIM(2), IJ_SM(2), IJ_TP(2), IJ_AC(2)
dimension IJ_OFF(2)
!dimension ISCR(2)
dimension ICGRP(16), ISGRP(16)
dimension H(*)

! Add or subtract for spindensity
if (IAB == 1) then
  XAB = 1.0d0
else
  XAB = -1.0d0
end if
! Local arrays
NTEST = 0
if (NTEST >= 1000) then
  write(6,*)
  write(6,*) ' ================'
  write(6,*) ' GSBBD1 in action'
  write(6,*) ' ================'
  write(6,*)
  write(6,*) ' Occupation of active left strings'
  call IWRTMA(ISEL,1,NGAS,1,NGAS)
  write(6,*) ' Occupation of active Right strings'
  call IWRTMA(ICEL,1,NGAS,1,NGAS)
  write(6,*) ' ISCSM, ICCSM = ',ISCSM,ICCSM

  write(6,*) ' GSBBD1, sclfac ',SCLFAC
end if

IFRST = 1
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
IJSM = STSTSX(ISCSM,ICCSM)
if (NTEST >= 1000) write(6,*) ' ISCSM,ICCSM IJSM ',ISCSM,ICCSM,IJSM
if (IJSM == 0) goto 1001
do IJTP=1,NSXTP
  ITYP = ITP(IJTP)
  JTYP = JTP(IJTP)
  if (NTEST >= 1000) write(6,*) ' ITYP JTYP ',ITYP,JTYP
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
  !ISCR(1) = ITYP
  !ISCR(2) = JTYP
  !IJ_TP(1) = ISCR(IJ_REO(1))
  !IJ_TP(2) = ISCR(IJ_REO(2))

  do ISM=1,NSMOB
    ! new i and j so new intermediate strings
    KFRST = 1

    JSM = ADSXA(ISM,IJSM)
    if (JSM == 0) goto 800
    if (NTEST >= 1000) write(6,*) ' ISM JSM ',ISM,JSM
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

    if (NTEST >= 2000) write(6,*) ' NIORB NJORB ',NIORB,NJORB
    if ((NIORB == 0) .or. (NJORB == 0)) goto 800

    ! For operator connecting to |Ka> and |Ja> i.e. operator 2
    SCLFACS = SCLFAC*SIGNIJ
    if (NTEST >= 1000) write(6,*) ' IJ_SM,IJ_TP,IJ_AC',IJ_SM(2),IJ_TP(2),IJ_AC(2)
    call ADAST_GAS(IJ_SM(2),IJ_TP(2),NGAS,ICGRP,ICCSM,I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFACS,IJ_AC(1))
    ! For operator connecting |Ka> and |Ia>, i.e. operator 1
    if (NKASTR == 0) goto 800
    ONE = 1.0d0
    call ADAST_GAS(IJ_SM(1),IJ_TP(1),NGAS,ISGRP,ISCSM,I2,XI2S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJ_AC(1))
    if (NKASTR == 0) goto 800
    ! Compress list to common nonvanishing elements
    IDOCOMP = 1
    if (IDOCOMP == 1) then
      call COMPRS2LST(I1,XI1S,IJ_DIM(2),I2,XI2S,IJ_DIM(1),NKASTR,NKAEFF)
    else
      NKAEFF = NKASTR
    end if
    !write(6,*) ' NKAEFF NKASTR',NKAEFF,NKASTR

    ! Loop over partitionings of N-1 strings
    KBOT = 1-MAXK
    KTOP = 0
700 continue
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

      if (NTEST >= 1000) then
        write(6,*) ' CSCR and SSCR'
        call WRTMAT(CSCR,IJ_DIM(2),NKI,IJ_DIM(2),NKI)
        call WRTMAT(SSCR,IJ_DIM(1),NKI,IJ_DIM(1),NKI)
      end if

      ! And then the hard work
      NKI = LKABTC*NIBTC
      FACTORC = 0.0d0
      FACTORAB = 1.0d0
      call MATML7(RHO1S,SSCR,CSCR,IJ_DIM(1),IJ_DIM(2),NKI,IJ_DIM(1),NKI,IJ_DIM(2),FACTORC,FACTORAB,1)

      if (NTEST >= 100) then
        write(6,*) ' Block to one-body density'
        call WRTMAT(RHO1S,IJ_DIM(1),IJ_DIM(2),IJ_DIM(1),IJ_DIM(2))
      end if
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
    if (KEND == 0) goto 700
    ! End of loop over I partitionings
800 continue
  end do
  ! (end of loop over symmetries)
end do
1001 continue

!stop ' enforced stop in RSBBD1'

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(SXSTST)
  call Unused_integer_array(ITSOB)
  call Unused_real_array(H)
  call Unused_integer(IUSE_PH)
end if

end subroutine GSBBD1_LUCIA
