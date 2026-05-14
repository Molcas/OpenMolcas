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
subroutine GSBBD2B(RHO2,RHO2S,RHO2A,IASM,IATP,IBSM,IBTP,NIA,NIB,JASM,JATP,JBSM,JBTP,NJA,NJB,IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC, &
                   SB,CB,MXPNGAS,NOBPTS,IOBPTS,MAXK,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,NSMOB,IUSEAB,CJRES,SIRES,NORB,SCLFAC, &
                   S2_TERM1,IPACK)
! SUBROUTINE GSBBD2B --> 52
!
! alpha-beta contribution to two-particle density matrix
! from given c-block and s-block.
!
! S2_TERM1 = - <L!a+i alpha a+jbeta a i beta a j alpha !R>
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
! NSMOB  : Number of symmetries of orbitals
! MAXK   : Largest number of inner resolution strings treated at simult.
! IPACK  : Should we pack the density?
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
! I1, XI1S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! I2, XI2S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! X : Space for block of two-electron integrals
!
! Jeppe Olsen, Fall of 1996

use Symmetry_Info, only: Mul
use Para_Info, only: MyRank, nProcs
use lucia_data, only: LOFFI
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: RHO2(*), RHO2S(*), RHO2A(*), S2_TERM1
integer(kind=iwp), intent(in) :: IASM, IATP, IBSM, IBTP, NIA, NIB, JASM, JATP, JBSM, JBTP, NJA, NJB, IAGRP, IBGRP, NGAS, IAOC(*), &
                                 IBOC(*), JAOC(*), JBOC(*), MXPNGAS, NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*), MAXK, NSMOB, IUSEAB, NORB
integer(kind=iwp), intent(_OUT_) :: I1(*), I2(*), I3(*), I4(*)
real(kind=wp), intent(in) :: SB(*), CB(*), SCLFAC
real(kind=wp), intent(_OUT_) :: XI1S(*), XI2S(*), XI3S(*), XI4S(*), X(*), CJRES(*), SIRES(*)
logical(kind=iwp), intent(in) :: IPACK
integer(kind=iwp) :: I, IDOCOMP, II, IJ, IJSM, IJTYP, IKABTC, IKORD, IOFF, ISM, ITP(20), ITYP, J, JI, JJ, JOFF, JSM, JTP(20), &
                     JTYP, KABOT, KATOP, KLSM, KLTYP, KOFF, KSM, KTP(20), KTYP, LKABTC, LOFF, LSM, LTP(20), LTYP, NI, NIJ, NIJTYP, &
                     NJ, NK, NKABTC, NKABTCSZ, NKAEFF, NKASTR, NKBSTR, NKLTYP, NL
real(kind=wp), allocatable :: OFFI(:)

#ifdef _DEBUGPRINT_
write(u6,*) ' ================'
write(u6,*) ' GSBBD2B speaking'
write(u6,*) ' ================'
!write(u6,*) ' NJAS NJB = ',NJA,NJB
!write(u6,*) ' IAGRP IBGRP = ',IAGRP,IBGRP
!write(u6,*) ' MXPNGAS = ',MXPNGAS
!write(u6,*) ' NSMOB = ',NSMOB
#endif

! Symmetry of allowed excitations
IJSM = Mul(IASM,JASM)
KLSM = Mul(IBSM,JBSM)
if ((IJSM == 0) .or. (KLSM == 0)) return
#ifdef _DEBUGPRINT_
write(u6,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
write(u6,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
#endif
! Types of SX that connects the two strings
call SXTYP_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC)
call SXTYP_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC)
if ((NIJTYP == 0) .or. (NKLTYP == 0)) return
! Repeated allocation/deallocation inside ADSTN_GAS has been
! outerlooped to here. OFFI added to call parameters of
! ADSTN_GAS. PAM March 2006.
call mma_allocate(OFFI,lOFFI,Label='OFFI')
OFFI(:) = Zero

do IJTYP=1,NIJTYP
  ITYP = ITP(IJTYP)
  JTYP = JTP(IJTYP)
  do ISM=1,NSMOB
    JSM = Mul(ISM,IJSM)
    if (JSM == 0) cycle
#   ifdef _DEBUGPRINT_
    write(u6,*) ' ISM JSM ',ISM,JSM
#   endif
    IOFF = IOBPTS(ITYP,ISM)
    JOFF = IOBPTS(JTYP,JSM)
    NI = NOBPTS(ITYP,ISM)
    NJ = NOBPTS(JTYP,JSM)
    if ((NI == 0) .or. (NJ == 0)) cycle
    ! Generate annihilation mappings for all Ka strings
    ! a+j!ka> = +/-/0 * !Ja>
    call ADSTN_GAS(OFFI,JSM,JTYP,JATP,JASM,IAGRP,I1,XI1S,NKASTR,SCLFAC)
    if (NKASTR == 0) cycle
    ! a+i!ka> = +/-/0 * !Ia>
    call ADSTN_GAS(OFFI,ISM,ITYP,IATP,IASM,IAGRP,I3,XI3S,NKASTR,One)
    if (NKASTR == 0) cycle
    ! Compress list to common nonvanishing elements
    IDOCOMP = 1
    if (IDOCOMP == 1) then
      !    COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
      call COMPRS2LST(I1,XI1S,NJ,I3,XI3S,NI,NKASTR,NKAEFF)
    else
      NKAEFF = NKASTR
    end if

    ! Loop over batches of KA strings
    NKABTC = 0
    do
      NKABTC = NKABTC+NPROCS
      NKABTCSZ = max(NKAEFF-1,0)/NKABTC+1
      if (NKABTCSZ <= MAXK) exit
    end do

    do IKABTC=1+MYRANK,NKABTC,NPROCS
      KABOT = (IKABTC-1)*NKABTCSZ+1
      KATOP = min(KABOT+NKABTCSZ-1,NKAEFF)
      LKABTC = KATOP-KABOT+1
      if (LKABTC <= 0) exit
      ! Obtain C(ka,J,JB) for Ka in batch
      do JJ=1,NJ
        call GET_CKAJJB(CB,NJ,NJA,CJRES,LKABTC,NJB,JJ,I1(KABOT+(JJ-1)*NKASTR),XI1S(KABOT+(JJ-1)*NKASTR))
      end do
      ! Obtain S(ka,i,Ib) for Ka in batch
      do II=1,NI
        call GET_CKAJJB(SB,NI,NIA,SIRES,LKABTC,NIB,II,I3(KABOT+(II-1)*NKASTR),XI3S(KABOT+(II-1)*NKASTR))
      end do

      do KLTYP=1,NKLTYP
        KTYP = KTP(KLTYP)
        LTYP = LTP(KLTYP)

        do KSM=1,NSMOB
          LSM = Mul(KSM,KLSM)
          if (LSM == 0) cycle
          KOFF = IOBPTS(KTYP,KSM)
          LOFF = IOBPTS(LTYP,LSM)
          NK = NOBPTS(KTYP,KSM)
          NL = NOBPTS(LTYP,LSM)
          ! If IUSEAB is used, only terms with i >= k will be generated so
          IKORD = 0
          if ((IUSEAB == 1) .and. (ISM > KSM)) cycle
          if ((IUSEAB == 1) .and. (ISM == KSM) .and. (ITYP < KTYP)) cycle
          if ((IUSEAB == 1) .and. (ISM == KSM) .and. (ITYP == KTYP)) IKORD = 1

          if ((NK == 0) .or. (NL == 0)) cycle
          ! Obtain all connections a+l!Kb> = +/-/0!Jb>
          call ADSTN_GAS(OFFI,LSM,LTYP,JBTP,JBSM,IBGRP,I2,XI2S,NKBSTR,One)
          if (NKBSTR == 0) cycle
          ! Obtain all connections a+k!Kb> = +/-/0!Ib>
          call ADSTN_GAS(OFFI,KSM,KTYP,IBTP,IBSM,IBGRP,I4,XI4S,NKBSTR,One)
          if (NKBSTR == 0) cycle

          ! Update two-electron density matrix
          ! Rho2b(ij,kl) =  Sum(ka)S(Ka,i,Ib)<Ib!Eb(kl)!Jb>C(Ka,j,Jb)

          X(1:NI*NJ*NK*NL) = Zero

          call ABTOR2(SIRES,CJRES,LKABTC,NKBSTR,X,NI,NJ,NK,NL,NKBSTR,I4,XI4S,I2,XI2S,IKORD)
          ! contributions to Rho2(ij,kl) has been obtained, scatter out
          !call wrtmat(x,ni*nj,nk*nl,ni*nj,nk*nl)
          ! Contribution to S2
          if ((KTYP == JTYP) .and. (KSM == JSM) .and. (ITYP == LTYP) .and. (ISM == LSM)) then
            do I=1,NI
              do J=1,NJ
                IJ = (J-1)*NI+I
                JI = (I-1)*NJ+J
                NIJ = NI*NJ
                S2_TERM1 = S2_TERM1-X((JI-1)*NIJ+IJ)
              end do
            end do
          end if

          call ADTOR2(RHO2,RHO2S,RHO2A,X,2,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB,IPACK)

          !write(u6,*) ' updated density matrix B ','norb = ',norb
          !write(u6,*) ' offset ','IOFF,JOFF,KOFF,LOFF',IOFF,JOFF,KOFF,LOFF
          !call prsym(rho2s,nTri_Elem(NORB))

        end do
      end do
    end do
    ! End of loop over partitioning of alpha strings
  end do
end do
! This 'flush' outerlooped here. Was previously inside ADSTN_GAS.
call mma_deallocate(OFFI)

end subroutine GSBBD2B
