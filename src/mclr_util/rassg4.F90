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
!               1996, Anders Bernhardsson                              *
!***********************************************************************

subroutine RASSG4(C,S,CB,SB,C2,ICOCOC,ISOCOC,ICSM,ISSM,ICBLTP,ISBLTP,NSSOA,NSSOB,NAEL,IAGRP,NBEL,IBGRP,NOCTPA,NOCTPB,NSM,NTSOB, &
                  IBTSOB,ITSOB,MAXK,MAXI,LC,LS,XINT,CSCR,SSCR,IAEL1,IAEL3,IBEL1,IBEL3,IDC,ISOOSC,NSOOSC,ISOOSE,NSOOSE,ICOOSC, &
                  NCOOSC,ICOOSE,NCOOSE,IASOOS,IACOOS,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,IDOH2,ISTRFL,PS,LUC,LUHC,IST,CJRES,SIRES, &
                  NOPART,TimeDep)
! LOOP OVER SIGMA AND C VECTOR
!
! Jeppe Olsen, Winter of 1991
! small modifications by eaw 96
!
! =====
! Input
! =====
!
! ICOCOC : Allowed type combinations for C
! ISOCOC : Allowed type combinations for S(igma)
! ICS    : Symmetry for C
! ISS    : Symmetry for S
! ICBLTP : Block types for C
! ISBLTP : Block types for S
!
! H     : Active one-body Hamiltonian with core contributions
! C     : CI vector
! CB    : Array able to hold largest STT block of C
! NSSOA : Number of strings per type and symmetry for alpha strings
! NAEL  : Number of active alpha electrons
! NSSOB : Number of strings per type and symmetry for beta strings
! NBEL  : Number of active beta electrons
! NTSOB : Number of orbitals per type and symmetry
! ITSOB : Orbitals of given type and symmetry
! IBTSOB: Offset for ITSOB
!
! MAXK  : Largest number of N-2,N-1 strings treated simultaneously
! MAXI  : Max number of N strings treated simultaneously
!
! LC : Length of scratch array for C
! LS : Length of scratch array for S
! XINT : Scratch array for integrals
! CSCR : Scratch array for C vector (space for res. matr.)
! SSCR : Scratch array for S vector (space for res. matr.)
!
! The C and S vectors are accessed through routines that
! either fetches/disposes symmetry blocks or
! Symmetry-occupation-occupation blocks
!
! IST :  = 1 => Singlet operator
!        = 2 => Triplet operator\
!
! IDOH2 : = 1 => both one and two particle parts
!         = 0 => only one-electron operator
! A triplet one electron operator is defined as E(aa)-E(bb)
! A triplet two-electron operator is defined as (E(aa)+E(bb))(E(aa)-E(bb))

use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: C(*), PS
real(kind=wp), intent(inout) :: S(*)
integer(kind=iwp), intent(in) :: NOCTPA, NOCTPB, ICOCOC(NOCTPA,NOCTPB), ISOCOC(NOCTPA,NOCTPB), ICSM, ISSM, NSM, ICBLTP(NSM), &
                                 ISBLTP(NSM), NSSOA(NOCTPA,NSM), NSSOB(NOCTPB,NSM), NAEL, IAGRP, NBEL, IBGRP, NTSOB(3,NSM), &
                                 IBTSOB(3,NSM), ITSOB(*), MAXK, MAXI, LC, LS, IAEL1(*), IAEL3(*), IBEL1(*), IBEL3(*), IDC, IDOH2, &
                                 ISTRFL(*), LUC, LUHC, IST, NOPART
real(kind=wp), intent(_OUT_) :: CB(*), SB(*), C2(*), XINT(*), CSCR(*), SSCR(*), XI1S(MAXK,*), XI2S(MAXK,*), XI3S(MAXK,*), &
                                XI4S(MAXK,*), CJRES(*), SIRES(*)
integer(kind=iwp), intent(out) :: ISOOSC(NOCTPA,NOCTPB,NSM), NSOOSC(NOCTPA,NOCTPB,NSM), ISOOSE(NOCTPA,NOCTPB,NSM), &
                                  NSOOSE(NOCTPA,NOCTPB,NSM), ICOOSC(NOCTPA,NOCTPB,NSM), NCOOSC(NOCTPA,NOCTPB,NSM), &
                                  ICOOSE(NOCTPA,NOCTPB,NSM), NCOOSE(NOCTPA,NOCTPB,NSM), IASOOS(NOCTPA,NOCTPB,NSM), &
                                  IACOOS(NOCTPA,NOCTPB,NSM)
integer(kind=iwp), intent(_OUT_) :: I1(MAXK,*), I2(MAXK,*), I3(MAXK,*), I4(MAXK,*)
logical(kind=iwp), intent(in) :: TimeDep
integer(kind=iwp) :: DUM(1), I1ASM, I1BSM, I1TA, I1TB, IASM, IATP, IBSM, IBTP, IC1SM, IC1TA, IC1TB, ICBLK, ICBSM, ICENSM, ICENTA, &
                     ICENTB, ICOFF, ICSTSM, ICSTTA, ICSTTB, IFINIC, IFRSTC, IFRSTS, IOFF, IPERM, IS1SM, IS1TA, IS1TB, ISBLK, &
                     ISENSM, ISENTA, ISENTB, ISFINI, ISOFF, ISSTSM, ISSTTA, ISSTTB, JASM, JATP, JBSM, JBTP, LASM(4), LATP(4), &
                     LBLK, LBSM(4), LBTP(4), LCOL, LLASM, LLATP, LLBSM, LLBTP, LROW, LSBLK, LSGN(5), LTRP(5), NCBLK, NCCMBC, &
                     NCCMBE, NIA, NIB, NJA, NJB, NLLA, NLLB, NONEWC, NONEWS, NPERM, NSBLK, NSCMBC, NSCMBE
real(kind=wp) :: PL, XNORM2
real(kind=wp), external :: DDot_

PL = Zero
! ================================
! 1 : Arrays for accessing C and S
! ================================

!***********************************************************************

! Sigma, compact form
call ZOOS(ISSM,ISBLTP,NSM,ISOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,IDC,ISOOSC,NSOOSC,NSCMBC,0)
! Sigma with expanded diagonal blocks
call ZOOS(ISSM,ISBLTP,NSM,ISOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,IDC,ISOOSE,NSOOSE,NSCMBE,1)
! C, compact form
call ZOOS(ICSM,ICBLTP,NSM,ISOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,IDC,ICOOSC,NCOOSC,NCCMBC,0)
! C, Full determinant form
call ZOOS(ICSM,ICBLTP,NSM,ISOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,1,ICOOSE,NCOOSE,NCCMBE,1)

!************************************************************************

! Initialize loop over batches of sigma blocks
ISENSM = 1
ISENTA = 1
ISENTB = 1
IFRSTS = 1
! Loop over batches over sigma blocks
if (LUHC > 0) rewind(LUHC)

outer: do
  ! Next batch of sigma blocks
  ISSTSM = ISENSM
  ISSTTA = ISENTA
  ISSTTB = ISENTB
  call INCOOS(IDC,ISBLTP,NSOOSE,NOCTPA,NOCTPB,ISSTSM,ISSTTA,ISSTTB,NSM,ISENSM,ISENTA,ISENTB,IASOOS,LS,ISFINI,NSBLK,IFRSTS,ISOCOC)
  if ((NSBLK == 0) .and. (ISFINI /= 0)) exit
  IFRSTS = 0
  ! Initialize sigma blocks
  IS1SM = ISSTSM
  IS1TA = ISSTTA
  IS1TB = ISSTTB
  LSBLK = 0
  do ISBLK=1,NSBLK
    LSBLK = LSBLK+NSOOSE(IS1TA,IS1TB,IS1SM)
    if (ISBLK /= NSBLK) call NXTBLK_MCLR(IS1TA,IS1TB,IS1SM,NOCTPA,NOCTPB,NSM,ISBLTP,IDC,NONEWS,ISOCOC)
  end do
  SB(1:LSBLK) = Zero
  ! Initialize loop over blocks over C vector
  ICENSM = 1
  ICENTA = 1
  ICENTB = 1
  if (LUC > 0) rewind(LUC)
  ! Loop over blocks of C vector
  IFRSTC = 1
  do
    ICSTSM = ICENSM
    ICSTTA = ICENTA
    ICSTTB = ICENTB

    call INCOOS(IDC,ICBLTP,NCOOSE,NOCTPA,NOCTPB,ICSTSM,ICSTTA,ICSTTB,NSM,ICENSM,ICENTA,ICENTB,IACOOS,LC,IFINIC,NCBLK,IFRSTC,ICOCOC)

    ! If no more C blocks go to next batch of sigma blocks
    if ((NCBLK == 0) .and. (IFINIC /= 0)) cycle outer
    IFRSTC = 0
    ! Read C blocks into core
    IC1SM = ICSTSM ! Symmetry alpha
    IC1TA = ICSTTA ! Type alpha string
    IC1TB = ICSTTB ! Type Beta string
    ICOFF = 1
    do ICBLK=1,NCBLK
      ICBSM = Mul(IC1SM,ICSM) ! Symmetry Beta string
      ! C CI vector in
      ! CB Block of CI vector out
      if (ICOCOC(IC1TA,IC1TB) == 1) &
        call GSTTBL_MCLR(C,CB(ICOFF),IC1TA,IC1SM,IC1TB,ICBSM,NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ICOOSC,IDC,PL,LUC,C2)

      ICOFF = ICOFF+NCOOSE(IC1TA,IC1TB,IC1SM)
      if (ICBLK /= NCBLK) call NXTBLK_MCLR(IC1TA,IC1TB,IC1SM,NOCTPA,NOCTPB,NSM,ICBLTP,IDC,NONEWC,ICOCOC)
    end do

    !*******************************************************************

    ! Loop over sigma and C blocks in core and obtain  contribution from
    ! given C block to given S block
    ISOFF = 1
    IASM = ISSTSM
    IATP = ISSTTA
    IBTP = ISSTTB
    do ISBLK=1,NSBLK
      IBSM = Mul(IASM,ISSM)
      NIA = NSSOA(IATP,IASM)
      NIB = NSSOB(IBTP,IBSM)
      if ((NIA /= 0) .and. (NIB /= 0)) then
        JASM = ICSTSM
        JATP = ICSTTA
        JBTP = ICSTTB
        ICOFF = 1
        do ICBLK=1,NCBLK
          JBSM = Mul(JASM,ICSM)
          NJA = NSSOA(JATP,JASM)
          NJB = NSSOB(JBTP,JBSM)
          XNORM2 = DDot_(NJA*NJB,CB(ICOFF),1,CB(ICOFF),1)
          if ((NIA /= 0) .and. (NIB /= 0) .and. (NJA /= 0) .and. (NJB /= 0) .and. (ISOCOC(IATP,IBTP) == 1) .and. &
              (ICOCOC(JATP,JBTP) == 1) .and. (XNORM2 /= Zero)) then
            ! Other symmetry blocks that can be obtained from this block
            !write(u6,*) 'Other symmetry blocks that can be obtained'
            !call xflush(u6)

            call PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
            do IPERM=1,NPERM
              LLASM = LASM(IPERM)
              LLBSM = LBSM(IPERM)
              LLATP = LATP(IPERM)
              LLBTP = LBTP(IPERM)
              NLLA = NSSOA(LLATP,LLASM)
              NLLB = NSSOB(LLBTP,LLBSM)
              if (LTRP(IPERM) == 1) then
                LROW = NSSOA(LATP(IPERM-1),LASM(IPERM-1))
                LCOL = NSSOB(LBTP(IPERM-1),LBSM(IPERM-1))
                call TRNSPS(LROW,LCOL,CB(ICOFF),C2)
                CB(iCOFF:iCOFF+LROW*LCOL-1) = C2(1:LROW*LCOL)
              end if
              if (LSGN(IPERM) == -1) CB(iCOFF:iCOFF+LROW*LCOL-1) = -CB(iCOFF:iCOFF+LROW*LCOL-1)

              ! Generation of contribution to sigma block
              ! from given CI block

              !write(u6,*) 'TimeDep in rassg4',TimeDep
              !call xflush(u6)
              if (TimeDep) then
                !write(u6,*) 'I call rssbcbn_td'
                call RSSBCBN_td(IASM,IATP,IBSM,IBTP,LLASM,LLATP,LLBSM,LLBTP,IAEL1(IATP),IAEL3(IATP),IBEL1(IBTP),IBEL3(IBTP), &
                                IAEL1(LLATP),IAEL3(LLATP),IBEL1(LLBTP),IBEL3(LLBTP),NAEL,NBEL,IAGRP,IBGRP,SB(ISOFF),CB(ICOFF), &
                                IDOH2,NTSOB,IBTSOB,ITSOB,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,C2,NSM,NIA,NIB, &
                                NLLA,NLLB,IST,CJRES,SIRES,NOPART,TimeDep)

                !call RECPRT('SSCR in rassg4',' ',SSCR,5,1) !yma
                !call xflush(u6)
              else
                call RSSBCBN_MCLR(IASM,IATP,IBSM,IBTP,LLASM,LLATP,LLBSM,LLBTP,IAEL1(IATP),IAEL3(IATP),IBEL1(IBTP),IBEL3(IBTP), &
                                  IAEL1(LLATP),IAEL3(LLATP),IBEL1(LLBTP),IBEL3(LLBTP),NAEL,NBEL,IAGRP,IBGRP,SB(ISOFF),CB(ICOFF), &
                                  IDOH2,NTSOB,IBTSOB,ITSOB,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,C2,NSM,NIA, &
                                  NIB,NLLA,NLLB,IST,CJRES,SIRES,NOPART,TimeDep)
              end if

            end do
            ! Transpose or scale to restore order ??
            if (LTRP(NPERM+1) == 1) then
              call TRNSPS(NJB,NJA,CB(ICOFF),C2)
              CB(ICOFF:ICOFF+NJA*NJB-1) = C2(1:NJA*NJB)
            end if
            if (LSGN(NPERM+1) == -1) CB(ICOFF:ICOFF+NJA*NJB-1) = -CB(ICOFF:ICOFF+NJA*NJB-1)

          end if
          ICOFF = ICOFF+NCOOSE(JATP,JBTP,JASM)
          ! NeXt C block
          if (ICBLK /= NCBLK) call NXTBLK_MCLR(JATP,JBTP,JASM,NOCTPA,NOCTPB,NSM,ICBLTP,IDC,NONEWC,ICOCOC)
        end do
        ! End of loop over C blocks in Batch
        ! NeXt S block
      end if
      ISOFF = ISOFF+NSOOSE(IATP,IBTP,IASM)
      if (ISBLK /= NSBLK) call NXTBLK_MCLR(IATP,IBTP,IASM,NOCTPA,NOCTPB,NSM,ISBLTP,IDC,NONEWS,ISOCOC)
    end do
    !*******************************************************************
    ! End of loop over S blocks in batch
    ! End of loop over batches of C blocks

    if (IFINIC /= 0) exit
  end do
  ! Transfer S block to permanent storage

  !call RECPRT('SB in rassg4',' ',S,5,1)

  I1ASM = ISSTSM
  I1TA = ISSTTA
  I1TB = ISSTTB
  IOFF = 1
  do ISBLK=1,NSBLK
    I1BSM = Mul(I1ASM,ISSM)
    if (ISOCOC(I1TA,I1TB) == 1) &
      call PSTTBL_MCLR(S,SB(IOFF),I1TA,I1ASM,I1TB,I1BSM,NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ISOOSC,2,IDC,LUHC,C2)
    IOFF = IOFF+NSOOSE(I1TA,I1TB,I1ASM)
    if (ISBLK /= NSBLK) call NXTBLK_MCLR(I1TA,I1TB,I1ASM,NOCTPA,NOCTPB,NSM,ISBLTP,IDC,NONEWS,ISOCOC)
  end do
  if (ISFINI /= 0) exit outer
end do outer
! End of loop over batches of S blocks
!***********************************************************************
if (LUHC > 0) then
  DUM(1) = -1
  call ITODS(DUM,1,LBLK,LUHC)
end if

end subroutine RASSG4
