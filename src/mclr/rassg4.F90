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

subroutine RASSG4(C,S,CB,SB,C2,ICOCOC,ISOCOC,ICSMOS,ISSMOS,ICBLTP,ISBLTP,NORB1,NORB2,NORB3,NACOB,NSSOA,ISSOA,NSSOB,ISSOB,NAEL, &
                  IAGRP,NBEL,IBGRP,NOCTPA,NOCTPB,NSMST,NSMOB,NSMSX,NSMDX,NTSOB,IBTSOB,ITSOB,MAXIJ,MAXK,MAXI,ICSMOD,IINMOD,LI,LC, &
                  LS,XINT,CSCR,SSCR,SXSTSM,IAEL1,IAEL3,IBEL1,IBEL3,IDC,ISOOSC,NSOOSC,ISOOSE,NSOOSE,ICOOSC,NCOOSC,ICOOSE,NCOOSE, &
                  IASOOS,IACOOS,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,IDOH2,ISTRFL,PS,LUC,LUHC,IST,CJRES,SIRES,NOPARt,TimeDep)
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
! ICSMOS : Symmetry array for C
! ISSMOS : Symmetry array for S
! ICBLTP : Block types for C
! ISBLTP : Block types for S
!
! NORB1(2,3) : Number of orbitals in RAS1(2,3)
! NACOB : Number of active orbitals
! H     : Active one-body Hamiltonian with core contributions
! C     : CI vector
! CB    : Array able to hold largest STT block of C
! NSSOA : Number of strings per type and symmetry for alpha strings
! ISSOA : Offset for strings if given type and symmetry, alpha strings
! NAEL  : Number of active alpha electrons
! NSSOB : Number of strings per type and symmetry for beta strings
! ISSOB : Offset for strings if given type and symmetry, beta strings
! NBEL  : Number of active beta electrons
! NTSOB : Number of orbitals per type and symmetry
! ITSOB : Orbitals of given type and symmetry
! IBTSOB: Offset for ITSOB
!
! MAXIJ : Largest allowed number of orbital pairs treated simultaneously
! MAXK  : Largest number of N-2,N-1 strings treated simultaneously
! MAXI  : Max number of N strings treated simultaneously
!
! ICSMOD : 1 => Single symmetry blocks of C and S are treated
!               simultaneously
! ICSMOD : 2 => Single symmetry-occ-occ blocks of C and S are treated
!               simultaneously
! IINMOD :
!
! LI : Length of scratch array for integrals
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

use Constants, only: Zero
use DetDim, only: MXPORB, MXPOBS

implicit none
! General input
integer NAEL, IAGRP, NBEL, IBGRP, NOCTPA, NOCTPB
integer NSMST, NSMOB, NSMSX, NSMDX
integer MAXIJ, MAXK, MAXI, ICSMOD, IINMOD, LI, LC, LS
integer ICOCOC(NOCTPA,NOCTPB), ISOCOC(NOCTPA,NOCTPB)
integer ICSMOS(NSMST), ISSMOS(NSMST)
integer ICBLTP(NSMST), ISBLTP(NSMST)
integer NORB1, NORB2, NORB3, NACOB
integer NSSOA(NOCTPA,nsmst), ISSOA(NOCTPA,nsmst)
integer NSSOB(NOCTPB,nsmst), ISSOB(NOCTPB,nsmst)
integer NTSOB(3,NSMOB), IBTSOB(3,NSMOB), ITSOB(mxporb)
integer SXSTSM(NSMSX,NSMST)
integer IAEL1(*), IAEL3(*)
integer IBEL1(*), IBEL3(*)
integer IDC
integer IDOH2
real*8 PS
integer LUC, LUHC, IST, NOPART
logical TimeDep
! Scratch
real*8 C(*), S(*)
real*8 SB(*), CB(*), C2(*)
real*8 XINT(*), CSCR(*), SSCR(*)
integer ISOOSC(NOCTPA,NOCTPB,NSMST), NSOOSC(NOCTPA,NOCTPB,NSMST)
integer ISOOSE(NOCTPA,NOCTPB,NSMST), NSOOSE(NOCTPA,NOCTPB,NSMST)
integer ICOOSC(NOCTPA,NOCTPB,NSMST), NCOOSC(NOCTPA,NOCTPB,NSMST)
integer ICOOSE(NOCTPA,NOCTPB,NSMST), NCOOSE(NOCTPA,NOCTPB,NSMST)
integer IASOOS(NOCTPA,NOCTPB,NSMST), IACOOS(NOCTPA,NOCTPB,NSMST)
integer I1(MAXK,*), I2(MAXK,*), I3(MAXK,*), I4(MAXK,*)
real*8 XI1S(MAXK,*), XI2S(MAXK,*), XI3S(MAXK,*), XI4S(MAXK,*)
integer ISTRFL(*)
real*8 CJRES(*), SIRES(*)
! Local variables
integer LASM(4), LBSM(4), LATP(4), LBTP(4), LSGN(5), LTRP(5)
real*8 PL
integer ISENSM, ISENTA, ISENTB, IFRSTS, ISSTSM, ISSTTA, NSBLK, ISFINI, IS1SM, IS1TA, IS1TB, LSBLK, ISBLK, ICENSM, ICENTA, ICENTB, &
        IFRSTC, ICSTSM, ICSTTA, ICSTTB, NCBLK, IFINIC, IC1SM, IC1TA, IC1TB, ICOFF, ICBLK, ICBSM, ISOFF, IASM, IBSM, IATP, IBTP, &
        NIA, NIB, JASM, JATP, NJA, NJB, IPERM, NPERM, LLASM, LLBSM, LLATP, LLBTP, NLLA, NLLB, LROW, LCOL, I1ASM, I1BSM, I1TA, &
        I1TB, IOFF, LBLK, NCCMBC, NCCMBE, NONEWC, NONEWS, NSCMBC, NSCMBE, ISSTTB, JBTP, JBSM, DUM(1)
real*8 XNORM2
real*8, external :: DDot_

PL = Zero
! ================================
! 1 : Arrays for accessing C and S
! ================================

!***********************************************************************

! Sigma, compact form
call ZOOS(ISSMOS,ISBLTP,NSMST,ISOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,IDC,ISOOSC,NSOOSC,NSCMBC,0)
! Sigma with expanded diagonal blocks
call ZOOS(ISSMOS,ISBLTP,NSMST,ISOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,IDC,ISOOSE,NSOOSE,NSCMBE,1)
! C, compact form
call ZOOS(ICSMOS,ICBLTP,NSMST,ISOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,IDC,ICOOSC,NCOOSC,NCCMBC,0)
! C, Full determinant form
call ZOOS(ICSMOS,ICBLTP,NSMST,ISOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,1,ICOOSE,NCOOSE,NCCMBE,1)

!************************************************************************

! Initialize loop over batches of sigma blocks
ISENSM = 1
ISENTA = 1
ISENTB = 1
IFRSTS = 1
! Loop over batches over sigma blocks
if (LUHC > 0) rewind LUHC

10001 continue
! Next batch of sigma blocks
ISSTSM = ISENSM
ISSTTA = ISENTA
ISSTTB = ISENTB
call INCOOS(IDC,ISBLTP,NSOOSE,NOCTPA,NOCTPB,ISSTSM,ISSTTA,ISSTTB,NSMST,ISENSM,ISENTA,ISENTB,IASOOS,LS,ISFINI,NSBLK,IFRSTS,ISOCOC)
if ((NSBLK == 0) .and. (ISFINI /= 0)) goto 10002
IFRSTS = 0
! Initialize sigma blocks
IS1SM = ISSTSM
IS1TA = ISSTTA
IS1TB = ISSTTB
LSBLK = 0
do ISBLK=1,NSBLK
  LSBLK = LSBLK+NSOOSE(IS1TA,IS1TB,IS1SM)
  if (ISBLK /= NSBLK) call NXTBLK_MCLR(IS1TA,IS1TB,IS1SM,NOCTPA,NOCTPB,NSMST,ISBLTP,IDC,NONEWS,ISOCOC)
end do
SB(1:LSBLK) = Zero
! Initialize loop over blocks over C vector
ICENSM = 1
ICENTA = 1
ICENTB = 1
if (LUC > 0) rewind LUC
! Loop over blocks of C vector
IFRSTC = 1
9001 continue
ICSTSM = ICENSM
ICSTTA = ICENTA
ICSTTB = ICENTB

call INCOOS(IDC,ICBLTP,NCOOSE,NOCTPA,NOCTPB,ICSTSM,ICSTTA,ICSTTB,NSMST,ICENSM,ICENTA,ICENTB,IACOOS,LC,IFINIC,NCBLK,IFRSTC,ICOCOC)

! If no more C blocks goto next batch of sigma blocks
if ((NCBLK == 0) .and. (IFINIC /= 0)) goto 10001
IFRSTC = 0
! Read C blocks into core
IC1SM = ICSTSM ! Symmetry alpha
IC1TA = ICSTTA ! Type alpha string
IC1TB = ICSTTB ! Type Beta string
ICOFF = 1
do ICBLK=1,NCBLK
  ICBSM = ICSMOS(IC1SM) ! Symmetry Beta string
  ! C CI vector in
  ! CB Block of CI vector out
  if (ICOCOC(IC1TA,IC1TB) == 1) &
    call GSTTBL_MCLR(C,CB(ICOFF),IC1TA,IC1SM,IC1TB,ICBSM,ICOCOC,NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ICOOSC,IDC,PL,LUC,C2)

  ICOFF = ICOFF+NCOOSE(IC1TA,IC1TB,IC1SM)
  if (ICBLK /= NCBLK) call NXTBLK_MCLR(IC1TA,IC1TB,IC1SM,NOCTPA,NOCTPB,NSMST,ICBLTP,IDC,NONEWC,ICOCOC)
end do

!***********************************************************************

! Loop over sigma and C blocks in core and obtain  contribution from
! given C block to given S block
ISOFF = 1
IASM = ISSTSM
IATP = ISSTTA
IBTP = ISSTTB
do ISBLK=1,NSBLK
  call xflush(6)
  IBSM = ISSMOS(IASM)
  NIA = NSSOA(IATP,IASM)
  NIB = NSSOB(IBTP,IBSM)
  if ((NIA /= 0) .and. (NIB /= 0)) then
    JASM = ICSTSM
    JATP = ICSTTA
    JBTP = ICSTTB
    ICOFF = 1
    do ICBLK=1,NCBLK
      call xflush(6)
      JBSM = ICSMOS(JASM)
      NJA = NSSOA(JATP,JASM)
      NJB = NSSOB(JBTP,JBSM)
      XNORM2 = DDot_(NJA*NJB,CB(ICOFF),1,CB(ICOFF),1)
      if ((NIA /= 0) .and. (NIB /= 0) .and. (NJA /= 0) .and. (NJB /= 0) .and. (ISOCOC(IATP,IBTP) == 1) .and. &
          (ICOCOC(JATP,JBTP) == 1) .and. (XNORM2 /= 0.0d0)) then
        ! Other symmetry blocks that can be obtained from this block
        !write(6,*) 'Other symmetry blocks that can be obtained'
        !call xflush(6)

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
            call dcopy_(LROW*LCOL,C2,1,CB(iCOFF),1)
          end if
          if (LSGN(IPERM) == -1) call DSCAL_(LROW*LCOL,-1.0d0,CB(ICOFF),1)

          ! Generation of contribution to sigma block
          ! from given CI block

          !write(6,*) 'TimeDep in rassg4',TimeDep
          !call xflush(6)
          if (TimeDep) then
            !write(6,*) 'I call rssbcbn_td'
            call RSSBCBN_td(IASM,IATP,IBSM,IBTP,LLASM,LLATP,LLBSM,LLBTP,IAEL1(IATP),IAEL3(IATP),IBEL1(IBTP),IBEL3(IBTP), &
                            IAEL1(LLATP),IAEL3(LLATP),IBEL1(LLBTP),IBEL3(LLBTP),NAEL,NBEL,IAGRP,IBGRP,SB(ISOFF),CB(ICOFF),IDOH2, &
                            NTSOB,IBTSOB,ITSOB,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,C2,NSMOB,NSMST,NSMSX, &
                            NSMDX,NIA,NIB,NLLA,NLLB,MXPOBS,IST,CJRES,SIRES,NOPART,TimeDep)

            !call RECPRT('SSCR in rassg4',' ',SSCR,5,1) !yma
            !call xflush(6)
          else
            call RSSBCBN_MCLR(IASM,IATP,IBSM,IBTP,LLASM,LLATP,LLBSM,LLBTP,IAEL1(IATP),IAEL3(IATP),IBEL1(IBTP),IBEL3(IBTP), &
                              IAEL1(LLATP),IAEL3(LLATP),IBEL1(LLBTP),IBEL3(LLBTP),NAEL,NBEL,IAGRP,IBGRP,SB(ISOFF),CB(ICOFF),IDOH2, &
                              NTSOB,IBTSOB,ITSOB,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,C2,NSMOB,NSMST,NSMSX, &
                              NSMDX,NIA,NIB,NLLA,NLLB,MXPOBS,IST,CJRES,SIRES,NOPART,TimeDep)
          end if

        end do
        ! Transpose or scale to restore order ??
        if (LTRP(NPERM+1) == 1) then
          call TRNSPS(NJB,NJA,CB(ICOFF),C2)
          call dcopy_(NJA*NJB,C2,1,CB(ICOFF),1)
        end if
        if (LSGN(NPERM+1) == -1) call DSCAL_(NJA*NJB,-1.0d0,CB(ICOFF),1)

      end if
      ICOFF = ICOFF+NCOOSE(JATP,JBTP,JASM)
      ! NeXt C block
      if (ICBLK /= NCBLK) call NXTBLK_MCLR(JATP,JBTP,JASM,NOCTPA,NOCTPB,NSMST,ICBLTP,IDC,NONEWC,ICOCOC)
    end do
    ! End of loop over C blocks in Batch
    ! NeXt S block
  end if
  ISOFF = ISOFF+NSOOSE(IATP,IBTP,IASM)
  if (ISBLK /= NSBLK) call NXTBLK_MCLR(IATP,IBTP,IASM,NOCTPA,NOCTPB,NSMST,ISBLTP,IDC,NONEWS,ISOCOC)
end do
!***********************************************************************
! End of loop over S blocks in batch
! End of loop over batches of C blocks

if (IFINIC == 0) goto 9001
! Transfer S block to permanent storage

!call RECPRT('SB in rassg4',' ',S,5,1)

I1ASM = ISSTSM
I1TA = ISSTTA
I1TB = ISSTTB
IOFF = 1
do ISBLK=1,NSBLK
  I1BSM = ISSMOS(I1ASM)
  if (ISOCOC(I1TA,I1TB) == 1) &
    call PSTTBL_MCLR(S,SB(IOFF),I1TA,I1ASM,I1TB,I1BSM,ISOCOC,NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ISOOSC,2,IDC,LUHC,C2)
  IOFF = IOFF+NSOOSE(I1TA,I1TB,I1ASM)
  if (ISBLK /= NSBLK) call NXTBLK_MCLR(I1TA,I1TB,I1ASM,NOCTPA,NOCTPB,NSMST,ISBLTP,IDC,NONEWS,ISOCOC)
end do
if (ISFINI == 0) goto 10001
! End of loop over batches of S blocks
10002 continue
!***********************************************************************
if (LUHC > 0) then
  DUM(1) = -1
  call ITODS(DUM,1,LBLK,LUHC)
end if

! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(NORB1)
  call Unused_integer(NORB2)
  call Unused_integer(NORB3)
  call Unused_integer(NACOB)
  call Unused_integer_array(ISSOA)
  call Unused_integer_array(ISSOB)
  call Unused_integer(MAXIJ)
  call Unused_integer(ICSMOD)
  call Unused_integer(IINMOD)
  call Unused_integer(LI)
  call Unused_integer_array(SXSTSM)
end if

end subroutine RASSG4
