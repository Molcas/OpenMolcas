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
! Copyright (C) 1991,1995, Jeppe Olsen                                 *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine GASDN2_MCLR(I12,RHO1,RHO2,R,L,CB,SB,C2,ICOCOC,ISOCOC,ICSM,ISSM,ICBLTP,ISBLTP,NACOB,NSSOA,NSSOB,NAEL,IAGRP,NBEL,IBGRP, &
                       IOCTPA,IOCTPB,NOCTPA,NOCTPB,NSM,MXPNGAS,NOBPTS,IOBPTS,MAXK,MAXI,LC,LS,CSCR,SSCR,NGAS,NELFSPGPA,NELFSPGPB, &
                       IDC,ISOOSC,NSOOSC,ISOOSE,NSOOSE,ICOOSC,NCOOSC,ICOOSE,NCOOSE,IASOOS,IACOOS,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S, &
                       X,RHO1S,LUL,LUR,PSL,PSR,ieaw,n1,n2)
! Jeppe Olsen, Winter of 1991
! GAS modifications, August 1995
!
! =====
! Input
! =====
!
! I12    : = 1 => calculate one-electron density matrix
!          = 2 => calculate one-and two-electron density matrix
! RHO1   : Initial one-electron density matrix
! RHO2   : Initial two-electron density matrix
!
! ICOCOC : Allowed type combinations for C
! ISOCOC : Allowed type combinations for S(igma)
! ICS    : Symmetry for C
! ISS    : Symmetry for S
! ICBLTP : Block types for C
! ISBLTP : Block types for S
!
! NACOB : Number of active orbitals
! NSSOA : Number of strings per type and symmetry for alpha strings
! NAEL  : Number of active alpha electrons
! NSSOB : Number of strings per type and symmetry for beta strings
! NBEL  : Number of active beta electrons
!
! MAXIJ : Largest allowed number of orbital pairs treated simultaneously
! MAXK  : Largest number of N-2,N-1 strings treated simultaneously
! MAXI  : Max number of N strings treated simultaneously
!
! LC : Length of scratch array for C
! LS : Length of scratch array for S
! RHO1S: Scratch array for one body
! CSCR : Scratch array for C vector
! SSCR : Scratch array for S vector
!
! The L and R vectors are accessed through routines that
! either fetches/disposes symmetry blocks or
! Symmetry-occupation-occupation blocks

use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: I12, NOCTPA, NOCTPB, ICOCOC(NOCTPA,NOCTPB), ISOCOC(NOCTPA,NOCTPB), ICSM, ISSM, ICBLTP(*), &
                                 ISBLTP(*), NACOB, NSM, NSSOA(NSM,NOCTPA), NSSOB(NSM,NOCTPB), NAEL, IAGRP, NBEL, IBGRP, IOCTPA, &
                                 IOCTPB, MXPNGAS, NOBPTS(*), IOBPTS(*), MAXK, MAXI, LC, LS, NELFSPGPA(3,*), NELFSPGPB(3,*), IDC, &
                                 LUL, LUR, ieaw, n1, n2
real(kind=wp), intent(inout) :: RHO1(*), RHO2(*)
real(kind=wp), intent(in) :: R(*), L(*), PSL, PSR
real(kind=wp), intent(out) :: CB(*), SB(*), C2(*), CSCR(*), SSCR(*), XI1S(*), XI2S(*), XI3S(*), XI4S(*), X(*), RHO1S(*)
integer(kind=iwp), intent(inout) :: NGAS
integer(kind=iwp), intent(out) :: ISOOSC(NOCTPA,NOCTPB,NSM), NSOOSC(NOCTPA,NOCTPB,NSM), ISOOSE(NOCTPA,NOCTPB,NSM), &
                                  NSOOSE(NOCTPA,NOCTPB,NSM), ICOOSC(NOCTPA,NOCTPB,NSM), NCOOSC(NOCTPA,NOCTPB,NSM), &
                                  ICOOSE(NOCTPA,NOCTPB,NSM), NCOOSE(NOCTPA,NOCTPB,NSM), IASOOS(NOCTPA,NOCTPB,NSM), &
                                  IACOOS(NOCTPA,NOCTPB,NSM)
integer(kind=iwp), intent(_OUT_) :: I1(*), I2(*), I3(*), I4(*)
integer(kind=iwp) :: IASM, IATP, IBSM, IBTP, IC1SM, IC1TA, IC1TB, ICBLK, ICBSM, ICENSM, ICENTA, ICENTB, ICOFF, ICSTSM, ICSTTA, &
                     ICSTTB, IFINIC, IFRSTC, IFRSTS, IIASM, IIATP, IIBSM, IIBTP, ILPERM, IRPERM, IS1SM, IS1TA, IS1TB, ISBLK, &
                     ISBSM, ISENSM, ISENTA, ISENTB, ISFINI, ISOFF, ISSTSM, ISSTTA, ISSTTB, ISTRFL(1), JASM, JATP, JBSM, JBTP, &
                     JJASM, JJATP, JJBSM, JJBTP, LASM(4), LATP(4), LBSM(4), LBTP(4), LCOL, LROW, LSGN(5), LTRP(5), NCBLK, NCCMBC, &
                     NCCMBE, NIA, NIB, NIIA, NIIB, NJA, NJB, NJJA, NJJB, NLPERM, NONEWC, NONEWS, NRPERM, NSBLK, NSCMBC, NSCMBE, &
                     RASM(4), RATP(4), RBSM(4), RBTP(4), RSGN(5), RTRP(5)
real(kind=wp) :: PLL, PLR, XNORM2
real(kind=wp), external :: dDot_

PLL = Zero
PLR = Zero

! ================================
! 1 : Arrays for accessing L and R
! ================================
! L Compact form
call ZOOS(ISSM,ISBLTP,NSM,ISOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,IDC,ISOOSC,NSOOSC,NSCMBC,0)
! L Full determinant form
call ZOOS(ISSM,ISBLTP,NSM,ISOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,1,ISOOSE,NSOOSE,NSCMBE,1)
! R compact form
call ZOOS(ICSM,ICBLTP,NSM,ICOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,IDC,ICOOSC,NCOOSC,NCCMBC,0)
! R Full determinant form
call ZOOS(ICSM,ICBLTP,NSM,ICOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,1,ICOOSE,NCOOSE,NCCMBE,1)
! Initialize loop over batches of L blocks
ISENSM = 1
ISENTA = 1
ISENTB = 1
IFRSTS = 1
! Loop over batches over L blocks
outer: do
  ! Next batch of L blocks
  ISSTSM = ISENSM
  ISSTTA = ISENTA
  ISSTTB = ISENTB
  call INCOOS(IDC,ISBLTP,NSOOSE,NOCTPA,NOCTPB,ISSTSM,ISSTTA,ISSTTB,NSM,ISENSM,ISENTA,ISENTB,IASOOS,LS,ISFINI,NSBLK,IFRSTS,ISOCOC)
  if ((NSBLK == 0) .and. (ISFINI /= 0)) exit outer
  IFRSTS = 0
  ! Obtain L blocks
  IS1SM = ISSTSM
  IS1TA = ISSTTA
  IS1TB = ISSTTB
  ISOFF = 1
  do ISBLK=1,NSBLK
    ISBSM = Mul(IS1SM,ISSM)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' ISBLK ISOFF ',ISBLK,ISOFF
#   endif
    if (ISOCOC(IS1TA,IS1TB) == 1) &
      call GSTTBL_MCLR(L,SB(ISOFF),IS1TA,IS1SM,IS1TB,ISBSM,NOCTPA,NOCTPB,NSSOA,NSSOB,PSL,ISOOSC,IDC,PLL,LUL,C2)
    ISOFF = ISOFF+NSOOSE(IS1TA,IS1TB,IS1SM)
    if (ISBLK /= NSBLK) call NXTBLK_MCLR(IS1TA,IS1TB,IS1SM,NOCTPA,NOCTPB,NSM,ISBLTP,IDC,NONEWS,ISOCOC)
  end do
  ! Initialize loop over blocks over L vector
  ICENSM = 1
  ICENTA = 1
  ICENTB = 1
  ! Loop over blocks of R vector
  IFRSTC = 1
  do
#   ifdef _DEBUGPRINT_
    write(u6,*) ' >>> next batch of R blocks'
#   endif
    ICSTSM = ICENSM
    ICSTTA = ICENTA
    ICSTTB = ICENTB

    call INCOOS(IDC,ICBLTP,NCOOSE,NOCTPA,NOCTPB,ICSTSM,ICSTTA,ICSTTB,NSM,ICENSM,ICENTA,ICENTB,IACOOS,LC,IFINIC,NCBLK,IFRSTC,ICOCOC)
    ! If no more R blocks go to next batch of L blocks
    if ((NCBLK == 0) .and. (IFINIC /= 0)) cycle outer
    IFRSTC = 0
    ! Read L blocks into core
    IC1SM = ICSTSM
    IC1TA = ICSTTA
    IC1TB = ICSTTB
    ICOFF = 1
    do ICBLK=1,NCBLK
      ICBSM = Mul(IC1SM,ICSM)
      if (ICOCOC(IC1TA,IC1TB) == 1) &
        call GSTTBL_MCLR(R,CB(ICOFF),IC1TA,IC1SM,IC1TB,ICBSM,NOCTPA,NOCTPB,NSSOA,NSSOB,PSR,ICOOSC,IDC,PLR,LUR,C2)
      ICOFF = ICOFF+NCOOSE(IC1TA,IC1TB,IC1SM)
      if (ICBLK /= NCBLK) call NXTBLK_MCLR(IC1TA,IC1TB,IC1SM,NOCTPA,NOCTPB,NSM,ICBLTP,IDC,NONEWC,ICOCOC)
    end do
    ! Loop over L and R blocks in core and obtain  contribution from
    ! given L and R blocks
    ISOFF = 1
    IASM = ISSTSM
    IATP = ISSTTA
    IBTP = ISSTTB
    do ISBLK=1,NSBLK
      IBSM = Mul(IASM,ISSM)
      NIA = NSSOA(IASM,IATP)
      NIB = NSSOB(IBSM,IBTP)
      ! Possible permutations of L blocks
      call PRMBLK(IDC,ISTRFL,IASM,IBSM,IATP,IBTP,PSL,PLR,LATP,LBTP,LASM,LBSM,LSGN,LTRP,NLPERM)
      do ILPERM=1,NLPERM
        IIASM = LASM(ILPERM)
        IIBSM = LBSM(ILPERM)
        IIATP = LATP(ILPERM)
        IIBTP = LBTP(ILPERM)
        NIIA = NSSOA(IIASM,IIATP)
        NIIB = NSSOB(IIBSM,IIBTP)

        if (LTRP(ILPERM) == 1) then
          LROW = NSSOA(LASM(ILPERM-1),LATP(ILPERM-1))
          LCOL = NSSOB(LBSM(ILPERM-1),LBTP(ILPERM-1))
          call TRPMT3(SB(ISOFF),LROW,LCOL,C2)
          SB(ISOFF:ISOFF+LROW*LCOL-1) = C2(1:LROW*LCOL)
        end if
        if (LSGN(ILPERM) == -1) SB(ISOFF:ISOFF+NIA*NIB-1) = -SB(ISOFF:ISOFF+NIA*NIB-1)

        JASM = ICSTSM
        JATP = ICSTTA
        JBTP = ICSTTB
        ICOFF = 1
        do ICBLK=1,NCBLK
          JBSM = Mul(JASM,ICSM)
          NJA = NSSOA(JASM,JATP)
          NJB = NSSOB(JBSM,JBTP)
          XNORM2 = dDot_(NJA*NJB,CB(ICOFF),1,CB(ICOFF),1)
          if ((NIA*NIB*NJA*NJB /= 0) .and. (ISOCOC(IATP,IBTP) == 1) .and. (ICOCOC(JATP,JBTP) == 1) .and. (XNORM2 /= Zero)) then
            ! Possible permutations of this block
            call PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PSR,PLR,RATP,RBTP,RASM,RBSM,RSGN,RTRP,NRPERM)
            do IRPERM=1,NRPERM
              if (RTRP(IRPERM) == 1) then
                LROW = NSSOA(RASM(IRPERM-1),RATP(IRPERM-1))
                LCOL = NSSOB(RBSM(IRPERM-1),RBTP(IRPERM-1))
                call TRPMT3(CB(ICOFF),LROW,LCOL,C2)
                CB(ICOFF:ICOFF+LROW*LCOL-1) = C2(1:LROW*LCOL)
              end if
              if (RSGN(IRPERM) == -1) CB(ICOFF:ICOFF+NJA*NJB-1) = -CB(ICOFF:ICOFF+NJA*NJB-1)
              JJASM = RASM(IRPERM)
              JJBSM = RBSM(IRPERM)
              JJATP = RATP(IRPERM)
              JJBTP = RBTP(IRPERM)
              NJJA = NSSOA(JJASM,JJATP)
              NJJB = NSSOB(JJBSM,JJBTP)
              call GSDNBB2_MCLR(I12,RHO1,RHO2,IIASM,IIATP,IIBSM,IIBTP,JJASM,JJATP,JJBSM,JJBTP,NGAS,NELFSPGPA(:,IOCTPA-1+IIATP), &
                                NELFSPGPB(:,IOCTPB-1+IIBTP),NELFSPGPA(:,IOCTPA-1+JJATP),NELFSPGPB(:,IOCTPB-1+JJBTP),NAEL,NBEL, &
                                IAGRP,IBGRP,SB(ISOFF),CB(ICOFF),C2,MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,I3, &
                                XI3S,I4,XI4S,X,NSM,NIIA,NIIB,NJJA,NJJB,NACOB,RHO1S,ieaw,n1,n2)
            end do
            ! Transpose or scale R block to restore order ??
            if (RTRP(NRPERM+1) == 1) then
              call TRPMT3(CB(ICOFF),NJB,NJA,C2)
              CB(ICOFF:ICOFF+NJA*NJB-1) = C2(1:NJA*NJB)
            end if
            if (RSGN(NRPERM+1) == -1) CB(ICOFF:ICOFF+NJA*NJB-1) = -CB(ICOFF:ICOFF+NJA*NJB-1)

          end if
          ICOFF = ICOFF+NCOOSE(JATP,JBTP,JASM)
          ! Next C block
          if (ICBLK /= NCBLK) call NXTBLK_MCLR(JATP,JBTP,JASM,NOCTPA,NOCTPB,NSM,ICBLTP,IDC,NONEWC,ICOCOC)
        end do
        ! End of loop over R blocks in Batch
      end do
      ! Transpose or scale L block to restore order ??
      if (LTRP(NLPERM+1) == 1) then
        call TRPMT3(SB(ISOFF),NIB,NIA,C2)
        SB(ISOFF:ISOFF+NIA*NIB-1) = C2(1:NIA*NIB)
      end if
      if (LSGN(NLPERM+1) == -1) SB(ISOFF:ISOFF+NIA*NIB-1) = -SB(ISOFF:ISOFF+NIA*NIB-1)
      ! Next L block
      ISOFF = ISOFF+NSOOSE(IATP,IBTP,IASM)
      if (ISBLK /= NSBLK) call NXTBLK_MCLR(IATP,IBTP,IASM,NOCTPA,NOCTPB,NSM,ISBLTP,IDC,NONEWS,ISOCOC)
    end do
    ! End of loop over L blocks in batch
    ! End of loop over batches of R blocks
    if (IFINIC /= 0) exit
  end do
  if (ISFINI /= 0) exit outer
end do outer
! End of loop over batches of L blocks

end subroutine GASDN2_MCLR
