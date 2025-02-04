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
! Copyright (C) 1991,1995,1997-1999, Jeppe Olsen                       *
!***********************************************************************

subroutine GASDN2_LUCIA(I12,RHO1,RHO2,RHO2S,RHO2A,L,R,CB,SB,C2,ICOCOC,ISOCOC,ICSMOS,ISSMOS,ICBLTP,ISBLTP,NACOB,NSSOA,NSSOB,NAEL, &
                        IAGRP,NBEL,IBGRP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,NSMST,NSMOB,MXPNGAS,NOBPTS,IOBPTS,MAXK,MAXI,CSCR,SSCR,STSTSX, &
                        SXDXSX,ADSXA,NGAS,NELFSPGP,IDC,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,MXPOBS,IPRNT,RHO1S,LUL,LUR,PSL,PSR, &
                        NBATCHL,LBATL,I1BATL,IBLOCKL,NBATCHR,LBATR,I1BATR,IBLOCKR,ICONSPA,ICONSPB,SCLFAC_L,SCLFAC_R,S2_TERM1, &
                        IPHGAS,IDOSRHO1,SRHO1,IPACK)
! SUBROUTINE GASDN2_LUCIA --> 89
!
! Jeppe Olsen, Winter of 1991
! GAS modificatios, August 1995
!
! Table driven, June 97
!
! Last revision : Jan. 98 (IUSE_PH,IPHGAS added)
!                 Jan. 99 (IDOSRHO1,SRHO1 added)
!
! =====
! Input
! =====
!
! I12    : = 1 => calculate one-electrondensity matrix
!          = 2 => calculate one-and two- electrondensity matrix
! RHO1   : Initial one-electron density matrix
! RHO2   : Initial two-electron density matrix
! RHO2S  : Initial symmetric two-electron density matrix
! RHO2A  : Initial anti-symmetric two-electron density matrix
!
! ICOCOC : Allowed type combinations for C
! ISOCOC : Allowed type combinations for S(igma)
! ICSMOS : Symmetry array for C
! ISSMOS : Symmetry array for S
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
! IPACK : Logical: If true calculate symmetry packed densities (i.e.
!         RHO2S and RHO2A instead of RHO2
!
! RHO1S: Scratch array for one body
! CSCR : Scratch array for C vector
! SSCR : Scratch array for S vector
!
! The L and R vectors are accessed through routines that
! either fetches/disposes symmetry blocks or
! Symmetry-occupation-occupation blocks

use lucia_data, only: IDISK
use Constants, only: Zero, One
use Definitions, only: u6

implicit none
integer I12, NACOB, NAEL, IAGRP, NBEL, IBGRP, IOCTPA, IOCTPB, NOCTPA, NOCTPB, NSMST, NSMOB, MXPNGAS, MAXK, MAXI, NGAS, IDC, &
        MXPOBS, IPRNT, LUL, LUR, NBATCHL, NBATCHR, IDOSRHO1
real*8 PSL, PSR, S2_TERM1
! General input
integer ICOCOC(NOCTPA,NOCTPB), ISOCOC(NOCTPA,NOCTPB)
integer ICSMOS(NSMST), ISSMOS(NSMST)
integer ICBLTP(*), ISBLTP(*)
integer NSSOA(NSMST,NOCTPA)
integer NSSOB(NSMST,NOCTPB)
integer STSTSX(NSMST,NSMST)
integer ADSXA(MXPOBS,2*MXPOBS)
integer SXDXSX(2*MXPOBS,4*MXPOBS)
integer NOBPTS(MXPNGAS,NSMOB), IOBPTS(MXPNGAS,NSMOB)
integer NELFSPGP(MXPNGAS,*)
logical IPACK
integer IPHGAS(*)
real*8 SRHO1(*)
! Info on batches and blocks
integer LBATL(NBATCHL), I1BATL(NBATCHL), IBLOCKL(8,*)
integer LBATR(NBATCHR), I1BATR(NBATCHR), IBLOCKR(8,*)
! Interaction between supergroups
integer ICONSPA(NOCTPA,NOCTPA), ICONSPB(NOCTPB,NOCTPB)
! Scratch
real*8 SB(*), CB(*), C2(*)
real*8 CSCR(*), SSCR(*)
integer I1(*), I2(*), I3(*), I4(*)
real*8 XI1S(*), XI2S(*), XI3S(*), XI4S(*)
real*8 X(*)
real*8 RHO1S(*)
real*8 SCLFAC_L(*), SCLFAC_R(*)
integer ICOOSC(NOCTPA,NOCTPB), ISOOSC(NOCTPA,NOCTPB)
integer LASM(4), LBSM(4), LATP(4), LBTP(4), LSGN(5), LTRP(5)
integer RASM(4), RBSM(4), RATP(4), RBTP(4), RSGN(5), RTRP(5)
real*8 L(*), R(*)
! Output
real*8 RHO1(*), RHO2(*), RHO2S(*), RHO2A(*)
integer ISTRFL(1), LBL(1), IDUMMY(1)
integer INTERACT, NTEST, IBATCHL, NBLKL, IIL, IL, IATP, IBTP, IASM, IBSM, IOFF, ISCALE, IBATCHR, NBLKR, IIR, IR, JATP, JBTP, JASM, &
        JBSM, JOFF, IPERM, NPERM, IIASM, IIBSM, IIATP, IIBTP, IAEXC, IBEXC, IABEXC, NIA, NIB, ILPERM, NLPERM, NIIA, NIIB, LROW, &
        LCOL, NJA, NJB, NRPERM
real*8 PS, PL, PLR, FACTOR, SCLFAC, PCL

! Some dummy initializations
INTERACT = 0 ! jwk-cleanup

NTEST = 0
NTEST = max(NTEST,IPRNT)
if (NTEST >= 20) then
  write(u6,*) ' ================='
  write(u6,*) ' GASDN2 speaking :'
  write(u6,*) ' ================='
  write(u6,*)
  write(u6,*) ' NACOB,MAXK,NGAS,IDC,MXPOBS',NACOB,MAXK,NGAS,IDC,MXPOBS
  write(u6,*) ' LUL, LUR ',LUL,LUR
end if
if (NTEST >= 100) then
  write(u6,*) ' Initial L vector'
  if (LUL == 0) then
    call WRTRS2(L,ISSMOS,ISBLTP,ISOCOC,NOCTPA,NOCTPB,NSSOA,NSSOB,NSMST)
  else
    call WRTVCD(L,LUL,1,-1)
  end if
  write(u6,*) ' Initial R vector'
  if (LUR == 0) then
    call WRTRS2(R,ICSMOS,ICBLTP,ICOCOC,NOCTPA,NOCTPB,NSSOA,NSSOB,NSMST)
  else
    call WRTVCD(R,LUR,1,-1)
  end if
end if
! Loop over batches over L blocks
if (LUL /= 0) IDISK(LUL) = 0
do IBATCHL=1,NBATCHL
  ! Obtain L blocks
  NBLKL = LBATL(IBATCHL)
  if (NTEST >= 200) write(u6,*) ' Left batch, number of blocks',IBATCHL,NBLKL
  do IIL=1,NBLKL
    IL = I1BATL(IBATCHL)-1+IIL
    IATP = IBLOCKL(1,IL)
    IBTP = IBLOCKL(2,IL)
    IASM = IBLOCKL(3,IL)
    IBSM = IBLOCKL(4,IL)
    IOFF = IBLOCKL(5,IL)
    if (NTEST >= 200) write(u6,*) 'IATP IBTP IASM IBSM',IATP,IBTP,IASM,IBSM
    ISCALE = 0
    if (NTEST >= 200) write(u6,*) 'IOFF ',IOFF
    call GSTTBL(L,SB(IOFF),IATP,IASM,IBTP,IBSM,NOCTPA,NOCTPB,NSSOA,NSSOB,PSL,ISOOSC,IDC,PSL,LUL,C2,NSMST,ISCALE,SCLFAC_L(IL))
  end do
  ! Loop over batches  of R vector
  if (LUR /= 0) IDISK(LUR) = 0
  do IBATCHR=1,NBATCHR
    ! Read R blocks into core
    NBLKR = LBATR(IBATCHR)
    if (NTEST >= 200) write(u6,*) ' Right batch, number of blocks',IBATCHR,NBLKR
    do IIR=1,NBLKR
      IR = I1BATR(IBATCHR)-1+IIR
      JATP = IBLOCKR(1,IR)
      JBTP = IBLOCKR(2,IR)
      JASM = IBLOCKR(3,IR)
      JBSM = IBLOCKR(4,IR)
      JOFF = IBLOCKR(5,IR)
      if (NTEST >= 200) write(u6,*) ' JATP JBTP JASM JBSM ',JATP,JBTP,JASM,JBSM
      ! Read R blocks into core

      ! Only blocks interacting with current batch of L are read in
      ! Loop over L  blocks in batch
      do IIL=1,NBLKL
        IL = I1BATL(IBATCHL)-1+IIL
        IATP = IBLOCKL(1,IL)
        IBTP = IBLOCKL(2,IL)
        IASM = IBLOCKL(3,IL)
        IBSM = IBLOCKL(4,IL)
        ! Well, permutations of L blocks
        PS = One
        PL = One
        call PRMBLK(IDC,ISTRFL,IASM,IBSM,IATP,IBTP,PS,PL,LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
        do IPERM=1,NPERM
          IIASM = LASM(IPERM)
          IIBSM = LBSM(IPERM)
          IIATP = LATP(IPERM)
          IIBTP = LBTP(IPERM)

          IAEXC = ICONSPA(IIATP,JATP)
          IBEXC = ICONSPB(IIBTP,JBTP)
          if ((IAEXC == 0) .and. (IIASM /= JASM)) IAEXC = 1
          if ((IBEXC == 0) .and. (IIBSM /= JBSM)) IBEXC = 1
          IABEXC = IAEXC+IBEXC
          if (IABEXC <= I12) INTERACT = 1
        end do
      end do
      ! End of checking whether C-block is needed
      ISCALE = 0
      if (INTERACT == 1) then
        ISCALE = 0
        call GSTTBL(R,CB(JOFF),JATP,JASM,JBTP,JBSM,NOCTPA,NOCTPB,NSSOA,NSSOB,PSR,ICOOSC,IDC,PCL,LUR,C2,NSMST,ISCALE,SCLFAC_R(IR))
      else
        !write(u6,*) ' TTSS for C block skipped'
        !call IWRTMA(IBLOCKR(1,IR),4,1,4,1)
        call IDAFILE(LUR,2,LBL,1,IDISK(LUR))
        call IDAFILE(LUR,2,IDUMMY,1,IDISK(LUR))
        call SKPRCD2(LBL(1),-1,LUR)
        SCLFAC_R(IR) = Zero
      end if

      if (NTEST >= 100) then
        if (INTERACT == 1) then
          write(u6,*) ' TTSS for C block read in'
          call IWRTMA(IBLOCKR(1,IR),4,1,4,1)
        else
          write(u6,*) ' TTSS for C block skipped'
          call IWRTMA(IBLOCKR(1,IR),4,1,4,1)
        end if
      end if
    end do

    ! Loop over L and R blocks in batches and obtain  contribution from
    ! given L and R blocks
    PLR = One
    do IIL=1,NBLKL
      IL = I1BATL(IBATCHL)-1+IIL
      if (SCLFAC_L(IL) /= Zero) then
        IATP = IBLOCKL(1,IL)
        IBTP = IBLOCKL(2,IL)
        IASM = IBLOCKL(3,IL)
        IBSM = IBLOCKL(4,IL)
        IOFF = IBLOCKL(5,IL)

        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
        ! Possible permutations of L blocks
        call PRMBLK(IDC,ISTRFL,IASM,IBSM,IATP,IBTP,PSL,PLR,LATP,LBTP,LASM,LBSM,LSGN,LTRP,NLPERM)
        do ILPERM=1,NLPERM
          !write(u6,*) ' Loop 9999 ILPERM = ',ILPERM
          IIASM = LASM(ILPERM)
          IIBSM = LBSM(ILPERM)
          IIATP = LATP(ILPERM)
          IIBTP = LBTP(ILPERM)
          NIIA = NSSOA(IIASM,IIATP)
          NIIB = NSSOB(IIBSM,IIBTP)

          if (LTRP(ILPERM) == 1) then
            LROW = NSSOA(LASM(ILPERM-1),LATP(ILPERM-1))
            LCOL = NSSOB(LBSM(ILPERM-1),LBTP(ILPERM-1))
            call TRPMT3(SB(IOFF),LROW,LCOL,C2)
            call COPVEC(C2,SB(IOFF),LROW*LCOL)
          end if
          if (LSGN(ILPERM) == -1) call SCALVE(SB(IOFF),-One,NIA*NIB)

          do IIR=1,NBLKR
            IR = I1BATR(IBATCHR)-1+IIR
            if (SCLFAC_R(IR) /= Zero) then
              JATP = IBLOCKR(1,IR)
              JBTP = IBLOCKR(2,IR)
              JASM = IBLOCKR(3,IR)
              JBSM = IBLOCKR(4,IR)
              JOFF = IBLOCKR(5,IR)

              NJA = NSSOA(JASM,JATP)
              NJB = NSSOB(JBSM,JBTP)

              IAEXC = ICONSPA(JATP,IIATP)
              IBEXC = ICONSPB(JBTP,IIBTP)

              if ((IAEXC == 0) .and. (JASM /= IIASM)) IAEXC = 1
              if ((IBEXC == 0) .and. (JBSM /= IIBSM)) IBEXC = 1
              IABEXC = IAEXC+IBEXC

              if (IABEXC <= I12) then
                INTERACT = 1
              else
                INTERACT = 0
              end if

              if (INTERACT == 1) then
                ! Possible permutations of this block
                call PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PSR,PLR,RATP,RBTP,RASM,RBSM,RSGN,RTRP,NRPERM)
                ! Well, spin permutations are simple to handle
                ! if there are two terms just calculate and and multiply with
                ! 1+PSL*PSR
                if (NRPERM == 1) then
                  FACTOR = One
                else
                  FACTOR = One+PSL*PSR
                end if
                SCLFAC = FACTOR*SCLFAC_L(IL)*SCLFAC_R(IR)
                if ((INTERACT == 1) .and. (SCLFAC /= Zero)) then
                  if (NTEST >= 20) then
                    write(u6,*) ' RSDNBB will be called for'
                    write(u6,*) ' L block :'
                    write(u6,'(A,5I5)') ' IIASM IIBSM IIATP IIBTP',IIASM,IIBSM,IIATP,IIBTP
                    write(u6,*) ' R  block :'
                    write(u6,'(A,5I5)') ' JASM JBSM JATP JBTP',JASM,JBSM,JATP,JBTP
                    write(u6,*) ' IOFF,JOFF ',IOFF,JOFF
                    write(u6,*) ' SCLFAC = ',SCLFAC
                  end if

                  call GSDNBB2_LUCIA(I12,RHO1,RHO2,RHO2S,RHO2A,IIASM,IIATP,IIBSM,IIBTP,JASM,JATP,JBSM,JBTP,NGAS, &
                                     NELFSPGP(1,IOCTPA-1+IIATP),NELFSPGP(1,IOCTPB-1+IIBTP),NELFSPGP(1,IOCTPA-1+JATP), &
                                     NELFSPGP(1,IOCTPB-1+JBTP),NAEL,NBEL,IAGRP,IBGRP,SB(IOFF),CB(JOFF),C2,ADSXA,STSTSX,SXDXSX, &
                                     MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,NSMOB,NSMST,NIIA, &
                                     NIIB,NJA,NJB,MXPOBS,IPRNT,NACOB,RHO1S,SCLFAC,S2_TERM1,IPHGAS,IDOSRHO1,SRHO1,IPACK)

                  ! CALL GSDNBB2_LUCIA --> 66

                  if (NTEST >= 500) then
                    write(u6,*) ' Updated rho1'
                    call wrtmat(rho1,nacob,nacob,nacob,nacob)
                    write(u6,*) ' Updated srho1'
                    call wrtmat(srho1,nacob,nacob,nacob,nacob)
                  end if

                end if
              end if
            end if
          end do
          ! End of loop over R blocks in Batch
        end do
        ! Transpose or scale L block to restore order ??
        if (LTRP(NLPERM+1) == 1) then
          call TRPMT3(SB(IOFF),NIB,NIA,C2)
          call COPVEC(C2,SB(IOFF),NIA*NIB)
        end if
        if (LSGN(NLPERM+1) == -1) call SCALVE(SB(IOFF),-One,NIA*NIB)

      end if
    end do
    ! End of loop over L blocks in batch
  end do
  ! End of loop over batches of R blocks
end do
! End of loop over batches of L blocks

end subroutine GASDN2_LUCIA
