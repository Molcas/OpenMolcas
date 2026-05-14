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
! Copyright (C) 1993, Jeppe Olsen                                      *
!***********************************************************************

subroutine H0CSF(H0,IPQCSF,IPQCNF,MXP1DM,MXP2DM,MXQDM,DTOC,IPRODT,ICONF,IREFSM,NACTOB,SCR,NCONF,NEL,NAEL,NBEL,IPWAY,NP1CSF,NP1CNF, &
                 NP2CSF,NP2CNF,NQCSF,NQCNF,NPQCSF,NPQCNF,DIAG,DIAGCN,INTSPC,ICOMBI,PSSIGN)
! Obtain H0 subspace defined by the three parameters
! MXP1DM,MXP2DM,MXQDM and obtain
! explicit representation of hamilton matrix in subspace
!
! The H0 space consist of three subspaces
! P1,P2 AND Q
! The H0 matrix can be pictured as
!
!              P1    P2        Q
!             ***************************
!             *    *     *              *
!         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
!             ***************************         is used in this block
!         P2  *    *     *              *
!             * Ex *  Ex *     Diag     *    Diag : Diagonal
!             ************              *           appriximation used
!             *    *      *             *
!             *    *        *           *
!             * Ex *  Diag    *         *
!         Q   *    *            *       *
!             *    *              *     *
!             *    *                *   *
!             *    *                  * *
!             ***************************
!
! The exact Hamiltonian is therefore calculated in subspace P1 and P2
! but only the interaction between a larger space Q and subspace P1
! is calculated exactly
!
! Lucia Version, September 1993, Jeppe Olsen
!
! ===========
! ARGUMENTS :
! ===========
!   H0 : combined block to contain :
!      PHP : hamilton matrix in subspace P1+P2 output)
!      PHQ : Hamiltonian matrix block with row in P1 and column in Q
!      QHQ : Diagonal approximation to matrix in Q-Q space
!   IPQCSF : CSFs defining subspace (output)
!   IPQCNF : Configurations defining subspace (output)
!   MXP1DM : Largest allowed dimension of subspace P1 (Input)
!   MXP2DM : Largest allowed dimension of subspace P2 (Input)
!   MXQDM  : Largest allowed dimension of subspace Q  (Input)
!   DTOC : Transformation matrix between CSFs and DETs (input)
!   IPRODT : Prototype determinants (input)
!   ICONF : List of configurations  (input)
!   IREFSM : symmetry of considered CI space (input)
!   NACTOB : Number of active orbitals (input)
!   SCR    : Scratch array of length ????
!   NCONF : Number of configurations of symmetry IREFSM
!   IPWAY : Defines way of choosing Primary space
!    = 1  : use the first configurations (until at most
!           MXPDIM CSFs have been included)
!   NP1CNF : Number of primary configurations obtained in P1 (output)
!   NP1CSF : Number of primary CSFs obtained in P1 (OUTPUT)
!   NP2CNF : Number of primary configurations obtained in P2 (output)
!   NP2CSF : Number of primary CSFs obtained in P2 (OUTPUT)
!   NQCNF  : Number of primary configurations obtained in Q (output)
!   NQCSF  : Number of primary CSFs obtained in Q (OUTPUT)
!
!   DIAG   : Hamilton diagonal over CSFs (INPUT)
!   DIAGCN : space for diagonal over configurations
!   INTSPC : Internal space number of actual expansion
!   ICOMBI : = 0 => no spin combinations
!            = 1 =>    spin combinations
!   PSSIGN : Spin combination sign
!
! =========================================
! Jeppe Olsen, Spring of 90, from PHPCSF
! =========================================
!  Lucia Version, September 1993
! =========================================

use Index_Functions, only: nTri_Elem
use MCLR_Data, only: NCNATS, NCPCNT, NTYP
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: H0(*), SCR(*), DIAGCN(*)
integer(kind=iwp), intent(_OUT_) :: IPQCSF(*), IPQCNF(*)
integer(kind=iwp), intent(in) :: MXP1DM, MXP2DM, MXQDM, IPRODT(*), ICONF(*), IREFSM, NACTOB, NCONF, NEL, NAEL, NBEL, IPWAY, &
                                 INTSPC, ICOMBI
real(kind=wp), intent(in) :: DTOC(*), DIAG(*), PSSIGN
integer(kind=iwp), intent(out) :: NP1CSF, NP1CNF, NP2CSF, NP2CNF, NQCSF, NQCNF, NPQCSF, NPQCNF
integer(kind=iwp) :: i, ICNF, ICSFMN, ICSFOF, IDEG, IDGCNF, IDGCSF, IDGVL, IFINIT, IICNF, IICSF, IMIN, ITYP, KLCONF, KLPHP, KLPHQ, &
                     MXPQDM, NCSFMN, NDGVL, NIRREP, NJCNF, NPCNF, NPCSF
real(kind=wp) :: DIAVAL, XMAX, XMIN
integer(kind=iwp), allocatable :: ISCR(:)
real(kind=wp), external :: FNDMNX

! 1 : Obtain primary subspace

NP1CSF = 0
NP1CNF = 0
NP2CSF = 0
NP2CNF = 0
NPCSF = 0
NQCSF = 0
NQCNF = 0
NPQCSF = 0
NPQCNF = 0
ICSFMN = 0 ! dummy initialize
NCSFMN = 0 ! dummy initialize
MXPQDM = MXP1DM+MXP2DM+MXQDM

! ======================================================================
! 1 :                 Generate initial subspace
! ======================================================================

if (IPWAY == 1) then

  ! Just use the first CSFs as subspace

  ICNF = 0
  outer: do ITYP=1,NTYP
    NJCNF = NCNATS(ITYP,IREFSM)
    NIRREP = NCPCNT(ITYP)
    do IICNF=1,NJCNF
      ICNF = ICNF+1
      if (NP1CSF+NIRREP <= MXP1DM) then
        NPQCNF = NPQCNF+1
        NP1CNF = NP1CNF+1
        IPQCNF(NPQCNF) = ICNF
        IPQCSF(NPQCSF+1:NPQCSF+NIRREP) = [(IICSF,IICSF=NPQCSF+1,NPQCSF+NIRREP)]
        NPQCSF = NPQCSF+NIRREP
        NP1CSF = NP1CSF+NIRREP
      else if (NP2CSF+NIRREP <= MXP2DM) then
        NPQCNF = NPQCNF+1
        NP2CNF = NP2CNF+1
        IPQCNF(NPQCNF) = ICNF
        IPQCSF(NPQCSF+1:NPQCSF+NIRREP) = [(IICSF,IICSF=NPQCSF+1,NPQCSF+NIRREP)]
        NPQCSF = NPQCSF+NIRREP
        NP2CSF = NP2CSF+NIRREP
      else if (NQCSF+NIRREP <= MXQDM) then
        NPQCNF = NPQCNF+1
        NQCNF = NQCNF+1
        IPQCNF(NPQCNF) = ICNF
        IPQCSF(NPQCSF+1:NPQCSF+NIRREP) = [(IICSF,IICSF=NPQCSF+1,NPQCSF+NIRREP)]
        NPQCSF = NPQCSF+NIRREP
        NQCSF = NQCSF+NIRREP
      else
        exit outer
      end if
    end do
  end do outer

else if (IPWAY == 2) then
  ! Obtain lowest CSFs

  ! 1 :  local Diagonal elements over configurations

  IICSF = 1
  IICNF = 1

  call mma_allocate(ISCR,MXPQDM+NEL,Label='ISCR')
  KLCONF = MXPQDM+1

  do ITYP=1,NTYP
    NJCNF = NCNATS(ITYP,IREFSM)
    NIRREP = NCPCNT(ITYP)
    do ICNF=1,NJCNF
      DIAGCN(IICNF) = DIAG(IICSF)
      IICNF = IICNF+1
      IICSF = IICSF+NIRREP
    end do
  end do
  ! Largest element
  XMAX = FNDMNX(DIAGCN,NCONF,2)
  ! loop over lowest configurations

  IFINIT = 0

  do
    XMIN = XMAX+One
    IMIN = 0

    IICNF = 1
    ICSFOF = 1
    do ITYP=1,NTYP
      NIRREP = NCPCNT(ITYP)
      do ICNF=1,NCNATS(ITYP,IREFSM)
        if (DIAGCN(IICNF) < XMIN) then
          XMIN = DIAGCN(IICNF)
          IMIN = IICNF
          ICSFMN = ICSFOF
          NCSFMN = NIRREP
        end if

        IICNF = IICNF+1
        ICSFOF = ICSFOF+NIRREP
      end do
    end do

    ! Next lowest element has been found

    if (NPQCSF+NCSFMN <= MXPQDM) then
      ! 1  add new configuration
      NPQCNF = NPQCNF+1
      IPQCNF(NPQCNF) = IMIN
      SCR(NPQCNF) = XMIN
      IPQCSF(NPQCSF+1:NPQCSF+NCSFMN) = [(i,i=ICSFMN,ICSFMN+NCSFMN-1)]
      NPQCSF = NPQCSF+NCSFMN

      ! Mask
      DIAGCN(IMIN) = XMAX+One
    else
      IFINIT = 1
      ! 2  No space for this configuration, remove previous
      !    configurations with the same diagonal value
      IICNF = NPQCNF+1
      do
        IICNF = IICNF-1
        DIAVAL = DIAGCN(IPQCNF(IICNF))
        if (abs(DIAVAL-XMIN) > 1.0e-10_wp) exit
        NPQCNF = NPQCNF-1
        call GETCNF_MCLR(ISCR(KLCONF),ITYP,IPQCNF(IICNF),ICONF,IREFSM,NEL)
        NPQCSF = NPQCSF-NCPCNT(ITYP)
      end do
    end if
    if ((IFINIT /= 0) .or. (NPQCNF >= NCONF)) exit
  end do
  ! NPQCSF has now been collected, obtain P1,P2 and Q space
  ! so that degenerate configurations are not  in
  ! different  subspaces

  ! Arrange selected configurations in degenerate pairs

  call DEGVEC(SCR,NPQCNF,NDGVL,ISCR)
  ! Number of configurations in P1,P2,Q
  ICNF = 0
  do IDEG=1,NDGVL
    ! Number of CSFs in this group of degenerate values
    IDGVL = ISCR(IDEG)
    IDGCSF = 0
    do IDGCNF=1,IDGVL
      ICNF = ICNF+1
      call GETCNF_MCLR(ISCR(KLCONF),ITYP,IPQCNF(ICNF),ICONF,IREFSM,NEL)
      IDGCSF = IDGCSF+NCPCNT(ITYP)
    end do
    if ((NP1CSF+IDGCSF <= MXP1DM) .and. (NP2CSF+NQCSF == 0)) then
      ! Add to P1
      NP1CSF = NP1CSF+IDGCSF
      NP1CNF = NP1CNF+IDGVL
    else if ((NP2CSF+IDGCSF <= MXP2DM) .and. (NQCSF == 0)) then
      ! Add to P2
      NP2CSF = NP2CSF+IDGCSF
      NP2CNF = NP2CNF+IDGVL
    else if (NQCSF+IDGCSF <= MXQDM) then
      ! Add to Q
      NQCSF = NQCSF+IDGCSF
      NQCNF = NQCNF+IDGVL
    else
      ! No space for configuration so
      exit
    end if
  end do

  NPCSF = NP1CSF+NP2CSF
  NPCNF = NP1CNF+NP2CNF

  NPQCSF = NPCSF+NQCSF
  NPQCNF = NPCNF+NQCNF

  call mma_deallocate(ISCR)

end if
! End if for IPWAY = 2

! This is not beautiful, but necessary
!MXP1DM = NP1CSF
!MXP2DM = NP2CSF
!MXQDM = NQCSF

! ============================================================
! 2          Construct Hamiltonian matrix in subspace
! ============================================================

! Do not add core energy to subspace Hamiltonian, add to eigenvalues
! Pointers in H0
KLPHP = 1
KLPHQ = KLPHP+nTri_Elem(NPCSF)

! PHP matrix

call CNHCNM(H0(KLPHP),1,IPQCNF,NPCNF,IPQCNF,NPCNF,NPCSF,DIAGCN,ICONF,NEL,IREFSM,NAEL,NBEL,NACTOB,IPRODT,DTOC,INTSPC,ICOMBI,PSSIGN)

! PHQ matrix

call CNHCNM(H0(KLPHQ),0,IPQCNF,NP1CNF,IPQCNF(1+NPCNF),NQCNF,NP1CSF,DIAGCN,ICONF,NEL,IREFSM,NAEL,NBEL,NACTOB,IPRODT,DTOC,INTSPC, &
            ICOMBI,PSSIGN)

end subroutine H0CSF
