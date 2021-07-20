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
! Copyright (C) 1989, Jeppe Olsen                                      *
!               1989, Markus P. Fuelscher                              *
!               Giovanni Li Manni                                      *
!***********************************************************************

subroutine DiagOrd(PHPCSF,PHPCNF,IPORDCSF,IPORDCNF,MXPDIM,condition,iter,DTOC,IPRODT,ICONF,IREFSM,ONEBOD,ECORE,NACTOB,SCR,NCONF, &
                   NEL,NAEL,NBEL,TUVX,NTEST,ExFac,IREOTS)
! Obtain primary subspace and obtain
! explicit representation of Hamiltonian matrix in subspace
!
! ARGUMENTS :
! ===========
! PHPCSF  : Diagonal Hamiltonian matrix elements in CSF (output)
! PHPCNF  : Diagonal Hamiltonian matrix elements in CNF (output)
! IPORDCSF: index array containing energetic order in CSF (output)
! IPORDCNF: index array containing energetic order in CNF (output)
! IPCSF   : CSF's defining subspace (input)
! IPCNF   : Configurations defining subspace (input)
! MXPDIM  : Largest allowed dimension of subspace (Input)
! DTOC    : Transformation matrix between CSF's and DET's (input)
! IPRODT  : Prototype determinants (input)
! ICONF   : List of configurations  (input)
! IREFSM  : symmetry of considered CI space (input)
! Onebod  : one body hamilton matrix in rectangular form (input)
! ECORE   : Core energy (input)
! NACTOB  : Number of active orbitals (input)
! SCR     : Scratch array of length ????
! NCONF   : Number of configurations of symmetry IREFSM
! NPCNF   : = NCONF
! NPCSF   : = MXPDIM
! TUVX    : Two-electron integrals (MO space)
!
! IREOTS : Type => symmetry reordering array
!
!*********** Author: GLMJ ****************
!   history:
!     Jeppe Olsen, Summer of '89
!     adapted to DETRAS by M.P. Fuelscher, October 1989

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Constants, only: One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: MXPDIM, IPORDCSF(MXPDIM), iter, IPRODT(*), ICONF(*), IREFSM, NACTOB, NCONF, NEL, NAEL, NBEL, &
                                 IREOTS(NACTOB)
real(kind=wp), intent(out) :: PHPCSF(MXPDIM), PHPCNF(NCONF)
integer(kind=iwp), intent(out) :: IPORDCNF(NCONF)
real(kind=wp), intent(in) :: condition, DTOC(*), ONEBOD(NACTOB,NACTOB), ECORE, TUVX(*), ExFac
real(kind=wp), intent(_OUT_) :: SCR(*)
integer(kind=iwp), intent(inout) :: NTEST
integer(kind=iwp) :: ICSFMN, IICNF, IICSF, IILACT, IILB, ILRI, ILTYP, IMIN, iTmpDimBlockA, iTmpDimBlockACNF, KLCNFO, KLCONF, &
                     KLCSFO, KLFREE, KLPHPS, KLSCRS, MXCSFC, NCSFL, NCSFMN, NIRREP, NJCNF, NPCNF, NPCSF
real(kind=wp) :: Acc, RefSplit, XMAX, XMIN
real(kind=wp), external :: FNDMNX
#include "spinfo.fh"
#include "splitcas.fh"

call DIAGORD_INTERNAL(SCR)

! This is to allow type punning without an explicit interface
contains

subroutine DIAGORD_INTERNAL(SCR)

  real(kind=wp), target :: SCR(*)
  integer(kind=iwp), pointer :: iSCR(:)
  integer(kind=iwp) :: i, ICNF, ICNL, IIL, ITYP

  ICSFMN = 0 ! dummy initialize to avoid warning

  ! construct the diagonal array out of the Hamiltonian matrix

  MXCSFC = 0
  do ITYP=1,NTYP
    MXCSFC = max(MXCSFC,NCSFTP(ITYP))
  end do

  KLFREE = 1
  KLCONF = KLFREE
  KLFREE = KLFREE+NEL
  KLFREE = KLFREE+NEL
  KLPHPS = KLFREE
  KLFREE = KLFREE+MXCSFC*MXCSFC
  KLSCRS = KLFREE
  KLFREE = KLFREE+MXCSFC
  KLCSFO = KLFREE
  KLFREE = KLFREE+MXPDIM
  KLCNFO = KLFREE
  KLFREE = KLFREE+NCONF

  IILB = 1
  do ICNL=1,NCONF
    !write(u6,*) 'IILB',IILB
    call c_f_pointer(c_loc(SCR(KLCONF)),iSCR,[1])
    call GETCNF_LUCIA(iSCR,ILTYP,ICNL,ICONF,IREFSM,NEL)
    nullify(iSCR)
    !call GETCNF_LUCIA(iSCR,ILTYP,IPCNF(ICNL),ICONF,IREFSM,NEL)
    NCSFL = NCSFTP(ILTYP)
    !write(u6,*) 'NCSFL = ',NCSFL
    !call xflush(u6)
    call c_f_pointer(c_loc(SCR(KLCONF)),iSCR,[1])
    call CNHCN(iSCR,ILTYP,iSCR,ILTYP,SCR(KLPHPS),SCR(KLFREE),NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,NTEST,ExFac,IREOTS)
    nullify(iSCR)
    do IIL=1,NCSFL
      IILACT = IILB-1+IIL
      ILRI = IIL*IIL
      PHPCSF(IILACT) = SCR(KLPHPS-1+ILRI)
      !write(u6,*) 'IILACT, ILRI = ',IILACT,ILRI
      !write(u6,*) 'PHPCSF(IILACT)',PHPCSF(IILACT)
      !call xflush(u6)
      SCR(KLSCRS-1+IIL) = SCR(KLPHPS-1+ILRI)
    end do
    XMAX = -FNDMNX(SCR(KLSCRS),NCSFL,2)
    PHPCNF(ICNL) = XMAX
    !write(u6,*) 'PHPCNF(IILACT)',PHPCNF(ICNL)
    !call xflush(u6)
    IILB = IILB+NCSFL
  end do
  !*********************** ORDERING ***************************
  Acc = 1.0e-13_wp
  XMAX = FNDMNX(PHPCNF,NCONF,2)
  RefSplit = -XMAX
  NPCSF = 0
  NPCNF = 0
  ! To avoid warnings:
  iTmpDimBlockACNF = 0
  iTmpDimBlockA = 0

  !IFINIT = 0
  do
    XMIN = XMAX+One
    IMIN = 0
    IICNF = 1
    IICSF = 1
    do ITYP=1,NTYP
      NJCNF = NCNFTP(ITYP,IREFSM)
      NIRREP = NCSFTP(ITYP)
      do ICNF=1,NJCNF
        if (PHPCNF(IICNF)+Acc < XMIN) then
          XMIN = PHPCNF(IICNF)
          IMIN = IICNF
          !write(u6,*) 'IMIN         :',IMIN
          ICSFMN = IICSF
          NCSFMN = NIRREP
        end if
        IICNF = IICNF+1
        IICSF = IICSF+NIRREP
      end do
    end do
    !write(u6,*) 'XMIN         :',XMIN
    !write(u6,*) 'NPCSF+NCSFMN :',NPCSF+NCSFMN
    !if ((NPCSF+NCSFMN) <= MXPDIM) then
    SCR(KLCNFO+NPCNF) = PHPCNF(IMIN)
    NPCNF = NPCNF+1
    IPORDCNF(NPCNF) = IMIN
    !write(u6,*)'NPCNF, IPORDCNF(NPCNF) :', NPCNF,IPORDCNF(NPCNF)
    call ISTVC2(IPORDCSF(NPCSF+1),ICSFMN-1,1,NCSFMN)
    !write(u6,*) 'IPORDCSF(NPCSF) :',(IPORDCSF(NPCSF+i),i=1,NCSFMN)
    !call IVCPRT('IPORDCSF(NPCSF) :',' ', IPORDCSF(NPCSF+1),NCSFMN)
    NPCSF = NPCSF+NCSFMN

    if (iter == 1) then
      if (EnerSplit) then
        if ((PHPCNF(IMIN)-RefSplit) <= condition) then
          iDimBlockACNF = NPCNF
          iDimBlockA = NPCSF
        end if
      else if (PerSplit) then
        if ((real(NPCSF)/real(MXPDIM))*1.0e2_wp <= condition+Acc) then
          iDimBlockACNF = NPCNF
          iDimBlockA = NPCSF
        end if
      end if
    else
      if (NPCNF <= iDimBlockACNF) then
        iTmpDimBlockACNF = NPCNF
        iTmpDimBlockA = NPCSF
      end if
    end if

    PHPCNF(IMIN) = XMAX+One
    !else
    !  ! No space for this configuration, remove previous
    !  ! configurations with the same diagonal value
    !  IFINIT = 1
    !  IICNF = NPCNF+1
    !  do
    !    IICNF = IICNF-1
    !    DIAVAL = PHPCNF(IPORDCNF(IICNF))
    !    if (abs(DIAVAL-XMIN) > 1.0e-10_wp) exit
    !    NPCNF = NPCNF-1
    !    call GETCNF_LUCIA(SCR(KLFREE),ITYP,IPCNF(IICNF),ICONF,IREFSM,NEL)
    !    call GETCNF_LUCIA(PHPCNF(NCONF+1),ITYP,IPCNF(IICNF),ICONF,IREFSM,NEL)
    !    NPCSF = NPCSF-NCSFTP(ITYP)
    !  end do
    !end if
    !if ((IFINIT /= 0) .or. (NPCNF >= NCONF)) exit
    if (NPCNF >= NCONF) exit
  end do

  if (iter /= 1) then
    iDimBlockACNF = iTmpDimBlockACNF
    iDimBlockA = iTmpDimBlockA
    !write(u6,*) 'iDimBlockACNF, iDimBlockACSF',iDimBlockACNF,iDimBlockA
  end if

  do i=0,MXPDIM-1
    SCR(KLCSFO+i) = PHPCSF(IPORDCSF(i+1))
  end do

  call dcopy_(NCONF,SCR(KLCNFO),1,PHPCNF,1)
  call dcopy_(MXPDIM,SCR(KLCSFO),1,PHPCSF,1)

  !write(u6,*) ' OUTPUT from DiagOrd'
  !write(u6,*) ' =================='
  !write(u6,*) ' Number of Configurations in primary subspace ',NCONF
  !write(u6,*) ' Number of CSFs in primary subspace ',MXPDIM
  !write(u6,*) ' Configurations included :'
  !call IWRTMA(IPORDCNF,1,NCONF,1,NCONF)
  !write(u6,*) ' CSFs included :'
  !CALL IWRTMA(IPORDCSF,1,MXPDIM,1,MXPDIM)
  !write(u6,*) 'Ordered Diagonal array in CSF basis :'
  !do i=1,MXPDIM
  !  write(u6,*) PHPCSF(IPORDCSF(i))
  !  write(u6,*) PHPCSF(i)
  !end do
  !write(u6,*) 'Ordered Diagonal array in CNF basis :'
  !do i=1,NCONF
  !  write(u6,*) PHPCNF(i)
  !end do

  return

end subroutine DIAGORD_INTERNAL

end subroutine DiagOrd
