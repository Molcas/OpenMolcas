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
!***********************************************************************

subroutine PHPCSF(PHP,IPCSF,IPCNF,MXPDIM,DTOC,IPRODT,ICONF,IREFSM,ONEBOD,ECORE,NACTOB,SCR,NCONF,NEL,NAEL,NBEL,NPCSF,NPCNF,DIAG, &
                  TUVX,NTEST,ExFac,IREOTS)
! Obtain primary subspace and obtain
! explicit representation of hamilton matrix in subspace
!
! ARGUMENTS :
! ===========
! PHP    : hamilton matrix in subspace (output)
! IPCSF  : CSF's defining subspace (output)
! IPCNF  : Configurations defining subspace (output)
! MXPDIM : Largest allowed dimension of subspace (Input)
! DTOC   : Transformation matrix between CSF's and DET's (input)
! IPRODT : Prototype determinants (input)
! ICONF  : List of configurations  (input)
! IREFSM : symmetry of considered CI space (input)
! Onebod : one body hamilton matrix in rectangular form (input)
! ECORE  : Core energy (input)
! NACTOB : Number of active orbitals (input)
! SCR    : Scratch array of length ????
! NCONF  : Number of configurations of symmetry IREFSM
! NPCNF  : Number of primary configurations obtained (output)
! NPCSF  : Number of primary CSF's obtained  (output)
! TUVX   : Two-electron integrals (MO space)
! DIAG   : Hamilton diagonal over CSF's
!
! IREOTS : Type => symmetry reordering array
!
! Jeppe Olsen , Summer of '89
! adapted to DETRAS by M.P. Fuelscher, Oktober 1989

implicit real*8(A-H,O-Z)
#include "spinfo.fh"
#include "warnings.h"
dimension PHP(*), IPCSF(*), IPCNF(*), DIAG(*)
dimension DTOC(*), IPRODT(*), ICONF(*)
dimension ONEBOD(*), SCR(*)
dimension TUVX(*), IREOTS(*)

call PHPCSF_INTERNAL(SCR)

! This is to allow type punning without an explicit interface
contains

subroutine PHPCSF_INTERNAL(SCR)

  use iso_c_binding

  real*8, target :: SCR(*)
  integer, pointer :: iSCRl(:), iSCRr(:)

  ! Assumed machine accuray (give and take)

  Acc = 1.0D-13

  ! construct the diagonal of the Hamilton matrix in CNF basis
  ICSFMN = 0
  IICNF = 1
  IICSF = 1
  do ITYP=1,NTYP
    NJCNF = NCNFTP(ITYP,IREFSM)
    NIRREP = NCSFTP(ITYP)
    do ICNF=1,NJCNF
      SCR(IICNF) = DIAG(IICSF)
      IICNF = IICNF+1
      IICSF = IICSF+NIRREP
    end do
  end do
  if (NTEST >= 30) call RecPrt('SCR',' ',SCR,1,IICNF-1)

  ! find the elements of the subspace
  ! (MXPDIM elements smallest in energy)

  XMAX = FNDMNX(SCR,NCONF,2)
  NPCSF = 0
  NPCNF = 0
  IFINIT = 0
400 continue
  XMIN = XMAX+1.0d0
  IMIN = 0
  IICNF = 1
  IICSF = 1
  do ITYP=1,NTYP
    NJCNF = NCNFTP(ITYP,IREFSM)
    NIRREP = NCSFTP(ITYP)
    do ICNF=1,NJCNF
      if (SCR(IICNF)+Acc < XMIN) then
        XMIN = SCR(IICNF)
        IMIN = IICNF
        ICSFMN = IICSF
        NCSFMN = NIRREP
      end if
      IICNF = IICNF+1
      IICSF = IICSF+NIRREP
    end do
  end do
  if ((NPCSF+NCSFMN) <= MXPDIM) then
    ! add new configuration
    NPCNF = NPCNF+1
    IPCNF(NPCNF) = IMIN
    call ISTVC2(IPCSF(NPCSF+1),ICSFMN-1,1,NCSFMN)
    NPCSF = NPCSF+NCSFMN
    SCR(IMIN) = XMAX+1.0d0
  else
    ! No space for this configuration , remove previous
    ! configurations with the same diagonal value
    ! PAM 2011: Do not remove any previous configurations.
    IFINIT = 1
    !IICNF = NPCNF+1
    !600 continue
    !IICNF = IICNF-1
    !DIAVAL = SCR(IPCNF(IICNF))
    !if (abs(DIAVAL-XMIN) <= 1.0D-10) then
    !  NPCNF = NPCNF-1
    !  call GETCNF_LUCIA(SCR(NCONF+1),ITYP,IPCNF(IICNF),ICONF,IREFSM,NEL)
    !  NPCSF = NPCSF-NCSFTP(ITYP)
    !  goto 600
    !end if
  end if
  if ((IFINIT == 0) .and. (NPCNF < NCONF)) goto 400

  if (NPCNF == 0) then
    call WarningMessage(2,'Making explicit Hamiltonian failed.')
    write(6,*) ' An unforeseen catastrophic failure occurred'
    write(6,*) ' in the CI solver. The size of the explicit'
    write(6,*) ' part of the CI Hamiltonian matrix was not'
    write(6,*) ' sufficient. Suggested fix: Change the size'
    write(6,*) ' by adding ''SDAV=XXXXX'' to the rasscf input.'
    write(6,*) ' XXXXX is some integer at least ',MXPDIM+NCSFMN
    write(6,*) ' Sorry about this. Consider telling the'
    write(6,*) ' Molcas group about this failure.'
    call Quit(_RC_INTERNAL_ERROR_)
  end if

  if (NTEST >= 30) then
    write(6,*) ' Output from PHPCSF '
    write(6,*) ' ================== '
    write(6,*) ' Number of Configurations in primary subspace ',NPCNF
    write(6,*) ' Number of CSFs in primary subspace ',NPCSF
    write(6,*) ' Configurations included : '
    call IWRTMA(IPCNF,1,NPCNF,1,NPCNF)
    write(6,*) ' CSFs included : '
    call IWRTMA(IPCSF,1,NPCSF,1,NPCSF)
  end if

  ! construct Hamiltonian matrix in subspace

  MXCSFC = 0
  do ITYP=1,NTYP
    MXCSFC = max(MXCSFC,NCSFTP(ITYP))
  end do

  KLFREE = 1
  KLCONF = KLFREE
  KLFREE = KLFREE+NEL
  KRCONF = KLFREE
  KLFREE = KLFREE+NEL
  KLPHPS = KLFREE
  KLFREE = KLFREE+MXCSFC*MXCSFC

  IILB = 1
  do ICNL=1,NPCNF
    call c_f_pointer(c_loc(SCR(KLCONF)),iSCRl,[1])
    call GETCNF_LUCIA(iSCRl,ILTYP,IPCNF(ICNL),ICONF,IREFSM,NEL)
    nullify(iSCRl)
    NCSFL = NCSFTP(ILTYP)
    IIRB = 1
    do ICNR=1,ICNL
      call c_f_pointer(c_loc(SCR(KRCONF)),iSCRr,[1])
      call GETCNF_LUCIA(iSCRr,IRTYP,IPCNF(ICNR),ICONF,IREFSM,NEL)
      nullify(iSCRr)
      NCSFR = NCSFTP(IRTYP)
      call c_f_pointer(c_loc(SCR(KLCONF)),iSCRl,[1])
      call c_f_pointer(c_loc(SCR(KRCONF)),iSCRr,[1])
      call CNHCN(iSCRl,ILTYP,iSCRr,IRTYP,SCR(KLPHPS),SCR(KLFREE),NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,NTEST,ExFac,IREOTS)
      nullify(iSCRl,iSCRr)
      do IIL=1,NCSFL
        if (IILB == IIRB) then
          IIRMAX = IIL
        else
          IIRMAX = NCSFR
        end if
        do IIR=1,IIRMAX
          IIRACT = IIRB-1+IIR
          IILACT = IILB-1+IIL
          ILRI = (IIR-1)*NCSFL+IIL
          ILRO = ((IILACT*IILACT-IILACT)/2)+IIRACT
          PHP(ILRO) = SCR(KLPHPS-1+ILRI)
        end do
      end do
      IIRB = IIRB+NCSFR
    end do
    IILB = IILB+NCSFL
  end do

  return

end subroutine PHPCSF_INTERNAL

end subroutine PHPCSF
