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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Readin_vb()
!***********************************************************************
!                                                                      *
!     Read the input                                                   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use gas_data, only: IGSOCCX, NGAS, NGSSH
use rasscf_global, only: iRlxRoot, iROOT, iZROT, NAC, NACPAR, NACPR2, NIN, NO2M, NORBT, NROOTS, NSEC, NTIT, NTOT3, NTOT4, Title
use jobiph_j, only: ispin_j, lsym_j, nactel_j, ndel_j, nelec3_j, nfro_j, nhole1_j, nish_j, nrs1_j, nrs2_j, nrs3_j, title_j
use general_data, only: INVEC, ISPIN, NACTEL, NASH, NBAS, NDEL, NDELT, NELEC3, NFRO, NFROT, NHOLE1, NISH, NORB, NRS1, NRS1T, NRS2, &
                        NRS2T, NRS3, NRS3T, NSSH, NSYM, NTOT, NTOT1, NTOT2, NTOTSP, STSYM, STSYM
use Molcas, only: MxSym
use RASDim, only: MxTit
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IGAS, II, ISYM, ITU, J, NAO, NGSSH_HI, NGSSH_LO, NT, NU

!---  set INVEC -> get MOs from JOBIPH file ---------------------------*
INVEC = 3
!---  Initialize ------------------------------------------------------*
do j=1,18
  Title(j) = ' '
end do

!---  Read input from standard input ----------------------------------*
!---  process TITLE    command ----------------------------------------*
ntit = 0
do ii=1,mxtit
  if (len_trim(title_j(ii)) > 0) then
    ntit = ntit+1
    title(ntit) = title_j(ii)
  end if
end do
!---  process NACT command --------------------------------------------*
nactel = nactel_j
nhole1 = nhole1_j
nelec3 = nelec3_j
!---  process SPIN command --------------------------------------------*
ispin = ispin_j
!---  process SYMM command --------------------------------------------*
stsym = lsym_j
!---  process FROZ command --------------------------------------------*
nfro(:) = nfro_j(:)
!---  process INAC command --------------------------------------------*
nish(:) = nish_j(:)
!---  process RAS1 command --------------------------------------------*
nrs1(:) = nrs1_j(:)
!---  process RAS2 command --------------------------------------------*
nrs2(:) = nrs2_j(:)
!---  process RAS3 command --------------------------------------------*
nrs3(:) = nrs3_j(:)
!---  process DELE command --------------------------------------------*
ndel(:) = ndel_j(:)

if ((nroots > 1) .and. (irlxroot == 0)) iRlxRoot = iroot(nroots)
if (nroots == 1) iRlxRoot = 0
!---  complete orbital specifications ---------------------------------*
do iSym=1,mxsym
  NASH(ISYM) = NRS1(ISYM)+NRS2(ISYM)+NRS3(ISYM)
  NORB(ISYM) = NBAS(ISYM)-NFRO(ISYM)-NDEL(ISYM)
  NSSH(ISYM) = NORB(ISYM)-NISH(ISYM)-NASH(ISYM)
end do

! SVC: convert CAS/RAS to general GAS description here, then we only
! need to copy it for lucia later, which always uses GAS description.
NGSSH(1,1:NSYM) = NRS1(1:NSYM)
NGSSH(2,1:NSYM) = NRS2(1:NSYM)
NGSSH(3,1:NSYM) = NRS3(1:NSYM)
IGSOCCX(1,1) = max(2*sum(NRS1(1:NSYM))-NHOLE1,0)
IGSOCCX(1,2) = 2*sum(NRS1(1:NSYM))
IGSOCCX(2,1) = NACTEL-NELEC3
IGSOCCX(2,2) = NACTEL
IGSOCCX(3,1) = NACTEL
IGSOCCX(3,2) = NACTEL

!---  Compute IZROT. IZROT is a matrix (lower triangular over the -----*
!     active space), which specifies which t,u rotations should be
!     avoided, since the orbitals belong to the same RAS space.
!     This is the only way the RAS concept is explicitly used in the
!     SX section of the program.
ITU = 0
do ISYM=1,mxsym
  NAO = NASH(ISYM)
  if (NAO > 1) then
    do NT=2,NAO
      do NU=1,NT-1
        ITU = ITU+1
        IZROT(ITU) = 0
        !SVC: check if NU<NT are included in the same gas space
        NGSSH_LO = 0
        do IGAS=1,NGAS
          NGSSH_HI = NGSSH_LO+NGSSH(IGAS,ISYM)
          if ((NU > NGSSH_LO) .and. (NT <= NGSSH_HI)) IZROT(ITU) = 1
          NGSSH_LO = NGSSH_HI
        end do
      end do
    end do
  end if
end do
!---  complete the input processing -----------------------------------*
NTOT = 0
NTOT1 = 0
NTOT2 = 0
NO2M = 0
NIN = 0
NAC = 0
NDELT = 0
NFROT = 0
NSEC = 0
NORBT = 0
NTOT3 = 0
NTOTSP = 0
NTOT4 = 0
NRS1T = 0
NRS2T = 0
NRS3T = 0
do ISYM=1,NSYM
  NTOT = NTOT+NBAS(ISYM)
  NTOT1 = NTOT1+(NBAS(ISYM)*(NBAS(ISYM)+1))/2
  NTOT2 = NTOT2+NBAS(ISYM)**2
  NO2M = max(NO2M,NBAS(ISYM)**2)
  NRS1T = NRS1T+NRS1(ISYM)
  NRS2T = NRS2T+NRS2(ISYM)
  NRS3T = NRS3T+NRS3(ISYM)
  NFROT = NFROT+NFRO(ISYM)
  NIN = NIN+NISH(ISYM)
  NAC = NAC+NASH(ISYM)
  NDELT = NDELT+NDEL(ISYM)
  NSEC = NSEC+NSSH(ISYM)
  NORBT = NORBT+NORB(ISYM)
  NTOT3 = NTOT3+(NORB(ISYM)+NORB(ISYM)**2)/2
  NTOTSP = NTOTSP+(NASH(ISYM)*(NASH(ISYM)+1)/2)
  NTOT4 = NTOT4+NORB(ISYM)**2
end do
NACPAR = (NAC+NAC**2)/2
NACPR2 = (NACPAR+NACPAR**2)/2

call Put_iArray('nIsh',nIsh,nSym)
call Put_iArray('nAsh',nAsh,nSym)

end subroutine Readin_vb
