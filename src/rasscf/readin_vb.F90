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

use Index_Functions, only: nTri_Elem
use gas_data, only: IGSOCCX, NGAS, NGSSH
use rasscf_global, only: iRlxRoot, iROOT, iZROT, NAC, NACPAR, NACPR2, NIN, NO2M, NORBT, NROOTS, NSEC, NTIT, NTOT3, NTOT4, Title
use jobiph_j, only: ispin_j, lsym_j, nactel_j, ndel_j, nelec3_j, nfro_j, nhole1_j, nish_j, nrs1_j, nrs2_j, nrs3_j, title_j
use general_data, only: INVEC, ISPIN, NACTEL, NASH, NBAS, NDEL, NDELT, NELEC3, NFRO, NFROT, NHOLE1, NISH, NORB, NRS1, NRS1T, NRS2, &
                        NRS2T, NRS3, NRS3T, NSSH, NSYM, NTOT, NTOT1, NTOT2, NTOTSP, STSYM, STSYM
use Molcas, only: MxSym
use RASDim, only: MxTit
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IGAS, II, ISYM, ITU, NAO, NGSSH_HI, NGSSH_LO, NT, NU

!---  set INVEC -> get MOs from JOBIPH file ---------------------------*
INVEC = 3
!---  Initialize ------------------------------------------------------*
Title(:) = ' '

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
NASH(1:mxsym) = NRS1(1:mxsym)+NRS2(1:mxsym)+NRS3(1:mxsym)
NORB(1:mxsym) = NBAS(1:mxsym)-NFRO(1:mxsym)-NDEL(1:mxsym)
NSSH(1:mxsym) = NORB(1:mxsym)-NISH(1:mxsym)-NASH(1:mxsym)

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
NTOT = sum(NBAS(1:NSYM))
NTOT1 = sum(nTri_Elem(NBAS(1:NSYM)))
NTOT2 = sum(NBAS(1:NSYM)**2)
NO2M = maxval(NBAS(1:NSYM)**2)
NIN = sum(NISH(1:NSYM))
NAC = sum(NASH(1:NSYM))
NDELT = sum(NDEL(1:NSYM))
NFROT = sum(NFRO(1:NSYM))
NSEC = sum(NSSH(1:NSYM))
NORBT = sum(NORB(1:NSYM))
NTOT3 = sum(nTri_Elem(NORB(1:NSYM)))
NTOTSP = sum(nTri_Elem(NASH(1:NSYM)))
NTOT4 = sum(NORB(1:NSYM)**2)
NRS1T = sum(NRS1(1:NSYM))
NRS2T = sum(NRS2(1:NSYM))
NRS3T = sum(NRS3(1:NSYM))
NACPAR = nTri_Elem(NAC)
NACPR2 = nTri_Elem(NACPAR)

call Put_iArray('nIsh',nIsh,nSym)
call Put_iArray('nAsh',nAsh,nSym)

end subroutine Readin_vb
