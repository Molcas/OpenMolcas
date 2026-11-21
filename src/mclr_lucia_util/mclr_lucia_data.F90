!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

!> The data of MCLR's copy of Lucia, the CI library of MCLR behind the CI
!> interface (CI_Operations, backend MCLR_Lucia_CI). Only the Lucia
!> routines of mclr_util use it; MCLR itself (src/mclr) does not.
!>
!> Operator in action: I12 (2: one- and two-electron integrals), IST
!>  (1: singlet, 2: triplet operator), Square (integrals square, i.e. a
!>  non-symmetric operator), KAIN1/KINT2 (the integrals), pINT1 (offsets of
!>  the one-electron symmetry blocks), IREFSM (reference symmetry).
!> MS2 (2*MS), IDC (MS combinations), PSSIGN (their PS factor).
!> CSF <-> SD tables: i1, iAllo, iAnders, lconf, lldet; the file LUCSF2SD
!>  (FNCSF2SD).
!> CI space: NICISP, IASTFI/IBSTFI (alpha/beta string types), IACTI,
!>  MNR1IC/MXR1IC/MNR3IC/MXR3IC (RAS1/RAS3 electrons), NELCI/NAELCI/NBELCI
!>  (electrons), XISPSM (determinants per symmetry and CI space), MXSB,
!>  MXSOOB (largest blocks).
!> Run control: ICISTR, NOCSF, IDIAG, MAXI, MAXK, NOPART.
!> Orbitals: NACOB, NOCOB, NTOOB, NORB1..3, IREOTS (type -> symmetry
!>  order), ISMFTO, IBSO, NTSOB, IBTSOB, ITSOB, NOBPT, NOBPTS.
!> Spin: MULTSP, MS2P, MINOP, MAXOP, NTYP, NDPCNT, NCPCNT, NDTASM, NCSASM,
!>  NCNASM; configuration types of the CSF basis (intcsf): NTYP types of
!>  MINOP..MINOP+NTYP-1 open shells, NCNFTP(type,sym) configurations and
!>  NCSFTP(type) CSFs per configuration of each type.
module MCLR_Lucia_Data

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: MXINKA = 200, MXPCSM = 8, MXPCTP = 30, MXPICI = 30, MXPNGAS = 3, MXPOBS = 20, MXPORB = 500, &
                                MXPR4T = 10
integer(kind=iwp) :: i1, I12, IACTI(MXPICI), iAllo = 0, iAnders, IASTFI(MXPICI), IBSO(MXPOBS), IBSTFI(MXPICI), IBTSOB(3,MXPOBS), &
                     ICISTR, IDC, IDIAG, IREFSM, IREOTS(MXPORB), ISMFTO(MXPORB), IST, ITSOB(MXPORB), lconf, lldet, LUCSF2SD = 20, &
                     MAXI, MAXK, MAXOP, MINOP, MNR1IC(MXPICI), MNR3IC(MXPICI), MS2, MS2P, MULTSP, MXR1IC(MXPICI), MXR3IC(MXPICI), &
                     MXSB, MXSOOB, NACOB, NAELCI(MXPICI), NBELCI(MXPICI), NCNASM(MXPCSM), NCNFTP(MXPCTP,MXPCSM), NCPCNT(MXPCTP), &
                     NCSASM(MXPCSM), NCSFTP(MXPCTP), NDPCNT(MXPCTP), NDTASM(MXPCSM), NELCI(MXPICI), NICISP, NOBPT(MXPR4T+6), &
                     NOBPTS(MXPR4T+6,MXPOBS), NOCOB, NOCSF, NOPART, NORB1, NORB2, NORB3, NTOOB, NTSOB(3,MXPOBS), NTYP
real(kind=wp) :: PSSIGN, XISPSM(MXPCSM,MXPICI) = Zero
logical(kind=iwp) :: square
! TD_Ints: the two-body integrals of the run are in the square-pair layout
! (time-dependent response), set when the library is initialised
logical(kind=iwp) :: TD_Ints = .false.
! The CI space, the copies of the arguments of MCLR_Lucia_Init (the copy
! of Lucia uses no data of the programs): symmetries, orbitals, inactive
! and active orbitals and RAS spaces per symmetry, holes in RAS1 and
! electrons in RAS3, active electrons, multiplicity, symmetry of the CI
! state, CSFs per symmetry
integer(kind=iwp) :: L_iSpin = 0, L_nActEl = 0, L_nAsh(8) = 0, L_nCSF(8) = 0, L_nElec3 = 0, L_nHole1 = 0, L_nIsh(8) = 0, &
                     L_nOrb(8) = 0, L_nRs1(8) = 0, L_nRs2(8) = 0, L_nRs3(8) = 0, L_nSym = 0, L_STSym = 0
character(len=8) :: FNCSF2SD = 'REORD'
integer(kind=iwp), allocatable :: pINT1(:)
real(kind=wp), pointer :: KAIN1(:), KINT2(:), KINT2A(:)

public :: FNCSF2SD, i1, I12, IACTI, iAllo, iAnders, IASTFI, IBSO, IBSTFI, IBTSOB, ICISTR, IDC, IDIAG, IREFSM, IREOTS, ISMFTO, &
          IST, ITSOB, KAIN1, KINT2, KINT2A, lconf, lldet, LUCSF2SD, MAXI, MAXK, MAXOP, MINOP, MNR1IC, MNR3IC, MS2, MS2P, MULTSP, &
          MXINKA, MXPCSM, MXPNGAS, MXR1IC, MXR3IC, MXSB, MXSOOB, NACOB, NAELCI, NBELCI, NCNASM, NCNFTP, NCPCNT, NCSASM, NCSFTP, &
          NDPCNT, NDTASM, NELCI, NICISP, NOBPT, NOBPTS, NOCOB, NOCSF, NOPART, NORB1, NORB2, NORB3, NTOOB, NTSOB, NTYP, pINT1, &
          PSSIGN, square, TD_Ints, XISPSM
public :: L_iSpin, L_nActEl, L_nAsh, L_nCSF, L_nElec3, L_nHole1, L_nIsh, L_nOrb, L_nRs1, L_nRs2, L_nRs3, L_nSym, L_STSym

end module MCLR_Lucia_Data
