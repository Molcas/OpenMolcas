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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine ORBCTL(CMO,NCMO,TORB,NTORB,FIFA,nFIFA,FIMO,nFIMO)
! Calculate transformation matrix to PT2 orbitals, defined as those
! that have standard Fock matrix FIFA diagonal within inactive,
! active, and secondary subblocks.

use fciqmc_interface, only: DoFCIQMC
use Printlevel, only: DEBUG, VERBOSE
use caspt2_global, only: iPrGlb
use caspt2_module, only: bName, EPS, nBas, nBasT, nDel, nFro, nOrb, nSym, OutFmt, PrOrb, ThrEne, ThrOcc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Five
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NCMO, NTORB, nFIFA, nFIMO
real(kind=wp), intent(inout) :: CMO(NCMO)
real(kind=wp), intent(out) :: TORB(NTORB)
real(kind=wp), intent(inout) :: FIFA(nFIFA), FIMO(nFIMO)
integer(kind=iwp) :: I1, I2, ISYM
real(kind=wp) :: OCC_DUM(1)
real(kind=wp), allocatable :: OrbE(:)

! Determine PT2 orbitals, and transform CI coeffs.
if (IPRGLB >= DEBUG) write(u6,*) ' ORBCTL calling MKRPTORB...'
! The CMO coefficient array is changed by orthonormal transformations of
! each of the inactive,Ras1,Ras2,Ras3,and secondary (i.e. virtual) orbitals
! in each symmetry. The transformation matrices are stored as a sequence of
! square matrices in TORB. The transformation is such that each of the
! diagonal blocks of the Fock matrix FIFA is diagonalized.
! MKRPTORB will at the same time transform each of the CI arrays on file
! such that, with new CMO vectors, they still represent the original
! wave function.

! The CI arrays are on file with unit number LUCIEX. There is NSTATE
! CI arrays, stored sequentially. The original set starts at disk address
! IDCIEX(1), the transformed ones are written after IDTCEX(1).

call MKRPTORB(FIFA,nFIFA,TORB,nTORB,CMO,NCMO)

if (IPRGLB >= DEBUG) write(u6,*) ' ORBCTL back from MKRPTORB.'

! Use the transformation matrices to change the FIMO and FIFA arrays:
if (.not. DoFCIQMC) then
  call TRANSFOCK(TORB,nTORB,FIMO,size(FIMO),1)
  call TRANSFOCK(TORB,nTORB,FIFA,nFIFA,1)

  if (IPRGLB >= DEBUG) write(u6,*) ' ORBCTL back from TRANSFOCK.'
end if

if (IPRGLB >= VERBOSE) then
  ! Print new orbitals. First, form array of orbital energies.
  call mma_allocate(ORBE,NBAST,Label='ORBE')
  I1 = 1
  I2 = 1
  do ISYM=1,NSYM
    if (NFRO(ISYM) > 0) then
      ORBE(I2:I2+NFRO(ISYM)-1) = Zero
      I2 = I2+NFRO(ISYM)
    end if
    if (NORB(ISYM) > 0) then
      ORBE(I2:I2+NORB(ISYM)-1) = EPS(I1:I1+NORB(ISYM)-1)
      I1 = I1+NORB(ISYM)
      I2 = I2+NORB(ISYM)
    end if
    if (NDEL(ISYM) > 0) then
      ORBE(I2:I2+NDEL(ISYM)-1) = Zero
      I2 = I2+NDEL(ISYM)
    end if
  end do
  ! Then call utility routine PRIMO.
  write(u6,*) ' The internal wave function representation has been changed to use quasi-canonical orbitals:'
  write(u6,*) ' those which diagonalize the Fock matrix within inactive-inactive,'
  write(u6,*) ' active-active and virtual-virtual submatrices.'
  if (.not. PRORB) then
    write(u6,*) ' On user''s request, the quasi-canonical orbitals'
    write(u6,*) ' will not be printed.'
    call mma_deallocate(ORBE)
    return
  end if

  ! Print orbitals. Different options:
  if (OUTFMT == 'LONG    ') then
    THRENE = Two**31
    THROCC = -Two**31
  else if (OUTFMT == 'DEFAULT ') then
    THRENE = Five
    THROCC = 5.0e-4_wp
  end if
  call PRIMO(' Quasi-canonical orbitals',.false.,.true.,THROCC,THRENE,NSYM,NBAS,NBAS,BNAME,ORBE,OCC_DUM,CMO,-1)
  call mma_deallocate(ORBE)
end if

end subroutine ORBCTL
