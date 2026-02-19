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

subroutine CREIPH()
! RASSCF program: version IBM 3090: input section
!
! Purpose: to set up an result file (JOBIPH = 15) where the most
!          important results of the calculation are collected as
!          an interphase to other programs. JOBIPH is a direct
!          acces file where the first record is a table of content.
!          The address to this record is zero. The table of content
!          has presently the length 15 and is stored in
!          COMMON/INTAUX/ as IADR15(15).
! Subsections of JOBIPH:
! IADR15(1):  List of input data written here with WR_RASSCF_Info
! IADR15(2):  Molecular orbitals first written here and then updated
!             at the end of SXCTL in each MCSCF itration.
!             To be used as starting orbitals for restarts|
!             Modified finally by NEWORB to pseudo-nat orbitals
!             Average occupation numbers appended.
! IADR15(3):  density matrices for active orbitals
!             (4 matrices per root)
!             D  : one-body density matrix
!             DS : spin density matrix
!             P  : symmetric two-body density matrix
!             PA : antsymmetric two-body density matrix
! IADR15(4):  The CI-vectors corresponding to the LROOTS lowest
!             roots written in CICTL after each iteration.
! IADR15(5):  The Fock matrix for the occupied orbitals.
!             written in FOCKOC at the end of the calculation.
! IADR15(6):  Energies given in array ENER(10,50) in COMMON/RELAUX/.
! IADR15(7):  Convergence parameters in array CONV(6,50) in RELAUX.
! IADR15(8):  Not used
! IADR15(9):  Canonical MO's constructed and written in FCKPT2 at
!             the end of the calculation.
! IADR15(10): Inactive Fock matrix for CASPT2 in MO basis
!             Note: Frozen and deleted orbitals not included,
!             Modified to include the CASPT2 matrix FP 900608 (BOR)
! IADR15(11): The diagonal of the CASPT Fock matrix.
! IADR15(12): Natural orbitals for the final wave functions.
!             Note: These orbitals do not correspond to the same
!             CI wave function and should only be used to compute
!             properties.
!             One set of orbitals is stored for each root.
!             Each set of NO's is followed by a record
!             containing the corresponding occupation numbers.
! IADR15(13): One-body spin density matrices (one for each root)
! IADR15(14): Spin orbitals for the final wave functions.
!             Note: These orbitals do not correspond to the same
!             CI wave function and should only be used to compute
!             spin properties.
!             One set of orbitals for each root in an average
!             calculation. Stored sequentially (891219).
!             Each set of NO's is followed by a record
!             containing the corresponding occupation numbers.
! IADR15(15): First unwritten disk address, if positive;
!             If = -1, it means a new JOBIPH structure is used,
!             same as before, except:
! IADR15(16): First unwritten disk address (new layout),
! IADR15(17): Effective Hamiltonian, written here by MS-CASPT2,
!             LROOTS**2 elements, Fortran standard layout.
! IADR15(18): Table that translates between levels and absolute
!             orbital indices, setup by setsxci.
! IADR15(19)--IADR15(30): Presently unused.
! ********** IBM 3090 MOLCAS Release 90 02 22 **********

use sxci, only: IDXCI, IDXSX
use rasscf_global, only: BName, header, IADR15, IPT2, iRoot, lRoots, NACPAR, NACPR2, nOrbT, nRoots, NTOT3, POTNUC, Title, Weight
use general_data, only: ISPIN, JOBIPH, NACTEL, NASH, NBAS, NCONF, NDEL, NELEC3, NFRO, NHOLE1, NISH, NRS1, NRS2, NRS3, NSYM, NTOT, &
                        NTOT2, STSYM
use Molcas, only: LenIn, MxAct, MxOrb, MxRoot, MxSym
use RASDim, only: MxIter, MxTit
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: I, IAD15, NFOCK
real(kind=wp) :: Dum(1)
real(kind=wp), allocatable :: HEFF(:,:)

IADR15(1:15) = 0
! Dummy write table of contents.
! New layout scheme; length is 30 integers.
IAD15 = 0
call IDAFILE(JOBIPH,1,IADR15,30,IAD15)
IADR15(1) = IAD15

! DUMMY WRITE THE REMAINING RECORDS TO OBTAIN ADDRESSES

call WR_RASSCF_Info(JOBIPH,1,iAD15,NACTEL,ISPIN,NSYM,STSYM,NFRO,NISH,NASH,NDEL,NBAS,MxSym,BName,(LenIn+8)*mxOrb,NCONF,HEADER,144, &
                    TITLE,4*18*mxTit,POTNUC,LROOTS,NROOTS,IROOT,MxRoot,NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2,WEIGHT)
IADR15(2) = IAD15

call DDAFILE(JOBIPH,0,DUM,NTOT2,IAD15)
call DDAFILE(JOBIPH,0,DUM,NTOT,IAD15)
IADR15(3) = IAD15
do i=1,lRoots
  call DDafIle(JOBIPH,0,Dum,NACPAR,IAD15)
  call DDafIle(JOBIPH,0,Dum,NACPAR,IAD15)
  call DDafIle(JOBIPH,0,Dum,NACPR2,IAD15)
  call DDafIle(JOBIPH,0,Dum,NACPR2,IAD15)
end do
IADR15(4) = IAD15
do i=1,lRoots
  call DDafIle(JOBIPH,0,Dum,nConf,IAD15)
end do
IADR15(5) = IAD15
NFOCK = sum((NISH(1:NSYM)+NASH(1:NSYM))**2)
call DDAFILE(JOBIPH,0,DUM,NFOCK,IAD15)
IADR15(6) = IAD15
call DDAFILE(JOBIPH,0,DUM,mxRoot*mxIter,IAD15)
IADR15(7) = IAD15
call DDAFILE(JOBIPH,0,DUM,6*mxIter,IAD15)
IADR15(8) = IAD15
IADR15(9) = IAD15
call DDAFILE(JOBIPH,0,DUM,NTOT2,IAD15)
IADR15(10) = IAD15
call DDAFILE(JOBIPH,0,DUM,NTOT3,IAD15)
call DDAFILE(JOBIPH,0,DUM,NTOT3,IAD15)
IADR15(11) = IAD15
call DDAFILE(JOBIPH,0,DUM,NORBT,IAD15)
IADR15(12) = IAD15
do i=1,lRoots
  call DDAFILE(JOBIPH,0,DUM,NTOT2,IAD15)
  call DDAFILE(JOBIPH,0,DUM,NTOT,IAD15)
end do
IADR15(13) = IAD15
IADR15(14) = IAD15
do i=1,lroots
  call DDAFILE(JOBIPH,0,DUM,NTOT2,IAD15)
  call DDAFILE(JOBIPH,0,DUM,NTOT,IAD15)
end do
! New layout scheme:
IADR15(15) = -1
call mma_allocate(HEFF,LROOTS,LROOTS,Label='HEFF')
call unitmat(HEFF,LROOTS)
HEFF(:,:) = 1.0e12_wp*HEFF(:,:)
IADR15(17) = IAD15
call DDAFILE(JOBIPH,1,HEFF,LROOTS**2,IAD15)
call mma_deallocate(HEFF)
IADR15(18) = IAD15
!SVC: translates levels to orbital index (L2ACT in caspt2)
call IDAFILE(JOBIPH,1,IDXSX,mxAct,IAD15)
!SVC: translates orbital index to levels (LEVEL in caspt2)
call IDAFILE(JOBIPH,1,IDXCI,mxAct,IAD15)
IADR15(19:30) = 0
! First unused disk address at IADR(16):
IADR15(16) = IAD15

! Write table of content: IADR15

IAD15 = 0
call IDAFILE(JOBIPH,1,IADR15,30,IAD15)

end subroutine CREIPH
