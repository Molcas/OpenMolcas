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
! Copyright (C) 1997, Per Ake Malmqvist                                *
!***********************************************************************

subroutine CREIPH_CASPT2(Heff,Ueff,U0,nState)
! Normal operation: A new file, 'JOBMIX', will be created, with the
! CMO's and CI arrays of the JOBIPH, except that the CI arrays have
! been modified. They are now linear combinations of the original ones,
! using coefficients taken from the eigenvectors of the effective
! Hamiltonian.
! Also, replace the original CASSCF energies with CASPT2 or MS-CASPT2
! energies.

use fciqmc_interface, only: DoFCIQMC
use PrintLevel, only: USUAL
use REFWFN, only: IADR15, REFWFN_FILENAME
use Molcas, only: LenIn, MxAct, MxLev, MxOrb, MxRoot
use RASDim, only: MxIter, MxTit
use sguga_states, only: SGS
use caspt2_global, only: CMO, CMO_Internal, iPrGlb, NCMO, Weight
use general_data, only: ISPIN, NACTEL, NASH,  NELE3=>nElec3, NHOLE1, STSYM
use caspt2_module, only: BNAME, CITHR, DOCUMULANT, ENERGY, HEADER, IFMIX, IFMSCOUP, IFQCAN, IFRMS, IFXMS, IROOT, ISCF, &
                         LROOTS, MSTATE, MSTATE, MXCI, NBAS, NBSQT, NCONF, NDEL, NFRO, NISH, NRAS1, &
                         NRAS2, NRAS3, NROOTS, NSYM, POTNUC, TITLE
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Nstate
real(kind=wp), intent(in) :: Heff(Nstate,Nstate), Ueff(Nstate,Nstate), U0(Nstate,Nstate)
integer(kind=iwp) :: I, IAD15, ID, IDISK, IDR, IDW, IISTATE, ISNUM, ISTATE, J, JOBIPH, JOBMIX, JSNUM, MROOTS, NIDIST, NOLDE
integer(kind=iwp), allocatable :: IDIST(:), JROOT(:), xL2Act(:), xLevel(:)
real(kind=wp), allocatable :: CI1(:), CI2(:), EFFCP(:), OLDE(:), Weight_(:)
integer(kind=iwp), parameter :: jstate=1

! Not called, if .not. IFMIX, then only the new CI coefficients are
! printed, no JOBMIX file is created.
if (IFMSCOUP .and. (ISCF == 0)) then
  if (.not. IFMIX) then
    if (IPRGLB >= USUAL) call PRINT_CI_MIX(Ueff)
    return
  end if
end if
if (DOCUMULANT .or. (.not. IFMIX)) return

if (DoFCIQMC) return

if (IFMSCOUP) then
  if (IPRGLB >= USUAL) then
    write(u6,*) ' THE ORIGINAL CI ARRAYS ARE NOW MIXED AS LINEAR'
    write(u6,*) ' COMBINATIONS, GIVEN BY THE EIGENVECTORS.'
  end if
end if

if (IPRGLB >= USUAL) then
  write(u6,*) ' A NEW JOBIPH FILE NAMED ''JOBMIX'' IS PREPARED.'
  write(u6,'(A)') repeat('*',80)
end if

! Note that JOBIPH file will contain all the RASSCF CI vectors
! plus a CASPT2 effective Hamiltonian for the selected states.
! The effective Hamiltonian for the states not included in the
! CASPT2 treatment will be diagonal with RASSCF energies!

! The JOBMIX will contain the (possibly mixed) CI vectors,
! with CASPT2 energies for the selected states and zero energy
! for states not included in the CASPT2 treatment
! If NoMulti was specified, the original state indexing is
! maintained, otherwise the new states are just 1, 2, 3...

call mma_allocate(CI1,MXCI,Label='CI1')
call mma_allocate(CI2,MXCI,Label='CI2')
JOBIPH = 15
call DANAME(JOBIPH,refwfn_filename)
JOBMIX = 11
call DANAME(JOBMIX,'JOBMIX')
! IADR15 is already known (it is the table of contents of the
! JOBIPH file). When copying/modifying selected data from JOBIPH
! to JOBMIX, we use the same TOC array, IADR15.
IAD15 = 0
call IDAFILE(JOBIPH,2,IADR15,30,IAD15)
IAD15 = 0
call IDAFILE(JOBMIX,1,IADR15,30,IAD15)
IAD15 = IADR15(1)
! Modify root index in case of MS
call mma_allocate(JROOT,MXROOT,LABEL='JROOT')
if (IFMSCOUP) then
  JROOT(1:NSTATE) = [(ISTATE,ISTATE=1,NSTATE)]
  JROOT(NSTATE+1:) = 0
  MROOTS = NSTATE
else
  JROOT(:) = IROOT(:)
  MROOTS = NROOTS
end if
call mma_allocate(Weight_,MxRoot,Label='Weight_')
WEIGHT_(1:NROOTS) = WEIGHT(1:NROOTS)
WEIGHT_(NROOTS+1:) = Zero
call WR_RASSCF_INFO(JOBMIX,1,iAd15,NACTEL,ISPIN,NSYM,STSYM,NFRO,NISH,NASH,NDEL,NBAS,8,BNAME,(LenIn+8)*MXORB,NCONF,HEADER,144, &
                    TITLE,4*18*MXTIT,POTNUC,LROOTS,MROOTS,JROOT,MXROOT,NRAS1,NRAS2,NRAS3,NHOLE1,NELE3,IFQCAN,Weight_)
call mma_deallocate(Weight_)
call mma_deallocate(JROOT)
! Copy MO coefficients from JOBIPH to JOBMIX
NCMO = NBSQT
call mma_allocate(CMO_Internal,NCMO,Label='CMO_Internal')
CMO => CMO_Internal
IAD15 = IADR15(9)
call DDAFILE(JOBIPH,2,CMO,NCMO,IAD15)
IAD15 = IADR15(9)
call DDAFILE(JOBMIX,1,CMO,NCMO,IAD15)
! If IFQCAN == 0, there is also an additional CMO set:
if (IFQCAN == 0) then
  IAD15 = IADR15(2)
  call DDAFILE(JOBIPH,2,CMO,NCMO,IAD15)
  IAD15 = IADR15(2)
  call DDAFILE(JOBMIX,1,CMO,NCMO,IAD15)
end if
call mma_deallocate(CMO_Internal)
nullify(CMO)
! Copy all CI coefficients
IDR = IADR15(4)
IDW = IADR15(4)
do I=1,LROOTS
  call DDAFILE(JOBIPH,2,CI1,NCONF,IDR)
  call DDAFILE(JOBMIX,1,CI1,NCONF,IDW)
end do
! Replace old energy array with (MS-)CASPT2 energy values:
NOLDE = MXROOT*MXITER
call mma_allocate(OLDE,NOLDE,LABEL='OLDE')
OLDE(:) = Zero
if (IFMSCOUP) then
  OLDE(1:NSTATE) = ENERGY(1:NSTATE)
else
  do ISTATE=1,NSTATE
    OLDE(MSTATE(ISTATE)) = ENERGY(ISTATE)
  end do
end if
IAD15 = IADR15(6)
call DDAFILE(JOBMIX,1,OLDE,NOLDE,IAD15)
call mma_deallocate(OLDE)
IAD15 = IADR15(18)
!SVC: translates levels to orbital index
!Copy to local array since L2Act and Level are protected.
call mma_allocate(xL2Act,MxLev,Label='xL2Act')
call mma_allocate(xLevel,MxLev,Label='xLevel')
XL2Act(:) = SGS(jstate)%L2Act(:)
call IDAFILE(JOBMIX,1,xL2ACT,mxAct,IAD15)
!SVC: translates orbital index to levels
XLevel(:) = SGS(jstate)%Level(:)
call IDAFILE(JOBMIX,1,xLEVEL,mxAct,IAD15)
call mma_deallocate(xL2Act)
call mma_deallocate(xLevel)

! PAM07: Eliminate unsafe IPOSFILE calls, use instead dummy i/o operations
! to find disk addresses to CI arrays:
NIDIST = maxval(MSTATE(1:NSTATE))
call mma_allocate(IDIST,NIDIST,Label='IDIST')
ID = IADR15(4)
do JSNUM=1,NIDIST
  IDIST(JSNUM) = ID
  ! This dummy operation does nothing, merely updates file pointer ID
  call DDAFILE(JOBIPH,0,CI1,NCONF,ID)
end do
! PAM07: Now IDIST() is used, instead of IPOSFILE, below!

! PAM05: Now CREIPH is called also in the NOMULT=1 case, to allow making
! a JOBMIX file also when NOMULT was ordered. Then the energies
! will be state-specific, of course. But no mixing of CI vectors.
if (IFMSCOUP) then
  ! Also write effective Hamiltonian on Jobiph file:
  call mma_allocate(EFFCP,LROOTS**2,Label='EFFCP')
  ! Read the effective Hamiltonian on JobIph file:
  IAD15 = IADR15(17)
  call DDAFILE(JOBIPH,2,EFFCP,LROOTS**2,IAD15)
  ! Replace the relevant elements:
  do I=1,NSTATE
    ISNUM = MSTATE(I)
    do J=1,NSTATE
      JSNUM = MSTATE(J)
      EFFCP(ISNUM+LROOTS*(JSNUM-1)) = HEFF(I,J)
    end do
  end do
  ! Write the present effective Hamiltonian:
  IAD15 = IADR15(17)
  call DDAFILE(JOBIPH,1,EFFCP,LROOTS**2,IAD15)
  ! Write a diagonal Hamiltonian in the JOBMIX:
  IAD15 = IADR15(17)
  EFFCP(:) = Zero
  do ISTATE=1,NSTATE
    EFFCP((ISTATE-1)*LROOTS+ISTATE) = ENERGY(ISTATE)
  end do
  call DDAFILE(JOBMIX,1,EFFCP,LROOTS**2,IAD15)
  call mma_deallocate(EFFCP)
  ! Now 'mix' those states that were treated in the multi-state CASPT2
  if (IPRGLB >= USUAL) then
    write(u6,*)
    call CollapseOutput(1,'Mixed CI coefficients:')
  end if
  do ISTATE=1,NSTATE
    CI2(:) = Zero
    do IISTATE=1,NSTATE
      JSNUM = MSTATE(IISTATE)
      IDISK = IDIST(JSNUM)
      call DDAFILE(JOBIPH,2,CI1,NCONF,IDISK)
      CI2(1:NCONF) = CI2(1:NCONF)+Ueff(IISTATE,ISTATE)*CI1(1:NCONF)
    end do
    if (ISCF == 0) then
      if (IPRGLB >= USUAL) then
        write(u6,'(1x,a,i3)') ' The CI coefficients for the MIXED state nr. ',ISTATE
        call PRWF_CP2(STSYM,NCONF,CI2,CITHR)
      end if
    end if
    IDISK = IDIST(ISTATE)
    call DDAFILE(JOBMIX,1,CI2,NCONF,IDISK)
  end do
  if (IPRGLB >= USUAL) then
    call CollapseOutput(0,'Mixed CI coefficients:')
    write(u6,*)
  end if
else if (IFXMS .or. IFRMS) then
  ! In case of XMS/XDW/RMS and NOMUL, the CI vectors are replaced by the
  ! rotated zeroth-order states (they should have been printed earlier,
  ! in grpini)
  do ISTATE=1,NSTATE
    ISNUM = MSTATE(ISTATE)
    CI2(:) = Zero
    do IISTATE=1,NSTATE
      JSNUM = MSTATE(IISTATE)
      IDISK = IDIST(JSNUM)
      call DDAFILE(JOBIPH,2,CI1,NCONF,IDISK)
      CI2(1:NCONF) = CI2(1:NCONF)+U0(IISTATE,ISTATE)*CI1(1:NCONF)
    end do
    IDISK = IDIST(ISNUM)
    call DDAFILE(JOBMIX,1,CI2,NCONF,IDISK)
  end do
end if

call mma_deallocate(IDIST)
call mma_deallocate(CI1)
call mma_deallocate(CI2)

call DACLOS(JOBIPH)
call DACLOS(JOBMIX)

end subroutine CREIPH_CASPT2
