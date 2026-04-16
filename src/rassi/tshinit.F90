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

subroutine TSHinit(Energy)

use rasdef, only: NRAS, NRASEL, NRS1, NRS1T, NRS2, NRS3, NRSPRT
use rassi_aux, only: ipglob
use rassi_global_arrays, only: JBNUM, LROOT, PART
use gugx, only: CIStruct, EXStruct, SGStruct
use Cntrl, only: ChkHop, IRREP, ISTATE1, ISTATE2, LSYM1, LSYM2, MLTPLT, NACTE, nCI1, nCI2, NELE3, NHOLE1, NSTATE, RASTYP
use Symmetry_Info, only: nIrrep
use rassi_data, only: NDEL, NFRO, NISH, NSSH
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Energy(nState)
integer(kind=iwp) :: iRlxRoot, JOB1, JOB2, MPLET1, MPLET2, NACTE1, NACTE2, NELE31, NELE32, NHOL11, NHOL12
real(kind=wp) :: EDIFF
logical(kind=iwp) :: LOWROOT, UPROOT
type(SGStruct) :: SGS(2)
type(CIStruct) :: CIS(2)
type(EXStruct) :: EXS(2)
character(len=8) :: WFTYP1, WFTYP2
real(kind=wp), allocatable :: CI1(:), CI2(:)
real(kind=wp), parameter ::  Ethr = 0.03_wp

! Print a banner

if (IPGLOB >= 2) then
  write(u6,*)
  write(u6,*)
  write(u6,'(6X,A)') repeat('*',100)
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A,36X,A,37X,A)') '*',' Surface hopping section ','*'
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A)') repeat('*',100)
  write(u6,*)
  write(u6,*)

  write(u6,'(6X,A)') 'Surface hopping section'
  write(u6,'(6X,A)') '-----------------------'
end if

! Get the current state and print it's energy

call Get_iScalar('Relax CASSCF root',iRlxRoot)
if (IPGLOB >= 2) then
  write(u6,'(6X,A,I14)') 'The current state is:',iRlxRoot
  write(u6,'(6X,A,6X,ES15.6,A,/)') 'Its energy is:',ENERGY(iRlxRoot),' a.u.'
end if

! Get wave function parameters for current state

ISTATE1 = iRlxRoot
JOB1 = JBNUM(ISTATE1)
NACTE1 = NACTE(JOB1)
MPLET1 = MLTPLT(JOB1)
LSYM1 = IRREP(JOB1)
NHOL11 = NHOLE1(JOB1)
NELE31 = NELE3(JOB1)
WFTYP1 = RASTYP(JOB1)
SGS(1)%IFRAS = 1
SGS(2)%IFRAS = 1

! Set the variables for the wave function of the current state

if (WFTYP1 == 'GENERAL') then
  NRSPRT = 3
  NRAS(:,1) = NRS1(:)
  NRAS(:,2) = NRS2(:)
  NRAS(:,3) = NRS3(:)
  NRASEL(1) = 2*NRS1T-NHOL11
  NRASEL(2) = NACTE1-NELE31
  NRASEL(3) = NACTE1
  call SGINIT(nIrrep,NACTE1,MPLET1,SGS(1),CIS(1))
  if (IPGLOB > 4) then
    write(u6,*) 'Split-graph structure for JOB1=',JOB1
    call SGPRINT(SGS(1))
  end if
  call CXINIT(SGS(1),CIS(1),EXS(1))
  ! CI sizes, as function of symmetry, are now known.
  NCI1 = CIS(1)%NCSF(LSYM1)
else
  ! The only other cases are HISPIN, CLOSED or EMPTY.
  ! NOTE: The HISPIN case is suspected to be buggy. Not used now.
  NCI1 = 1
end if
call mma_allocate(CI1,NCI1,Label='CI1')

! Check for the possibility to hop to a lower root

if ((ISTATE1-1) >= 1) then
  LOWROOT = .true.
  if (IPGLOB >= 2) write(u6,'(6X,A,I2)') 'There is a lower root, which is: ',LROOT(ISTATE1-1)
else
  LOWROOT = .false.
  if (IPGLOB >= 2) write(u6,'(6X,A)') 'There is no lower root'
end if

! Check for the possibility to hop to an upper root

if ((ISTATE1+1) <= NSTATE) then
  UPROOT = .true.
  if (IPGLOB >= 2) write(u6,'(6X,A,I2,/)') 'There is an upper root, which is: ',LROOT(ISTATE1+1)
else
  UPROOT = .false.
  if (IPGLOB >= 2) write(u6,'(6X,A,/)') 'There is no upper root'
end if

! Check surface hopping to a root lower than the current one

if (LOWROOT) then
  ISTATE2 = ISTATE1-1
  if (IPGLOB >= 2) then
    write(u6,'(6X,A,I3)') 'The lower state is:',ISTATE2
    write(u6,'(6X,A,6X,ES15.6,A)') 'Its energy is:',ENERGY(ISTATE2),' a.u.'
    write(u6,'(6X,A,12X,ES15.6,A,/)') 'Ediff = ',ENERGY(iRlxRoot)-ENERGY(ISTATE2),' a.u.'
  end if
  !--------------------------------------------------------------------*
  ! Adapted from RASSI subroutines GTMCTL and READCI                   *

  ! Get wave function parameters for ISTATE2
  JOB2 = JBNUM(ISTATE2)
  NACTE2 = NACTE(JOB2)
  MPLET2 = MLTPLT(JOB2)
  LSYM2 = IRREP(JOB2)
  NHOL12 = NHOLE1(JOB2)
  NELE32 = NELE3(JOB2)
  WFTYP2 = RASTYP(JOB2)
  call NEWPRTTAB(nIrrep,NFRO,NISH,NRS1,NRS2,NRS3,NSSH,NDEL)
  if (IPGLOB >= 4) call PRPRTTAB(PART)
  ! For the second wave function
  if (WFTYP2 == 'GENERAL') then
    NRSPRT = 3
    NRAS(:,1) = NRS1(:)
    NRAS(:,2) = NRS2(:)
    NRAS(:,3) = NRS3(:)
    NRASEL(1) = 2*NRS1T-NHOL12
    NRASEL(2) = NACTE2-NELE32
    NRASEL(3) = NACTE2
    call SGINIT(nIrrep,NACTE2,MPLET2,SGS(2),CIS(2))
    if (IPGLOB > 4) then
      write(u6,*) 'Split-graph structure for JOB2=',JOB2
      call SGPRINT(SGS(2))
    end if
    call CXINIT(SGS(2),CIS(2),EXS(2))
    ! CI sizes, as function of symmetry, are now known.
    NCI2 = CIS(2)%NCSF(LSYM2)
  else
    ! Presently, the only other cases are HISPIN, CLOSED or EMPTY.
    NCI2 = 1
  end if
  call mma_allocate(CI2,NCI2,Label='CI2')
  ! Check if the Energy gap is smaller than the threshold.
  Ediff = abs(ENERGY(ISTATE2)-ENERGY(ISTATE1))
  if (Ediff <= Ethr) then
    ChkHop = .true.
  else
    ChkHop = .false.
  end if
  if (IPGLOB >= 2) then
    write(u6,'(6X,A,I8,4X,A,I8)') 'ISTATE1=',iRlxRoot,'ISTATE2=',ISTATE2
    write(u6,'(6X,A,I11,4X,A,I11)') 'NCI1=',NCI1,'NCI2=',NCI2
  end if
# ifdef _DEBUGPRINT_
  write(u6,*) ' TSHinit calls TSHop.'
# endif
  call TSHop(CI1,CI2)
# ifdef _DEBUGPRINT_
  write(u6,*) ' TSHinit back from TSHop.'
# endif
  if (WFTYP2 == 'GENERAL') call MkGUGA_Free(SGS(2),CIS(2),EXS(2))
  call mma_deallocate(CI2)
  call mma_deallocate(PART)
end if

! Check surface hopping to a root higher than the current one

if (UPROOT) then
  ISTATE2 = ISTATE1+1
  if (IPGLOB >= 2) then
    write(u6,'(6X,A,I3)') 'The upper state is:',ISTATE2
    write(u6,'(6X,A,6X,ES15.6,A)') 'Its energy is:',ENERGY(ISTATE2),' a.u.'
    write(u6,'(6X,A,12X,ES15.6,A,/)') 'Ediff = ',ENERGY(ISTATE2)-ENERGY(iRlxRoot),' a.u.'
  end if
  !--------------------------------------------------------------------*
  ! Adapted from RASSI subroutines GTMCTL and READCI                   *

  ! Get wave function parameters for ISTATE2
  JOB2 = JBNUM(ISTATE2)
  NACTE2 = NACTE(JOB2)
  MPLET2 = MLTPLT(JOB2)
  LSYM2 = IRREP(JOB2)
  NHOL12 = NHOLE1(JOB2)
  NELE32 = NELE3(JOB2)
  WFTYP2 = RASTYP(JOB2)
  call NEWPRTTAB(nIrrep,NFRO,NISH,NRS1,NRS2,NRS3,NSSH,NDEL)
  if (IPGLOB >= 4) call PRPRTTAB(PART)
  ! For the second wave function
  if (WFTYP2 == 'GENERAL') then
    NRSPRT = 3
    NRAS(:,1) = NRS1(:)
    NRAS(:,2) = NRS2(:)
    NRAS(:,3) = NRS3(:)
    NRASEL(1) = 2*NRS1T-NHOL12
    NRASEL(2) = NACTE2-NELE32
    NRASEL(3) = NACTE2
    call SGINIT(nIrrep,NACTE2,MPLET2,SGS(2),CIS(2))
    if (IPGLOB > 4) then
      write(u6,*) 'Split-graph structure for JOB2=',JOB2
      call SGPRINT(SGS(2))
    end if
    call CXINIT(SGS(2),CIS(2),EXS(2))
    ! CI sizes, as function of symmetry, are now known.
    NCI2 = CIS(2)%NCSF(LSYM2)
  else
    ! Presently, the only other cases are HISPIN, CLOSED or EMPTY.
    NCI2 = 1
  end if
  call mma_allocate(CI2,NCI2,Label='CI2')
  ! Check if the Energy gap is smaller than the threshold.
  Ediff = abs(ENERGY(ISTATE2)-ENERGY(ISTATE1))
  if (Ediff <= Ethr) then
    ChkHop = .true.
  else
    ChkHop = .false.
  end if
  if (IPGLOB >= 2) then
    write(u6,'(6X,A,I8,4X,A,I8)') 'ISTATE1=',iRlxRoot,'ISTATE2=',ISTATE2
    write(u6,'(6X,A,I11,4X,A,I11)') 'NCI1=',NCI1,'NCI2=',NCI2
  end if
# ifdef _DEBUGPRINT_
  write(u6,*) ' TSHinit calls TSHop.'
# endif
  call TSHop(CI1,CI2)
# ifdef _DEBUGPRINT_
  write(u6,*) ' TSHinit back from TSHop.'
# endif
  if (WFTYP2 == 'GENERAL') call MkGUGA_Free(SGS(2),CIS(2),EXS(2))
  call mma_deallocate(CI2)
  call mma_deallocate(PART)
end if

if (WFTYP1 == 'GENERAL') then
  call MkGUGA_Free(SGS(1),CIS(1),EXS(1))
end if
call mma_deallocate(CI1)

end subroutine TSHinit
