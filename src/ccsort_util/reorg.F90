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
!               1995,1996, Pavel Neogrady                              *
!***********************************************************************

subroutine REORG(run_triples,IRETURN)
!----------------------------------------------------------------------*
!     1994  PER-AAKE MALMQUIST                                         *
!     DEPARTMENT OF THEORETICAL CHEMISTRY                              *
!     UNIVERSITY OF LUND, SWEDEN                                       *
!                                                                      *
! modified by P.N. (Dec. 1995)                                         *
! all allocation memory routines removed by P.N. (8.03.1996)           *
!----------------------------------------------------------------------*
!
! FILES USED:
!     TRAINT    2 electron MO INTEGRALS
!     JOBIPH    THE JOB-INTERFACE FILE AS PRODUCED BY THE RASSCF
!               PROGRAM
! INPUT
!     AT PRESENT:  THE INPUT FILE IS SEARCHED FOR THE STRING
!     '&REORG '. Input is only TITLE
!
! LIMITATIONS
!     Like in CASPT2
!
!***********************************************************************

use ccsort_global, only: clopkey, Escf, fullprint, IADR15, JOBIPH, LROOT, LUINTM, mbas, NASH, NDEL, ndelr, NFRO, nfror, NISH, noa, &
                         nob, NORB, NSSH, NSYM, nva, nvb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(out) :: run_triples
integer(kind=iwp), intent(out) :: IRETURN
#include "rasdim.fh"
integer(kind=iwp) :: i, iad15, ij, j, lad15, NOIPSB(106), ntot2, ntot3
logical(kind=iwp) :: run_sort
real(kind=wp), allocatable :: Ene(:,:), EPS(:), EPSRAS(:), FI(:), FIRAS(:), FOKA(:), FOKB(:)
integer(kind=iwp), external :: iPrintLevel

fullprint = 0
if (iPrintLevel(-1) <= 0) fullprint = -1
call mma_allocate(FIRAS,mbas*mbas,Label='FIRAS')
call mma_allocate(FI,mbas*mbas,Label='FI')

! READ AND ECHO INPUT DATA, READ JOBIPH, PRINT INPUT DATA,
call RDINPPN(run_triples,run_sort)
if (fullprint >= 0) call PRINPPN()
call CHKINP_CCSORT()

if (run_sort) then

  ! read FI from JOBIPH
  ntot3 = 0
  ntot2 = 0
  do i=1,nsym
    ntot3 = ntot3+(norb(i)*(norb(i)+1))/2
    ntot2 = ntot2+norb(i)
  end do

  ! pick the total energy from the JOBIPH file

  call mma_allocate(Ene,mxRoot,mxIter)

  iad15 = iadr15(6)
  lad15 = mxroot*mxiter
  call dDaFile(JOBIPH,2,Ene,lad15,iad15)
  EScf = Zero
  i = 1

  ! take the last non-zero energy stored

  do while ((Ene(LROOT,i) /= Zero) .and. (i <= mxIter))
    Escf = Ene(LROOT,i)
    i = i+1
  end do
  call mma_deallocate(Ene)
  if (fullprint >= 0) then
    write(u6,*)
    write(u6,'(6X,A,F16.8)') 'SCF energy:',Escf
    write(u6,'(6X,A)') '-----------'
    write(u6,*)
  end if

  ! get fi from previous RASSCF

  iad15 = iadr15(10)
  call ddafile(JOBIPH,2,firas(1),ntot3,iad15)

  ! get eps from previous RASSCF

  call mma_allocate(eps,mbas,label='eps')
  call mma_allocate(epsras,ntot2,label='epsras')

  iad15 = iadr15(11)
  call ddafile(JOBIPH,2,epsras,ntot2,iad15)

  ! reduce fi,eps and update n's
  call mod1(nsym,nfro,nish,nssh,ndel,norb,nfror,ndelr,firas,fi,epsras,eps)

  call mma_deallocate(epsras)

  ! def diagonal Fok for closed shell

  if (clopkey == 2) call mod2(nsym,nish,nash,norb,fi,eps)

  ! define noa,nob,nva,nvb

  noa(1:nsym) = nish(1:nsym)+nash(1:nsym)
  nob(1:nsym) = nish(1:nsym)
  nva(1:nsym) = nssh(1:nsym)
  nvb(1:nsym) = nssh(1:nsym)+nash(1:nsym)

  if (nsym < 8) then
    noa(nsym+1:8) = 0
    nob(nsym+1:8) = 0
    nva(nsym+1:8) = 0
    nvb(nsym+1:8) = 0
  end if

  if (fullprint > 1) then
    write(u6,*)
    write(u6,'(6X,A)') 'Diagonal Fock matrix elements and orbital energies:'
    write(u6,'(6X,A)') '---------------------------------------------------'
    write(u6,*)
    write(u6,'(6X,A)') '----------------------------------------'
    write(u6,'(6X,A)') '   i      F(i,i)           eps(i)       '
    write(u6,'(6X,A)') '----------------------------------------'
    ij = 0
    do i=1,norb(1)
      do j=1,i
        ij = ij+1
        if (i == j) write(u6,'(6X,I4,2F18.10)') i,fi(ij),eps(i)
      end do
    end do
    write(u6,'(6X,A)') '----------------------------------------'
    write(u6,*)
  end if

  ! prepare address (stupid)

  !FUE The unit number of the transformed two electron integrals
  !FUE must be 40, 50, 60, 70, 80 or 90. Any other number will
  !FUE not be compatible with the I/O driver in MOLCAS.

  LUINTM = 40
  !JR  call DANAME(LUINTM,'TRAINT')
  call DANAME_MF(LUINTM,'TRAINT')
  call mkaddress(NOIPSB)

  ! open TRAINT and call action

  call mma_allocate(FOKA,mbas*(mbas+1)/2,label='FOKA')
  call mma_allocate(FOKB,mbas*(mbas+1)/2,label='FOKB')

  call action_ccsort(FOKA,FOKB,fi,eps)
  call mma_deallocate(FOKA)
  call mma_deallocate(FOKB)
  call mma_deallocate(eps)

  ! close files

  call daclos(luintm)
  call daclos(jobiph)

else
  ! case, when SORT was skipped
  write(u6,*) ' SORT part was skipped'
  write(u6,*) ' Input parameters are from last actual run of SORT'
end if

ireturn = 0

call mma_deallocate(FIRAS)
call mma_deallocate(FI)

return

end subroutine REORG
