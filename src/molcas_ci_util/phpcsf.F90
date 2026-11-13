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
!               2026, Roland Lindh                                     *
!***********************************************************************

subroutine PHPCSF(PHP,IPCSF,MXPDIM,IREFSM,ECORE,NPCSF,DIAG,TU,nTU,TUVX,nTUVX,NTEST)
! Obtain primary subspace and obtain
! explicit representation of hamilton matrix in subspace
!
! The routine is generic with respect to the CI library: the layout of the
! CSF basis by configuration and the explicit Hamiltonian are obtained
! through the CI interface (Mk_Cnf_Layout, Mk_IHJ).
!
! ARGUMENTS :
! ===========
! PHP    : hamilton matrix in subspace, lower triangle (output)
! IPCSF  : CSF's defining subspace (output)
! MXPDIM : Largest allowed dimension of subspace (input)
! IREFSM : symmetry of considered CI space (input)
! ECORE  : Core energy (input)
! NPCSF  : Number of primary CSF's obtained  (output)
! DIAG   : Hamilton diagonal over CSF's (input)
! TU     : One-electron integrals, triangular packing (input)
! TUVX   : Two-electron integrals (MO space), triangular packing (input)
! NTEST  : Print level (input)
!
! Jeppe Olsen, Summer of '89
! adapted to DETRAS by M.P. Fuelscher, October 1989
! made generic with respect to the CI library, R. Lindh, 2026

use CI_Interfaces, only: Mk_Cnf_Layout, Mk_IHJ
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: MXPDIM, IREFSM, nTU, nTUVX
real(kind=wp), intent(out) :: PHP(nTri_Elem(MXPDIM))
integer(kind=iwp), intent(out) :: IPCSF(MXPDIM), NPCSF
real(kind=wp), intent(in) :: ECORE, DIAG(*), TU(nTU), TUVX(nTUVX)
integer(kind=iwp), intent(inout) :: NTEST
integer(kind=iwp) :: i, ICNF, ICSFMN, IFINIT, IICSF, IMIN, NCONF, NCSFMN, NPCNF
real(kind=wp) :: XMAX, XMIN
integer(kind=iwp), allocatable :: nCSFCnf(:), IPCNF(:), iCSF1(:)
real(kind=wp), allocatable :: CnfDia(:)
real(kind=wp), parameter :: Acc = 1.0e-13_wp ! Assumed machine accuray (give and take)
real(kind=wp), external :: FNDMNX

! Layout of the CSF basis by configuration, and the diagonal value and
! first CSF of each configuration.
call Mk_Cnf_Layout(IREFSM,nCSFCnf)
NCONF = size(nCSFCnf)
call mma_allocate(CnfDia,NCONF,Label='CnfDia')
call mma_allocate(iCSF1,NCONF,Label='iCSF1')
call mma_allocate(IPCNF,NCONF,Label='IPCNF')
IICSF = 1
do ICNF=1,NCONF
  CnfDia(ICNF) = DIAG(IICSF)
  iCSF1(ICNF) = IICSF
  IICSF = IICSF+nCSFCnf(ICNF)
end do
if (NTEST >= 30) call RecPrt('CnfDia',' ',CnfDia,1,NCONF)

! find the elements of the subspace
! (MXPDIM elements smallest in energy, whole configurations)

XMAX = FNDMNX(CnfDia,NCONF,2)
NPCSF = 0
NPCNF = 0
IFINIT = 0
ICSFMN = 0
NCSFMN = 0
do
  XMIN = XMAX+One
  IMIN = 0
  do ICNF=1,NCONF
    if (CnfDia(ICNF)+Acc < XMIN) then
      XMIN = CnfDia(ICNF)
      IMIN = ICNF
      ICSFMN = iCSF1(ICNF)
      NCSFMN = nCSFCnf(ICNF)
    end if
  end do
  if ((NPCSF+NCSFMN) <= MXPDIM) then
    ! add new configuration
    NPCNF = NPCNF+1
    IPCNF(NPCNF) = IMIN
    IPCSF(NPCSF+1:NPCSF+NCSFMN) = [(i,i=ICSFMN,ICSFMN+NCSFMN-1)]
    NPCSF = NPCSF+NCSFMN
    CnfDia(IMIN) = XMAX+One
  else
    ! No space for this configuration.
    ! PAM 2011: Do not remove any previous configurations.
    IFINIT = 1
  end if
  if ((IFINIT /= 0) .or. (NPCNF >= NCONF)) exit
end do

if (NPCNF == 0) then
  call WarningMessage(2,'Making explicit Hamiltonian failed.')
  write(u6,*) ' An unforeseen catastrophic failure occurred'
  write(u6,*) ' in the CI solver. The size of the explicit'
  write(u6,*) ' part of the CI Hamiltonian matrix was not'
  write(u6,*) ' sufficient. Suggested fix: Change the size'
  write(u6,*) ' by adding ''SDAV=XXXXX'' to the rasscf input.'
  write(u6,*) ' XXXXX is some integer at least ',MXPDIM+NCSFMN
  write(u6,*) ' Sorry about this. Consider telling the'
  write(u6,*) ' Molcas group about this failure.'
  call Abend()
end if

if (NTEST >= 30) then
  write(u6,*) ' Output from PHPCSF'
  write(u6,*) ' =================='
  write(u6,*) ' Number of Configurations in primary subspace ',NPCNF
  write(u6,*) ' Number of CSFs in primary subspace ',NPCSF
  write(u6,*) ' Configurations included :'
  call IWRTMA(IPCNF,1,NPCNF,1,NPCNF)
  write(u6,*) ' CSFs included :'
  call IWRTMA(IPCSF,1,NPCSF,1,NPCSF)
end if

! construct Hamiltonian matrix in subspace

call Mk_IHJ(PHP,nTri_Elem(MXPDIM),IPCNF,NPCNF,IREFSM,ECORE,TU,nTU,TUVX,nTUVX,NTEST)

call mma_deallocate(IPCNF)
call mma_deallocate(iCSF1)
call mma_deallocate(CnfDia)
call mma_deallocate(nCSFCnf)

end subroutine PHPCSF
