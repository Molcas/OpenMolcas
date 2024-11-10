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
! Copyright (C) 2018, Andrew M. Sand                                   *
!               2019, Thais R. Scott                                   *
!               2021, Jie J. Bao                                       *
!***********************************************************************
!                                                                      *
! 2021, Jie J. Bao - created file                                      *
! 2024, Matthew R. Hennefarth - upgraded to F90                        *
!***********************************************************************

Subroutine SaveFock_PDFT(CMO,hcore,coul,D1Act,NQ,p2d,state)
  use definitions,only:iwp,wp,u6
  use constants,only:zero
  use printlevel,only:debug
  use mcpdft_output,only:iPrLoc
  use stdalloc,only:mma_allocate,mma_deallocate
  use wadr,only:fockocc
  use rasscf_global,only:nFint,ISTORP,ntot4
  use mspdftgrad,only:F1MS,F2MS,FocMS

! Notes: Two references will be referred to in the comments.
! Ref1:  Sand, et al. JCTC, 2018, 14,  126.
! Ref2: Scott, et al. JCP,  2020, 153, 014106.
  Implicit None

  real(kind=wp) :: CMO(*),D1Act(*),hcore(*),coul(*)
  real(kind=wp) :: P2D(*)

  integer(kind=iwp) :: NQ,state

#include "rasdim.fh"
#include "general.fh"

  integer(kind=iwp) :: I,ISA,iprlev
  real(kind=wp),allocatable :: ONTOPT(:),ONTOPO(:)
  real(kind=wp),allocatable :: FA_V(:),FI_V(:)
  real(kind=wp),allocatable :: fa_t(:),Q(:),focki(:),dm2(:)
  real(kind=wp) :: fock(ntot4)

  fock(:) = zero

  write(u6,'(2X,A)') 'Calculating potentials for analytic gradients for MS-PDFT'

  IPRLEV = IPRLOC(3)

  call mma_allocate(focki,ntot1,label='focki')
  call ao2mo_1particle(cmo,hcore(:ntot1)+coul(:ntot1),focki,nsym,nbas,norb,nfro)

! loading one-electron potential and two-electron potential
! Used as F1 and F2 in equations 58 and 59 in Ref1.
  Call mma_allocate(ONTOPT,nfint,Label='OnTopT')
  Call mma_allocate(ONTOPO,ntot1,Label='OnTopO')
  Call Get_dArray('ONTOPT',OnTopT,NFINT)
  Call Get_dArray('ONTOPO',OnTopO,NTOT1)

  If(IPRLEV >= DEBUG) THEN
    write(u6,*) 'One-electron potentials'
    do i = 1,ntot1
      write(u6,*) OnTopO(i)
    enddo
    write(u6,*) 'Two-electron potentials'
    DO i = 1,nfint
      write(u6,*) OnTopT(i)
    ENDDO
  ENDIF

  CALL mma_allocate(FI_V,Ntot1,Label='FI_V')
  Call Get_dArray('FI_V',FI_V,NTOT1)

  F1MS(:,state) = focki(:)+fi_v(:)+ontopo(:)

  Call Get_TUVX(OnTopT,f2ms(:,state))

! ____________________________________________________________
! This next part is to generate the MC-PDFT generalized fock operator.

! The corrections (from the potentials) to FI and FA are built in the NQ
! part of the code, for efficiency's sake.  It still needs to be
! debugged.
  CALL mma_allocate(FA_V,Ntot1,Label='FA_V')
  Call Get_dArray('FA_V',FA_V,NTOT1)

  IF(IPRLEV >= DEBUG) THEN
    write(u6,*) "extra terms to update FI"
    DO i = 1,ntot1
      write(u6,*) FI_V(i)
    ENDDO
    write(u6,*) "extra terms to update FA"
    DO i = 1,ntot1
      write(u6,*) FA_V(i)
    ENDDO
    CALL mma_allocate(FA_t,Ntot1,Label='FA_t')
    fa_t(:) = ontopo(:)+fi_v(:)+fa_v(:)
    write(u6,*) "Total F additions:"
    Call TriPrt(' ','(5G18.10)',FA_T,norb(1))
    CALL mma_deallocate(FA_t)
  ENDIF

  fi_v(:) = focki(:)+OnTopO(:)+FI_V(:)

  Call mma_deallocate(OnTopO)

!Reordering of the two-body density matrix.
  IF(ISTORP(NSYM+1) > 0) THEN
    call mma_allocate(dm2,ISTORP(NSYM+1),label="dm2")
    CALL PMAT_RASSCF(P2d,dm2)
  else
    call mma_allocate(dm2,1,label="dm2")
  ENDIF

!Must add to existing FOCK operator (occ/act). FOCK is not empty.
  CALL mma_allocate(Q,NQ,Label='Q') ! q-matrix(1symmblock)
  CALL FOCK_update(FOCK,fi_v,fa_v,D1Act,dm2,Q,OnTopT,CMO)
  Call mma_deallocate(Q)
  Call mma_deallocate(OnTopT)
  call mma_deallocate(focki)
  call mma_deallocate(dm2)
  CALL mma_deallocate(FI_V)
  CALL mma_deallocate(FA_V)

  focms(:,state) = fockocc(:)
  IF(IPRLEV >= DEBUG) THEN
    write(u6,*) 'FOCC_OCC'
    call wrtmat(fockocc,1,ntot1,1,ntot1)
    write(u6,*) 'DONE WITH NEW FOCK OPERATOR'
  ENDIF

  iSA = 1
  Call Put_iScalar('SA ready',iSA)

EndSubroutine SaveFock_PDFT
