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
Subroutine SaveFock_PDFT(CMO,FockI,FockA,D1Act,Fock,P,NQ,PUVX,p2d)
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Jan. 04, 2021, created this file.               *
! ****************************************************************

  use definitions,only:iwp,wp,u6
  use constants,only:zero
  use mspdft,only:iIntS
  use printlevel,only:debug
  use mcpdft_output,only:iPrLoc
  use stdalloc,only:mma_allocate,mma_deallocate
  use wadr,only:FockOcc
  use rasscf_global,only:NACPR2,nFint,ISTORP
  use mspdftgrad,only:F1MS,F2MS,FocMS

! Notes: Two references will be referred to in the comments.
! Ref1:  Sand, et al. JCTC, 2018, 14,  126.
! Ref2: Scott, et al. JCP,  2020, 153, 014106.
  Implicit None

  real(kind=wp) :: CMO(*),FockI(*),FockA(*),D1Act(*),Fock(*)
  real(kind=wp) :: P(*),PUVX(*),P2D(*)

  integer(kind=iwp) :: NQ

#include "rasdim.fh"
#include "general.fh"
#include "timers.fh"
#include "SysDef.fh"

  integer(kind=iwp) :: i_off1,isym,I,IBAS,ISA,J,iprlev
  real(kind=wp),allocatable :: ONTOPT(:),ONTOPO(:),FOne(:)
  real(kind=wp),allocatable :: FA_V(:),FI_V(:),TUVX(:)
  real(kind=wp),allocatable :: fa_t(:),Q(:)

  write(u6,'(2X,A)') 'Calculating potentials for analytic gradients for MS-PDFT'

  IPRLEV = IPRLOC(3)

! loading one-electron potential and two-electron potential
! Used as F1 and F2 in equations 58 and 59 in Ref1.

  Call mma_allocate(ONTOPT,nfint,Label='OnTopT')
  Call mma_allocate(ONTOPO,ntot1,Label='OnTopO')
  OnTopT(:) = zero
  OnTopO(:) = zero

  Call Get_dArray('ONTOPT',OnTopT,NFINT)
  Call Get_dArray('ONTOPO',OnTopO,NTOT1)

  If(IPRLEV >= DEBUG) THEN
    write(u6,*) 'One-electron potentials'
    do i = 1,ntot1
      write(u6,*) OnTopO(i)
    enddo
    write(u6,*) 'Two-electron potentials'
    DO i = 1,nfint
      if(abs(puvx(i)) >= 1d-10) then
        write(u6,*) OnTopT(i),puvx(i)
      else
        write(u6,*) OnTopT(i),0.0d0
      endif
    ENDDO
  ENDIF

  CALL mma_allocate(FI_V,Ntot1,Label='FI_V')
  Call Get_dArray('FI_V',FI_V,NTOT1)

  !Focka = fi_v+OnTopO
  Call daxpy_(ntot1,1.0d0,FI_V,1,Focka,1)
  Call daxpy_(ntot1,1.0d0,OnTopO,1,Focka,1)

  Call mma_allocate(Fone,NTOT1,Label='FOne')
  FOne(:) = zero

  i_off1 = 1
  DO iSym = 1,nSym
    iBas = nBas(iSym)
    !FI + FA + V_oe
    Do i = 1,iBas
      do j = 1,i
        Fone(i_off1) = Fone(i_off1)+FockA(i_off1)
        i_off1 = i_off1+1
      enddo
    EndDo
  ENDDO

  IF(IPRLEV >= DEBUG) then
    write(u6,*) 'F1 to send'
    DO i = 1,NTOT1
      write(u6,*) Fone(i)
    ENDDO
  ENDIF

  F1MS(:,iIntS) = Fone(:)
  Call mma_deallocate(FOne)

  Call mma_allocate(TUVX,NACPR2,Label='TUVX')
  TUVX(:) = zero
  Call Get_TUVX(OnTopT,TUVX)

  CALL DCopy_(NACPR2,TUVX,1,F2MS(:,iIntS),1)
  Call mma_deallocate(TUVX)

! ____________________________________________________________
! This next part is to generate the MC-PDFT generalized fock operator.

  focka(:ntot1) = zero
  focki(:ntot1) = zero

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
    Call dcopy_(ntot1,[0.0d0],0,FA_t,1)
    Call DaXpY_(NTOT1,1.0D0,OnTopO,1,FA_t,1)
    Call Daxpy_(NTOT1,1.0D0,FI_V,1,FA_t,1)
    Call Daxpy_(NTOT1,1.0D0,FA_V,1,FA_t,1)
    write(u6,*) "Total F additions:"
    Call TriPrt(' ','(5G18.10)',FA_T,norb(1))
    CALL mma_deallocate(FA_t)
  ENDIF

  FockI(:ntot1) = OnTopO+FI_V
  FockA(:ntot1) = FA_V

  IF(IPRLEV >= DEBUG) THEN
    write(u6,*) "new FI"
    Call TriPrt(' ','(5G18.10)',FockI,norb(1))
    write(u6,*) "new FA"
    Call TriPrt(' ','(5G18.10)',FockA,norb(1))
  ENDIF

  CALL mma_deallocate(FI_V)
  CALL mma_deallocate(FA_V)
  Call mma_deallocate(OnTopO)

!Reordering of the two-body density matrix.
  IF(ISTORP(NSYM+1) > 0) THEN
    CALL DCOPY_(ISTORP(NSYM+1),[0.0D0],0,P,1)
    CALL PMAT_RASSCF(P2d,P)
  ENDIF

!Must add to existing FOCK operator (occ/act). FOCK is not empty.
  CALL mma_allocate(Q,NQ,Label='Q') ! q-matrix(1symmblock)
  CALL FOCK_update(FOCK,FockI,FockA,D1Act,P,Q,OnTopT,CMO)
  Call mma_deallocate(Q)
  Call mma_deallocate(OnTopT)

  CALL DCopy_(nTot1,FockOcc,1,FocMS(:,iIntS),1)
  IF(IPRLEV >= DEBUG) THEN
    write(u6,*) 'FOCC_OCC'
    call wrtmat(FockOcc,1,ntot1,1,ntot1)
    write(u6,*) 'DONE WITH NEW FOCK OPERATOR'
  ENDIF

  iSA = 1
  Call Put_iScalar('SA ready',iSA)

EndSubroutine SaveFock_PDFT
