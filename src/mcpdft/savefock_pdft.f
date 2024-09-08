************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2018, Andrew M. Sand                                   *
*               2019, Thais R. Scott                                   *
*               2021, Jie J. Bao                                       *
************************************************************************
      Subroutine SaveFock_PDFT(CMO,FockI,FockA,iD1Act,LFock,
     &                         LP,NQ,LQ,LPUVX,ip2d,istate)
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Jan. 04, 2021, created this file.               *
* ****************************************************************

      use mspdft, only: F1MS, F2MS, FocMS, iIntS
      use printlevel, only: debug
      use mcpdft_output, only: lf, iPrLoc
      use rctfld_module
      use stdalloc, only: mma_allocate, mma_deallocate
      use wadr, only: BM, FockOcc

* Notes: Two references will be referred to in the comments.
* Ref1:  Sand, et al. JCTC, 2018, 14,  126.
* Ref2: Scott, et al. JCP,  2020, 153, 014106.
      Implicit Real*8 (A-H,O-Z)
      Real*8 CMO(*), FockI(*), FockA(*)
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
#include "pamint.fh"
#include "timers.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"

      INTEGER iD1Act,LP,NQ,LQ,LPUVX,ip2d,istate,LFock


******Auxiliary Variables
      INTEGER i_off1,isym
      INTEGER ittTUVX,ifat
      INTEGER iPrLev
      External isfreeunit
      CHARACTER(len=64) FILENAME
      CHARACTER(len=8) STATENAME
      Real*8, Allocatable:: ONTOPT(:), ONTOPO(:), FOne(:)
      Real*8, Allocatable:: FA_V(:), FI_V(:)


      write(lf,'(2X,A)')
     &'Calculating potentials for analytic gradients for MS-PDFT'

      IPRLEV=IPRLOC(3)

      write(StateName,'(I3)') ISTATE
      write (FILENAME,fmt='(a,a)')
     &"TmpFock",trim(adjustl(STATENAME))


******ioading one-electron potential and two-electron potential
******Used as F1 and F2 in equations 58 and 59 in Ref1.


      Call mma_allocate(ONTOPT,nfint,Label='OnTopT')
      Call mma_allocate(ONTOPO,ntot1,Label='OnTopO')
      OnTopT(:)=0.0D0
      OnTopO(:)=0.0D0

      Call Get_dArray('ONTOPT',OnTopT,NFINT)
      Call Get_dArray('ONTOPO',OnTopO,NTOT1)


      If (IPRLEV.ge.DEBUG ) THEN
        write(lf,*) 'One-electron potentials'
        do i=1,ntot1
          write(lf,*) OnTopO(i)
        end do
        write(lf,*) 'Two-electron potentials'
        DO i=1,nfint
         if (abs(work(lpuvx-1+i)).ge.1d-10)then
           write(lf,*) OnTopT(i),work(lpuvx-1+i)
         else
           write(lf,*) OnTopT(i),0.0d0
         end if
        END DO
       END IF

      Call mma_allocate(Fone,NTOT1,Label='FOne')
      FOne(:)=0.0D0

      CALL mma_allocate(FI_V,Ntot1,Label='FI_V')
      Call Get_dArray('FI_V',FI_V,NTOT1)

*     Focka=fi_v+OnTopO
      Call daxpy_(ntot1,1.0d0,FI_V,1,Focka,1)
      Call daxpy_(ntot1,1.0d0,OnTopO,1,Focka,1)

      i_off1=1

*     F1=fiv+tmploeotp
      DO iSym = 1,nSym
       iBas = nBas(iSym)
       !FI + FA + V_oe
       Do i=1,iBas
        do j=1,i
         Fone(i_off1) = Fone(i_off1) + FockA(i_off1)
         i_off1 = i_off1 + 1
        end do
       End Do
      END DO

      IF ( IPRLEV.ge.DEBUG ) then
       write(lf,*) 'F1 to send'
       DO i=1,NTOT1
         write(lf,*) Fone(i)
       END DO
      END IF

      CALL DCopy_(nTot1,Fone,1,F1MS(:,iIntS),1)
      Call GetMem('ttTUVX','Allo','Real',ittTUVX,NACPR2)
      CALL DCOPY_(nacpr2,[0.0D0],0,WORK(ittTUVX),1)
      Call Get_TUVX(OnTopT,Work(ittTUVX))

      CALL DCopy_(NACPR2,Work(ittTUVX),1,F2MS(:,iIntS),1)
      Call mma_deallocate(FOne)
      Call GetMem('ttTUVX','Free','Real',ittTUVX,NACPR2)

!____________________________________________________________
!This next part is to generate the MC-PDFT generalized fock operator.

      CALL DCOPY_(ntot1,[0.0D0],0,FockA,1)
      CALL DCOPY_(ntot1,[0.0D0],0,Focki,1)

!The corrections (from the potentials) to FI and FA are built in the NQ
!part of the code, for efficiency's sake.  It still needs to be
!debugged.
      CALL mma_allocate(FA_V,Ntot1,Label='FA_V')
      Call Get_dArray('FA_V',FA_V,NTOT1)

      IF ( IPRLEV.GE.DEBUG ) THEN
       write(lf,*) "extra terms to update FI"
       DO i=1,ntot1
        write(lf,*) FI_V(i)
       END DO
       write(lf,*) "extra terms to update FA"
       DO i=1,ntot1
        write(lf,*) FA_V(i)
       END DO
       CALL GETMEM('FA_t','ALLO','REAL',ifat,Ntot1)
       Call dcopy_(ntot1,[0.0d0],0,work(ifat),1)
       Call DaXpY_(NTOT1,1.0D0,OnTopO,1,Work(ifat),1)
       Call Daxpy_(NTOT1,1.0D0,FI_V,1,Work(ifat),1)
       Call Daxpy_(NTOT1,1.0D0,FA_V,1,Work(ifat),1)
       write(lf,*) "Total F additions:"
       Call TriPrt(' ','(5G18.10)',Work(ifat),norb(1))
       CALL GETMEM('FA_t','free','REAL',ifat,Ntot1)
      END IF

      Call DaXpY_(NTOT1,1.0D0,OnTopO,1,FockI,1)
      Call Daxpy_(NTOT1,1.0D0,FI_V,1,FockI,1)
      Call Daxpy_(NTOT1,1.0D0,FA_V,1,FockA,1)

      IF ( IPRLEV.GE.DEBUG ) THEN
       write(lf,*) "new FI"
       Call TriPrt(' ','(5G18.10)',FockI,norb(1))
       write(lf,*) "new FA"
       Call TriPrt(' ','(5G18.10)',FockA,norb(1))
      END IF

      CALL mma_deallocate(FI_V)
      CALL mma_deallocate(FA_V)

!Reordering of the two-body density matrix.
      IF(ISTORP(NSYM+1).GT.0) THEN
       CALL DCOPY_(ISTORP(NSYM+1),[0.0D0],0,WORK(LP),1)
       CALL PMAT_RASSCF(Work(iP2d),WORK(LP))
      END IF
!Must add to existing FOCK operator (occ/act). FOCK is not empty.
      CALL mma_allocate(BM,NSXS,Label='BM')
      CALL GETMEM('SXLQ','ALLO','REAL',LQ,NQ) ! q-matrix(1symmblock)
      CALL FOCK_update(WORK(LFOCK),BM,FockI,
     &     FockA,Work(iD1Act),WORK(LP),
     &     WORK(LQ),OnTopT,IFINAL,CMO)

      CALL DCopy_(nTot1,FockOcc,1,FocMS(:,iIntS),1)
      IF ( IPRLEV.GE.DEBUG ) THEN
       write(lf,*) 'FOCC_OCC'
       call wrtmat(FockOcc,1,ntot1,1,ntot1)
       write(lf,*) 'DONE WITH NEW FOCK OPERATOR'
      END IF

      Call mma_deallocate(BM)
      CALL GETMEM('SXLQ','Free','REAL',LQ,NQ) ! q-matrix(1symmblock)
      Call mma_deallocate(OnTopO)
      Call mma_deallocate(OnTopT)

      iSA = 1
      Call Put_iScalar('SA ready',iSA)

      End Subroutine SaveFock_PDFT
