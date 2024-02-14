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
      Subroutine SaveFock_PDFT(CMO,IFockI,IFockA,iD1Act,LFock,
     &                         LP,NQ,LQ,LPUVX,ip2d,istate)
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Jan. 04, 2021, created this file.               *
* ****************************************************************

      Use KSDFT_Info, Only: ifav, ifiv
      use mspdft, only: iF1MS, iF2MS, iFocMS, iIntS
      use mcpdft_output, only: debug, lf, iPrLoc
      use rctfld_module
      use stdalloc, only: mma_allocate, mma_deallocate
      use wadr, only: BM, FockOcc

* Notes: Two references will be referred to in the comments.
* Ref1:  Sand, et al. JCTC, 2018, 14,  126.
* Ref2: Scott, et al. JCP,  2020, 153, 014106.
      Implicit Real*8 (A-H,O-Z)
      Real*8 CMO(*)
#include "rasdim.fh"
#include "general.fh"
#include "input_ras_mcpdft.fh"
#include "rasscf.fh"
#include "WrkSpc.fh"
#include "pamint.fh"
#include "timers.fh"
#include "SysDef.fh"
!      Logical TraOnly



      INTEGER IFockI,IFockA,iD1Act,LP,NQ,LQ,LPUVX,ip2d,istate,LFock


******Auxiliary Variables
      INTEGER i_off1,ifone,isym,ipTmpLTEOTP,ipTmpLOEOTP
      INTEGER ittTUVX,ifat
      INTEGER iPrLev
      External isfreeunit
      CHARACTER(len=64) FILENAME
      CHARACTER(len=8) STATENAME


      write(lf,'(2X,A)')
     &'Calculating potentials for analytic gradients for MS-PDFT'

      IPRLEV=IPRLOC(3)

      write(StateName,'(I3)') ISTATE
      write (FILENAME,fmt='(a,a)')
     &"TmpFock",trim(adjustl(STATENAME))


******ioading one-electron potential and two-electron potential
******Used as F1 and F2 in equations 58 and 59 in Ref1.


      Call GetMem('ONTOPT','ALLO','Real',ipTmpLTEOTP,nfint)
      Call GetMem('ONTOPO','ALLO','Real',ipTmpLOEOTP,ntot1)
      Call FZero(Work(iptmplteotp),Nfint)
      Call FZero(Work(iptmploeotp),ntot1)

      Call Get_dArray('ONTOPT',work(ipTmpLTEOTP),NFINT)
      Call Get_dArray('ONTOPO',work(ipTmpLOEOTP),NTOT1)


      If (IPRLEV.ge.DEBUG ) THEN
        write(lf,*) 'One-electron potentials'
        do i=1,ntot1
          write(lf,*) Work(iptmploeotp-1+i)
        end do
        write(lf,*) 'Two-electron potentials'
        DO i=1,nfint
         if (abs(work(lpuvx-1+i)).ge.1d-10)then
           write(lf,*) Work(iptmplteotp-1+i),work(lpuvx-1+i)
         else
           write(lf,*)Work(iptmplteotp-1+i),0.0d0
         end if
        END DO
       END IF

      Call GetMem('F_ONE','ALLO','Real',iFone,NTOT1)
      CALL DCOPY_(NTOT1,[0.0D0],0,WORK(iFone),1)

      CALL GETMEM('FI_V','ALLO','REAL',ifiv,Ntot1)
      Call Get_dArray('FI_V',work(ifiv),NTOT1)

*     Focka=fiv+tmploeotp
      Call daxpy_(ntot1,1.0d0,Work(ifiv),1,Work(iFocka),1)
      Call daxpy_(ntot1,1.0d0,Work(iptmploeotp),1,Work(iFocka),1)

      i_off1=0

*     F1=fiv+tmploeotp
      DO iSym = 1,nSym
       iBas = nBas(iSym)
       !FI + FA + V_oe
       Do i=1,iBas
        do j=1,i
         Work(iFone+i_off1) = Work(ifone+i_off1) +
     &   Work(ifocka+i_off1)
         i_off1 = i_off1 + 1
        end do
       End Do
      END DO

      IF ( IPRLEV.ge.DEBUG ) then
       write(lf,*) 'F1 to send'
       DO i=1,NTOT1
         write(lf,*) work(iFone-1+i)
       END DO
      END IF

      CALL DCopy_(nTot1,Work(iFone),1,WORK(iF1MS+(iIntS-1)*nTot1),1)
      Call GetMem('ttTUVX','Allo','Real',ittTUVX,NACPR2)
      CALL DCOPY_(nacpr2,[0.0D0],0,WORK(ittTUVX),1)
      Call Get_TUVX(Work(ipTmpLTEOTP),Work(ittTUVX))

      CALL DCopy_(NACPR2,Work(ittTUVX),1,WORK(iF2MS+(iIntS-1)*NACPR2),1)
      Call GetMem('F_ONE','Free','Real',iFone,NTOT1)
      Call GetMem('ttTUVX','Free','Real',ittTUVX,NACPR2)

!____________________________________________________________
!This next part is to generate the MC-PDFT generalized fock operator.

      CALL DCOPY_(ntot1,[0.0D0],0,WORK(iFocka),1)
      CALL DCOPY_(ntot1,[0.0D0],0,WORK(iFocki),1)

!The corrections (from the potentials) to FI and FA are built in the NQ
!part of the code, for efficiency's sake.  It still needs to be
!debugged.
      CALL GETMEM('FA_V','ALLO','REAL',ifav,Ntot1)
      Call Get_dArray('FA_V',work(ifav),NTOT1)

      IF ( IPRLEV.GE.DEBUG ) THEN
       write(lf,*) "extra terms to update FI"
       DO i=1,ntot1
        write(lf,*) Work(ifiv-1+i)
       END DO
       write(lf,*) "extra terms to update FA"
       DO i=1,ntot1
        write(lf,*) Work(ifav-1+i)
       END DO
       CALL GETMEM('FA_t','ALLO','REAL',ifat,Ntot1)
       Call dcopy_(ntot1,[0.0d0],0,work(ifat),1)
       Call DaXpY_(NTOT1,1.0D0,Work(ipTmpLOEOTP),1,Work(ifat),1)
       Call Daxpy_(NTOT1,1.0D0,Work(ifiv),1,Work(ifat),1)
       Call Daxpy_(NTOT1,1.0D0,Work(ifav),1,Work(ifat),1)
       write(lf,*) "Total F additions:"
       Call TriPrt(' ','(5G18.10)',Work(ifat),norb(1))
       CALL GETMEM('FA_t','free','REAL',ifat,Ntot1)
      END IF

      Call DaXpY_(NTOT1,1.0D0,Work(ipTmpLOEOTP),1,Work(ifocki),1)
      Call Daxpy_(NTOT1,1.0D0,Work(ifiv),1,Work(ifocki),1)
      Call Daxpy_(NTOT1,1.0D0,Work(ifav),1,Work(ifocka),1)

      IF ( IPRLEV.GE.DEBUG ) THEN
       write(lf,*) "new FI"
       Call TriPrt(' ','(5G18.10)',Work(ifocki),norb(1))
       write(lf,*) "new FA"
       Call TriPrt(' ','(5G18.10)',Work(ifocka),norb(1))
      END IF

      CALL GETMEM('FI_V','Free','REAL',ifiv,Ntot1)
      CALL GETMEM('FA_V','Free','REAL',ifav,Ntot1)

!Reordering of the two-body density matrix.
      IF(ISTORP(NSYM+1).GT.0) THEN
       CALL DCOPY_(ISTORP(NSYM+1),[0.0D0],0,WORK(LP),1)
       CALL PMAT_RASSCF(Work(iP2d),WORK(LP))
      END IF
!Must add to existing FOCK operator (occ/act). FOCK is not empty.
      CALL mma_allocate(BM,NSXS,Label='BM')
      CALL GETMEM('SXLQ','ALLO','REAL',LQ,NQ) ! q-matrix(1symmblock)
      CALL FOCK_update(WORK(LFOCK),BM,Work(iFockI),
     &     Work(iFockA),Work(iD1Act),WORK(LP),
     &     WORK(LQ),WORK(ipTmpLTEOTP),IFINAL,CMO)

      CALL DCopy_(nTot1,FockOcc,1,WORK(iFocMS+(iIntS-1)*nTot1),1)
      IF ( IPRLEV.GE.DEBUG ) THEN
       write(lf,*) 'FOCC_OCC'
       call wrtmat(FockOcc,1,ntot1,1,ntot1)
       write(lf,*) 'DONE WITH NEW FOCK OPERATOR'
      END IF

      Call mma_deallocate(BM)
      CALL GETMEM('SXLQ','Free','REAL',LQ,NQ) ! q-matrix(1symmblock)
      Call GetMem('ONTOPO','FREE','Real',ipTmpLOEOTP,ntot1)
      Call GetMem('ONTOPT','FREE','Real',ipTmpLTEOTP,nfint)

      iSA = 1
      Call Put_iScalar('SA ready',iSA)

      RETURN
      End Subroutine
