************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE CHO_CAS_DRV(rc,CMO,DI,FI,DA1,FA,DA2,TraOnly)
      Implicit real*8 (a-h,o-z)

      Integer   rc
      Integer   ISTAQ(8),Nscreen
      Real*8    DA1(*),DI(*),DA2(*),FI(*),FA(*),CMO(*)
      Integer   nForb(8),nIorb(8),nAorb(8),nChM(8),nChI(8)
      Integer   ipDSA2(8,8,8),nnA(8,8),ipAorb(2),ipKLT(2)
      Logical   DoActive,DoQmat,TraOnly,DoLocK,Deco,DoCholesky
      Integer   ALGO
      Real*8    dmpk
      Logical   DeAllocte_CVA

      COMMON /CHOTODO /DoActive,DoQmat,ipQmat
      COMMON /CHOPMAT / ipPL
      Common /CHLCAS / DoCholesky,ALGO
      Common /CHOLK / DoLocK,Deco,dmpk,Nscreen

      Character*11 SECNAM
      Parameter (SECNAM = 'CHO_CAS_DRV')

#include "rasdim.fh"
#include "wadr.fh"
#include "general.fh"
#include "rasscf.fh"
#include "WrkSpc.fh"
C ************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
C  **************************************************


      rc=0

      DeAllocte_CVA=.False.

      IF (TraOnly) THEN
c
c --- It only performs the MO transformation of FI and FA
c -------------------------------------------------------
c
*     transform FI from AO to MO basis  (LT-storage)
      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nOrb(iSym)
        iFro = nFro(iSym)
        Call GetMem('Scr1','Allo','Real',iTmp1,iBas*iBas)
        Call GetMem('Scr2','Allo','Real',iTmp2,iOrb*iBas)
        Call Square(FI(iOff1),Work(iTmp1),1,iBas,iBas)
C        Call MXMA(Work(iTmp1),1,iBas,
C     &            CMO(iOff2+(iFro*iBas)),1,iBas,
C     &            Work(iTmp2),1,iBas,
C     &            iBas,iBas,iOrb)
        Call DGEMM_('N','N',iBas,iOrb,iBas,1.0d0,Work(iTmp1),
     &               iBas,CMO(iOff2+(iFro*iBas)),max(iBas,iBas),
     &               0.0d0,Work(iTmp2),iBas)
        Call MXMT(Work(iTmp2),iBas,1,
     &            CMO(iOff2+(iFro*iBas)),1,iBas,
     &            FI(iOff3),
     &            iOrb,iBas)
        Call GetMem('Scr2','Free','Real',iTmp2,iOrb*iBas)
        Call GetMem('Scr1','Free','Real',iTmp1,iBas*iBas)
        iOff1 = iOff1 + (iBas*iBas+iBas)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
      End Do

*     transform FA from AO to MO basis  (LT-storage)
      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nOrb(iSym)
        iFro = nFro(iSym)
        Call GetMem('Scr1','Allo','Real',iTmp1,iBas*iBas)
        Call GetMem('Scr2','Allo','Real',iTmp2,iOrb*iBas)
        Call Square(FA(iOff1),Work(iTmp1),1,iBas,iBas)
        Call DGEMM_('N','N',iBas,iOrb,iBas,1.0d0,Work(iTmp1),
     &               iBas,CMO(iOff2+(iFro*iBas)),max(iBas,iBas),
     &               0.0d0,Work(iTmp2),iBas)

        Call MXMT(Work(iTmp2),iBas,1,
     &            CMO(iOff2+(iFro*iBas)),1,iBas,
     &            FA(iOff3),
     &            iOrb,iBas)
        Call GetMem('Scr2','Free','Real',iTmp2,iOrb*iBas)
        Call GetMem('Scr1','Free','Real',iTmp1,iBas*iBas)
        iOff1 = iOff1 + (iBas*iBas+iBas)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
      End Do

c**************************************************************************


      ELSE

c --- It only computes FI and FA  in AO-basis and returns
c --- the active integrals (tw|xy)
c --- If specified in input, the routine also computes
c --- the auxiliary Q-matrix stored as Q(av), where a is an AO index
c --- and v refers to the active orbitals only

      Do iSym=1,nSym
         nForb(iSym) = nFro(iSym)
         nIorb(iSym) = nIsh(iSym)
         nAorb(iSym) = nAsh(iSym)
      End Do


C --- Build the packed densities from the Squared ones
      Call Getmem('DILT','Allo','Real',ipDILT,NTot1)
      Call Getmem('DALT','Allo','Real',ipDALT,NTot1)

      Call Fold(nSym,nBas,DI,Work(ipDILT))
      Call Fold(nSym,nBas,DA1,Work(ipDALT))

      FactXI = -1.0D0

!AMS - should this be set differently for ExFac.ne.1?
!      FactXI = 0-ExFac

      If (Deco) Then

         FactXI = -0.5D0

!AMS - should this be set differently for ExFac.ne.1?
!         FactXI = 0-(ExFac*.5d0)

c --- decompose the Inactive density on request
         CALL GETMEM('choIn','allo','real',ipIna,NTot2)
         CALL GETMEM('ddec','allo','real',ipddec,NTot2)
         call dcopy_(NTot2,DI(1),1,Work(ipddec),1)

         ipInc = ipIna

         Thr = 1.0d-12
         ipd = ipddec
         ipV = ipIna
         incs=0
         Do i=1,nSym
            if((nForb(i)+nIorb(i)).gt.0)then
             CALL CD_InCore(Work(ipd),nBas(i),Work(ipV),nBas(i),
     &                     NumV,Thr,rc)
             If (rc.ne.0) Then
              write(6,*)SECNAM//': ill-defined dens decomp for Inact'
              write(6,*) 'rc value produced = ', rc
              Call qtrace()
              Call abend()
             EndIf
             nChI(i) = NumV
             if ( NumV .ne. nIsh(i)+nForb(i) ) then
               write(6,*)'Warning! The number of occupied from the deco'
     &                 //'mposition of the Inactive density matrix is ',
     &                   numV,' in symm. ',i
               write(6,*)'Expected value = ',nIsh(i)+nForb(i)
               incs=incs+1
               Ymax=0.0d0
               do ja=1,nBas(i)
                  jaa=ipd-1+nBas(i)*(ja-1)+ja
                  Ymax=Max(Ymax,Work(jaa))
               end do
               write(6,*)'Max diagonal of the density in symm. ',i,' is'
     &                   //' equal to ',Ymax
             endif
            else
             nChI(i) = 0
            endif
            ipd = ipd + nBas(i)**2
            ipV = ipV + nBas(i)**2
         End Do

         If (incs.gt.0 .and. DoLocK) then
            dmpk_old = dmpk
            dmpk = 1.0d-2*dmpk
            write(6,*)'LK-damping decreased from ',dmpk_old,' to ',dmpk
         EndIf

         CALL GETMEM('ddec','free','real',ipddec,NTot2)

c --- to get the right input arguments for CHO_FCAS_AO and CHO_FMCSCF
         If(.not.DoLocK)Then
           Do i=1,nSym
              nForb(i) = 0
              nIorb(i) = nChI(i)
           End Do
         EndIf


      Else

        ipInc = ip_of_Work(CMO(1))

        Do i=1,nSym
           nChI(i) = nForb(i)+nIorb(i)
        End Do

      EndIf

c --- Various offsets
c --------------------
        ISTAQ(1)=0
      DO ISYM=2,NSYM
        NB=NBAS(ISYM-1)
        NP=NCHI(ISYM-1)+NAORB(ISYM-1)
        NP2=NB*NP
        ISTAQ(ISYM)=ISTAQ(ISYM-1)+NP2 ! Reordered MOs coeff
      END DO

C --- Reordering of the MOs coefficients to fit cholesky needs
      If(.not.DoLocK)Then

        NTaq = ISTAQ(nSym) + nBas(nSym)*(nChI(nSym)+nAorb(nSym))

        Call Getmem('rMOs','Allo','Real',ipPorb,NTaq)

        nOcs=0
        ioff1=0
        Do iSym=1,nSym
           do ikk=1,nChI(iSym)
              ioff2=ioff1+nBas(iSym)*(ikk-1)
              call dcopy_(nBas(iSym),Work(ipInc+ioff2),1,
     &        Work(ipPorb+ISTAQ(iSym)+ikk-1),nChI(iSym))
           end do
           ioff2=ioff1+nBas(iSym)*(nForb(iSym)+nIorb(iSym))
           do ikk=1,nAorb(iSym)
              call dcopy_(nBas(iSym),CMO(1+ioff2+nBas(iSym)*(ikk-1)),1,
     &        Work(ipPorb+ISTAQ(iSym)+nChI(iSym)*nBas(iSym)+ikk-1),
     &             nAorb(iSym))
           end do
           ioff1=ioff1+nBas(iSym)**2
c           call recprt('ReorIMOs','',Work(ipPorb+ISTAQ(iSym)),
c     &                  nChI(iSym),nBas(iSym))
c           call recprt('ReorAMOs','',Work(ipPorb+ISTAQ(iSym)+nChI(iSym)*
c     &                  nBas(iSym)),
c     &                  nAorb(iSym),nBas(iSym))

            nOcs = nOcs + nAorb(iSym)**2

        End Do

      Else

C *** Only the active orbitals MO coeff need reordering
           nVB=0
           Do iSym=1,nSym
            nVB = nVB + nAorb(iSym)*nBas(iSym)
           End Do
           Call GetMem('Cva','Allo','Real',ipAorb(1),nVB)
           DeAllocte_CVA=.True.

           ioff1 = 0
           ioff3 = 0
           Do iSym=1,nSym
            ioff2 = ioff1 + nBas(iSym)*(nForb(iSym)+nIorb(iSym))
            do ikk=1,nAorb(iSym)
               call dcopy_(nBas(iSym),CMO(1+ioff2+nBas(iSym)*(ikk-1)),1,
     &                 Work(ipAorb(1)+ioff3+ikk-1),nAorb(iSym))
            end do
            ioff1 = ioff1 + nBas(iSym)**2
            ioff3 = ioff3 + nAorb(iSym)*nBas(iSym)
           End Do

      EndIf

C --- Optional Section for Q-matrix evaluation
C --- Reorder the 2-el density matrix to fit cholesky needs
      If (DoQmat.and.ALGO.ne.1) Then

         Call set_nnA(nSym,nAorb,nnA)

         nPmat=0   ! P[vw],xy
         Do iSymXY=1,nSym
            Do iSymy=1,nSym
               iSymx=MulD2h(iSymXY,iSymy)
               if (iSymx.le.iSymy) then
               Do iSymw=1,nSym
                  iSymv=MulD2h(iSymXY,iSymw)
                  nPmat = nPmat
     &                  + nAorb(iSymv)*nAorb(iSymw)*nnA(iSymx,iSymy)
               End Do
               endif
            End do
         End Do

         Call Getmem('P-mat','Allo','Real',ipPmat,nPmat)
         Call Getmem('P-scal','Allo','Real',ipPL,NACPR2)
         Call CHO_Pmat(DA2,Work(ipPmat))

         Call Fzero(Work(ipPmat),nPmat)

c         ipDA2 = ip_of_Work(DA2(1))
         ipDA2 = ipPL

         Call Reord_Pmat(ipDA2,ipPmat,ipDSA2)

         Call Getmem('P-scal','Free','Real',ipPL,NACPR2)

      EndIf



      If (DoActive) Then
C ---  Decompose the active density  -----------------------------

#ifdef _DEBUGPRINT_
       koff=0
       do i=1,nSym
          CALL CD_TESTER(rc,ipDALT+koff,nBas(i),.true.)
          write(6,*) 'DALT for sym=', i
          CALL TRIPRT('DALT',' ',Work(ipDALT+koff),nBas(i))
          koff = koff + nBas(i)*(nBas(i)+1)/2
       end do
#endif

        CALL GETMEM('choMOs','allo','real',ipVec,NTot2) !cholesk MOs
        CALL GETMEM('ddec','allo','real',ipddec,NTot2)
        call dcopy_(NTot2,DA1(1),1,Work(ipddec),1)

        Thr = 1.0d-12
        ipd = ipddec
        ipV = ipVec
        nTvec = 0
        Do i=1,nSym
           if(nAorb(i).gt.0)then
! NOTE(Giovanni): CD will proceed with approx. decompos for QMC
!                 This will avoid warnings for negative-definit
             call CD_InCore(Work(ipd),nBas(i),Work(ipV),nBas(i),
     &                      NumV,Thr,rc)
             If (rc.ne.0) Then
                write(6,*)SECNAM//': ill-defined dens decomp for active'
                write(6,*) 'rc value produced = ', rc
                Call qtrace()
                Call abend()
             EndIf
             nChM(i) = NumV
           else
             nChM(i) = 0
           endif
           ipd = ipd + nBas(i)**2
           ipV = ipV + nBas(i)**2
           nTvec = nTvec + nChM(i)*nBas(i)
        End Do

        CALL GETMEM('ddec','free','real',ipddec,NTot2)

      Else

        ipVec = -999999999  ! avoid compiler warnings

        Do i=1,nSym
           nChM(i) = 0
        End Do

      EndIf

      If (.not.DoLocK .and. DoActive) Then
c --- reorder "Cholesky MOs" to Cva storage

        CALL GETMEM('chM','allo','real',ipChM,nTvec) !cholesky reord MOs
        ioff1=0
        ioff2=0
         Do iSym=1,nSym
           If (nBas(iSym)*nChM(iSym).ne.0) Then
               do ikk=1,nChM(iSym)
                  koff = ioff1 + nBas(iSym)*(ikk-1)
                  call dcopy_(nBas(iSym),Work(ipVec+koff),1,
     &                 Work(ipChM+ioff2+ikk-1),nChM(iSym))
               end do
           EndIf
           ioff1=ioff1+nBas(iSym)**2
           ioff2=ioff2+nChM(iSym)*nBas(iSym)
         End Do

      Else

         ipChM = -9999999  ! avoid compiler warnings

      EndIf
C ----------------------------------------------------------------

      Call Fzero(FI(1),nTot1) ! LT-storage
      Call Fzero(FA(1),nTot1) ! LT-storage

      ipFI = ip_of_Work(FI(1))
      ipFA = ip_of_Work(FA(1))


      IF (ALGO.eq.1 .and. .not. DoLocK) THEN

         ipInt = lpwxy   ! (PU|VX) integrals are computed
         ipCM = ip_of_work(CMO(1))  ! MOs coeff. in C(a,p) storage

         CALL CHO_FMCSCF(rc,ipFA,ipFI,nForb,nIorb,nAorb,FactXI,
     &        ipPorb,ipDILT,ipDALT,DoActive,ipChM,nChM,ipInt,ExFac)


      ELSEIF (ALGO.eq.1 .and. DoLocK) THEN

         ipInt = lpwxy   ! (PU|VX) integrals are computed
         ipCM = ip_of_work(CMO(1))  ! MOs coeff. in C(a,p) storage
         ipAorb(2) = ipVec  ! decomposed active density

         Call Getmem('KILT','Allo','Real',ipKLT(1),NTot1)
         Call Fzero(Work(ipKLT(1)),NTot1)

         If (DoActive) Then
           Call Getmem('KALT','Allo','Real',ipKLT(2),NTot1)
           Call Fzero(Work(ipKLT(2)),NTot1)
         EndIf

         CALL CHO_LK_CASSCF(ipDILT,ipDALT,ipFI,ipFA,ipKLT,ipInc,ipInt,
     &                      FactXI,nChI,nAorb,nChM,ipAorb,DoActive,
     &                      nScreen,dmpK,abs(CBLBM),ExFac)

         If(DoActive) Call Getmem('KALT','Free','Real',ipKLT(2),NTot1)
         Call Getmem('KILT','Free','Real',ipKLT(1),NTot1)


      ELSEIF (ALGO.eq.2) THEN

         ipInt = LTUVX   ! (TU|VX) integrals only are computed

         CALL CHO_FCAS_AO(rc,ipFA,ipFI,ipQmat,nForb,nIorb,nAorb,FactXI,
     &    ipPorb,ipDILT,ipDALT,ipDSA2,DoActive,DoQmat,ipChM,nChM,ipInt,
     &    ExFac)

*  Synchronization of the Fock matrices
         Call GaDsum(Work(ipFI),NTot1)
         Call GaDsum(Work(ipFA),NTot1)
*  Synchronization of the Q-matrix
         lQ=0
         Do i=1,nSym
            lQ = lQ + nBas(i)*nAorb(i)
         End Do
         Call GaDsum(Work(ipQmat),lQ)
*  Synchronization of the (TU|VX) integrals
         Call GaDsum(Work(LTUVX),NACPR2)


      ELSE

         write(6,*)SECNAM//': wrong input parameter. ALGO= ',ALGO
         rc=55
         Return

      ENDIF


      If(.not.DoLocK .and. DoActive)Then

        CALL GETMEM('chM','free','real',ipChM,nTvec) !cholesky reord MOs

      EndIf

      If (DeAllocte_CVA) Call GetMem('Cva','free','Real',ipAorb(1),nVB)


      If(DoActive) CALL GETMEM('choMOs','free','real',ipVec,NTot2)

      If (DoQmat.and.ALGO.ne.1) Then
         Call Getmem('P-mat','Free','Real',ipPmat,NPmat)
      EndIf

      If(.not.DoLocK)Then
        Call Getmem('rMOs','Free','Real',ipPorb,NTaq)
      EndIf

      If (Deco) CALL GETMEM('choIn','free','real',ipIna,NTot2)

      Call Getmem('DALT','Free','Real',ipDALT,NTot1)
      Call Getmem('DILT','Free','Real',ipDILT,NTot1)



      ENDIF


      Return
      END
