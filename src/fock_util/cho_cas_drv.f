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
      SUBROUTINE CHO_CAS_DRV(rc,W_CMO,DI,FI,DA1,FA,DA2,TraOnly)
      Use Data_Structures, only: DSBA_Type, Allocate_DSBA,
     &                           Deallocate_DSBA
      Implicit real*8 (a-h,o-z)

      Integer   rc
      Real*8    DA1(*),DI(*),DA2(*),FI(*),FA(*),W_CMO(*)
      Integer   nForb(8),nIorb(8),nAorb(8),nChM(8),nChI(8)
      Integer   nnA(8,8)
      Logical   TraOnly

#include "real.fh"
#include "chotodo.fh"
#include "chopmat.fh"
#include "chlcas.fh"
#include "cholk.fh"

      Character(LEN=11), Parameter:: SECNAM = 'CHO_CAS_DRV'

#include "rasdim.fh"
#include "wadr.fh"
#include "general.fh"
#include "rasscf.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Type (DSBA_Type) CVa(2), POrb(3), Ddec, ChoIn,
     &                 DLT(2), FLT(2), MSQ, FLT_MO(2)

      Real*8, Allocatable:: Tmp1(:), Tmp2(:)
      Real*8, Allocatable:: PMat(:), PL(:)
C ************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
C  **************************************************

      rc=0


      Call Allocate_DSBA(FLT(1),nBas,nBas,nSym,Case='TRI',Ref=FI)
      Call Allocate_DSBA(FLT(2),nBas,nBas,nSym,Case='TRI',Ref=FA)

      IF (TraOnly) THEN
*
*        Let us construct a second pointer structure based on norb

         Call Allocate_DSBA(FLT_MO(1),nOrb,nOrb,nSym,Case='TRI',Ref=FI)
         Call Allocate_DSBA(FLT_MO(2),nOrb,nOrb,nSym,Case='TRI',Ref=FA)
c
c ------ It only performs the MO transformation of FI and FA
c ----------------------------------------------------------
c
*        transform FI from AO to MO basis  (LT-storage)
      iOff2 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nOrb(iSym)
        iFro = nFro(iSym)
        Call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
        Call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp1')
        Call Square(FLT(1)%SB(iSym)%A1,Tmp1,1,iBas,iBas)
        Call DGEMM_('N','N',iBas,iOrb,iBas,
     &              1.0d0,Tmp1,iBas,
     &                    W_CMO(iOff2+(iFro*iBas)),max(iBas,iBas),
     &              0.0d0,Tmp2,iBas)
        Call MXMT(Tmp2,iBas,1,
     &            W_CMO(iOff2+(iFro*iBas)),1,iBas,
     &            FLT_MO(1)%SB(iSym)%A1,
     &            iOrb,iBas)
        Call mma_deallocate(Tmp2)
        Call mma_deallocate(Tmp1)
        iOff2 = iOff2 + iBas*iBas
      End Do

*     transform FA from AO to MO basis  (LT-storage)
      iOff2 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nOrb(iSym)
        iFro = nFro(iSym)
        Call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
        Call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp1')
        Call Square(FLT(2)%SB(iSym)%A1,Tmp1,1,iBas,iBas)
        Call DGEMM_('N','N',iBas,iOrb,iBas,
     &               1.0d0,Tmp1,iBas,
     &                     W_CMO(iOff2+(iFro*iBas)),max(iBas,iBas),
     &               0.0d0,Tmp2,iBas)

        Call MXMT(Tmp2,iBas,1,
     &            W_CMO(iOff2+(iFro*iBas)),1,iBas,
     &            FLT_MO(2)%SB(iSym)%A1,
     &            iOrb,iBas)
        Call mma_deallocate(Tmp2)
        Call mma_deallocate(Tmp1)
        iOff2 = iOff2 + iBas*iBas
      End Do

         Call deallocate_DSBA(FLT_MO(1))
         Call deallocate_DSBA(FLT_MO(2))

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
      Call Allocate_DSBA(DLT(1),nBas,nBas,nSym,Case='TRI')
      Call Allocate_DSBA(DLT(2),nBas,nBas,nSym,Case='TRI')

      Call Fold(nSym,nBas,DI, DLT(1)%A0)
      Call Fold(nSym,nBas,DA1,DLT(2)%A0)

      FactXI = -1.0D0

!AMS - should this be set differently for ExFac.ne.1?
!      FactXI = 0-ExFac

      If (Deco) Then

         FactXI = -0.5D0

!AMS - should this be set differently for ExFac.ne.1?
!         FactXI = 0-(ExFac*.5d0)

c --- decompose the Inactive density on request
         Call Allocate_DSBA(ChoIn,nBas,nBas,nSym)
         Call Allocate_DSBA(DDec,nBas,nBas,nSym)
         call dcopy_(NTot2,DI(1),1,DDec%A0,1)

         Call Allocate_DSBA(MSQ,nBas,nBas,nSym,Ref=ChoIn%A0)

         Thr = 1.0d-12
         incs=0
         Do i=1,nSym
            if((nForb(i)+nIorb(i)).gt.0)then
             CALL CD_InCore(DDec%SB(i)%A2,nBas(i),ChoIn%SB(i)%A2,
     &                      nBas(i),NumV,Thr,rc)
             If (rc.ne.0) Then
              write(6,*)SECNAM//': ill-defined dens decomp for Inact'
              write(6,*) 'rc value produced = ', rc
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
                  Ymax=Max(Ymax,DDec%SB(i)%A2(ja,ja))
               end do
               write(6,*)'Max diagonal of the density in symm. ',i,' is'
     &                   //' equal to ',Ymax
             endif
            else
             nChI(i) = 0
            endif
         End Do

         If (incs.gt.0 .and. DoLocK) then
            dmpk_old = dmpk
            dmpk = 1.0d-2*dmpk
            write(6,*)'LK-damping decreased from ',dmpk_old,' to ',dmpk
         EndIf

         Call Deallocate_DSBA(DDEc)

c --- to get the right input arguments for CHO_FCAS_AO and CHO_FMCSCF
         If(.not.DoLocK)Then
           Do i=1,nSym
              nForb(i) = 0
              nIorb(i) = nChI(i)
           End Do
         EndIf


      Else

         Call Allocate_DSBA(MSQ,nBas,nBas,nSym,Ref=W_CMO)

        Do i=1,nSym
           nChI(i) = nForb(i)+nIorb(i)
        End Do

      EndIf

      ipInc = ip_of_Work(MSQ%A0(1))

C --- Reordering of the MOs coefficients to fit cholesky needs
      If(.not.DoLocK)Then

        Call Allocate_DSBA(POrb(1),nChI,nBas,nSym)
        Call Allocate_DSBA(POrb(3),nAOrb,nBas,nSym)

        nOcs=0
        ioff1=0
        Do iSym=1,nSym

           do ikk=1,nChI(iSym)
              ioff2=ioff1+nBas(iSym)*(ikk-1)
              POrb(1)%SB(iSym)%A2(ikk,:) =
     &           Work(ipInc+ioff2 : ipInc+ioff2 -1 + nBas(iSym))
           end do

           ioff2=ioff1+nBas(iSym)*(nForb(iSym)+nIorb(iSym))
           do ikk=1,nAorb(iSym)
              POrb(3)%SB(iSym)%A2(ikk,:) =
     &             W_CMO( ioff2+nBas(iSym)*(ikk-1) + 1 :
     &                  ioff2+nBas(iSym)*(ikk-1) + nBas(iSym))
           end do
           ioff1=ioff1+nBas(iSym)**2
           nOcs = nOcs + nAorb(iSym)**2

        End Do

      Else

C *** Only the active orbitals MO coeff need reordering
           Call Allocate_DSBA(CVa(1),nAorb,nBas,nSym)

           ioff1 = 0
           Do iSym=1,nSym
            ioff2 = ioff1 + nBas(iSym)*(nForb(iSym)+nIorb(iSym))
            do ikk=1,nAorb(iSym)
               ioff = ioff2+nBas(iSym)*(ikk-1)
               CVa(1)%SB(iSym)%A2(ikk,:) =
     &           W_CMO(ioff +  1 : ioff + nBas(iSym))
            end do
            ioff1 = ioff1 + nBas(iSym)**2
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

         Call mma_allocate(PMat,nPMat,Label='PMat')
         ipPmat = ip_of_Work(PMat(1))
         Call mma_allocate(PL,NACPR2,Label='PL')
         ipPL = ip_of_Work(PL(1))
         Call CHO_Pmat(DA2,Pmat)

         PMat(:)=Zero

         Call Reord_Pmat(ipPL,ipPmat)

         Call mma_deallocate(PL)

      EndIf



      If (DoActive) Then
C ---  Decompose the active density  -----------------------------

#ifdef _DEBUGPRINT_
       koff=0
       do i=1,nSym
          CALL CD_TESTER(rc,DALT(1+koff),nBas(i),.true.)
          write(6,*) 'DALT for sym=', i
          CALL TRIPRT('DALT',' ',DALT(1+koff),nBas(i))
          koff = koff + nBas(i)*(nBas(i)+1)/2
       end do
#endif

        Call Allocate_DSBA(CVa(2),nBas,nBas,nSym)
        Call Allocate_DSBA(DDec,nBas,nBas,nSym)
        call dcopy_(NTot2,DA1(1),1,DDec%A0,1)

        Thr = 1.0d-12
        Do i=1,nSym
           if(nAorb(i).gt.0)then
! NOTE(Giovanni): CD will proceed with approx. decompos for QMC
!                 This will avoid warnings for negative-definit
             call CD_InCore(DDec%SB(i)%A2,nBas(i),
     &                      CVa(2)%SB(i)%A2,nBas(i),
     &                      NumV,Thr,rc)
             If (rc.ne.0) Then
                write(6,*)SECNAM//': ill-defined dens decomp for active'
                write(6,*) 'rc value produced = ', rc
                Call abend()
             EndIf
             nChM(i) = NumV
           else
             nChM(i) = 0
           endif
        End Do

        Call Deallocate_DSBA(DDec)

      Else

        ! Dummy allocation
        Call Allocate_DSBA(CVa(2),[1],[1],1)
        nChM(:) = 0

      EndIf

      If (.not.DoLocK .and. DoActive) Then

c --- reorder "Cholesky MOs" to Cva storage

        Call Allocate_DSBA(POrb(2),nChM,nBas,nSym)
        Do iSym=1,nSym
           If (nBas(iSym)*nChM(iSym).ne.0) Then
               do ikk=1,nChM(iSym)
                  POrb(2)%SB(iSym)%A2(ikk,:) =
     &               CVa(2)%SB(iSym)%A2(:,ikk)
               end do
           EndIf
        End Do

      Else

        Call Allocate_DSBA(POrb(2),[1],[1],1)

      EndIf
C ----------------------------------------------------------------
      FLT(1)%A0(:)=Zero
      FLT(2)%A0(:)=Zero

      ipInt = lpwxy   ! (PU|VX) integrals are computed
      ipCM  = ip_of_work(W_CMO(1))  ! MOs coeff. in C(a,p) storage

      IF (ALGO.eq.1 .and. .not. DoLocK) THEN

         CALL CHO_FMCSCF(rc,FLT,nForb,nIorb,nAorb,FactXI,
     &                   DLT,DoActive,POrb,nChM,ipInt,ExFac)

      ELSEIF (ALGO.eq.1 .and. DoLocK) THEN

         CALL CHO_LK_CASSCF(DLT,FLT,MSQ,ipInt,
     &                      FactXI,nChI,nAorb,nChM,CVa,DoActive,
     &                      nScreen,dmpK,abs(CBLBM),ExFac)

      ELSE

         write(6,*)SECNAM//': wrong input parameter. ALGO= ',ALGO
         rc=55
         Return

      ENDIF

      Call Deallocate_DSBA(POrb(3))
      Call Deallocate_DSBA(POrb(2))
      Call Deallocate_DSBA(POrb(1))
      Call Deallocate_DSBA(CVa(1))
      Call Deallocate_DSBA(CVa(2))

      If (DoQmat.and.ALGO.ne.1) Call mma_deallocate(PMat)
      If (Deco) Call Deallocate_DSBA(ChoIn)

      Call Deallocate_DSBA(DLT(2))
      Call Deallocate_DSBA(DLT(1))

      Call deallocate_DSBA(MSQ)

      ENDIF

      Call deallocate_DSBA(FLT(2))
      Call deallocate_DSBA(FLT(1))

      Return
      END
