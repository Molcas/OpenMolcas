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

      Type (DSBA_Type) CVa(2), POrb(3), Ddec, ChoIn, CMO,
     &                 DLT(2), FLT(2), MSQ, FLT_MO(2)

      Real*8, Allocatable:: Tmp1(:), Tmp2(:)
*                                                                      *
************************************************************************
*                                                                      *
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
*                                                                      *
************************************************************************
*                                                                      *
      tmp=DA2(1)
      rc=0

      Call Allocate_DSBA(FLT(1),nBas,nBas,nSym,Case='TRI',Ref=FI)
      Call Allocate_DSBA(FLT(2),nBas,nBas,nSym,Case='TRI',Ref=FA)
      Call Allocate_DSBA(CMO,nBas,nBas,nSym,Ref=W_CMO)

*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      IF (TraOnly) THEN
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*        Let us construct a second pointer structure based on nOrb

         Call Allocate_DSBA(FLT_MO(1),nOrb,nOrb,nSym,Case='TRI',Ref=FI)
         Call Allocate_DSBA(FLT_MO(2),nOrb,nOrb,nSym,Case='TRI',Ref=FA)
c
c ------ It only performs the MO transformation of FI and FA
c ----------------------------------------------------------
c
*        transform FI/FA from AO to MO basis  (LT-storage)
         Do i = 1, 2
            Do iSym = 1,nSym
               iBas = nBas(iSym)
               iOrb = nOrb(iSym)
               iFro = nFro(iSym)
               Call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
               Call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp1')
               Call Square(FLT(i)%SB(iSym)%A1,Tmp1,1,iBas,iBas)
               Call DGEMM_('N','N',iBas,iOrb,iBas,
     &                     1.0d0,Tmp1,iBas,
     &                           CMO%SB(iSym)%A1(1+iFro*iBas),iBas,
     &                     0.0d0,Tmp2,iBas)
               Call MXMT(Tmp2,iBas,1,
     &                   CMO%SB(iSym)%A1(1+iFro*iBas),iBas,
     &                   FLT_MO(i)%SB(iSym)%A1,
     &                   iOrb,iBas)
               Call mma_deallocate(Tmp2)
               Call mma_deallocate(Tmp1)
            End Do
         End Do

         Call deallocate_DSBA(FLT_MO(1))
         Call deallocate_DSBA(FLT_MO(2))
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      ELSE
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
c --- It only computes FI and FA  in AO-basis and returns
c --- the active integrals (tw|xy)
c --- If specified in input, the routine also computes
c --- the auxiliary Q-matrix stored as Q(av), where a is an AO index
c --- and v refers to the active orbitals only

      nForb(:) = nFro(:)
      nIorb(:) = nIsh(:)
      nAorb(:) = nAsh(:)


C --- Build the packed densities from the Squared ones
      Call Allocate_DSBA(DLT(1),nBas,nBas,nSym,Case='TRI')
      Call Allocate_DSBA(DLT(2),nBas,nBas,nSym,Case='TRI')

      Call Fold(nSym,nBas,DI, DLT(1)%A0)
      Call Fold(nSym,nBas,DA1,DLT(2)%A0)

      FactXI = -1.0D0

!AMS - should this be set differently for ExFac.ne.1?
!     FactXI = 0-ExFac

      If (Deco) Then

         FactXI = -0.5D0

!AMS - should this be set differently for ExFac.ne.1?
!        FactXI = 0-(ExFac*.5d0)

c --- decompose the Inactive density on request
         Call Allocate_DSBA(ChoIn,nBas,nBas,nSym)
         Call Allocate_DSBA(DDec,nBas,nBas,nSym)
         DDec%A0(1:NTot2)=DI(1:NTot2)

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
         If (.not.DoLocK) Then
            nForb(:) = 0
            nIorb(:) = nChI(:)
         EndIf


      Else

         Call Allocate_DSBA(MSQ,nBas,nBas,nSym,Ref=W_CMO)

         nChI(:) = nForb(:)+nIorb(:)

      EndIf

      ipInc = ip_of_Work(MSQ%A0(1))

C --- Reordering of the MOs coefficients to fit cholesky needs
      If (.not.DoLocK) Then

         Call Allocate_DSBA(POrb(1),nChI,nBas,nSym)
         Call Allocate_DSBA(POrb(3),nAOrb,nBas,nSym)

*        nOcs=0
         ioff1=0
         Do iSym=1,nSym

            do ikk=1,nChI(iSym)
               ioff2=ioff1+nBas(iSym)*(ikk-1)
               POrb(1)%SB(iSym)%A2(ikk,:) =
     &           MSQ%SB(iSym)%A2(:,ikk)
            end do

            ioff2=ioff1+nBas(iSym)*(nForb(iSym)+nIorb(iSym))
            do ikk=1,nAorb(iSym)
               jkk = nForb(iSym) + nIorb(iSym) + ikk
               POrb(3)%SB(iSym)%A2(ikk,:) =
     &           CMO%SB(iSym)%A2(:,jkk)
            end do
            ioff1=ioff1+nBas(iSym)**2
*           nOcs = nOcs + nAorb(iSym)**2

         End Do

      Else

C *** Only the active orbitals MO coeff need reordering
         Call Allocate_DSBA(CVa(1),nAorb,nBas,nSym)

         ioff1 = 0
         Do iSym=1,nSym
            ioff2 = ioff1 + nBas(iSym)*(nForb(iSym)+nIorb(iSym))
            do ikk=1,nAorb(iSym)
               jkk = nForb(iSym) + nIorb(iSym) + ikk
               ioff = ioff2+nBas(iSym)*(ikk-1)
               CVa(1)%SB(iSym)%A2(ikk,:) =
     &           CMO%SB(iSym)%A2(:,jkk)
            end do
            ioff1 = ioff1 + nBas(iSym)**2
         End Do

      EndIf

      If (DoActive) Then
C -----  Decompose the active density  -----------------------------

#ifdef _DEBUGPRINT_
         do i=1,nSym
            CALL CD_TESTER(rc,DLT(2)%SB(i)%A1,nBas(i),.true.)
            write(6,*) 'DALT for sym=', i
            CALL TRIPRT('DALT',' ',DLT(2)%SB(i)%A1,nBas(i))
         end do
#endif

         Call Allocate_DSBA(CVa(2),nBas,nBas,nSym)
         Call Allocate_DSBA(DDec,nBas,nBas,nSym)
         DDec%A0(1:NTot2)=DA1(1:NTot2)

         Thr = 1.0d-12
         Do i=1,nSym
            if (nAorb(i).gt.0)then
! NOTE(Giovanni): CD will proceed with approx. decompos for QMC
!                 This will avoid warnings for negative-definit
               Call CD_InCore(DDec%SB(i)%A2,nBas(i),
     &                        CVa(2)%SB(i)%A2,nBas(i),
     &                        NumV,Thr,rc)
               If (rc.ne.0) Then
                  write(6,*)SECNAM
     &                      //': ill-defined dens decomp for active'
                  write(6,*) 'rc value produced = ', rc
                  Call abend()
               End If
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

      If (Deco) Call Deallocate_DSBA(ChoIn)

      Call Deallocate_DSBA(DLT(2))
      Call Deallocate_DSBA(DLT(1))

      Call deallocate_DSBA(MSQ)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      ENDIF
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Call deallocate_DSBA(CMO)
      Call deallocate_DSBA(FLT(2))
      Call deallocate_DSBA(FLT(1))

      Return
      END
