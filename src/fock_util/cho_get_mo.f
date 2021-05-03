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
      SUBROUTINE CHO_get_MO(iOK,nDen,nSym,nBas,nIsh,CM,ISLT,ISK,MSQ)
      use Data_Structures, only: DSBA_Type, Allocate_DSBA,
     &                           Deallocate_DSBA
      Implicit Real*8 (a-h,o-z)
      Integer  iOK, nDen, nSym
      Integer  nBas(nSym), nIsh(nSym), ISLT(nSym), ISK(nSym)
      Type (DSBA_Type)  CM(nDen), MSQ(nDen), SMat

#include "stdalloc.fh"

      Real*8, Allocatable:: SXMat(:)
      Real*8, Allocatable, Target:: Dmat0(:)
      Real*8, Pointer:: Dmat(:,:)
************************************************************************

      irc=0
      ikc=0

      i = isk(1)
      i = islt(1)
      nBm=nBas(1)
      Do iSym=2,nSym
         nBm=Max(nBm,nBas(iSym))
      End Do
      Call mma_allocate(Dmat0,nBm**2,Label='Dmat')

      iSym=1
      Do while (iSym .le. nSym)

        DMat(1:nBas(iSym),1:nBas(iSym)) => DMat0(1:nBas(iSym)**2)

        If (nBas(iSym).gt.0 .and. nIsh(iSym).gt.0) then

C --- Inactive D(a,b) = sum_i C(a,i)*C(b,i)

         Call DGEMM_('N','T',nBas(iSym),nBas(iSym),nIsh(iSym),
     &                      1.0d0,CM(1)%SB(iSYm)%A2,nBas(iSym),
     &                            CM(1)%SB(iSYm)%A2,nBas(iSym),
     &                      0.0d0,DMat,nBas(iSym))

         Ymax=0.0d0
         do ja=1,nBas(iSym)
            Ymax=Max(Ymax,DMat(ja,ja))
         end do
         Thr=1.0d-13*Ymax

         CALL CD_InCore(DMat,nBas(iSym),MSQ(1)%SB(iSym)%A2,nBas(iSym),
     &                  NumV,Thr,irc)

         If (NumV.ne.nIsh(iSym)) ikc=1

        EndIf

        If (irc.ne.0 .or. ikc.ne.0) iSym=nSym

        iSym=iSym+1

        DMat=>Null()

      End Do


      If (nDen.eq.2 .and. irc.eq.0 .and. ikc.eq.0) Then

         Call Allocate_DSBA(SMat,nBas,nBas,nSym,Case='TRI')
         Call mma_allocate(SXMat,nBm**2,Label='SXMat')

*        Read overlap integrals (LT-storage) and get Square-storage
         iRc=-1
         iOpt=2
         iComp=1
         iSyLbl=1
         Call RdOne(iRc,iOpt,'Mltpl  0',iComp,SMat%A0,iSyLbl)

*        Compute  X_b[a] = C_b U_a   where  U_a = C_a^T S X_a
*        ----------------------------------------------------
         Do i=1,nSym

           DMat(1:nBas(iSym),1:nBas(iSym)) => DMat0(1:nBas(iSym)**2)

           If (nBas(i).gt.0 .and. nIsh(i).gt.0) then

              CALL SQUARE(SMat%SB(i)%A1,DMat,1,NBas(i),NBas(i))

              call DGEMM_('N','N',nBas(i),nIsh(i),nBas(i),
     &                     1.0d0,DMat,nBas(i),
     &                           MSQ(1)%SB(i)%A2,nBas(i),
     &                     0.0d0,SXMat,nBas(i))

              DMat(1:nIsh(iSym),1:nIsh(iSym)) => DMat0(1:nIsh(iSym)**2)

              call DGEMM_('T','N',nIsh(i),nIsh(i),nBas(i),
     &                     1.0d0,CM(1)%SB(i)%A2,nBas(i),
     &                           SXMat,nBas(i),
     &                     0.0d0,DMat,nIsh(i))

c           write(6,*) ' U_a = C_a^T S X_a   for symmetry block: ',i
c           call cho_output(DMat,1,nIsh(i),1,nIsh(i),
c     &                               nIsh(i),nIsh(i),1,6)

              call DGEMM_('N','N',nBas(i),nIsh(i),nIsh(i),
     &                     1.0d0,CM(2)%SB(i)%A2,nBas(i),
     &                           DMat,nIsh(i),
     &                     0.0d0,MSQ(2)%SB(i)%A2,nBas(i))

           EndIf

           DMat=>Null()

         End Do

         Call mma_deallocate(SXMat)
         Call Deallocate_DSBA(SMat)

      EndIf

      Call mma_deallocate(DMat0)

      iOK=0
      If (irc.ne.0 .or. ikc.ne.0) iOK=1

      Return
      End
