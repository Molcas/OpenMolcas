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
* Copyright (C) 2008, Francesco Aquilante                              *
************************************************************************
      SUBROUTINE CHO_LR_MOs(iOK,nDen,nSym,nBas,nIsh,ipCM,ISLT,ISK,ipMSQ)

************************************************************************
*
*   Computes Left (L) and Right (R) Cholesky vectors to represent a
*            non-symmetric density matrix
*
*            D = C[a] * C[b]^T = L * R^T
*
*   The code is generalized to treat density matrices with any number
*            (nDen) of distinct blocks. The corresponding Cholesky
*            vectors are returned in the arrays defined by the
*            pointers ipMSQ
*
*   A non-zero value for the return code iOK indicates that something
*            went wrong and therefore the returned arrays may contain
*            junk
*
*   Author: F. Aquilante    (October 2008)
*
************************************************************************

      Implicit Real*8 (a-h,o-z)
      Integer  iOK, nDen, nSym
      Integer  nBas(nSym), nIsh(nSym), ISLT(nSym), ISK(nSym)
      Integer  ipCM(nDen), ipMSQ(nDen)

#include "WrkSpc.fh"
#include "stdalloc.fh"

      Real*8, Allocatable, Target:: DMat0(:), PMat0(:)
      Real*8, Pointer:: DMat(:,:)=>Null()
      Real*8, Pointer:: PMat(:,:,:,:)=>Null()

************************************************************************

      irc=0
      ikc=0

      nBm=nBas(1)
      Do iSym=2,nSym
         nBm=Max(nBm,nBas(iSym))
      End Do
      Call mma_allocate(DMat0,nBm**2,Label='DMat0')
      If (nDen.gt.1) Then
         Call mma_allocate(PMat0,2*(nDen*nBm)**2,Label='PMat0')
         ipV = ip_of_Work(PMat0(1+(nDen*nBm)**2))
      Else
         ipV = ipMSQ(1)
      EndIf

      iSym=1
      Do while (iSym .le. nSym)

        DMat(1:nBas(iSym),1:nBas(iSym)) => DMat0(1:nBas(iSym)**2)
        If (nDen==1) Then
           PMat(1:nBas(iSym),1:1,1:nBas(iSym),1:1) =>
     &          DMat0(1:nBas(iSym)**2)
        Else
           PMat(1:nBas(iSym),1:nDen,1:nBas(iSym),1:nDen) =>
     &          PMat0(1:(nDen*nBas(iSym))**2)
        End If

        If (nBas(iSym).gt.0 .and. nIsh(iSym).gt.0) then

C --- Inactive P[kl](a,b) = sum_i C[k](a,i)*C[l](b,i)

         n2b=nDen*nBas(iSym)

         Do k=1,nDen

            kCM = ipCM(k) + ISK(iSym)

            Do l=1,nDen

               lCM = ipCM(l) + ISK(iSym)

               Call DGEMM_('N','T',nBas(iSym),nBas(iSym),nIsh(iSym),
     &                      1.0d0,Work(kCM),nBas(iSym),
     &                            Work(lCM),nBas(iSym),
     &                      0.0d0,DMat,nBas(iSym))

               Do ib=1,nBas(iSym)
                  PMat(:,k,ib,l) = DMat(:,ib)
               End Do

            End Do

         End Do

         Ymax=0.0d0
         do ja=1,n2b  ! Loop over the compound index
            k = (ja-1)/nBas(iSym) + 1
            ia = ja - (k-1)*nBas(iSym)
            Ymax=Max(Ymax,PMat(ia,k,ia,k))
         end do
         Thr=1.0d-13*Ymax

         If (nDen.eq.1) ipV = ipMSQ(1) + ISK(iSym)

         CALL CD_InCore(PMat,n2b,Work(ipV),n2b,NumV,Thr,irc)

         If (NumV.ne.nIsh(iSym)) ikc=1

         If (nDen.gt.1) Then
            Do jden=1,nDen
               iMSQ = ipMSQ(jden) + ISK(iSym)
               Do jv=1,NumV
                  ifr=ipV+n2b*(jv-1)+nBas(iSym)*(jden-1)
                  ito=iMSQ+nBas(iSym)*(jv-1)
                  call dcopy_(nBas(iSym),Work(ifr),1,Work(ito),1)
               End Do
            End Do
         End If


        EndIf

        If (irc.ne.0 .or. ikc.ne.0) iSym=nSym

        iSym=iSym+1

        DMat => Null()
        PMat => Null()

      End Do


      If (nDen.gt.1) Call mma_deallocate(PMat0)
      Call mma_deallocate(DMat0)

      iOK=0
      If (irc.ne.0 .or. ikc.ne.0) iOK=1

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(ISLT)
      End
