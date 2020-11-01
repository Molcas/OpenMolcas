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
      SubRoutine SetUp_MCLR(DSYM)
      use Arrays, only: pInt1
*
*   Setup pointers and size of metrixes (includes in Pointers.fh)
*
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"
* for the integrals needed in sigma gen
#include "glbbas_mclr.fh"
#include "Input.fh"
#include "WrkSpc.fh"
      integer dsym
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      ip=1
      Call ICopy(8**3,[0],0,ipMO,1)
*
      nn=0
      nbmx=0
      Do iS=1,nsym
         nA(iS)=nn
         nB(is)=nAsh(is)+nIsh(is)
         nn=nn+nAsh(iS)
         nbmx=Max(nbmx,nBas(iS))
      End Do
*
      Call Set_nbmx(nbmx)
*
      nna=nn
      n2Dens=itri(nnA**2,nnA**2)

      n1Dens=nnA*nnA
      ip=1
      Do kS=1,nSym
         Do lS=1,nsym
            klS=iEOr(lS-1,kS-1)+1
            Do jS=1,nSym
               Do iS=1,nSym
                  ijS=iEOr(iS-1,jS-1)+1
                  ipp=nOrb(iS)*nAsh(jS)*nAsh(kS)*nAsh(lS)
                  If (iEOr(ijS-1,klS-1)+1.eq.DSym.and.(ipp.ne.0)) Then
                     ipMO(kS,lS,jS)=ip
                     ip=ip+ipp
                  End If
               End Do
            End Do
         End Do
      End Do
*
      nmba=ip-1
      nDens=0
      nDensLT=0
      nCMO=0
      nDensC=0
      Call iCopy(8**2,[0],0,ipMat,1)
      Call iCopy(8**2,[0],0,ipMC,1)
      Call iCopy(8**2,[0],0,ipMatLT,1)
      Call iCopy(8,[0],0,ipCM,1)

      Do jS=1,nSym
         Do iS=1,js
            If (iEOr(iS-1,jS-1).eq.DSym-1) Then
               If (is.lt.js) Then
                  ipMatLT(jS,iS)=ipntso(js-1,is-1,2**(dsym-1),nBas)+1
                  nDensLT=nDensLT+nBas(iS)*nBas(jS)
                  iExt0=nOrb(is)-nIsh(is)
                  iExt1=nOrb(iS)-nRs1(iS)
                  iExt2=nOrb(iS)-nRs2(iS)
                  iExt3=nOrb(iS)-nRs3(iS)
                  iInt4=nOrb(js)-nish(js)-nAsh(js)
                  iExt4=nIsh(is)+nAsh(is)
                  ipMC(jS,iS)=nDensC+1
                  nDensC=nDensC+iExt0*nIsh(js)+
     &                          iExt1*nRs1(js)+
     &                          iExt2*nRs2(js)+
     &                          iExt3*nRs3(js)+
     &                          iExt4*iInt4
               End If
               If (is.eq.js) Then
                  i1=nOrb(is)-nish(is)-nAsh(is)
                  ipMatLT(jS,iS)=nDensLT+1
                  nDensLT=nDensLT+nBas(iS)*(nBas(is)+1)/2
                  ipMC(jS,iS)=nDensC+1
                  iint0=nOrb(is)-nIsh(is)
                  iint1=i1+nRs2(is)+nRs3(is)
                  iint2=i1+nRs3(is)
                  iint3=i1
                  nDensC=nDensC+iint0*nIsh(is)+
     &                          iint1*nRs1(is)+
     &                          iint2*nRs2(is)+
     &                          iint3*nRs3(is)
               End If
            End If
         End Do
         ipCM(jS)=nCMO+1
         nCMO=nCMO+nBas(js)**2
      End Do
*
      If (TimeDep) nDensC=2*nDensC
      ndens2=0
      matab=1
      Do jS=1,nSym
         Do iS=1,nSym
            If (iEOr(iS-1,jS-1).eq.DSym-1) Then
               ipMatba(is,js)=matab
               ipMat(jS,iS)=nDens2+1
               nDens2=nDens2+nBas(iS)*nBas(jS)
               matab=matab+nash(js)*nOrb(iS)
            End If
         End Do
      End Do
      ndens=ndens2
*
*    Pointers stored in glbbas, just for making LUCIA happy
*
*
*  To begin with we assume that we have permutation symmetry
*
      if(iMethod.eq.2) then
      iOff=1
      Do iiSym=1,nSym
         iOrb=nRs1(iiSym)+nRs2(iiSym)+nRs3(iiSym)
         Do jjSym=1,iiSym
            jOrb=nRs1(jjSym)+nRs2(jjSym)+nRs3(jjSym)
*
            If (iEOr(iiSym-1,jjSym-1)+1.eq.Dsym) Then
               pINT1(iiSym)=iOff
               If (iiSym.eq.jjSym) Then
                  iOff=iOff+iOrb*(iOrb+1)/2
               Else
                  iOff=iOff+iOrb*jOrb
               End If
            End If
*
         End Do
      End Do
*
      iOff=1
      Do iiSym=1,nSym
         iOrb=nRs1(iiSym)+nRs2(iiSym)+nRs3(iiSym)
         Do jjSym=1,iiSym
            jOrb=nRs1(jjSym)+nRs2(jjSym)+nRs3(jjSym)
*
            ijSym=iEOr(iiSym-1,jjSym-1)+1
            klSym=iEOr(ijSym-1,DSym-1)+1
*
            ijNum= iiSym*(iiSym+1)/2+jjSym
*
            If (iiSym.eq.jjSym) Then
               ijOrb=iOrb*(iOrb+1)/2
            Else
               ijOrb=iOrb*jOrb
            End If
            Do kkSym=1,nsym
               kOrb=nRs1(kkSym)+nRs2(kkSym)+nRs3(kkSym)
*
               llSym=iEor(klSym-1,kkSym-1)+1
               lOrb=nRs1(llSym)+nRs2(llSym)+nRs3(llSym)
*
               If (llSym.gt.kkSym) Goto 100
*
               klNum=kkSym*(kkSym+1)/2+llSym
               If (klNum.gt.ijNum) Goto 100
*
               If (kkSym.eq.llSym) Then
                  klOrb=kOrb*(kOrb+1)/2
               Else
                  klOrb=kOrb*lOrb
               End If
               ip=iiSym-1+nSym*((jjSym-1)+nSym*(kkSym-1))
               If (ijNum.eq.klNum) Then
                  iPlus=ijOrb*(ijOrb+1)/2
               Else
                  iPlus=ijOrb*klOrb
               End If
*
               If (iPlus.gt.0) iWork(KpINT2+ip)=iOff
*
               iOff=iOff+iPlus
*
 100           Continue
            End Do
         End Do
      End Do
      Endif
*

      Call Get_iArray('nFro',nFro,nSym)
      Do i = 1, nSym
         If (nFro(i).ne.0) Then
            Call WarningMessage(2,
     &               'MCLR module can not handle frozen orbitals!')
            Call Abend()
         End If
      End Do
      Call Put_iArray('nFroPT',nFro,nSym)
      Call Get_iArray('nDel',nDel,nSym)
      Do i = 1, nSym
         If (nDel(i).ne.0) Then
            Call WarningMessage(2,
     &               'MCLR module can not handle deleted orbitals!')
            Call Abend()
         End If
      End Do
      Call Put_iArray('nDelPT',nDel,nSym)

*     Call iWrtMa( iWork(KpINT2),64,8,64,8)
*
      Return
      End

      Subroutine Set_nbmx(nbmx_)
      Implicit Real*8 (a-h,o-z)
#include "rasdim.fh"
#include "caspt2.fh"
*
      nbmx=nbmx_
*
      Return
      End
