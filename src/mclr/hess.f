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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      SubRoutine Hess(FockC,FockX,rCon,Temp1,Temp2,Temp3,
     &                Temp4, idsym,jdisp,idisp)
*
*     Constructs the connection parts that is dependend on the first
*     derivative of the connection.
*
      Use Arrays, only: Hss, CMO
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"

#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"

      Real*8 Temp1(nDens2),Temp2(nDens2),Temp3(nDens2),
     &       FockC(nDens2),FockX(nDens2),rcon(nDens2),
     &       temp4(*)
      Character*8 Label
*
      idum=1
      call dcopy_(ndens2,[0.0d0],0,Temp3,1)
      If (iAnd(ntpert(idisp),2**3).eq.8) Then
       Do iS=1,nSym
        js=iEOr(is-1,idSym-1)+1
        nnj=nOrb(js)!nash(js)+nash(js)
        If (nOrb(is)*nOrb(js).ne.0)
     &  Call DGEMM_('N','N',
     &              nOrb(is),nnj,nnj,
     &              1.0d0,rCon(ipMat(is,js)),nOrb(is),
     &                    Work(ipF0SQMO+ipCM(jS)-1),nOrb(js),
     &              0.0d0,Temp3(ipMat(is,js)),nOrb(is))

       End Do
       Call DScal_(ndens2,-0.5d0,Temp3,1)
       Call DaXpY_(nDens2,0.5d0,FockC,1,Temp3,1)
       Call DaXpY_(nDens2,1.0d0,FockX,1,Temp3,1)
      Else
        call dcopy_(nDens2,FockX,1,Temp3,1)
      End If
*
*              xa     ca      xa
*     Temp3=Y=F + 1/2F  -1/2 S  F
*
      Len=0
      Do iSym=1,nSym
       Len=Len+lDisp(iSym)*(lDisp(iSym)+1)/2
      End Do
*
      nIn=0
      mdisp=0
      Do iS=1,idsym-1
       mdisp=mdisp+ldisp(is)
       nIn=nIn+lDisp(is)*(lDisp(is)+1)/2
      End Do
*
      Do 310 kDisp=1,ldisp(idsym)
         mDisp=mdisp+1
         If (iAnd(ntpert(mdisp),2**3).eq.0) Goto 310
         iRC=-1
         iOpt=0
         Label='OvrGrd'
         iOp=2**idSym
         Call dRdMck(iRC,iOpt,Label,DspVec(mDisp),Temp1,iop)
         If (iRc.ne.0) Then
            Write (6,*) 'Hess: Error reading MCKINT'
            Write (6,'(A,A)') 'Label=',Label
            Call Abend()
         End If
         ip=1
         Do iS=1,nSym
          Do jS=1,iS
           If (nOrb(is)*nOrb(jS).ne.0) Then
           If (iEOr(iS-1,jS-1).eq.idsym-1) Then
              If (is.eq.js) Then
                Call Square(Temp1(ip),Temp2(ipMat(iS,jS)),
     &                     1,nBas(is),nBas(is))
                ip=ip+nBas(is)*(nBas(iS)+1)/2
              Else
               If (nBas(is)*nBas(js).ne.0)
     &         call dcopy_(nBas(iS)*nBas(jS),
     &                   Temp1(ip),1,Temp2(ipMat(iS,jS)),1)
               ip=ip+ nBas(is)*nBas(js)
              End If
              If (nBas(is)*nBas(js).ne.0) Then
              Call DGEMM_('T','N',
     &                    nOrb(iS),nBAs(jS),nBas(iS),
     &                    1.0d0,CMO(ipCM(iS)),nBas(iS),
     &                    Temp2(ipMat(iS,jS)),nBas(iS),
     &                    0.0d0,Temp4,nOrb(is))
              call dcopy_(nBas(is)*nBas(js),[0.0d0],0,
     &                    temp2(ipMat(is,js)),1)
              Call DGEMM_('N','N',
     &                    nOrb(iS),nOrb(jS),nBas(jS),
     &                    1.0d0,Temp4,nOrb(iS),
     &                    CMO(ipCM(jS)),nBas(jS),
     &                    0.0d0,Temp2(ipMat(iS,jS)),nOrb(iS))
*    &                    nOrb(iS),nBas(jS),nB(jS))
              if (is.ne.js) Then
              call dcopy_(nBas(is)*nBas(js),[0.0d0],0,
     &                   temp2(ipMat(js,is)),1)
              Call DGEMM_('T','T',
     &                    nOrb(js),nOrb(iS),nBas(js),
     &                    1.0d0,CMO(ipCM(js)),nBas(js),
     &                    Temp4,nOrb(is),
     &                    0.0d0,Temp2(ipMat(js,is)),nOrb(js))
*    &                    nbas(js),nBas(js),nB(iS))
              End If
              End If
            End If
           End If
          End Do
         End Do
         Fact=1.0d0
         If (kDisp.eq.jDisp) Fact=2.0d0
         Indx=nIn+Max(kDisp,jDisp)*(Max(kDisp,jDisp)-1)/2+
     &             Min(kDisp,jDisp)
         Hss(Indx)=Hss(Indx)-fact*ddot_(nDens2,Temp2,1,Temp3,1)
 310  Continue
      Return
      End
