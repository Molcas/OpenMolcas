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
       SubRoutine Read2_ns(rMO1,rMO2,FockI,FockA,
     &                  Temp1,nTemp,Temp2,Temp3,Temp4,
     &                  DI13,DI24,DI,
     &                  DA13,DA24,DA,
     &                  rkappa,idsym,
     &                  Signa,Fact,jSpin,lfat,lfit,lMOt,CMO)
************************************************************************
*                                        ~     ~                       *
*   Monster routine for construction of Fock, MO                       *
*   Handles non (anti) symmetric orbital rotations                     *
*                                                                      *
*   Input: rkappa : Rotation matrix                                    *
*          idsym  : symmetry of perturbation                           *
*          DIR    : Inactive One electron density                      *
*          DAR    : Inactive One electron density                      *
*                                                                      *
*   Scrtch:Temp1, Temp2,Temp3,Temp4,DR,DL                              *
*                                                                      *
*   Output:                                                            *
*         FockI(A):Fock matrix (one index transformed integrals)       *
*                  ~~                                                  *
*         rMO1    (pj|kl)\     Added together this gives oneindex      *
*                     ~~  -->  transformed integrals. They are         *
*         rMO2    (pj|kl)/     separated to make it easy to go         *
*                              to spindependent perturbations          *
*                                                                      *
*   Remember Coulomb type integrals are used to construct              *
*   exchange part Fock matrix and exchange integrals to construct      *
*   Coulomb part.                                                      *
*                                                                      *
*   Sign  =  1                                                         *
*                                                                      *
*   Sign  = -1  {I,K}=KI+signIK                                        *
*                                                                      *
*   jspin =  0  Fock matrixes and MO's needed for singlet              *
*               perturbations                                          *
*                                                                      *
*   jspin =  1  Fock matrixes and MO's needed for triplet              *
*               perturbations                                          *
*                                                                      *
************************************************************************
      Implicit Real*8(a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
      Real*8 rkappa(*),FockA(*),FockI(*),
     &       Temp1(ntemp),Temp2(*),
     &       temp3(*),Temp4(*),
     &       rmo1(*),rmo2(*),
     &       CMO(*),DA(*),DI(*),
     &       DA24(*),DI24(*),
     &       DA13(*),DI13(*)
      Parameter ( half  = 0.5d0 )
      Parameter ( One   = 1.0d0 )
      Parameter ( Two   = 2.0d0 )
      Logical lFAt,lFIT,lmot,singlet
*                                                                      *
************************************************************************
*                                                                      *
*  (mn|pq)=sum(o) T  (on|pq) + sign*T  (mo|pq)+T (mn|oq) +sign*T  (mn|po)
*                  mo                no         po              qo
*
*
*   DL = sum(po) D  T   C     (13)
*     bj          ij pi  bp
*
*                     *
*   DR  = sum(qo) D  T  C     (24)
*     ib           ij jp bp
*                                                                      *
************************************************************************
*                                                                      *
      If (jspin.eq.1) Then
        Singlet=.false.
      Else If (jspin.eq.0) Then
        Singlet=.true.
      Else
        Singlet=.false.
        Write(6,*) 'Error jspin=/=1,0'
        Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*                         t
*     Construct ( C* kappa * D )   &   ( C * kappa ) and store it with the
*     general index as the first index and the occupied as the second.
*     The general index is transformed to AO index (contravariant).
*
*                                                                      *
************************************************************************
*                                                                      *
      Do iS=1,nSym
         If (nBas(iS).ne.0) Then
            Do jS=1,nSym
               If (iEOr(iS-1,jS-1)+1.eq.idsym.and.nB(jS).ne.0) Then
*
                  Call DGEMM_('N','N',
     &                        nBas(iS),nB(jS),nB(jS),
     &                        1.0d0,rkappa(ipMat(is,js)),nBas(iS),
     &                        DI(ipCM(js)),nBas(jS),
     &                        0.0d0,DI24(ipMat(iS,jS)),nBas(iS))
*
                  Call DGEMM_('T','N',
     &                        nBas(iS),nB(jS),nB(jS),
     &                        1.0d0,rkappa(ipMat(js,is)),nBas(jS),
     &                        DI(ipCM(js)),nBas(jS),
     &                        0.0d0,DI13(ipMat(iS,jS)),nBas(iS))
*
                  If (iMethod.eq.2) Then
*
                     Call DGEMM_('N','N',
     &                           nBas(iS),nBas(jS),nBas(jS),
     &                           1.0d0,rkappa(ipMat(is,js)),nBas(iS),
     &                           DA(ipCM(js)),nBas(jS),
     &                           0.0d0,DA24(ipMat(iS,jS)),nBas(iS))
*
                     Call DGEMM_('T','N',
     &                           nBas(iS),nBas(jS),nB(jS),
     &                           1.0d0,rKappa(ipMat(js,iS)),nBas(jS),
     &                           DA(ipCM(js)),nBas(js),
     &                           0.0d0,DA13(ipMat(iS,jS)),nBas(iS))
*
                  End If
*
               End If
            End Do
         End If
      End Do
*
      sign=1.0d0
      If (iMethod.eq.2) Call DScal_(ndens2,signa,DA24,1)
      Call DScal_(ndens2,signa,Di24,1)
*                                                                      *
************************************************************************
*                                                                      *
*   Read in integrals                                                  *
*                                                                      *
************************************************************************
*                                                                      *
      Do iS=1,nSym
       Do jS=1,iS
        ijS=iEOr(iS-1,jS-1)+1
        Do kS=1,nSym
         Do lS=1,ks
          If (nBas(iS)*nBas(js)*nBas(kS)*nBAs(ls).ne.0) Then
           If (iEOr(kS-1,lS-1).eq.ijS-1) Then
            Do iB=1,nB(iS)
             nnB=nB(jS)
             If (iS.eq.jS) nnB=iB
             Do jB=1,nnB
*
              Call COUL(lS,kS,iS,jS,iB,jB,Temp2,Temp3)
*
***************************************************************************
*                                                                         *
*       A   C   T   I   V   E      F   O   C   K   M   A   T   R   I   X  *
*                                                                         *
***************************************************************************
*                                       ~
*-----------------Fpj = sum(i,l) Dil (pl|ij) = sum(ilr) Dil Klr (pr|ij) =
*
*                                IJ  I
*                  Fqj=sum(rI)  I  DL
*                                qr  r
*
              If (lFAT) Then
*
                If (iEOr(ls-1,is-1)+1.eq.idsym) Then
                 ipD=ipMat(lS,iS)+nBas(lS)*(ib-1)
                 ipF=ipMat(kS,jS)+nBas(kS)*(jB-1)
                 Call dGeMV_('T',nBas(lS),nBas(kS),
     &                     -Fact*Sign*half,
     &                     Temp2,nBas(lS),DA24(ipD),1,
     &                     One,FockA(ipF),1)
                End If
*
*                                IJ  I
*                  FqJ=sum(rI)  I  DL
*                                rq  r
*
                If (kS.ne.ls.and.(iEOr(kS-1,is-1)+1.eq.iDSym)) Then
                 ipD=ipMat(kS,iS)+nBas(kS)*(ib-1)
                 ipF=ipMat(lS,jS)+nBas(lS)*(jB-1)
                 Call dGeMV_('N', nBas(lS),nBas(kS),-Fact*Sign*half,
     &                     Temp2,nBas(lS),DA24(ipD),1,
     &                     One,FockA(ipF),1)
                End If
*
*                                JI  J
*                  FqI=sum(rJ)  I  DL
*                                qr  r

                If (iS.ne.jS.or.iB.ne.jb) Then
                 If (iEOr(lS-1,js-1)+1.eq.iDSym) Then
                  ipD=ipMat(lS,jS)+nBas(lS)*(jb-1)
                  ipF=ipMat(kS,iS)+nBas(kS)*(iB-1)
                  Call dGeMV_('T',nBas(lS),nBas(kS),-Fact*Sign*half,
     &                      Temp2,nBas(lS),DA24(ipD),1,
     &                      One,FockA(ipF),1)
                 End If
*
*                                 JI  J
*                  Fqr=sum(rJ) = I  DL
*                                 rq  r
*
                 If (kS.ne.ls.and.(iEOr(kS-1,js-1)+1.eq.iDSym)) Then
                  ipD=ipMat(kS,jS)+nBas(kS)*(jb-1)
                  ipF=ipMat(lS,iS)+nBas(lS)*(iB-1)
                  Call dGeMV_('N', nBas(lS),nBas(kS),-Fact*Sign*half,
     &                     Temp2,nBas(lS),DA24(ipD),1,
     &                     One,FockA(ipF),1)
                 End If

                End If
              End If
******************************************************************************
*                                                                            *
*  I   N   A   C   T   I   V   E      F   O   C   K   M   A   T   R   I   X  *
*                                                                            *
******************************************************************************
              If (lFIT) Then
*                                IJ  I
*                  Fqj=sum(rI)  I  DL
*                                qr  r
*
                If (iEOr(ls-1,is-1)+1.eq.idsym) Then
                 ipD=ipMat(lS,iS)+nBas(lS)*(ib-1)
                 ipF=ipMat(kS,jS)+nBas(kS)*(jB-1)
                 Call dGeMV_('T',nBas(lS),nBas(kS),
     &                     -Fact*Sign*half,
     &                     Temp2,nBas(lS),DI24(ipD),1,
     &                     One,FockI(ipF),1)
                End If
*
*                                IJ  I
*                  FqJ=sum(rI)  I  DL
*                                rq  r
*
                If (kS.ne.ls.and.(iEOr(kS-1,is-1)+1.eq.iDSym)) Then
                 ipD=ipMat(kS,iS)+nBas(kS)*(ib-1)
                 ipF=ipMat(lS,jS)+nBas(lS)*(jB-1)
                 Call dGeMV_('N', nBas(lS),nBas(kS),-Fact*Sign*half,
     &                     Temp2,nBas(lS),DI24(ipD),1,
     &                     One,FockI(ipF),1)
                End If

*
*                                JI  J
*                  FqI=sum(rJ)  I  DL
*                                qr  r

                If (iS.ne.jS.or.iB.ne.jb) Then
                 If (iEOr(lS-1,js-1)+1.eq.iDSym) Then
                  ipD=ipMat(lS,jS)+nBas(lS)*(jb-1)
                  ipF=ipMat(kS,iS)+nBas(kS)*(iB-1)
                  Call dGeMV_('T',nBas(lS),nBas(kS),-Fact*Sign*half,
     &                      Temp2,nBas(lS),DI24(ipD),1,
     &                      One,FockI(ipF),1)
                 End If
*
*                                 JI  J
*                  Fqr=sum(rJ) = I  DL
*                                 rq  r
*
                 If (kS.ne.ls.and.(iEOr(kS-1,js-1)+1.eq.iDSym)) Then
                  ipD=ipMat(kS,jS)+nBas(kS)*(jb-1)
                  ipF=ipMat(lS,iS)+nBas(lS)*(iB-1)
                  Call dGeMV_('N', nBas(lS),nBas(kS),-Fact*Sign*half,
     &                     Temp2,nBas(lS),DI24(ipD),1,
     &                     One,FockI(ipF),1)
                 End If

                End If
              End If
*                                                                      *
************************************************************************
*                                                                      *
*            M O !!!
*             ~           ~
*            (pj|kl) &  (pj|kl)
*                                                                      *
************************************************************************
*                                                                      *
              If (LMOT.and.jB.gt.nIsh(js).and.ib.gt.nish(is).and.
     &            (nBas(kS)*nBas(lS).ne.0)) Then
               iib=ib-nIsh(is)
               jjb=jb-nIsh(js)
               Call DGETMO(Temp2,nbas(ls),nbas(ls),nBas(kS),
     &                     Temp4,nbas(ks))
*

               ipS=iEOr(kS-1,idsym-1)+1
                ip1=ipMO(ls,js,is)+
     &             nBas(ips)*nAsh(ls)*((iiB-1)*nAsh(jS)+jjb-1)
                ip2=ipMO(ips,js,is)+
     &             nBas(ls)*nAsh(ips)*((iiB-1)*nAsh(jS)+jjb-1)
                ip3=ipMO(ls,is,js)+
     &             nBas(ips)*nAsh(ls)*((jjB-1)*nAsh(iS)+iib-1)
                ip4=ipMO(ips,is,js)+
     &             nBas(ls)*nAsh(ips)*((jjB-1)*nAsh(iS)+iib-1)
                If (nBas(ips)*nBas(lS).gt.0)
     &          Call DGEMM_('N','N',
     &                      nBas(ips),nBas(ls),nBas(ks),
     &                      1.0d0,rKappa(ipMat(ips,ks)),nBas(ipS),
     &                      Temp4,nBas(ks),
     &                      0.0d0,Temp3,nBas(ips))
*               ~
*              (pl|ij)
*
                If (nBas(ipS)*nAsh(ls).gt.0)
     &          Call DGEADD2(Fact,
     &                       Temp3(nbas(ips)*nish(ls)+1),nbas(ips),'N',
     &                       rmo1(ip1),nbas(ips),'N',
     &                       rmo1(ip1),nbas(ips),nbas(ips),nash(ls))
                If (nBas(lS)*nAsh(ipS).gt.0)
     &          Call DGEADD2(Fact,
     &                       Temp3(nish(ips)+1),nbas(ips),'T',
     &                       rmo2(ip2),nbas(ls),'N',
     &                       rmo2(ip2),nbas(ls),nbas(ls),nash(ips))
*               ~
*              (pl|ji)
*
                If (is.ne.js.or.ib.ne.jb) Then
                 If (nBas(ipS)*nAsh(lS).gt.0)
     &           Call Dgeadd2(Fact,
     &                      Temp3(nbas(ips)*nish(ls)+1),nbas(ips),'N',
     &                      rmo1(ip3),nbas(ips),'N',
     &                      rmo1(ip3),nbas(ips),nbas(ips),nash(ls))
                 If (nBas(ls)*nAsh(ipS).gt.0)
     &           Call Dgeadd2(Fact,
     &                      Temp3(nish(ips)+1),nbas(ips),'T',
     &                      rmo2(ip4),nbas(ls),'N',
     &                      rmo2(ip4),nbas(ls),nbas(ls),nash(ips))
               End If
               If (ks.ne.ls) Then
                 ipS=iEOr(lS-1,idsym-1)+1
                  ip1=ipMO(ks,js,is)+
     &               nBas(ips)*nAsh(ks)*((iiB-1)*nAsh(jS)+jjb-1)
                  ip2=ipMO(ips,js,is)+
     &               nBas(ks)*nAsh(ips)*((iiB-1)*nAsh(jS)+jjb-1)
                  ip3=ipMO(ks,is,js)+
     &               nBas(ips)*nAsh(ks)*((jjB-1)*nAsh(iS)+iib-1)
                  ip4=ipMO(ips,is,js)+
     &               nBas(ks)*nAsh(ips)*((jjB-1)*nAsh(iS)+iib-1)
                  If (nBas(ipS)*nBas(kS).gt.0)
     &            Call DGEMM_('N','T',
     &                        nBas(ips),nBas(ks),nBas(ls),
     &                        1.0d0,rKappa(ipmat(ips,ls)),nBas(ips),
     &                        Temp4,nBas(ks),
     &                        0.0d0,Temp3,nBas(ips))
*               ~
*              (pk|ij)
*
                  If (nBas(ipS)*nAsh(ks).gt.0)
     &            Call Dgeadd2(Fact,
     &                        Temp3(nbas(ips)*nish(ks)+1),nbas(ips),'N',
     &                        rmo1(ip1),nbas(ips),'N',
     &                        rmo1(ip1),nbas(ips),nbas(ips),nash(ks))
                  If (nBas(ks)*nAsh(ips).gt.0)
     &            Call Dgeadd2(Fact,
     &                        Temp3(nish(ips)+1),nbas(ips),'T',
     &                        rmo2(ip2),nbas(ks),'N',
     &                        rmo2(ip2),nbas(ks),nbas(ks),nash(ips))

*               ~
*              (pk|ji)
*

                  If (is.ne.js.or.ib.ne.jb) Then
                   If (nBas(ipS)*nAsh(ks).gt.0)
     &             Call DGEADD2(Fact,
     &                        Temp3(nbas(ips)*nish(ks)+1),nbas(ips),'N',
     &                        rmo1(ip3),nbas(ips),'N',
     &                        rmo1(ip3),nbas(ips),nbas(ips),nash(ks))
                   If (nBas(ks)*nAsh(ips).gt.0)
     &             Call DGEADD2(Fact,
     &                        Temp3(nish(ips)+1),nbas(ips),'T',
     &                        rmo2(ip4),nbas(ks),'N',
     &                        rmo2(ip4),nbas(ks),nbas(ks),nash(ips))
                 End If
               End If ! kl

               ipS=iEOr(lS-1,idsym-1)+1
                ip1=ipMO(ips,js,is)+
     &             nBas(kS)*nAsh(ips)*((iiB-1)*nAsh(jS)+jjb-1)
                ip2=ipMO(ks,js,is)+
     &             nBas(ipS)*nAsh(ks)*((iiB-1)*nAsh(jS)+jjb-1)
                ip3=ipMO(ips,is,js)+
     &             nBas(kS)*nAsh(ips)*((jjB-1)*nAsh(iS)+iib-1)
                ip4=ipMO(ks,is,js)+
     &             nBas(ipS)*nAsh(ks)*((jjB-1)*nAsh(iS)+iib-1)
                If (nBas(ks)*nBas(ips).gt.0)
     &          Call DGEMM_('N','N',
     &                      nBas(ks),nbas(ips),nBas(ls),
     &                      1.0d0,Temp4,nBas(kS),
     &                      rkappa(ipMat(ls,ips)),nBas(ls),
     &                      0.0d0,Temp3,nBas(ks))
*                ~
*              (pl|ji)
*
                If (nBas(ks)*nAsh(ipS).gt.0)
     &          Call DGEADD2(Fact*Signa,
     &                     Temp3(1+nbas(ks)*nish(ips)),nbas(ks),'N',
     &                     rmo1(ip1),nbas(ks),'N',
     &                     rmo1(ip1),nbas(ks),nbas(ks),nash(ips))
                If (nBas(ips)*nAsh(ks).gt.0)
     &          Call DGEADD2(Fact*Signa,
     &                     Temp3(nish(ks)+1),nbas(ks),'T',
     &                     rmo2(ip2),nbas(ips),'N',
     &                     rmo2(ip2),nbas(ips),nbas(ips),nash(ks))
*                ~
*              (pl|ij)
*
                 If (is.ne.js.or.ib.ne.jb) Then
                   If (nBas(ks)*nAsh(ips).gt.0)
     &             Call DGEADD2(Fact*Signa,
     &                          Temp3(1+nbas(ks)*nish(ips)),
     &                          nbas(ks),'N',
     &                          rmo1(ip3),nbas(ks),'N',
     &                          rmo1(ip3),nbas(ks),nbas(ks),nash(ips))
                   If (nBas(ips)*nAsh(ks).gt.0)
     &             Call DGEADD2(Fact*Signa,
     &                          Temp3(nish(ks)+1),nbas(ks),'T',
     &                          rmo2(ip4),nbas(ips),'N',
     &                          rmo2(ip4),nbas(ips),nbas(ips),nash(ks))
               End If

               If (ks.ne.ls) Then
                 ipS=iEOr(kS-1,idsym-1)+1
                  ip1=ipMO(ips,js,is)+
     &               nBas(ls)*nAsh(ips)*((iiB-1)*nAsh(jS)+jjb-1)
                  ip2=ipMO(ls,js,is)+
     &               nBas(ips)*nAsh(ls)*((iiB-1)*nAsh(jS)+jjb-1)
                  ip3=ipMO(ips,is,js)+
     &               nBas(ls)*nAsh(ips)*((jjB-1)*nAsh(iS)+iib-1)
                  ip4=ipMO(ls,is,js)+
     &               nBas(ips)*nAsh(ls)*((jjB-1)*nAsh(iS)+iib-1)
                  If (nBas(ls)*nBas(ips).gt.0)
     &            Call DGEMM_('T','N',
     &                        nBas(ls),nBas(ips),nBas(ks),
     &                        1.0d0,Temp4,nBas(ks),
     &                        rKappa(ipMat(ks,ips)),nBas(ks),
     &                        0.0d0,Temp3,nBas(ls))
*                ~
*              (pk|ij)
*
                  If (nBas(ls)*nAsh(ips).gt.0)
     &            Call Dgeadd2(Fact*Signa,
     &                         Temp3(1+nbas(ls)*nish(ips)),nbas(ls),'N',
     &                         rmo1(ip1),nbas(ls),'N',
     &                         rmo1(ip1),nbas(ls),nbas(ls),nash(ips))
                  If (nBas(ips)*nAsh(ls).gt.0)
     &            Call Dgeadd2(Fact*Signa,
     &                         Temp3(1+nish(ls)),nbas(ls),'T',
     &                         rmo2(ip2),nbas(ips),'N',
     &                         rmo2(ip2),nbas(ips),nbas(ips),nash(ls))
*                ~
*              (pk|ji)
*
                  If (is.ne.js.or.ib.ne.jb) Then
                   If (nBas(ls)*nAsh(ips).gt.0)
     &             Call DGeAdd2(Fact*Signa,
     &                         Temp3(1+nbas(ls)*nish(ips)),nbas(ls),'N',
     &                         rmo1(ip3),nbas(ls),'N',
     &                         rmo1(ip3),nbas(ls),nbas(ls),nash(ips))
                   If (nBas(ips)*nAsh(ls).gt.0)
     &             Call DGeAdd2(Fact*Signa,
     &                         Temp3(1+nish(ls)),nbas(ls),'T',
     &                         rmo2(ip4),nbas(ips),'N',
     &                         rmo2(ip4),nbas(ips),nbas(ips),nash(ls))
                 End If
               End If ! (kl)
*A
*|
*--------NOT DEBUGGED?
               End If ! lmo
*
             End Do
            End Do
          End If
          End If
         End Do
        End Do
       End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Do iS=1,nSym
       Do jS=1,nSym
        ijS=iEOr(iS-1,jS-1)+1
        Do kS=1,nSym
         Do lS=1,nSym
          If (nBas(is)*nBAs(js)*nBas(ks)*nBas(ls).ne.0) Then
          If (iEOr(kS-1,lS-1)+1.eq.ijS) Then
           Do lB=1,nB(lS)
            Do jB=1,nB(jS)
            Call EXCH(is,js,ks,ls,jb,lb,Temp1,Temp2)
*
******************************************************************************
*                                                                            *
*  I   N   A   C   T   I   V   E      F   O   C   K   M   A   T   R   I   X  *
*                                                                            *
******************************************************************************
             If (lFIT) Then
*                       ~
*             Fpi=Dkl(pi|kl)=Dkl Kkr (pi|rl)=
*
*                    il   l
*             F =   I   DR
*              pi    pr   r
*
              If (singlet) Then
               If (iEOr(ks-1,ls-1)+1.eq.idsym) Then
                ipD=ipMat(kS,lS)+nBas(kS)*(lb-1)
                ipF=ipMat(iS,jS)+nBas(iS)*(jB-1)
                Call dGeMV_('N', nBas(iS),nBas(kS),Fact,
     &                     Temp1,nBas(iS),DI13(ipD),1,
     &                     One,FockI(ipF),1)
*
*                       ~
*             Fpi=Dkl(pi|kl)=Dkl Klr (pi|lr)=(Dkl Klr (pi|rl)
*
*                    il   l
*             F =   I   DL
*              pi    pr   r
*

                ipD=ipMat(kS,lS)+nBas(kS)*(lb-1)
                ipF=ipMat(iS,jS)+nBas(iS)*(jB-1)
                Call dGeMV_('N', nBas(iS),nBas(kS),Fact*Sign,
     &                     Temp1,nBas(iS),DI24(ipD),1,
     &                     One,FockI(ipF),1)
               End If
              End If
*                       ~
*             Fpl=Dkj(pj|kl)=Dkj Kkr (pj|rl)
*
*                         jl   j
*                        I   DR
*                         pr   r
              If (iEOr(ks-1,js-1)+1.eq.idsym) Then
               ipD=ipMat(kS,jS)+nBas(kS)*(jb-1)
               ipF=ipMat(iS,lS)+nBas(iS)*(lB-1)
               Call dGeMV_('N', nBas(iS),nBas(kS),-Fact*half,
     &                     Temp1,nBas(iS),DI13(ipD),1,
     &                     One,FockI(ipF),1)


              End If
             End If
******************************************************************************
*                                                                            *
*          A   C   T   I   V   E      F   O   C   K   M   A   T   R   I   X  *
*                                                                            *
******************************************************************************
             If (lFAT) Then
              If (singlet) Then
*                       ~
*             Fpi=Dkl(pi|kl)=Dkl Kkr (pi|rl)=
*
*                    il   l
*             F =   I   DR
*              pi    pr   r
*
               If (iEOr(ks-1,ls-1)+1.eq.idsym) Then
                ipD=ipMat(kS,lS)+nBas(kS)*(lb-1)
                ipF=ipMat(iS,jS)+nBas(iS)*(jB-1)
                Call dGeMV_('N', nBas(iS),nBas(kS),Fact,
     &                     Temp1,nBas(iS),DA13(ipD),1,
     &                     One,FockA(ipF),1)
*
*                       ~
*             Fpi=Dkl(pi|kl)=Dkl Krl (pi|rl)=
*
*                    il   l
*             F =   I   DL
*              pi    pr   r
*

                ipD=ipMat(kS,lS)+nBas(kS)*(lb-1)
                ipF=ipMat(iS,jS)+nBas(iS)*(jB-1)
                Call dGeMV_('N', nBas(iS),nBas(kS),Fact*Sign,
     &                     Temp1,nBas(iS),DA24(ipD),1,
     &                     One,FockA(ipF),1)
               End If
              End If
*                       ~
*             Fpl=Dkj(pj|kl)=Dkj Kkr (pj|rl)
*
*                         jl   j
*                        I   DR
*                         pr   r
              If (iEOr(kS-1,js-1).eq.idsym-1) Then
                ipD=ipMat(kS,jS)+nBas(kS)*(jb-1)
                ipF=ipMat(iS,lS)+nBas(iS)*(lB-1)
                Call dGeMV_('N', nBas(iS),nBas(kS),-Fact*half,
     &                     Temp1,nBas(iS),DA13(ipD),1,
     &                     One,FockA(ipF),1)

              End If
             End If
*                                                                      *
************************************************************************
*                                                                      *
*            M O !!!
*                ~           ~
*            (pj|kl) &  (pj|kl)
*                                                                      *
************************************************************************
*                                                                      *
             If (lmot.and.(nAsh(js)*nAsh(ls).ne.0).and.
     &          ( (jb.gt.nish(js)).and.(lB.gt.nish(ls))) ) Then
              call dcopy_(nBas(iS)*nBas(kS),Temp1,1,Temp3,1)
*                 ~            ~
*             (pj|kl) &   (pj|lk)
*
*              JL                      JL
*             I   k    (iJ|kL)      & I   k    (iJLk)
*              ir  kr                  ir  kr
              ips=iEOR(kS-1,iDsym-1)+1
              If (nBas(is)*nAsh(ipS).ne.0)
     &         Call DGEMM_('N','T',
     &                     nBas(iS),nAsh(ipS),nBas(kS),
     &                     1.0d0,Temp3,nBas(iS),
     &                     rKappa(ipMat(ips,ks)+nish(ips)),nBas(ips),
     &                     0.0d0,Temp4,nBas(iS))
              ija=jB-nIsh(jS)
              ila=lB-nIsh(lS)
*                 ~
*             (pj|kl)
*
              ip2=ipMO(js,ips,ls)+nBas(iS)*(ija-1)+
     &               nBas(is)*nAsh(js)*nAsh(ips)*(ilA-1)
              ip3=ipMO(js,ls,ips)+nBas(iS)*(ija-1)+
     &               nBas(is)*nAsh(js)*(ilA-1)

              ip1=1
              Do ipa=1,nAsh(ips)
               Call DaXpY_(nBas(iS),Fact,Temp4(ip1),1,rMO1(ip2),1)
               Call DaXpY_(nBas(iS),Fact,Temp4(ip1),1,rMO2(ip3),1)
               ip2=ip2+nBas(is)*nAsh(js)
               ip3=ip3+nBas(is)*nAsh(js)*nash(ls)
               ip1=ip1+nBas(is)
              End Do
*                  ~
*             (pj|lk)
*
              If (nBas(iS)*nAsh(ipS).ne.0)
     &         Call DGEMM_('N','N',
     &                     nBas(iS),nAsh(ipS),nBas(kS),
     &                     1.0d0,Temp3,nBas(iS),
     &                     rKappa(ipMat(ks,ips)+nbas(ks)*nish(ips)),
     &                     nBas(ks),
     &                     0.0d0,Temp4,nBas(iS))
              ip2=ipMO(js,ls,ips)+nBas(iS)*(ija-1)+
     &               nBas(is)*nAsh(js)*(ilA-1)
              ip3=ipMO(js,ips,ls)+nBas(iS)*(ija-1)+
     &               nBas(is)*nAsh(js)*nash(ips)*(ilA-1)
              ip1=1
              Do ipa=1,nAsh(ips)
               Call DaXpY_(nBas(iS),Fact*signa,
     &                    Temp4(ip1),1,rMO1(ip2),1)
               Call DaXpY_(nBas(iS),Fact*signa,
     &                    Temp4(ip1),1,rMO2(ip3),1)
               ip2=ip2+nBas(is)*nAsh(js)*nAsh(lS)
               ip3=ip3+nBas(is)*nAsh(js)
               ip1=ip1+nBas(is)
              End Do
             End If
*
            End Do
           End Do
          End If
          End If
         End Do
        End Do
       End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real_array(CMO)
      End If
      End
      Subroutine dgeAdd2(r,A,LDA,FORMA,B,LDB,FORMB,C,LDC,M,N)
C
C     MATRIX Addition FOR GENERAL MATRICES
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 FORMA,FORMB
      REAL*8 A(*), B(*), C(*)

      IF (FORMA.EQ.'N' .AND. FORMB.EQ.'N') THEN
         Do 10 iRow=0,m-1
         Do 11 iCol=0,n-1
            c(iRow+iCol*ldc+1)=r*a(iRow+iCol*lda+1)+b(iRow+iCol*ldb+1)
   11    Continue
   10    Continue
      ELSE IF (FORMA.EQ.'T' .AND. FORMB.EQ.'N') THEN
         Do 20 iRow=0,m-1
         Do 21 iCol=0,n-1
            c(iRow+iCol*ldc+1)=r*a(iCol+iRow*lda+1)+b(iRow+iCol*ldb+1)
   21    Continue
   20    Continue
      ELSE IF (FORMA.EQ.'N' .AND. FORMB.EQ.'T') THEN
         Do 30 iRow=0,m-1
         Do 31 iCol=0,n-1
            c(iRow+iCol*ldc+1)=r*a(iRow+iCol*lda+1)+b(iCol+iRow*ldb+1)
   31    Continue
   30    Continue
      ELSE IF (FORMA.EQ.'T' .AND. FORMB.EQ.'T') THEN
         Do 40 iRow=0,m-1
         Do 41 iCol=0,n-1
            c(iRow+iCol*ldc+1)=r*a(iCol+iRow*lda+1)+b(iCol+iRow*ldb+1)
   41    Continue
   40    Continue
      ELSE
         Write(6,*) FORMA,FORMB
         Call Abend()
      END IF
      RETURN
      END
