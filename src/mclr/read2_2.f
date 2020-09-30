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
       SubRoutine Read2_2(rMO1,rMO2,FockI,FockA,
     &                    Temp1,nTemp,Temp2,Temp3,Temp4,
     &                    DI13,DI24,DI,
     &                    DA13,DA24,DA,
     &                    rkappa,idsym,
     &                    Signa,Fact,jSpin,lfat,lfit,lMOt,CMO)
********************************************************************
*                                        ~     ~                   *
*   Monster routine for construction of Fock, MO                   *
*                                                                  *
*   Input: rkappa : Rotation matrix                                *
*          idsym  : symmetry of perturbation                       *
*          DIR    : Inactive One electron density                  *
*          DAR    : Inactive One electron density                  *
*                                                                  *
*   Scrtch:Temp1, Temp2,Temp3,Temp4,DR,DL                          *
*                                                                  *
*   Output:                                                        *
*         FockI(A):Fock matrix (one index transformed integrals)   *
*                  ~~                                              *
*         rMO1    (pj|kl)\     Added together this gives oneindex  *
*                     ~~  -->  transformed integrals. They are     *
*         rMO2    (pj|kl)/     separated to make it easy to go     *
*                              to spindependent perturbations      *
*                                                                  *
*   Remember Coulomb type integrals are used to construct          *
*   exchange part Fock matrix and exchange integrals to construct  *
*   Coulomb part.                                                  *
*                                                                  *
*   Sign  =  1                                                     *
*                                                                  *
*   Sign  = -1  {I,K}=KI+signIK                                    *
*                                                                  *
*   jspin =  0  Fock matrixes and MO's needed for singlet          *
*               perturbations                                      *
*                                                                  *
*   jspin =  1  Fock matrixes and MO's needed for triplet          *
*               perturbations                                      *
*                                                                  *
********************************************************************
      Implicit Real*8(a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
#include "WrkSpc.fh"
#include "intgrl.fh"
      Real*8 rkappa(nDens2),FockA(nDens2),FockI(nDens2),
     &       Temp1(ntemp),Temp2(nDens2),
     &       temp3(nDens2),Temp4(nDens2),
     &       rmo1(nMba),rmo2(nmba),
     &       CMO(nCMO),DA(nCMO),DI(nCMO),
     &       DA24(nDens2),DI24(nDens2),
     &       DA13(nDens2),DI13(nDens2)
      Parameter ( half  = 0.5d0 )
      Parameter ( One   = 1.0d0 )
      Parameter ( Two   = 2.0d0 )
      Logical lFAt,lFIT,lmot,singlet,triplet
*                                                                      *
************************************************************************
*                                                                      *
* (mn|pq)=sum(o) T  (on|pq) + sign*T  (mo|pq)+T (mn|oq) +sign*T  (mn|po)
*                 mo                no         po              qo
*
*
*   DL = sum(po) D  T   C     (13)
*     bj          ij pi  bp
*
*                     *
*   DR  = sum(qo) D  T  C     (24)
*     ib           ij jp bp
*
*
      If (jspin.eq.1) Then
        Triplet=.true.
        Singlet=.false.
      Else If (jspin.eq.0) Then
        Singlet=.true.
        Triplet=.false.
      Else
        Singlet=.false.
        Triplet=.false.
        Write(6,*) 'Error jspin=/=1,0'
        Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*                         t
*     Construct ( C* kappa * D )   &   ( C * kappa ) and store it with
*     the general index as the 1st index and the occupied as the 2nd.
*     The general index is transformed to AO index (contravariant).
*                                                                      *
************************************************************************
*                                                                      *
      Do iS=1,nSym
         If (nOrb(iS).ne.0) Then
            Do jS=1,nSym
               If (iEOr(iS-1,jS-1)+1.eq.idsym.and.nB(jS).ne.0) Then
*
                  Call DGEMM_('N','N',
     &                        nOrb(iS),nB(jS),nOrb(jS),
     &                        1.0d0,rkappa(ipMat(is,js)),nOrb(iS),
     &                        DI(ipCM(js)),nOrb(jS),
     &                        0.0d0,DI24(ipMat(iS,jS)),nOrb(iS))
*
                  Call DGEMM_('T','N',
     &                        nOrb(iS),nB(jS),nOrb(jS),
     &                        1.0d0,rkappa(ipMat(js,is)),nOrb(jS),
     &                        DI(ipCM(js)),nOrb(jS),
     &                        0.0d0,DI13(ipMat(iS,jS)),nOrb(iS))
*
                  If (iMethod.eq.2) Then
*
                     Call DGEMM_('N','N',
     &                           nOrb(iS),nB(jS),nOrb(jS),
     &                           1.0d0,rkappa(ipMat(is,js)),nOrb(iS),
     &                           DA(ipCM(js)),nOrb(jS),
     &                           0.0d0,DA24(ipMat(iS,jS)),nOrb(iS))
*
                     Call DGEMM_('T','N',
     &                           nOrb(iS),nB(jS),nOrb(jS),
     &                           1.0d0,rKappa(ipMat(js,iS)),nOrb(jS),
     &                           DA(ipCM(js)),nOrb(js),
     &                           0.0d0,DA13(ipMat(iS,jS)),nOrb(iS))
*
                  End If
               End If
            End Do
         End If
      End Do
*
      sign=1.0d0
      If (imethod.eq.2) Call DScal_(ndens2,signa,DA24,1)
      Call DScal_(ndens2,signa,Di24,1)
*                                                                      *
************************************************************************
*                                                                      *
*     Read in integrals                                                *
*                                                                      *
************************************************************************
*                                                                      *
      Do iS=1,nSym
       Do jS=1,iS
*
        If (nOrb(iS)*nOrb(jS).ne.0) Then
*
        ijS=iEOr(iS-1,jS-1)+1
        Do kS=1,nSym
         Do lS=1,ks
          If (nOrb(kS)*nOrb(lS).ne.0) Then
*
*
           If (iEOr(kS-1,lS-1).eq.ijS-1) Then
            Do iB=1,nB(iS)
             nnB=nB(jS)
             If (iS.eq.jS) nnB=iB
             Do jB=1,nnB
*
              Call COUL(lS,kS,IS,jS,IB,JB,Temp2,Temp4)

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
                 ipD=ipMat(lS,iS)+nOrb(lS)*(ib-1)
                 ipF=ipMat(kS,jS)+nOrb(kS)*(jB-1)
                 Call dGeMV_('T',nOrb(lS),nOrb(kS),
     &                     -Fact*Sign*half,
     &                     Temp2,nOrb(lS),DA24(ipD),1,
     &                     One,FockA(ipF),1)
                End If
*
*                                IJ  I
*                  FqJ=sum(rI)  I  DL
*                                rq  r
*
                If (kS.ne.ls.and.(iEOr(kS-1,is-1)+1.eq.iDSym)) Then
                 ipD=ipMat(kS,iS)+nOrb(kS)*(ib-1)
                 ipF=ipMat(lS,jS)+nOrb(lS)*(jB-1)
                 Call dGeMV_('N', nOrb(lS),nOrb(kS),-Fact*Sign*half,
     &                     Temp2,nOrb(lS),DA24(ipD),1,
     &                     One,FockA(ipF),1)
                End If
*
*                                JI  J
*                  FqI=sum(rJ)  I  DL
*                                qr  r

                If (iS.ne.jS.or.iB.ne.jb) Then
                 If (iEOr(lS-1,js-1)+1.eq.iDSym) Then
                  ipD=ipMat(lS,jS)+nOrb(lS)*(jb-1)
                  ipF=ipMat(kS,iS)+nOrb(kS)*(iB-1)
                  Call dGeMV_('T',nOrb(lS),nOrb(kS),-Fact*Sign*half,
     &                      Temp2,nOrb(lS),DA24(ipD),1,
     &                      One,FockA(ipF),1)
                 End If
*
*                                 JI  J
*                  Fqr=sum(rJ) = I  DL
*                                 rq  r
*
                 If (kS.ne.ls.and.(iEOr(kS-1,js-1)+1.eq.iDSym)) Then
                  ipD=ipMat(kS,jS)+nOrb(kS)*(jb-1)
                  ipF=ipMat(lS,iS)+nOrb(lS)*(iB-1)
                  Call dGeMV_('N', nOrb(lS),nOrb(kS),-Fact*Sign*half,
     &                     Temp2,nOrb(lS),DA24(ipD),1,
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
                 ipD=ipMat(lS,iS)+nOrb(lS)*(ib-1)
                 ipF=ipMat(kS,jS)+norb(kS)*(jB-1)
                 Call dGeMV_('T',nOrb(lS),nOrb(kS),
     &                     -Fact*Sign*half,
     &                     Temp2,nOrb(lS),DI24(ipD),1,
     &                     One,FockI(ipF),1)
                End If
*
*                                IJ  I
*                  FqJ=sum(rI)  I  DL
*                                rq  r
*
                If (kS.ne.ls.and.(iEOr(kS-1,is-1)+1.eq.iDSym)) Then
                 ipD=ipMat(kS,iS)+nOrb(kS)*(ib-1)
                 ipF=ipMat(lS,jS)+norb(lS)*(jB-1)
                 Call dGeMV_('N', nOrb(lS),nOrb(kS),-Fact*Sign*half,
     &                     Temp2,nOrb(lS),DI24(ipD),1,
     &                     One,FockI(ipF),1)
                End If

*
*                                JI  J
*                  FqI=sum(rJ)  I  DL
*                                qr  r

                If (iS.ne.jS.or.iB.ne.jb) Then
                 If (iEOr(lS-1,js-1)+1.eq.iDSym) Then
                  ipD=ipMat(lS,jS)+nOrb(lS)*(jb-1)
                  ipF=ipMat(kS,iS)+norb(kS)*(iB-1)
                  Call dGeMV_('T',nOrb(lS),nOrb(kS),-Fact*Sign*half,
     &                      Temp2,nOrb(lS),DI24(ipD),1,
     &                      One,FockI(ipF),1)
                 End If
*
*                                 JI  J
*                  Fqr=sum(rJ) = I  DL
*                                 rq  r
*
                 If (kS.ne.ls.and.(iEOr(kS-1,js-1)+1.eq.iDSym)) Then
                  ipD=ipMat(kS,jS)+nOrb(kS)*(jb-1)
                  ipF=ipMat(lS,iS)+norb(lS)*(iB-1)
                  Call dGeMV_('N', nOrb(lS),nOrb(kS),-Fact*Sign*half,
     &                     Temp2,nOrb(lS),DI24(ipD),1,
     &                     One,FockI(ipF),1)
                 End If

                End If
              End If
*
******************************************************************************
*
*             The rest of the Fock matrix is taken from {F,K}
*             and exchange type integrals
*
******************************************************************************
*************************************************************************
*
*            M O !!!
*             ~           ~
*            (pj|kl) &  (pj|kl)
*
************************************************************************
*--------NOT DEBUGGED?
*|
*V
              If (LMOT.and.jB.gt.nIsh(js).and.ib.gt.nish(is).and.
     &           nOrb(ls)*nOrb(ks).gt.0) Then
               iib=ib-nIsh(is)
               jjb=jb-nIsh(js)

               Call DGETMO(Temp2,nOrb(ls),nOrb(ls),nOrb(kS),
     &                     Temp4,nOrb(ks))

*
               ipS=iEOr(kS-1,idsym-1)+1
               If (nOrb(ips)*nAsh(lS).gt.0) Then
                ipI=norb(ks)*nIsh(ls)+1
                ip1=ipMO(ls,js,is)+
     &             nOrb(ips)*nAsh(ls)*((iiB-1)*nAsh(jS)+jjb-1)
                ip4=ipMO(ls,is,js)+
     &             nOrb(ips)*nAsh(ls)*((jjB-1)*nAsh(iS)+iib-1)
                Call DGEMM_('N','N',
     &                      nOrb(ips),nAsh(ls),nOrb(ks),
     &                      1.0d0,rKappa(ipMat(ips,ks)),nOrb(ipS),
     &                      Temp4(ipI),nOrb(ks),
     &                      0.0d0,Temp3,nOrb(ips))
*               ~
*              (pl|ij)
*
                call daxpy_(nOrb(ips)*nash(ls),
     &                     Fact,Temp3,1,rmo1(ip1),1)

*               ~
*              (pl|ji)
*
                If (is.ne.js.or.ib.ne.jb)
     &              call daxpy_(nOrb(ips)*nash(ls),
     &                         Fact,temp3,1,rmo1(ip4),1)
               End If ! nOrb(ips)*nAsh(lS)
*
               If (ks.ne.ls) Then
                 ipS=iEOr(lS-1,idsym-1)+1
                 If (nOrb(ips)*nAsh(ks).gt.0) Then
                  ip1=ipMO(ks,js,is)+
     &               nOrb(ips)*nAsh(ks)*((iiB-1)*nAsh(jS)+jjb-1)
                  ip4=ipMO(ks,is,js)+
     &               nOrb(ips)*nAsh(ks)*((jjB-1)*nAsh(iS)+iib-1)
                  Call DGEMM_('N','T',
     &                        nOrb(ips),nAsh(ks),nOrb(ls),
     &                        1.0d0,rKappa(ipmat(ips,ls)),nOrb(ips),
     &                        Temp4(nIsh(ks)+1),nOrb(ks),
     &                        0.0d0,Temp3,nOrb(ips))
*               ~
*              (pk|ij)
*
                  call daxpy_(nOrb(ips)*nash(ks),
     &                 Fact,temp3,1,rmo1(ip1),1)

*               ~
*              (pk|ji)
*
                  If (is.ne.js.or.ib.ne.jb)
     &             call daxpy_(nOrb(ips)*nash(ks),Fact,
     6                       temp3,1,rmo1(ip4),1)
                 End If !(nOrb(ips)
               End If ! kl

               ipS=iEOr(lS-1,idsym-1)+1
               If (nOrb(ks)*nAsh(ips).gt.0) Then
                ip2=ipMO(ips,js,is)+
     &             nOrb(kS)*nAsh(ips)*((iiB-1)*nAsh(jS)+jjb-1)
                ip3=ipMO(ips,is,js)+
     &             nOrb(kS)*nAsh(ips)*((jjB-1)*nAsh(iS)+iib-1)
                Call DGEMM_('N','N',
     &                      nOrb(ks),nAsh(ips),nOrb(ls),
     &                      1.0d0,Temp4,nOrb(kS),
     &                      rkappa(ipMat(ls,ips)+nOrb(ls)*nIsh(ips)),
     &                      nOrb(ls),
     &                      0.0d0,Temp3,nOrb(ks))
*                ~
*              (pl|ji)
*
                call daxpy_(nOrb(ks)*nash(ips),Fact*Signa,
     &                     temp3,1,rmo1(ip2),1)
*                ~
*              (pl|ij)
*
                 If (is.ne.js.or.ib.ne.jb)
     &             call daxpy_(nOrb(ks)*nash(ips),Fact*Signa,
     &                        temp3,1,rmo1(ip3),1)
               End If ! nOrb(ips)

               If (ks.ne.ls) Then
                 ipS=iEOr(kS-1,idsym-1)+1
                 If (nOrb(ls)*nAsh(ips).gt.0) Then
                  ip1=ipMO(ips,js,is)+
     &               nOrb(ls)*nAsh(ips)*((iiB-1)*nAsh(jS)+jjb-1)
                  ip4=ipMO(ips,is,js)+
     &               nOrb(ls)*nAsh(ips)*((jjB-1)*nAsh(iS)+iib-1)
                  Call DGEMM_('T','N',
     &                        nOrb(ls),nAsh(ips),nOrb(ks),
     &                        1.0d0,Temp4,nOrb(ks),
     &                        rKappa(ipMat(ks,ips)+nOrb(ks)*nIsh(ips)),
     &                        nOrb(ks),
     &                        0.0d0,Temp3,nOrb(ls))
*                ~
*              (pk|ij)
*
                  call daxpy_(nOrb(ls)*nash(ips),
     &                      Fact*Signa,temp3,1,rmo1(ip1),1)
*                ~
*              (pk|ji)
*
                  If (is.ne.js.or.ib.ne.jb)
     &             call daxpy_(nOrb(ls)*nash(ips),Fact*Signa,
     &                        temp3,1,rmo1(ip4),1)
                 End If !(nOrb)
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
*
        End If
*
       End Do
      End Do
*
********************************************************************
*
      Do iS=1,nSym
       Do jS=1,nSym
*
        If (nOrb(iS)*nOrb(jS).ne.0) Then
*
        ijS=iEOr(iS-1,jS-1)+1
        Do kS=1,nSym
         Do lS=1,nSym
*
          If (iEOr(kS-1,lS-1)+1.eq.ijS) Then
          If (nOrb(kS)*nOrb(lS).ne.0) Then
*
           Do lB=1,nB(lS)
            Do jB=1,nB(jS)
*
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
                ipD=ipMat(kS,lS)+nOrb(kS)*(lb-1)
                ipF=ipMat(iS,jS)+norb(iS)*(jB-1)
                Call dGeMV_('N', nOrb(iS),nOrb(kS),Fact,
     &                     Temp1,nOrb(iS),DI13(ipD),1,
     &                     One,FockI(ipF),1)
*
*                       ~
*             Fpi=Dkl(pi|kl)=Dkl Klr (pi|lr)=(Dkl Klr (pi|rl)
*
*                    il   l
*             F =   I   DL
*              pi    pr   r
*

                ipD=ipMat(kS,lS)+nOrb(kS)*(lb-1)
                ipF=ipMat(iS,jS)+nOrb(iS)*(jB-1)
                Call dGeMV_('N', nOrb(iS),nOrb(kS),Fact*Sign,
     &                     Temp1,nOrb(iS),DI24(ipD),1,
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
               ipD=ipMat(kS,jS)+nOrb(kS)*(jb-1)
               ipF=ipMat(iS,lS)+nOrb(iS)*(lB-1)
               Call dGeMV_('N', nOrb(iS),nOrb(kS),-Fact*half,
     &                     Temp1,nOrb(iS),DI13(ipD),1,
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
                ipD=ipMat(kS,lS)+nOrb(kS)*(lb-1)
                ipF=ipMat(iS,jS)+nOrb(iS)*(jB-1)
                Call dGeMV_('N', nOrb(iS),nOrb(kS),Fact,
     &                     Temp1,nOrb(iS),DA13(ipD),1,
     &                     One,FockA(ipF),1)
*
*                       ~
*             Fpi=Dkl(pi|kl)=Dkl Krl (pi|rl)=
*
*                    il   l
*             F =   I   DL
*              pi    pr   r
*

                ipD=ipMat(kS,lS)+nOrb(kS)*(lb-1)
                ipF=ipMat(iS,jS)+nOrb(iS)*(jB-1)
                Call dGeMV_('N', nOrb(iS),nOrb(kS),Fact*Sign,
     &                     Temp1,nOrb(iS),DA24(ipD),1,
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
                ipD=ipMat(kS,jS)+nOrb(kS)*(jb-1)
                ipF=ipMat(iS,lS)+nOrb(iS)*(lB-1)
                Call dGeMV_('N', nOrb(iS),nOrb(kS),-Fact*half,
     &                     Temp1,nOrb(iS),DA13(ipD),1,
     &                     One,FockA(ipF),1)

              End If
             End If
*
*************************************************************************
*
*            M O !!!
*                ~           ~
*            (pj|kl) &  (pj|kl)
*
************************************************************************
*
             If (lmot.and.(nAsh(js)*nAsh(ls).ne.0).and.
     &          ( (jb.gt.nish(js)).and.(lB.gt.nish(ls))) ) Then

              call dcopy_(nOrb(iS)*nOrb(kS),Temp1,1,Temp3,1)

*                 ~            ~
*             (pj|kl) &   (pj|lk)
*
*              JL                      JL
*             I   k    (iJ|kL)      & I   k    (iJLk)
*              ir  kr                  ir  kr
              ips=iEOR(kS-1,iDsym-1)+1
              If (nOrb(iS)*nAsh(ipS).ne.0)
     &         Call DGEMM_('N','T',
     &                     nOrb(iS),nAsh(ipS),nOrb(kS),
     &                     1.0d0,Temp3,nOrb(iS),
     &                     rKappa(ipMat(ips,ks)+nIsh(ips)),nOrb(ips),
     &                     0.0d0,Temp4,nOrb(iS))
              ija=jB-nIsh(jS)
              ila=lB-nIsh(lS)
*                 ~
*             (pj|kl)
*
              ip2=ipMO(js,ips,ls)+nOrb(iS)*(ija-1)+
     &               nOrb(is)*nAsh(js)*nAsh(ips)*(ilA-1)
              ip1=1
              Do ipa=1,nAsh(ips)
               Call DaXpY_(nOrb(iS),Fact,Temp4(ip1),1,rMO2(ip2),1)
               ip2=ip2+nOrb(is)*nAsh(js)
               ip1=ip1+nOrb(is)
                          End Do
*                  ~
*             (pj|lk)
*
              If (nOrb(is)*nAsh(ipS).ne.0)
     &         Call DGEMM_('N','N',
     &                     nOrb(iS),nAsh(ipS),nOrb(kS),
     &                     1.0d0,Temp3,nOrb(iS),
     &                     rKappa(ipMat(ks,ips)+nOrb(ks)*nIsh(ipS)),
     &                     nOrb(ks),
     &                     0.0d0,Temp4,nOrb(iS))
              ip3=ipMO(js,ls,ips)+nOrb(iS)*(ija-1)+
     &               nOrb(is)*nAsh(js)*(ilA-1)
              ip1=1
              Do ipa=1,nAsh(ips)
               Call DaXpY_(nOrb(iS),Fact*signa,
     &                    Temp4(ip1),1,rMO2(ip3),1)
               ip3=ip3+nOrb(is)*nAsh(js)*nAsh(lS)
               ip1=ip1+nOrb(is)
              End Do
             End If
*
            End Do
           End Do
          End If
          End If
         End Do
        End Do
        End If
       End Do
      End Do
*
*

      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real_array(CMO)
      End If
      End
