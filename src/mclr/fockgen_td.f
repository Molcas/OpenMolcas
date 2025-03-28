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
       SubRoutine FockGen_td(d_0,rDens1,rdens2,fock,idsym)
********************************************************************
*                                                                  *
*   Constructs active fockmatrix and Q matrix                      *
*                                                                  *
*   Input: rkappa : Rotation matrix                                *
*          idsym  : symmetry of perturbation                       *
*                                                                  *
*                                                                  *
*   Output:MO     :MO integrals                                    *
*          Fock   :Fock matrix (one index transformed integrals)   *
*          MOtilde:MO (one index transformed integrals)            *
*                                                                  *
********************************************************************
      use Arrays, only: FIMO
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: zero, half, two
      use MCLR_Data, only: nDens2, nNA, ipMat, ipCM, nA
      use input_mclr, only: nSym,nAsh,nIsh,nBas
      Implicit None
      Real*8 d_0
      Integer idSym
      Real*8 Fock(nDens2),
     &       rdens2(*),rDens1(nna,nna)

      Real*8, Allocatable:: MO(:), Scr(:), TQ(:)
      Integer n1, iS, n2, ipS, kS, jS, iA, iAA, jA, jAA, ipF, ipM, kA,
     &        ip1, ip2, ip3
      Real*8 rd, rd1, rd2
      Interface
         SubRoutine AddGrad2(rMat,fact)
         Implicit None
         Real*8 fact
         Real*8 rMat(*)
         End SubRoutine AddGrad2
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      Fock(:)=Zero
*
      n1=0
      Do iS = 1, nSym
         n1=Max(n1,nBas(iS))
      End Do
      n2=n1**2
      Call mma_allocate(MO,n2,Label='MO')
      Call mma_allocate(Scr,n2,Label='Scr')
*
      Do ipS=1,nSym
         Do kS=1,nSym
            Do iS=1,nSym
               jS=iEor(iEor(ipS-1,kS-1),iS-1)+1
*                                                                      *
************************************************************************
*                                                                      *
*              Coulomb term: F  =2(pk|ji)d
*                             kp          ij
*                                                                      *
************************************************************************
*                                                                      *
               If (iEOr(ipS-1,kS-1)+1.eq.iDsym .and.
     &             nBas(ipS)*nIsh(kS).gt.0           ) Then
                  Do iA=1,nAsh(iS)
                     iAA=iA+nIsh(iS)
                     Do jA=1,nAsh(jS)
                        jAA=jA+nIsh(jS)
*
                        Call Coul(ipS,kS,iS,jS,iAA,jAA,MO,Scr)
*
                        rD=rDens1(iA+nA(iS),jA+nA(jS))*Two
                        Call DaXpY_(nBas(ipS)*nIsh(kS),rd,
     &                             MO,1,Fock(ipMat(ipS,Ks)),1)
*
                     End Do
                  End Do
               End If
*                                                                      *
************************************************************************
*                                                                      *
*
*              Exchange term   F = -(pk|ji)d     (i>j)  OBS KAN VARA FEL
*                               pl          kj
*
               If (iEOr(ips-1,iS-1)+1.eq.iDsym .and.
     &             nBas(ipS).gt.0                   ) Then
                  Do iA=1,nIsh(iS)
                     ipF=ipMat(ipS,iS)+nBas(ipS)*(iA-1)
                     Do jA=1,nAsh(jS)
                        jAA=jA+nIsh(js)
*
                        Call Coul(ipS,kS,iS,jS,iA,jAA,MO,Scr)
*
                        ipM=1+nIsh(kS)*nBas(ipS)
                        Do kA=1,nAsh(ks)
*
*                          Two different densities for the exchange term.
*
                           rd1=rDens1(jA+nA(jS),kA+nA(ks))*Two
                           rd2=rDens1(kA+nA(kS),jA+nA(js))*Two
                           Call DaXpY_(nBas(ipS),-rd1*Half,
     &                                MO(ipM),1,Fock(ipF),1)
                           Call DaXpY_(nBas(ipS),rd2*Half,
     &                                MO(ipM),1,
     &                                Fock(ipMat(is,ips)+iA-1),nbas(is))
                           ipM = ipM + nBas(ipS)
                        End Do
                     End Do
                  End Do
               End If
*                                                                      *
************************************************************************
*                                                                      *
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*   KAN GAA FEL
*
*    $F_{pb}=F_{pa}^ID_{ba}$
*
      Do iS=1,nSym
         If (nBas(iS).gt.0) Then
            jS=iEOr(is-1,iDSym-1)+1
            Do iA=1,nAsh(is)
               Do jA=1,nAsh(js)
                  rd2=rDens1(iA+nA(iS),jA+nA(js))
                  rd1=rDens1(jA+nA(jS),iA+nA(is))
                  ip1=nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)
                  ip2=nBas(iS)*(nIsh(js)+jA-1) +ipmat(is,js)
                  ip3=nIsh(js)+jA-1 +ipmat(js,is)
                  Call DaxPy_(nBas(iS),Rd1,FIMO(ip1),1,Fock(ip2),1)
                  Call DaxPy_(nBAs(iS),-Rd2,FIMO(ip1),1,
     &                       Fock(ip3),nbas(js))
               End Do
            End Do
         End If
      End Do
c QB is calc here

      Call CreQADD(Fock,rdens2,idsym,MO,Scr,n2)

      Call mma_allocate(TQ,ndens2,Label='TQ')
      TQ(:)=0.0d0

c QA here
      Call CreQADD2(TQ,rdens2,idsym,MO,Scr,n2)

      Call mma_deallocate(Scr)
      Call mma_deallocate(MO)

      Do iS=1,nsym
         jS=ieor(is-1,idsym-1)+1
         If (nBas(iS)*nBas(jS).gt.0)
     &      Call DGeSub(Fock(ipMat(is,js)),nbas(is),'N',
     &                  TQ(ipMat(js,is)),nbas(js),'T',
     &                  Fock(ipmat(is,js)),nbas(is),
     &                  nbas(is),nbas(js))
      End Do
*
      If (idSym.eq.1) Call AddGrad2(Fock,d_0)
*
      Call DScal_(nDens2,2.0d0,Fock,1)
*
      Call mma_deallocate(TQ)
*                                                                      *
************************************************************************
*                                                                      *
      End SubRoutine FockGen_td
