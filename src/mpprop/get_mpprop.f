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
      Subroutine Get_MpProp(nPrim,nBas,nAtoms,nCenters,!nOcOb,
     &                      nMltPl,!oNum,nOrb,
!     &                      oCof,
     &                      ip_D_p,ECENTX,ECENTY,ECENTZ,LNearestAtom,
     &                      LFirstRun,LLumOrb)
      Implicit Real*8 (a-h,o-z)
#include "MpParam.fh"
#include "WrkSpc.fh"
#include "Address.fh"
#include "MolProp.fh"
!      Dimension oNum(nOrb)
!      Dimension oCof(nBas,nPrim)
      Dimension ECENTX(nPrim*(nPrim+1)/2)
      Dimension ECENTY(nPrim*(nPrim+1)/2)
      Dimension ECENTZ(nPrim*(nPrim+1)/2)
      Logical   LNearestAtom
      Logical   LFirstRun
      Logical   LLumOrb

*
      Call Allocate_iWork(ip_iCompMat,(nMltPl+1)**3)
      Call Allocate_Work(ip_Qexp,nPrim**2)
!      Call Allocate_Work(ip_DenMat,nPrim**2)
*
      Call Get_MpProp_(nPrim,nBas,nAtoms,nCenters,!nOcOb,
     &                 nMltPl,!oNum,nOrb,
!     &                 oCof,
     &                 ip_D_p,ECENTX,ECENTY,ECENTZ,LNearestAtom,
     &                 iWork(ip_iCompMat),Work(ip_Qexp),
     &                 LFirstRun,LLumOrb)
*
!      Call Free_Work(ip_DenMat)
      Call Free_Work(ip_Qexp)
      Call Free_iWork(ip_iCompMat)
*
      Return
      End
      Subroutine Get_MpProp_(nPrim,nBas,nAtoms,nCenters,!nOcOb,
     &                      nMltPl,!oNum,nOrb,
!     &                      oCof,
     &                      ip_D_p,ECENTX,ECENTY,ECENTZ,LNearestAtom,
     &                      iCompMat,Qexp,LFirstRun,LLumOrb)
      Implicit Real*8 (a-h,o-z)
#include "MpParam.fh"
#include "WrkSpc.fh"
#include "Address.fh"
#include "MolProp.fh"
!      Dimension oNum(nOrb)
!      Dimension oCof(nBas,nPrim)
      Dimension CorP(3),CorN(3)
      Dimension iCompMat(0:nMltPl,0:nMltPl,0:nMltPl)
      Dimension ECENTX(nPrim*(nPrim+1)/2)
      Dimension ECENTY(nPrim*(nPrim+1)/2)
      Dimension ECENTZ(nPrim*(nPrim+1)/2)
      Dimension Qexp(nPrim,nPrim)
!      Dimension DenMat(nPrim,nPrim)
      Logical   LNearestAtom
      Logical   LFirstRun
      Logical   LLumOrb
*                                                                      *
************************************************************************
*                                                                      *
* SOME KIND OF PROLOG
*
      iStdout = 6
      Write(iStdOut,*)
      If(LLumOrb) Then
        Write(iStdOut,*) ' Using the INPORB file'
        If(LFirstRun) Then
          Write(iStdOut,*) ' This is first run of Get_MpProp'
          Write(iStdOut,*)
          If(Method.eq.'UHF-SCF') Then
            Write(iStdOut,*) ' Number of alpha electrons', nOcOB
          Else
            Write(iStdOut,*) ' Occupation Number ', nOcOB
          EndIf
          Write(iStdOut,*) ' Number of Orbitals' , nOrb
        Else
          Write(iStdOut,*) ' Running Get_MpProp for the second time'
          Write(iStdOut,*) ' Number beta electrons ', nOcOB
          Write(iStdOut,*) ' Number of Orbitals' , nOrb
        EndIf
      Else
        Write(iStdOut,*)' Using the densities from a Molcas calculation'
      EndIf
      Write(iStdOut,*)
* Set the varible that knowes the component
      Do iMltpl = 0,nMltPl
         iComp=0
         Do np=iMltpl,0,-1
            Do nq=iMltpl-np,0,-1
               nl=iMltpl-np-nq
               iComp=iComp+1
               iCompMat(np,nq,nl)=iComp
            EndDo
         EndDo
      EndDo
* Translate dipoles to the center for the quadrupoles
      If(LFirstRun) Then
!      Do iMltpl = 0,1
! An error written by me DH
       Do i=1,nPrim
         Do j=1,i
           Do k=1,3
             Work(iWork(iMltPlAd(1)+k-1)+i*(i-1)/2+j-1)=
     &       Work(iWork(iMltPlAd(1)+k-1)+i*(i-1)/2+j-1)+
     &       (CordMltPl(k,0)-CordMltPl(k,2))*
     &       Work(iWork(iMltPlAd(0))+i*(i-1)/2+j-1)
           EndDo
         EndDo
       EndDo
       Do iMltpl = 0,1
         CordMltPl(1,iMltPl)=CordMltPl(1,2)
         CordMltPl(2,iMltPl)=CordMltPl(2,2)
         CordMltPl(3,iMltPl)=CordMltPl(3,2)
       EndDo
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
*     GET THE DENSITY MATRIX
*
!      If(LLumOrb) Then
!       Do K = 1 , nPrim
!        Do L = 1 , K
!         DenMat(K,L) = 0.0d0
!        EndDo
!       EndDo
!       Do K=1,nPrim
!        Do L=1,K
!         Do i=1,nOcOb
!          DenMat(K,L) = DenMat(K,L) + oNum(I)*oCof(I,K)*oCof(I,L)
!         EndDo
!         DenMat(L,K) = DenMat(K,L)
!        EndDo
!       EndDo
!      Else
!        Do K=1,nPrim
!          Do L=1,K
!            DenMat(L,K) = Work(ip_D_p+K*(K-1)/2+L-1)
!            If(K.ne.L) Then
!              DenMat(K,L) = DenMat(L,K)
!            EndIf
!          EndDo
!        EndDo
!      EndIf
*                                                                      *
************************************************************************
*                                                                      *
* GET THE INTERACTIONSITES
      Write(6,*) !Dummy added due to compiler bug on adils, intelfast.
      Do i=1,nPrim
         Do j=1,i
            Qexp(i,j)=work(iwork(iMltPlAd(0))+i*(i-1)/2+j-1)*
     &      Work(ip_D_p+i*(i-1)/2+j-1)
            Qexp(j,i)=Qexp(i,j)
         EndDo
      EndDo
*
      iA = -99
      Do nA=1,nAtoms

C
C NOW SUM UP THE OTHER MULTIPOLE MOMENTS
C
         Do iMltpl = 0,nMltPl
            iComp=0
            Do np=iMltpl,0,-1
               Do nq=iMltpl-np,0,-1
                  nl=iMltpl-np-nq
                  iComp=iComp+1
                  sum=0.0d0
                  Do ip=0,np
                     Call NoverP(np,ip,rnPoveriP)
                     If(np.eq.ip) Then
                        xfac=rnpoverip
                     Else
                        xfac=rnPoveriP*(CordMltPl(1,ip)
     &                  -COR(1,nA,nA))**(np-ip)
                     EndIf
                     Do iq=0,nq
                        Call NoverP(nq,iq,rnqoveriq)
                        If(nq.eq.iq) Then
                           yfac=rnqoveriq
                        Else
                           yfac=rnqoveriq*(CordMltPl(2,iq)
     &                     -COR(2,nA,nA))**(nq-iq)
                        EndIf
                        Do il=0,nl
                           Call NoverP(nl,il,rnloveril)
                           If(nl.eq.il) Then
                              zfac=rnloveril
                           Else
                              zfac=rnloveril*(CordMltPl(3,il)
     &                        -COR(3,nA,nA))**(nl-il)
                           EndIf
                           If(xfac.eq.0.0D0.or.yfac.eq.0.0D0.or.
     &                     zfac.eq.0.0D0) Goto 10
                           Do iPBas = 1,nAtomPBas(nA)
                              i = iAtPrTab(iPBas,nA)
                              Do jPBas = 1,nAtomPBas(nA)
                                 !!!! Check i>j
                                 If(iAtPrTab(jPBas,nA).gt.i) Then
                                    jj = i
                                    ii = iAtPrTab(jPBas,nA)
                                 Else
                                    jj = iAtPrTab(jPBas,nA)
                                    ii = i
                                 EndIf
                                 sum=sum+
     &                           xfac*yfac*zfac*
     &                           Work(ip_D_p+ii*(ii-1)/2+jj-1)*
     &                           Work(iWork(iMltPlAd(ip+iq+il)+
     &                           iCompMat(ip,iq,il)-1)+ii*(ii-1)/2+jj-1)
                              EndDo
                           EndDo
10                         Continue
                        EndDo
                     EndDo
                  EndDo
                  Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)
     &            +nA*(nA+1)/2-1)=
     &            Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)
!                                minus from the negative sign of the electron
     &            +nA*(nA+1)/2-1)-sum
                  Work(iAtBoMltPlAdCopy(iMltpl)+nCenters*(iComp-1)
     &            +nA*(nA+1)/2-1)=
     &            Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)
     &            +nA*(nA+1)/2-1)
                  Work(iAtMltPlAd(iMltPl)+nAtoms*(iComp-1)+nA-1)=
     &            Work(iAtMltPlAd(iMltPl)+nAtoms*(iComp-1)+nA-1)-sum
               EndDo
            EndDo
         EndDo
C
C THE NUCLEAR CHARGE IS ADDED ON
C
         If(LFirstRun) Then
            Work(iAtMltPlAd(0)+nA-1)=
     &      Work(iAtMltPlAd(0)+nA-1)+Work(iQnuc+nA-1)
            Work(iAtBoMltPlAd(0)+nA*(nA+1)/2-1)=
     &      Work(iAtBoMltPlAd(0)+nA*(nA+1)/2-1)+
     &      Work(iQnuc+nA-1)
            Work(iAtBoMltPlAdCopy(0)+nA*(nA+1)/2-1)=
     &      Work(iAtBoMltPlAdCopy(0)+nA*(nA+1)/2-1)+
     &      Work(iQnuc+nA-1)
         EndIf

         Do nB = 1, nA-1
C
C MULTIPOLE BETWEEN NA,NB
C
            Qp = 0.0
            Qn = 0.0
            Do i=1,3
               CorP(i)= 0.0
               CorN(i)= 0.0
            EndDo
            Do iPBas = 1,nAtomPBas(NA)
              i = iAtPrTab(iPBas,NA)
              Do jPBas = 1,nAtomPBas(NB)
                j = iAtPrTab(jPBas,NB)
                !!!! Check i>j
                If(iAtPrTab(jPBas,nB).gt.i) Then
                   jj = i
                   ii = iAtPrTab(jPBas,nB)
                Else
                   jj = iAtPrTab(jPBas,nB)
                   ii = i
                EndIf
                IF(Qexp(i,j).GT.(0.0D0)) THEN
                   Qp=Qp+Qexp(i,j)
                   CorP(1)=CorP(1)+QEXP(I,J)*ECENTX(ii*(ii-1)/2+jj)
                   CorP(2)=CorP(2)+QEXP(I,J)*ECENTY(ii*(ii-1)/2+jj)
                   CorP(3)=CorP(3)+QEXP(I,J)*ECENTZ(ii*(ii-1)/2+jj)
                ELSE
                   Qn=Qn+Qexp(i,j)
                   CorN(1)= CorN(1)+QEXP(I,J)*ECENTX(ii*(ii-1)/2+J)
                   CorN(2)= CorN(2)+QEXP(I,J)*ECENTY(ii*(ii-1)/2+J)
                   CorN(3)= CorN(3)+QEXP(I,J)*ECENTZ(ii*(ii-1)/2+J)
                ENDIF
              EndDo
            EndDo

            if(Qp.ne.0.0D00) then
               CorP(1)= CorP(1)/Qp
               CorP(2)= CorP(2)/Qp
               CorP(3)= CorP(3)/Qp
            endif
            if(Qn.ne.0.0D00) then
               CorN(1)= CorN(1)/Qn
               CorN(2)= CorN(2)/Qn
               CorN(3)= CorN(3)/Qn
            endif
            Qs=Qp-Qn
            if(Qs.ne.0.0D00) then
               FracP= Qp/Qs
               FracN=1.0D00 - FracP
            else
               FracP=0.0D00
               FracN=0.0D00
            endif
            Cor(1,nA,nB)=CorP(1)*FracP+CorN(1)*FracN
            Cor(2,nA,nB)=CorP(2)*FracP+CorN(2)*FracN
            Cor(3,nA,nB)=CorP(3)*FracP+CorN(3)*FracN

C
C CALCULATE THE WEIGTH OF EACH ATOMIC CENTER BY A SIMPLE GEOMETRIC RATIO
C
            Rwei = SQRT( (Cor(1,nA,nA)-Cor(1,nA,nB))**2
     &      +(Cor(2,nA,nA)-Cor(2,nA,nB))**2
     &      +(Cor(3,nA,nA)-Cor(3,nA,nB))**2 )
            Rtot = SQRT( (Cor(1,nA,nA)-Cor(1,nB,nB))**2
     &      +(Cor(2,nA,nA)-Cor(2,nB,nB))**2
     &      +(Cor(3,nA,nA)-Cor(3,nB,nB))**2 )
C
C FRACTION OF MULTIPOLE EXPANSION TO BE RELATED TO THE PAIR ATOMS
C
            FracB = Rwei/Rtot
            FracA = 1.0D00 - FracB
            Frac(nA,nB) = FracA

C Find the closest atom
            If(LNearestAtom.and.(.not.BondMat(nA,nB))) Then
               Smallest=1.0D200
               Do i=1,nAtoms
                  R=sqrt((Cor(1,nA,nB)-Cor(1,i,i))**2+
     &                    (Cor(2,nA,nB)-Cor(2,i,i))**2+
     &                    (Cor(3,nA,nB)-Cor(3,i,i))**2)
                  If(R.lt.Smallest) Then
                     Smallest=R
                     iA=i
                  EndIf
               EndDo
               RA=sqrt((Cor(1,nA,nB)-Cor(1,nA,nA))**2+
     &                  (Cor(2,nA,nB)-Cor(2,nA,nA))**2+
     &                  (Cor(3,nA,nB)-Cor(3,nA,nA))**2)
               RB=sqrt((Cor(1,nA,nB)-Cor(1,nB,nB))**2+
     &                  (Cor(2,nA,nB)-Cor(2,nB,nB))**2+
     &                  (Cor(3,nA,nB)-Cor(3,nB,nB))**2)
               If( ((iA.ne.nA) .and. (iA.ne.nB) ) .and.
     &             ( (Smallest.lt.RA) .and.  (Smallest.lt.RB) ) ) Then
                  iA=iA
                  FracA=1.0D0
                  FracB=0.0D0
                  Write(iStdOut,*)
                  Write(iStdOut,*)' Moving bond between the atoms  ',
     &            Labe(nA),Labe(nB)
                  Write(iStdOut,*)' to the atom                    ',
     &            Labe(iA)
                  Write(iStdOut,*)
               Else
                  iA=nA
               EndIf
            Else
               iA=nA
            EndIf

*                                                                      *
************************************************************************
*                                                                      *
*     TRANSFORM THE MULTIPOLES TO BONDS AND ATOMS
*
*Do multipoles
*
            Do iMltpl = 0,nMltPl
            iComp=0
            Do np=iMltpl,0,-1
               Do nq=iMltpl-np,0,-1
                  nl=iMltpl-np-nq
                  iComp=iComp+1
                  sum=0.0
                  Do ip=0,np
                     Call NoverP(np,ip,rnPoveriP)
                     If(np.eq.ip) Then
                        xfac=rnpoverip
                     Else
                        xfac=rnPoveriP*
     &                  (CordMltPl(1,ip)-COR(1,nA,nB))**(np-ip)
                     EndIf
                     Do iq=0,nq
                        Call NoverP(nq,iq,rnqoveriq)
                        If(nq.eq.iq) Then
                           yfac=rnqoveriq
                        Else
                           yfac=rnqoveriq*
     &                     (CordMltPl(2,iq)-COR(2,nA,nB))**(nq-iq)
                        EndIf
                        Do il=0,nl
                           Call NoverP(nl,il,rnloveril)
                           If(nl.eq.il) Then
                              zfac=rnloveril
                           Else
                              zfac=rnloveril*
     &                        (CordMltPl(3,il)-COR(3,nA,nB))**(nl-il)
                           EndIf
                           If(xfac.eq.0.0D0.or.yfac.eq.0.0D0.or.
     &                     zfac.eq.0.0D0) Goto 20
                           Do iPBas = 1,nAtomPBas(nA)
                              i = iAtPrTab(iPBas,nA)
                              Do jPBas = 1,nAtomPBas(nB)
                                 !!!! Check i>j
                                 If(iAtPrTab(jPBas,nB).gt.i) Then
                                    jj = i
                                    ii = iAtPrTab(jPBas,nB)
                                 Else
                                    jj = iAtPrTab(jPBas,nB)
                                    ii = i
                                 EndIf
                                 sum=sum+
     &                           xfac*yfac*zfac*2.0d0*
     &                           Work(ip_D_p+ii*(ii-1)/2+jj-1)*
     &                           Work(iWork(iMltPlAd(ip+iq+il)+
     &                           iCompMat(ip,iq,il)-1)+ii*(ii-1)/2+jj-1)
                              EndDo
                           EndDo
20                         Continue
                        EndDo
                     EndDo
                  EndDo
                  Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)
     &            +nA*(nA-1)/2+nB-1)=
     &            Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)
!                                   minus from the negative sign of the electron
     &            +nA*(nA-1)/2+nB-1)-sum
                  ! Copy the multipole arrays
                  Work(iAtBoMltPlAdCopy(iMltpl)+nCenters*(iComp-1)
     &            +nA*(nA-1)/2+nB-1)=
     &            Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)
     &            +nA*(nA-1)/2+nB-1)
               EndDo
            EndDo
            EndDo
*
* If one should move the bond to the closest atom iA would be equal to
* that otherwise iA should be equal to nA
*
           If( ((Method.eq.'UHF-SCF').and.(LFirstRun.eqv..False.)) .or.
     &         ((Method.ne.'UHF-SCF').and.(LFirstRun.eqv..True.)) ) Then
            Do iMltpl = 0,nMltPl
            iComp=0
            Do np=iMltpl,0,-1
               Do nq=iMltpl-np,0,-1
                  nl=iMltpl-np-nq
                  iComp=iComp+1
                  sum_a=0.0D0
                  sum_b=0.0D0
                  Do ip=0,np
                     Call NoverP(np,ip,rnPoveriP)
                     If(np.eq.ip) Then
                        xfac_a=rnPoveriP
                        xfac_b=rnPoveriP
                     Else
                        xfac_a=(COR(1,nA,nB)-COR(1,iA,iA))**(np-ip)*
     &                  rnPoveriP
                        xfac_b=(COR(1,nA,nB)-COR(1,nB,nB))**(np-ip)*
     &                  rnPoveriP
                     EndIf
                     Do iq=0,nq
                        Call NoverP(nq,iq,rnqoveriq)
                        If(nq.eq.iq) Then
                           yfac_a=rnqoveriq
                           yfac_b=rnqoveriq
                        Else
                           yfac_a=(COR(2,nA,nB)-COR(2,iA,iA))**(nq-iq)*
     &                     rnqoveriq
                           yfac_b=(COR(2,nA,nB)-COR(2,nB,nB))**(nq-iq)*
     &                     rnqoveriq
                        EndIf
                        Do il=0,nl
                           Call NoverP(nl,il,rnloveril)
                           If(nl.eq.il) Then
                              zfac_a=rnloveril
                              zfac_b=rnloveril
                           Else
                             zfac_a=(COR(3,nA,nB)-COR(3,iA,iA))**(nl-il)
     &                       *rnloveril
                             zfac_b=(COR(3,nA,nB)-COR(3,nB,nB))**(nl-il)
     &                       *rnloveril
                           EndIf
                           sum_a=sum_a+
     &                     xfac_a*yfac_a*zfac_a*FracA*
     &                     Work(iAtBoMltPlAd(ip+iq+il)+
     &                     nCenters*(iCompMat(ip,iq,il)-1)
     &                     +nA*(nA-1)/2+nB-1)
                           sum_b=sum_b+
     &                     xfac_b*yfac_b*zfac_b*FracB*
     &                     Work(iAtBoMltPlAd(ip+iq+il)+
     &                     nCenters*(iCompMat(ip,iq,il)-1)
     &                     +nA*(nA-1)/2+nB-1)
                        EndDo
                     EndDo
                  EndDo
                  If(BondMat(nA,nB)) Then
                    !If bonding
                    !Do atoms
                    Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+iA-1)=
     &              Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+iA-1)+sum_a
                    Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+nB-1)=
     &              Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+nB-1)+sum_b
                  Else
                    !If not bonding
                    !Do Atoms
                    Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+iA-1)=
     &              Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+iA-1)+sum_a
                    Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+nB-1)=
     &              Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+nB-1)+sum_b
                    !Do bonds
                    Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)
     &              +iA*(iA+1)/2-1)=
     &              Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)
     &              +iA*(iA+1)/2-1)+sum_a
                    Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)
     &              +nB*(nB+1)/2-1)=
     &              Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)
     &              +nB*(nB+1)/2-1)+sum_b
                  EndIf
               EndDo
            EndDo
            EndDo
           EndIf
         EndDo
      EndDo

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nBas)
      End
