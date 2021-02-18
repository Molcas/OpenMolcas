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
      Subroutine get_polar(nPrim,nBas,nAtoms,nCenters,NOCOB,
!EB     &OENE,ONUM,nOrb,OCOF,RCHC,LNearestAtom)
     &OENE,nOrb,OCOF,RCHC,LNearestAtom,LFirstRun)
      Implicit Real*8 (a-h,o-z)

#include "WrkSpc.fh"
#include "MpParam.fh"
#include "Address.fh"
#include "MolProp.fh"

      Dimension OENE(nOrb) ! EB, ONUM(nOrb)
      Dimension oCof(nBas,nPrim)
      Dimension Pol(6,nAtoms,nAtoms), Pd(3,nAtoms)
      Dimension RCHC(3,nBas)
      Logical   LNearestAtom
      Logical   LFirstRun

      iStdOut = 6
      iA = 0 ! Added by EB

      WRITE(iStdOut,*) '  '
      WRITE(iStdOut,*) ' CALCULATE THE POLARIZATION TENSOR '
      WRITE(iStdOut,*) '  '
C
C CONSTRUCT THE POLARIZATION CONTRIBUTION FROM EACH PAIR OF ATOMS
C
C NBOND assigns the array value where bondvaluearray starts
      Do nA=1,nAtoms
         Do nB=1,nAtoms
            Do i=1,6
               POL(i,nA,nB)=0.0D0
            EndDo
         EndDo
      EndDo
      Write(iStdOut,*)
      Write(iStdOut,*)'No occupied orbitals',NOCOB
      Write(iStdOut,*)'No orbitals',nBas
      Write(iStdOut,*)

      Do i = 1 , NOCOB
         Do j = NOCOB+1,nBas
C
C ORBITAL ENERGIES
C
            FOE  = 4.0d0/(OENE(J)-OENE(I))
C
C CALCULATE EXPANSION CENTER FOR THE TWO ORBITALS
C
            RIJX = 0.5d0*(RCHC(1,I)+RCHC(1,J))
            RIJY = 0.5d0*(RCHC(2,I)+RCHC(2,J))
            RIJZ = 0.5d0*(RCHC(3,I)+RCHC(3,J))
            Do nA = 1,nAtoms
               PAX = 0.0d0
               PAY = 0.0d0
               PAZ = 0.0d0
               Do iPBas = 1,nAtomPBas(NA)
                  K = iAtPrTab(iPBas,NA)
                  Do L =  1 , nPrim
                     If(L.gt.K)Then
                     KK = L
                     LL = K
                     Else
                     KK = K
                     LL = L
                     EndIf
                     PAX = PAX+OCOF(I,K)*OCOF(J,L)
     &               *(work(iwork(iMltPlAd(1))+kk*(kk-1)/2+ll-1)
     &               +work(iwork(iMltPlAd(0))+kk*(kk-1)/2+ll-1)
     &               *(CordMltPl(1,1)-RIJX))

                     PAY = PAY+OCOF(I,K)*OCOF(J,L)
     &               *(work(iwork(iMltPlAd(1)+1)+kk*(kk-1)/2+ll-1)
     &               +work(iwork(iMltPlAd(0))+kk*(kk-1)/2+ll-1)
     &               *(CordMltPl(2,1)-RIJY))

                     PAZ = PAZ + OCOF(I,K)*OCOF(J,L)
     &               *(work(iwork(iMltPlAd(1)+2)+kk*(kk-1)/2+ll-1)
     &               +work(iwork(iMltPlAd(0))+kk*(kk-1)/2+ll-1)
     &               *(CordMltPl(3,1)-RIJZ))
                  EndDo
               EndDo
               PD(1,nA)=PAX
               PD(2,nA)=PAY
               PD(3,nA)=PAZ
            EndDo
            Do nA=1,nAtoms
               Do nB=1,nAtoms
                  POL(1,nA,nB)=POL(1,nA,nB)+PD(1,nA)*PD(1,nB)*FOE
                  POL(2,nA,nB)=POL(2,nA,nB)+PD(1,nA)*PD(2,nB)*FOE
                  POL(3,nA,nB)=POL(3,nA,nB)+PD(1,nA)*PD(3,nB)*FOE
                  POL(4,nA,nB)=POL(4,nA,nB)+PD(2,nA)*PD(2,nB)*FOE
                  POL(5,nA,nB)=POL(5,nA,nB)+PD(2,nA)*PD(3,nB)*FOE
                  POL(6,nA,nB)=POL(6,nA,nB)+PD(3,nA)*PD(3,nB)*FOE
               EndDo
            EndDo
         EndDo
      EndDo
      Do nA = 1,nAtoms
         Do i=1,6
           Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA+1)/2-1)=POL(i,nA,nA)
           Work(iAtPolAd+nAtoms*(i-1)+nA-1)=POL(i,nA,nA)
         EndDo
         Do nB = 1,nA-1
            Do i=1,6
               Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1) =
     &         POL(i,nA,nB)+POL(i,nB,nA)
            EndDo
         EndDo
      EndDo
      Do nA = 1,nAtoms
         Do nB = 1,nA-1
            FracA = Frac(nA,nB)
            FracB = 1.0 - FracA
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
!                  Write(iStdOut,*)
!                  Write(iStdOut,*)
!     &            ' Moving polaizability between the atoms', nA, nB
!                  Write(iStdOut,*)' to the atom                  ', iA
!                  Write(iStdOut,*)
               Else
                  iA=nA
               EndIf
            Else
               iA=nA
            EndIf
            IF(BondMat(nA,nB)) THEN
               Do i=1,6
                  Work(iAtPolAd+nAtoms*(i-1)+iA-1)=
     &            Work(iAtPolAd+nAtoms*(i-1)+iA-1)+
     &            Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracA
                  Work(iAtPolAd+nAtoms*(i-1)+nB-1)=
     &            Work(iAtPolAd+nAtoms*(i-1)+nB-1)+
     &            Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracB
               EndDo
            ELSE
               Do i=1,6
                  Work(iAtPolAd+nAtoms*(i-1)+iA-1)=
     &            Work(iAtPolAd+nAtoms*(i-1)+iA-1)+
     &            Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracA
                  Work(iAtPolAd+nAtoms*(i-1)+nB-1)=
     &            Work(iAtPolAd+nAtoms*(i-1)+nB-1)+
     &            Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracB
                  Work(iAtBoPolAd+nCenters*(i-1)+iA*(iA+1)/2-1)=
     &            Work(iAtBoPolAd+nCenters*(i-1)+iA*(iA+1)/2-1)+
     &            Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracA
                  Work(iAtBoPolAd+nCenters*(i-1)+nB*(nB+1)/2-1)=
     &            Work(iAtBoPolAd+nCenters*(i-1)+nB*(nB+1)/2-1)+
     &            Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*FracB
               EndDo
            ENDIF
         EndDo
      EndDo
! If a UHF system is used, correct for the dubble counting of alpha and beta electrons
      If(LFirstRun.eqv..False.) Then
         FracA = 0.5d0
         Do nA = 1,nAtoms
            Do i=1,6
               Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA+1)/2-1)=
     &         Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA+1)/2-1)*FracA
               Work(iAtPolAd+nAtoms*(i-1)+nA-1)=
     &         Work(iAtPolAd+nAtoms*(i-1)+nA-1)*FracA
            EndDo
            Do nB = 1,nA-1
              If(BondMat(nA,nB).eqv..True.) THEN
                  Do i=1,6
                     Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1) =
     &               Work(iAtBoPolAd+nCenters*(i-1)+nA*(nA-1)/2+nB-1)*
     &               FracA
                  EndDo
               EndIf
            EndDo
         EndDo
      EndIf

      Return
      End
