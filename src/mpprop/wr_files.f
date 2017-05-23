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
*
* OBSERVE! If the output to the mpprop file that this subroutine generates
*   is modified, then notify the person responsible for qmstat, since that
*   program uses the mpprop outputfile.
!EB      Subroutine Wr_Files(nAtoms,nCenters,nMltPl,NORBI,NOCOB,OENE,iBond,
      Subroutine Wr_Files(nAtoms,nCenters,nMltPl,NORBI,NOCOB,NOCOB_b,
     &           OENE,OENE_b,LAllCenters)

      Implicit Real*8 (a-h,o-z)

#include "MpParam.fh"
#include "WrkSpc.fh"
#include "Address.fh"
#include "MolProp.fh"

      Dimension OENE(NOCOB)
      Dimension OENE_b(NOCOB_b)
!EB      Dimension iBond(2,nCenters)
      Real*8 MolPol(6)
      Character*8 Name
      Logical Exist
      Logical LAllCenters

      Do i=1,6
         MolPol(i)=0.0D0
      EndDo

      Lu=30
      Name='MPPROP'
      Call OpnFl(Name,Lu,Exist)
      Rewind(Lu)
*
      Write(Lu,6) '**************************************************'
      Write(Lu,6) '* Molecule'
      Write(Lu,6) Title
      Write(Lu,6) '* Method'
      Write(Lu,6) Method
      Write(Lu,6) '* Level of Multipoles and Polarizabilities'
      Write(Lu,3) nMltPl,1
*
C
C WE NOW HAVE THE MULTIPOLES IN TWO FORMS
C
C A) SUMMED ONTO ATOMIC TERMS ONLY
C B) SUMMED ONTO ATOMS AND BONDS
C
C EXPANSION TO BE USED FOR ELECTROSTATICS
C
C       CENTER
C       MULTIPOLE
C       POLARIZABILITY
C
      If(.not.LAllCenters) Then
      Write(Lu,'(A)') '* Atomic centers '
      Write(Lu,2) nAtoms
      Do i=1,nAtoms
         Write(Lu,3) iAtomType(i),iAtomPar(i),Cen_Lab(i*(i+1)/2)
         write(Lu,1) cor(1,i,i),cor(2,i,i),cor(3,i,i)
         Do iMltPl=0,nMltPl
            nComp=(iMltPl+1)*(iMltPl+2)/2
            Write(Lu,1)(Work(iAtMltPlAd(iMltPl)+nAtoms*(iComp-1)+i-1),
     &      iComp=1,nComp)
         EndDo
         Do iMltPl=0,2
            nComp=(iMltPl+1)*(iMltPl+2)/2
            Write(Lu,1)(Work(iAtBoMltPlAdCopy(iMltPl)+nCenters*
     &      (iComp-1)+i*(i-1)/2+i-1),iComp=1,nComp)
         EndDo
         Write(Lu,1) (Work(iAtPolAd+nAtoms*j+i-1),j=0,5)
         Do j=0,5
            MolPol(j+1)=MolPol(j+1)+Work(iAtPolAd+nAtoms*j+i-1)
         EndDo
      EndDo

      Else

      Write(Lu,6) '* All centers'
      WRITE(Lu,2) nCenters
      Do i=1,nAtoms
         WRITE(Lu,3) iAtomType(i),iAtomPar(i),Cen_Lab(i*(i+1)/2)
         WRITE(Lu,1) Cor(1,i,i),Cor(2,i,i),Cor(3,i,i)
         Do iMltPl=0,nMltPl
            nComp=(iMltPl+1)*(iMltPl+2)/2
            Write(Lu,1)(Work(iAtBoMltPlAd(iMltPl)+nCenters*
     &      (iComp-1)+i*(i-1)/2+i-1),iComp=1,nComp)
         EndDo
         Do iMltPl=0,2
            nComp=(iMltPl+1)*(iMltPl+2)/2
            Write(Lu,1)(Work(iAtBoMltPlAdCopy(iMltPl)+nCenters*
     &      (iComp-1)+i*(i-1)/2+i-1),iComp=1,nComp)
         EndDo
! Begin EB
! if "pointer" iAtBoPolAd is not associated with memory
! print zero
         If( iAtBoPolAd .Eq. 0 ) Then
            Write(Lu,1)0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0
         Else
! End EB
            Write(Lu,1)(Work(iAtBoPolAd+nCenters*k+i*(i-1)/2+i-1),k=0,5)
! Begin EB
         End If
! End EB
      EndDo
      Do i=1,nAtoms
         Do j=1,i-1
            If(BondMat(i,j)) Then
               WRITE(Lu,5)0,iBondPar(i),iAtomType(i),iAtomPar(i),
     &                     iAtomType(j),iAtomPar(j),Cen_Lab(i*(i-1)/2+j)
               WRITE(Lu,1) Cor(1,i,j),Cor(2,i,j),Cor(3,i,j)
               Do iMltPl=0,nMltPl
                  nComp=(iMltPl+1)*(iMltPl+2)/2
                  Write(Lu,1)(Work(iAtBoMltPlAd(iMltPl)+nCenters*
     &            (iComp-1)+i*(i-1)/2+j-1),iComp=1,nComp)
               EndDo
               Do iMltPl=0,2
                  nComp=(iMltPl+1)*(iMltPl+2)/2
                  Write(Lu,1)(Work(iAtBoMltPlAdCopy(iMltPl)+nCenters*
     &            (iComp-1)+i*(i-1)/2+j-1),iComp=1,nComp)
               EndDo
! Begin EB
! if "pointer" iAtBoPolAd is not associated with memory
! print zero
               If( iAtBoPolAd .Eq. 0 ) Then
                  Write(Lu,1)0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0
               Else
! End EB
                  Write(Lu,1)
     &                 (Work(iAtBoPolAd+nCenters*k+i*(i-1)/2+j-1),k=0,5)
! Begin EB
               End If
! End EB
            EndIf
         EndDo
      EndDo

      EndIf

      WRITE(Lu,6) '* Molecule properties '
      Do iMltPl=0,nMltPl
         nComp=(iMltPl+1)*(iMltPl+2)/2
         Write(Lu,1)(Work(iAtMltPlTotAd(iMltPl)+iComp-1),iComp=1,nComp)
      EndDo
      Write(Lu,1) (MolPol(i),i=1,6)
      WRITE(Lu,'(A)') '* Orbital information'
      WRITE(Lu,3) NORBI, NOCOB
      WRITE(Lu,1) EneV
      If(NOCOB.ne.0) WRITE(Lu,1) (OENE(I),I=1,NOCOB)
      If(Method.eq.'UHF-SCF') Then
         WRITE(Lu,3) NOCOB_b
         If(NOcOb_b.ne.0) WRITE(Lu,1) (OENE_b(I),I=1,NOCOB_b)
      EndIf

      Close(Lu)

1     format(3F20.10)
2     format(I5)
3     format(2I5,4X,2A)
!EB 4     format(I5,A1)
5     format(6I5,4X,A)
6     format(A)

      RETURN
      END
