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
      Subroutine wr_mpprop(nAtoms,nCenters,nMltPl,iPol)
      Implicit Real*8 (a-h,o-z)

#include "MpParam.fh"
#include "WrkSpc.fh"
#include "Address.fh"
#include "MolProp.fh"

      Parameter (mxComp=(mxMltPl+1)*(mxMltPl+2)/2)
      Character*16 MltPlLab, MltPlLabs(0:mxMltpl,mxComp) ! "0:" added by EB
      Character*16 String(0:16), PolString
      Real*8 MolPol(6)

      iStdOut = 6

      Do j=1,6
         MolPol(j)=0.0D0
      EndDo
      String(0)='Charge          '
      String(1)='Dipole          '
      String(2)='Quadrupole      '
      String(3)='Octupole        '
      String(4)='Hexadecapole    '
*      Do i=5,16
*         Write(String(i),'(I2,A)') i, '-th       Cartesian'
*      EndDo
      String(5)='5-th       Carte' !sian
      String(6)='6-th       Carte' !sian
      String(7)='7-th       Carte' !sian
      String(8)='8-th       Carte' !sian
      String(9)='9-th       Carte' !sian
      String(10)='10-th       Car' !tesian
      String(11)='11-th       Car' !tesian
      String(12)='12-th       Car' !tesian
      String(13)='13-th       Car' !tesian
      String(14)='14-th       Car' !tesian
      String(15)='15-th       Car' !tesian
      String(16)='16-th       Car' !tesian
      PolString='Polarizability  '

      If(16.lt.mxMltPl) Then
         Write(6,*) 'Increase length of MltPlLab'
         Call Abend()
      EndIf
      Do iMltPl=0,nMltPl
      ilab=0
      Do ix=iMltpl,0,-1
         Do iy=iMltpl-ix,0,-1
            iz=iMltpl-ix-iy
            ilab=ilab+1
            Do i=1,16
               MltPlLab(i:i)=' '
            EndDo
            Do izz=iMltPl,iMltPl-iz+1,-1
               MltPlLab(izz:izz)='Z'
            End Do
            Do iyy=iMltPl-iz,iMltPl-iz-iy+1,-1
               MltPlLab(iyy:iyy)='Y'
            End Do
                               !max(               ,1)   added by EB
            Do ixx=iMltPl-iz-iy,max(iMltPl-nMltPl+1,1),-1
               MltPlLab(ixx:ixx)='X'
            End Do
            MltPlLabs(iMltPl,ilab)=MltPlLab
         End Do
      End Do
      EndDo

      WRITE(iStdOut,*)
      WRITE(iStdOut,*) ' The name of the molecule will be', Title
      WRITE(iStdOut,*)
      WRITE(iStdOut,*)
     &' THE MULTIPOLES AND POLARIZABILITY FOR THE EXPANSION'
      WRITE(iStdOut,*)
     &' ***************************************************',
     &'*************'
      WRITE(iStdOut,*)
      WRITE(iStdOut,*) ' TOTAL NUMBER OF CENTERS   : ',nCenters
      WRITE(iStdOut,*) ' OF WHICH THERE ARE ATOMS  : ',nAtoms
      WRITE(iStdOut,*) ' AND BONDS IN THE MOLECULE : ',nCenters-nAtoms
      WRITE(iStdOut,*)

      Write(iStdOut,'(1x,a16,3a16)')'Coord           ',
     &(MltPlLabs(1,iComp)(1:1),iComp=1,3)
      Do iMltPl=0,nMltPl
         nComp=(iMltPl+1)*(iMltPl+2)/2
         Do iComp=1,nComp,6
            Write(6,'(1x,a,6a16)') String(iMltPl),
     &      (MltPlLabs(iMltPl,jComp)(1:iMltPl),
     &      jComp=iComp,min(iComp+5,nComp))
         EndDo
      EndDo
      If(iPol.gt.0) Then
      Write(iStdOut,'(1x,a,6a16)') PolString,(MltPlLabs(2,iComp)(1:2),
     &iComp=1,6)
      EndIf
      iCount=0
      WRITE(iStdOut,*)
      WRITE(iStdOut,*)
      WRITE(iStdOut,*)
      WRITE(iStdOut,*) 'Multipole expansion for Atoms and Bonds'
      WRITE(iStdOut,*) '***************************************'
      Do i = 1,nAtoms
         iCount=iCount+1
         Write(iStdOut,*)
         Write(iStdOut,'(I5,A8,I5,A3,A10)')
     &   iCount,' Center ',i*(i+1)/2,'   ',CEN_LAB(i*(i+1)/2)
         Write(iStdOut,*)
     &   '**************************************************'
         WRITE(iStdOut,'(1x,a16,3f16.8)')'Coord           ',
     &   Cor(1,I,I),Cor(2,I,I),Cor(3,I,I)
         Do iMltPl=0,nMltPl
         nComp=(iMltPl+1)*(iMltPl+2)/2
         Do iComp=1,nComp,6
         Write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl),
     &   (Work(iAtBoMltPlAd(iMltPl)+
     &   nCenters*(jComp-1)+i*(i+1)/2-1),
     &   jComp=iComp,min(iComp+5,nComp))
         EndDo
         EndDo
         If(iPol.gt.0) Then
         Write(iStdOut,'(1x,a16,6f16.8)')PolString,
     &   (Work(iAtBoPolAd+nCenters*(iComp-1)+i*(i+1)/2-1),
     &   iComp=1,6)
         EndIf
      EndDo
      Do i = 1,nAtoms
         Do j= 1,i-1
         If(BondMat(i,j)) Then
         iCount=iCount+1
         Write(iStdOut,*)
         Write(iStdOut,'(I5,A8,I5,A3,A10)')
     &   iCount,' Center ',i*(i-1)/2+j,'   ',CEN_LAB(i*(i-1)/2+j)
         Write(iStdOut,*)
     &   '**************************************************'
         WRITE(iStdOut,'(1x,a16,3f16.8)')'Coord           ',
     &   Cor(1,I,j),Cor(2,I,j),Cor(3,I,j)
         Do iMltPl=0,nMltPl
         nComp=(iMltPl+1)*(iMltPl+2)/2
         Do iComp=1,nComp,6
         Write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl),
     &   (Work(iAtBoMltPlAd(iMltPl)+
     &   nCenters*(jComp-1)+i*(i-1)/2+j-1),
     &   jComp=iComp,min(iComp+5,nComp))
         EndDo
         EndDo
         If(iPol.gt.0) Then
         Write(iStdOut,'(1x,a16,6f16.8)')PolString,
     &   (Work(iAtBoPolAd+nCenters*(iComp-1)+i*(i-1)/2+j-1),
     &   iComp=1,6)
         EndIf
         EndIf
         EndDo
      EndDo
      Write(iStdOut,*)
      Write(iStdOut,*)
      Write(iStdOut,*)
      Write(iStdOut,*) 'Multipole expansion for Atoms'
      Write(iStdOut,*) '*****************************'

      Do i = 1,nAtoms
         Write(iStdOut,*)
         Write(iStdOut,'(I5,A8,I5,A3,A10)')i,' Atom   ',i,'   ',
     &   CEN_LAB(i*(i+1)/2)
         Write(iStdOut,*)
     &   '**************************************************'
         WRITE(iStdOut,'(1x,a16,3f16.8)')'Coord           ',
     &   Cor(1,I,I),Cor(2,I,I),Cor(3,I,I)
         Do iMltPl=0,nMltPl
         nComp=(iMltPl+1)*(iMltPl+2)/2
         Do iComp=1,nComp,6
         Write(iStdOut,'(1x,a16,6f16.8)') String(iMltPl),
     &   (Work(iAtMltPlAd(iMltPl)+
     &   nAtoms*(jComp-1)+i-1),
     &   jComp=iComp,min(iComp+5,nComp))
         EndDo
         EndDo
         If(iPol.gt.0) Then
         Write(iStdOut,'(1x,a16,6f16.8)')PolString,
     &   (Work(iAtPolAd+nAtoms*(iComp-1)+i-1),
     &   iComp=1,6)
         Do j=1,6
            MolPol(j)=MolPol(j)+Work(iAtPolAd+nAtoms*(j-1)+i-1)
         EndDo
         EndIf
      EndDo
      Write(iStdOut,*)
      Write(iStdOut,*)
      Write(iStdOut,*)
      Write(iStdOut,*)
     &' SUMMED MULTIPOLES AND POLARIZABILITY FOR THE MOLECULE'
      Write(iStdOut,*)
     &' *****************************************************'
      Write(iStdOut,*)
      WRITE(iStdOut,'(1x,a16,3f16.8)')'Coord              ',
     &0.0D0,0.0D0,0.0D0
      Do iMltPl=0,nMltPl
        nComp=(iMltPl+1)*(iMltPl+2)/2
        Do iComp=1,nComp,6
        Write(iStdOut,'(1x,a16,6f16.8)')String(iMltPl),
     &  (Work(iAtMltPlTotAd(iMltPl)+jComp-1),
     &  jComp=iComp,min(iComp+5,nComp))
        EndDo
      EndDo
      If(iPol.gt.0)Write(iStdOut,'(1x,a16,6f16.8)')
     &PolString,(MolPol(j),j=1,6)
      Write(iStdOut,*)
      Write(iStdOut,*)
      WRITE(iStdOut,'(1x,a16,3f16.8)')'Coord              ',
     &0.0D0,0.0D0,0.0D0
      Do iMltPl=0,nMltPl
        nComp=(iMltPl+1)*(iMltPl+2)/2
        Do iComp=1,nComp,6
        Write(iStdOut,'(1x,a16,6f16.8)')String(iMltPl),
     &  (Work(iAtBoMltPlTotAd(iMltPl)+jComp-1),
     &  jComp=iComp,min(iComp+5,nComp))
        EndDo
      EndDo

!EB 10000 format(F16.6)
!EB 10003 format(3F16.6 )
      Return
      End
