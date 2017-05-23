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
      Subroutine Wr_Prop(nAtoms,nCenters,nBas,nMltPl,NOCOB,NOCOB_b,
     &           orbe,orbe_b,iPol,LAllCenters)

      IMPLICIT REAL*8 (A-H,O-Z)

      Character*8 MemLabel
      Dimension iCompMat(0:nMltPl,0:nMltPl,0:nMltPl)
      Dimension orbe(NOCOB)
      Dimension orbe_b(NOCOB_b)
      Logical   LAllCenters

#include "MpParam.fh"
#include "WrkSpc.fh"
#include "Address.fh"
#include "MolProp.fh"


      nTotCen=0
      Do i=1,nAtoms
         nTotCen=nTotCen+1
         WRITE(CEN_LAB(i*(i+1)/2),'(A)') Labe(i)
         Do j=1,i
            If(BondMat(i,j)) Then
              nTotCen=nTotCen+1
              WRITE(CEN_LAB(i*(i-1)/2+j),'(3A)')LABE(i),'- ',LABE(j)
            EndIf
         EndDo
      EndDo

      Do iMltpl = 0,nMltPl
         iComp=0
         nComp=(iMltPl+1)*(iMltPl+2)/2
         Write(MemLabel,'(A4,I4.4)') 'AtTo',iMltPl
         Call GetMem(MemLabel,'Allo','Real',iAtMltPlTotAd(iMltPl),nComp)
         Do i=0,nComp-1
            Work(iAtMltPlTotAd(iMltPl)+i)=0.0D0
         EndDo
         Write(MemLabel,'(A4,I4.4)') 'BoTo',iMltPl
         Call GetMem(MemLabel,'Allo','Real',iAtBoMltPlTotAd(iMltPl),
     &   nComp)
         Do i=0,nComp-1
            Work(iAtBoMltPlTotAd(iMltPl)+i)=0.0D0
         EndDo
         Do np=iMltpl,0,-1
            Do nq=iMltpl-np,0,-1
               nl=iMltpl-np-nq
               iComp=iComp+1
               iCompMat(np,nq,nl)=iComp
               Do nA=1,nAtoms
                  Do ip=0,np
                     Call NoverP(np,ip,rnPoveriP)
                     If(np.eq.ip) Then
                        xfac=rnpoverip
                     Else
                        xfac=rnPoveriP*(Cor(1,nA,nA)
     &                  )**(np-ip)
                     EndIf
                     Do iq=0,nq
                        Call NoverP(nq,iq,rnqoveriq)
                        If(nq.eq.iq) Then
                           yfac=rnqoveriq
                        Else
                           yfac=rnqoveriq*(Cor(2,nA,nA)
     &                      )**(nq-iq)
                        EndIf
                        Do il=0,nl
                           Call NoverP(nl,il,rnloveril)
                           If(nl.eq.il) Then
                              zfac=rnloveril
                           Else
                              zfac=rnloveril*(Cor(3,nA,nA)
     &                        )**(nl-il)
                           EndIf
                           fac=xfac*yfac*zfac*
     &                     Work(iAtMltPlAd(ip+iq+il)+
     &                     nAtoms*(iCompMat(ip,iq,il)-1)+
     &                     nA-1)
                           Work(iAtMltPlTotAd(iMltpl)+iComp-1)=
     &                     Work(iAtMltPlTotAd(iMltpl)+iComp-1)+fac
                        EndDo
                     EndDo
                  EndDo
                  Do nB=1,nA
                     If(nA.eq.nB.or.BondMat(nA,nB)) Then
                     Do ip=0,np
                        Call NoverP(np,ip,rnPoveriP)
                        If(np.eq.ip) Then
                           xfac=rnpoverip
                        Else
                           xfac=rnPoveriP*(Cor(1,nA,nB)
     &                     )**(np-ip)
                        EndIf
                        Do iq=0,nq
                           Call NoverP(nq,iq,rnqoveriq)
                           If(nq.eq.iq) Then
                              yfac=rnqoveriq
                           Else
                              yfac=rnqoveriq*(Cor(2,nA,nB)
     &                        )**(nq-iq)
                           EndIf
                           Do il=0,nl
                              Call NoverP(nl,il,rnloveril)
                              If(nl.eq.il) Then
                                 zfac=rnloveril
                              Else
                                 zfac=rnloveril*(Cor(3,nA,nB)
     &                           )**(nl-il)
                              EndIf
                              fac=xfac*yfac*zfac*
     &                        Work(iAtBoMltPlAd(ip+iq+il)+
     &                        nCenters*(iCompMat(ip,iq,il)-1)+
     &                        nA*(nA-1)/2+nB-1)
                              Work(iAtBoMltPlTotAd(iMltpl)+iComp-1)=
     &                        Work(iAtBoMltPlTotAd(iMltpl)+iComp-1)+fac
                           EndDo
                        EndDo
                     EndDo
                     EndIf
                  EndDo
               EndDo
            EndDo
         EndDo
      EndDo

      Call Wr_MpProp(nAtoms,nCenters,nMltPl,iPol)
!EB      Call Wr_Files(nAtoms,nCenters,nMltPl,nBas,NOCOB,orbe,iBond,
      Call Wr_Files(nAtoms,nCenters,nMltPl,nBas,NOCOB,NOCOB_b,
     &     orbe,orbe_b,LAllCenters)

      Do iMltpl = 0,nMltPl
         nComp=(iMltPl+1)*(iMltPl+2)/2
         Write(MemLabel,'(A4,I4.4)') 'AtTo',iMltPl
         Call GetMem(MemLabel,'Free','Real',iAtMltPlTotAd(iMltPl),
     &   iMltPl*nComp)
         Write(MemLabel,'(A4,I4.4)') 'BoTo',iMltPl
         Call GetMem(MemLabel,'Free','Real',iAtBoMltPlTotAd(iMltPl),
     &   iMltPl*nComp)
      EndDo

!EB 111   Format(A,3F15.5,6F10.3,10I5)
      Return
      END
