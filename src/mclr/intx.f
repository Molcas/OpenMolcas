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
      subRoutine INTX(FockI,Temp1,Temp2,Temp3,Temp4,Fock,
     &                rMo,loper,idisp,r)
*
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "disp_mclr.fh"
#include "WrkSpc.fh"
*
      Character*8 Label
      Real*8 FockI(nDens2),Temp2(ndens2),Temp3(nDens2),Temp4(ndens2),
     &       Temp1(nDens2),Fock(nDens2),rMO(*)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
************************************************************************
      If (nDens2.eq.0) Return
      jDisp=DspVec(idisp)
*
*           x
*     Read P
*
*     Remember Two electron contributions is saved in MO base
*     but one electron contributions is saved in AO base.
*
*
*
      If (iAnd(ntpert(idisp),2**2).eq.4) Then   ! 2 el contribution
       If (iAnd(ntpert(idisp),2**4).eq.2**4) Then  ! from mckinley
       If (iMethod.eq.2) Then
*
*-------------------------------------------------------------------*
*
*     RASSCF
*
*-------------------------------------------------------------------*
*
        Label='TOTAL'
        iop=2**loper
        irc=-1
        iopt=0
        Call dRdMck(iRC,iOpt,Label,jDisp,Fock,iop)
        If (iRc.ne.0) Then
           Write (6,*) 'IntX: Error reading MCKINT'
           Write (6,'(A,A)') 'Label=',Label
           Call QTrace
           Call Abend()
        End If
        Call ReLoad(Fock,loper+1,nbas,norb)
*
        Label='INACTIVE'
        iop=2**loper
        irc=-1
        iopt=0
        Call dRdMck(iRC,iOpt,Label,jDisp,Focki,iop)
        If (iRc.ne.0) Then
           Write (6,*) 'IntX: Error reading MCKINT'
           Write (6,'(A,A)') 'Label=',Label
           Call QTrace
           Call Abend()
        End If
        Call ReLoad(Focki,loper+1,nbas,norb)
*
        Label='MOPERT'
        iop=2**loper
        iopt=0
        irc=-1
        Call dRdMck(iRC,iOpt,Label,jDisp,rMO,iop)
        If (iRc.ne.0) Then
           Write (6,*) 'IntX: Error reading MCKINT'
           Write (6,'(A,A)') 'Label=',Label
           Call QTrace
           Call Abend()
        End If
       Else
*
*-------------------------------------------------------------------*
*
*     SCF
*
*-------------------------------------------------------------------*
*
       Label='TOTAL'
       iop=2**loper
       irc=-1
       iopt=0
       Call dRdMck(iRC,iOpt,Label,jDisp,Focki,iop)
       If (iRc.ne.0) Then
          Write (6,*) 'IntX: Error reading MCKINT'
          Write (6,'(A,A)') 'Label=',Label
          Call QTrace
          Call Abend()
       End If
       call dcopy_(ndens2,[0.0d0],0,fock,1)
       Do iS=1,nSym
        js=iEOR(is-1,loper)+1
        Call Dyax(nOrb(is)*nIsh(js),2.0d0,
     &     Focki(ipMat(is,js)),1,Fock(ipMat(is,js)),1)
       End Do
       End If
       End If
       End If
*
*-------------------------------------------------------------------*
*
* Two electron fock matrix done!
* Lets fix  the one electron matrixes
*
*-------------------------------------------------------------------*
*
      If (iAnd(nTPert(iDisp),2**1).eq.2) Then ! 1 el contribution
      iop=2**loper
      If (iAnd(ntpert(idisp),2**4).eq.2**4) Then  ! from mckinley
       Label='ONEGRD'
       iopt=0
       irc=-1
       Call dRdMck(iRC,iOpt,Label,jDisp,Temp1,iop)
       If (iRc.ne.0) Then
          Write (6,*) 'IntX: Error reading MCKINT'
          Write (6,'(A,A)') 'Label=',Label
          Call QTrace
          Call Abend()
       End If
      Else                                         ! or seward
       Label=SwLbl(idisp)
       iopt=0
       irc=-1
       Call RdOne(irc,iopt,Label,jDisp,temp1,iop)
       If (iRc.ne.0) Then
          Write (6,*) 'IntX: Error reading MCKINT'
          Write (6,'(A,A)') 'Label=',Label
          Call QTrace
          Call Abend()
       End If
       Call DSCAL_(ndens2,1.0d0,Temp1,1)
      End If
      End If
*
      call dcopy_(nDens,[0.0d0],0,Temp2,1)
      ip=1
      Do iS=1,nSym
       Do jS=1,is
        If (nBas(is)*nBas(js).ne.0) Then
        If (iEOr(iS-1,jS-1).eq.loper) Then
           If (is.eq.js) Then
             Call Square(Temp1(ipMatLt(iS,jS)),
     &                   Temp4,
     &                   1,nBas(is),nBas(is))
           Else
             call dcopy_(nBas(iS)*nBas(jS),
     &                  temp1(ipMatLT(iS,Js)),1,
     &                  Temp4,1)
           End If
           Call DGEMM_('T','N',
     &                 nOrb(iS),nBas(jS),nBAs(iS),
     &                 1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(is),
     &                 Temp4,nBas(iS),
     &                 0.0d0,Temp3,nOrb(iS))
           Call DGEMM_('N','N',
     &                 nOrb(is),nB(jS),nBas(jS),
     &                 1.0d0,Temp3,nOrb(iS),
     &                 Work(ipCMO+ipCM(jS)-1),nBas(jS),
     &                 0.0d0,Temp2(ipMat(iS,jS)),nOrb(iS))
           If (is.ne.js) Then
           Call DGEMM_('T','T',
     &                 nOrb(jS),nOrb(iS),nBAs(jS),
     &                 1.0d0,Work(ipCMO+ipCM(jS)-1),nBas(js),
     &                 Temp4,nBas(iS),
     &                 0.0d0,Temp3,nOrb(jS))
           Call DGEMM_('N','N',
     &                 nOrb(js),nB(iS),nBas(iS),
     &                 1.0d0,Temp3,nOrb(jS),
     &                 Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 0.0d0,Temp2(ipMat(jS,iS)),nOrb(jS))
          End If

        End If
        End If
       End Do
      End Do

*
       call dcopy_(ndens2,[0.0d0],0,Temp3,1)
       Do iS=1,nSym
        js=iEOR(is-1,loper)+1
        Do j=1,nAsh(is)+nish(is)
         Do i=1,nAsh(is)+nIsh(is)
          If (i.eq.j.and.i.le.nish(is).and.j.le.nish(is))
     &    Then
           rde=2.0d0
          Else If  (i.gt.nish(is).and.j.gt.nish(is)) Then
           rde=Work(ipG1-1+itri(i-nish(is)+nA(is),
     &           j-nIsh(is)+nA(is)))
          Else
           rde=0.0d0
          end if
          If (rde.ne.0.0d0)
     &    Call DaXpY_(nOrb(js),rDe,
     &               Temp2(ipMat(js,is)+(j-1)*nOrb(js)),1,
     &               Temp3(ipMat(js,is)+(i-1)*nOrb(js)),1)
         End Do
        End Do
       End Do
************************************************************************
*
*
*-------------------------------------------------------------------*
*
* One electron fock matrix done!
*
*-------------------------------------------------------------------*
*
      If (iAnd(ntpert(idisp),2**2).eq.4) Then
        call daxpy_(nDens2,1.0d0,Temp2,1,Focki,1)
        call daxpy_(nDens2,1.0d0,Temp3,1,Fock,1)
      Else
        call dcopy_(ndens2,Temp2,1,FockI,1)
        call dcopy_(ndens2,Temp3,1,Fock,1)
      End If

*
*
      return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(r)
      End
      Subroutine ReLoad(A,idsym,NBAS1,NBAS2)
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Real*8 A(*)
      Integer nbas2(nsym),nbas1(nsym)
      Call GetMem('A','ALLO','REAL',ipA,ndens2)
      Do iS=1,nsym
       js=ieor(is-1,idsym-1) +1
       Do j=0,Min(nbas2(js),nbas1(js))-1
        call dcopy_(Min(nbas1(is),nbas2(is)),
     &             A(ipMat(is,js)+j*nbas1(is)),1,
     &             Work(ipA-1+ipmat(is,js)+j*nbas2(is)),1)
       End Do
      End Do
      call dcopy_(ndens2,Work(ipA),1,A,1)
      Call GetMem('A','FREE','REAL',ipA,ndens2)
      Return
      End
