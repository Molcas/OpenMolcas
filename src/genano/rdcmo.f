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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
************************************************************************
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine RdCmo
      Implicit Real*8 (a-h,o-z)
#include "parm.fh"
#include "common.fh"
      Character*6 OneInt,NatOrb,RunFile
      Character*72 line
      Dimension Dummy(1),iDummy(1)
*----------------------------------------------------------------------*
      If(isUHF.eq.1) Then
         Call dCopy_(MxCmo, Cmo2,1, Cmo,1)
         Call dCopy_(MxBasX, Occ2,1, Occ,1)
         Return
      End If
*----------------------------------------------------------------------*
      iSymOne=1
      kSet=kSet+1
      OneInt='ONE'
      NatOrb='NAT'
      RunFile='RUN'
      Write(OneInt(4:6),'(i3.3)') kSet
      Write(NatOrb(4:6),'(i3.3)') kSet
      Write(RunFile(4:6),'(i3.3)') kSet
      Call NameRun(RunFile)
*----------------------------------------------------------------------*
      Write(6,*)
      Write(6,*) '--------------------------------------------------'
      Write(6,*)
      Write(6,'(a,i3,a,f7.3)') ' Adding density matrix',kSet,
     &                         ' with weight',wSet(kSet)
      Write(6,*)
      Write(6,*) 'Reading one-el. file: ',OneInt
      Lu_One=2
      Call OpnOne(irc,0,OneInt,Lu_One)
      Call get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      nDim=0
      Do iSym = 1, nSym
         nDim=nDim+nBas(iSym)
      End Do
      Call Get_cArray('Unique Basis Names',Name,(LENIN8)*nDim)
      Call ClsOne(irc,0)
      Write(6,'(a,i5)') ' nSym:',nSym
      Write(6,'(a,8i5)') ' nBas:',(nBas(i),i=1,nSym)
*     Write(6,'(a,1x,a)') (Name(1,i),Name(2,i),i=1,nDim)
*----------------------------------------------------------------------*
* hack to fix Rolands inconsistent labels
      Do i=1,nDim
         If(Name(i)(LENIN3:LENIN3).eq.'s') Name(i)(LENIN1:LENIN3)='01s'
         If(Name(i)(LENIN3:LENIN3).eq.'p') Name(i)(LENIN1:LENIN3)='02p'
         If(Name(i)(LENIN3:LENIN3).eq.'d') Name(i)(LENIN1:LENIN3)='03d'
         If(Name(i)(LENIN3:LENIN3).eq.'f') Name(i)(LENIN1:LENIN3)='04f'
         If(Name(i)(LENIN3:LENIN3).eq.'g') Name(i)(LENIN1:LENIN3)='05g'
         If(Name(i)(LENIN3:LENIN3).eq.'h') Name(i)(LENIN1:LENIN3)='06h'
         If(Name(i)(LENIN3:LENIN3).eq.'i') Name(i)(LENIN1:LENIN3)='07i'
         If(Name(i)(LENIN3:LENIN3).eq.'k') Name(i)(LENIN1:LENIN3)='08k'
      End Do
*     Write(6,'(a,1x,a)') (Name(1,i),Name(2,i),i=1,nDim)
*----------------------------------------------------------------------*
      IF(kSet.eq.1) Then
         Call Init_GenANO
      Else
         Call Check_genano
      End If
*----------------------------------------------------------------------*
      If(kSet.eq. kRfSet) Then
         Lu_One=2
         Call OpnOne(irc,0,OneInt,Lu_One)
         Call RdOne(irc,6,'Mltpl  0',1,Cmo,iSymOne)
         Call CpOvlp(Cmo,Ssym)
         Call ClsOne(irc,0)
      End If
*----------------------------------------------------------------------*
      Write(6,*)
      Write(6,*) 'Reading orbital file: ',NatOrb
      Lu_=17
      Call chk_vec_UHF(NatOrb,Lu_,isUHF)
      If(isUHF.eq.1) Then
         Call rdvec_(NatOrb,Lu_,'CO',1,nSym,nBas,nBas,
     &     Cmo,Cmo2, Occ, Occ2, Dummy, Dummy,
     &     iDummy,line,0,iErr,iWFtype)
         Write(6,'(a)') '***'
         Write(6,'(a)') '*** rdcmo: fix reading of eps for uhf!!!'
         Write(6,'(a)') '***'
      Else
         Call RdVec(NatOrb,Lu_,'COE',nSym,nBas,nBas,
     &      Cmo,occ, eps, iDummy, line,0, iErr)
      End If
      Write(6,*) 'Orbital set: ',line(:mylen(line))
*----------------------------------------------------------------------*
      If(lftdeg) Then
         indx=1
         Do 100 iSym=1,nSym
            Do 110 iOrb=1,nBas(iSym)
               occ(indx)=(1.0d0+1.0d-3/iOrb)*occ(indx)
               indx=indx+1
110         Continue
100      Continue
      End If
*----------------------------------------------------------------------*
      If(rydgen) Then
         Call RdVec(NatOrb,Lu_,'COE',nSym,nBas,nBas,
     &      Cmo,occ, eps, iDummy, line,0, iErr)
         eps0=1.0d6
         indx=1
         Do 200 iSym=1,nSym
            Do 210 iOrb=1,nBas(iSym)
               If(occ(indx).lt.1.0d-2) Then
                  eps0=Min(eps0,eps(indx))
               End If
               indx=indx+1
210         Continue
200      Continue
*        Write(6,'(a,f12.6)') 'eps0',eps0
         indx=1
         Do 300 iSym=1,nSym
            Do 310 iOrb=1,nBas(iSym)
               If(occ(indx).gt.1.0d-2) Then
                  occ(indx)=0.0d0
               Else If(eps(indx).lt.0.0d0) Then
                  eta=exp(6.9d0*(eps(indx)/eps0-1.0d0))
*                 Write(6,'(a,2f12.6)') 'eps/eta',eps(indx),eta
                  occ(indx)=eta
               Else
                  occ(indx)=0.0d0
               End If
               indx=indx+1
310         Continue
300      Continue
      End If
*----------------------------------------------------------------------*
      Call NameRun('#Pop')
      Return
*----------------------------------------------------------------------*
      Write(6,*) 'Error while reading vector file'
      Call Quit_OnUserError()
      End
