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
* Copyright (C) 2008, Bjorn O. Roos                                    *
*               2008, Valera Veryazov                                  *
************************************************************************
      subroutine expbas(ireturn)
************************************************************************
*                                                                      *
*     Objective: Expand MOs to larger basis set                        *
*                                                                      *
*     B. O. Roos, University of Lund, April 2008.                      *
*                                                                      *
************************************************************************
      Implicit real*8 (a-h,o-z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "info_expbas.fh"
      Dimension nBas1(mxsym),nBas2(mxsym),Occ1(maxbfn),Eorb1(maxbfn),
     &          Occ2(maxbfn),Eorb2(maxbfn)
      Integer indt1(maxbfn),indt2(maxbfn),Indtype(56)
      Character*(LENIN8) Bas1(maxbfn),Bas2(maxbfn)
      Character*80 VecTit
      Character*512 FName
      Logical Exist_1,Exist_2,okay
*----------------------------------------------------------------------*
*     Start program and say Hello                                      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*     Read information from Runfile 1                                  *
*----------------------------------------------------------------------*
      FName='RUNFIL1'
      iLen=mylen(FName)
      Call f_Inquire(FName(:iLen),Exist_1)
      If(.not.Exist_1) then
       Write (6,*) 'Error finding file '//FName(:iLen)
        Call QTrace()
        Call Abend()
      Endif
      Call namerun(FName(:iLen))
      Call get_iScalar('nSym',nSym1)
      Call Get_iArray('nBas',nBas1,nSym1)
      nDim1=0
      nTot1=0
      Do iSym=1,nSym1
         nDim1=nDim1+nBas1(iSym)
         nTot1=nTot1+nBas1(iSym)**2
      End Do
      Call Get_cArray('Unique Basis Names',Bas1,(LENIN8)*nDim1)
*----------------------------------------------------------------------*
*     Read information from Runfile 2                                  *
*----------------------------------------------------------------------*
      FName='RUNFIL2'
      iLen=mylen(FName)
      Call f_Inquire(FName(:iLen),Exist_2)
      If(.not.Exist_2) then
       Write (6,*) 'Error finding file '//FName(:iLen)
        Call QTrace()
        Call Abend()
      Endif
      Call namerun(FName(:iLen))
      Call get_iScalar('nSym',nSym2)
      Call Get_iArray('nBas',nBas2,nSym2)
      nDim2=0
      ntot2=0
      Do iSym=1,nSym2
         nDim2=nDim2+nBas2(iSym)
         ntot2=ntot2+nBas2(iSym)**2
      End Do
      Call Get_cArray('Unique Basis Names',Bas2,(LENIN8)*nDim2)
*----------------------------------------------------------------------*
*     Read MO coefficients from a formatted vector file                *
*----------------------------------------------------------------------*
      Call GetMem('CMO1','Allo','Real',ipCMO1,nTot1)
      Call GetMem('CMO2','Allo','Real',ipCMO2,nTot2)
      FName=EB_FileOrb
      If (mylen(FName).eq.0) FName='INPORB'
      iLen=mylen(FName)
      call f_Inquire (FName(:iLen),okay)
      If ( okay ) Then
        lOcc=0
        LuInpOrb=50
        Call RdVec(FName(:iLen),LuInpOrb,'COEI',
     &  nSym1,nBas1,nBas1,Work(ipCMO1),Occ1,Eorb1,indt1,VecTit,1,iErr)
      Else
        Write (6,*) 'RdCMO: Error finding MO file'
        Call QTrace()
        Call Abend()
      End If

*----------------------------------------------------------------------*
*     Print and check input information                                *
*----------------------------------------------------------------------*
*      Write(6,910) 'Start of option expand.'
  910 Format(/1x,A)
      Write(6,1000) Vectit(:mylen(Vectit))
 1000 Format(/1x,'Header on input orbitals file:'/A)
      Write(6,910) 'Information from input runfile'
      Write(6,920) 'Number of symmetries',nSym1
  920 Format(1x,A30,8i5)
      Write(6,920) 'Number of basis functions',(nBas1(i),i=1,nSym1)
C
      Write(6,910) 'Information from expanded basis set runfile'
      Write(6,920) 'Number of symmetries',nSym2
      Write(6,920) 'Number of basis functions',(nBas2(i),i=1,nSym2)
C     Check for inconsistensies:
      If(nSym1.ne.nSym2) then
       write(6,*) 'Symmetries are not equal. Stop here',nSym1, nSym2
       Call QTrace()
       Call Abend()
      Endif
      Do isym=1,nSym1
       if(nBas1(isym).gt.nBas2(isym)) then
        write(6,*) 'Second basis set must be larger than first'
        write(6,*) 'not fulfilled in sym',isym,'basis functions are',
     &  nBas1(isym),nBas2(isym)
        Call QTrace()
        Call Abend()
       Endif
      Enddo
*----------------------------------------------------------------------*
*     Build the new orbitals                                           *
*----------------------------------------------------------------------*
      ist1=0
      ist2=0
      ib1=1
      ib2=1
      Do isym=1,nsym1
       nb1=nBas1(isym)
       nb2=nBas2(isym)
       If(nb2.gt.0) then
        call expandbas(Bas1(ib1),nb1,Bas2(ib2),nb2,
     & Work(ist1+ipCMO1),Work(ist2+ipCMO2),occ1(ib1),eorb1(ib1),
     & indt1(ib1),occ2(ib2),eorb2(ib2),indt2(ib2))
        ist1=ist1+nb1**2
        ist2=ist2+nb2**2
        ib1=ib1+nb1
        ib2=ib2+nb2
       Endif
      Enddo
*----------------------------------------------------------------------*
*     Write the new orbitals in to the file EXPORB                     *
*----------------------------------------------------------------------*
C     First resort indt to standard
C
      Do i=1,56
       Indtype(i)=0
      Enddo
      ind=0
      ishift=0
      Do isym=1,nSym2
       nb2=nBas2(isym)
       if(nb2.ne.0) then
        Do ib2=1,nb2
         ind=ind+1
         Indtype(ishift+indt2(ind))=Indtype(ishift+indt2(ind))+1
        Enddo
       Endif
       ishift=ishift+7
      Enddo
C
      VecTit='Basis set expanded orbital file EXPORB'
      Lu_=60
      Call WRVEC('EXPORB',LU_,'COEI',nSym2,nBas2,nBas2,Work(ipCMO2),
     & occ2,eorb2,Indtype,VecTit)
      Write(6,*)'New orbitals have been built in file EXPORB'
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      Call GetMem('CMO1','Free','Real',ipCMO1,nTot2)
      Call GetMem('CMO2','Free','Real',ipCMO2,nTot2)
*
      ireturn=0
      return
      End
