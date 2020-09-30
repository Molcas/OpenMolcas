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
      SubRoutine InpOne
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Files_mclr.fh"
#include "glbbas_mclr.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
cnf
#include "real.fh"
#include "rctfld.fh"
      Logical Do_ESPF,First,Dff,Do_DFT,NonEq
cnf
      Character*8 Label
      Dimension idum(1)
*
      iRc=-1
      iOpt=1
      ndens2=0
      iisym=2**0
      Do iS=1,nSym
        nDens2=nDens2+Nbas(is)**2
      End Do
      Label='ONEHAM'
      Call iRdOne(iRc,iOpt,Label,1,idum,iisym)
      leng=idum(1)
      If (iRC.ne.0)  Then
         Write (6,*) 'InpOne: Error reading ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
      iisym=2**0
      iRc=-1
      iOpt=0
      Call GetMem('ONEHAM','ALLO','REAL',kint1,ndens2)
      kain1=kint1
      Call GetMem('HAM1','ALLO','REAL',itemp1,leng+10)
      Call GetMem('HAM2','ALLO','REAL',itemp2,ndens2)
      Call GetMem('HAM3','ALLO','REAL',itemp3,ndens2)
      Call RdOne(iRc,iOpt,Label,1,Work(itemp1),iisym)
      If (iRC.ne.0)  Then
         Write (6,*) 'InpOne: Error reading ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
cnf
*
*     Modify the one electron Hamiltonian for reaction
*     field and ESPF calculations
*
      Tot_Nuc_Charge=Zero
      Call GetMem('Nuclear charges','Allo','Real',ipNuc,nAtoms)
      Call Get_dArray('Effective nuclear Charge',Work(ipNuc),nAtoms)
      Do iNuc = 0, nAtoms-1
        Tot_Nuc_Charge = Tot_Nuc_Charge + Work(ipNuc+iNuc)
      End Do
      Call GetMem('Nuclear charges','Free','Real',ipNuc,nAtoms)
      Tot_El_Charge = Zero
      Do iSym = 1, nSym
         Tot_El_Charge = Tot_El_Charge
     &                 - Two*DBLE(nFro(iSym)+nIsh(iSym))
      End Do
      Tot_El_Charge = Tot_El_Charge - DBLE(nActEl)
      Tot_Charge = Tot_Nuc_Charge + Tot_El_Charge
      iCharge = Int(Tot_Charge)
      ERF = Zero
      ERFhi = Zero
      Call DecideOnESPF(Do_ESPF)
      If ( Do_ESPF .or. lRF) then
         If (lRF) Then
            Write(6,*) 'Sorry, MCLR+RF NYI'
            Call Quit_OnUserError()
         End If
*
*------ Scratch for one- and two-electron type contributions
*------ + variational density-matrix
*
         Call GetMem('Htmp','Allo','Real',iHtmp,leng)
         Call GetMem('Gtmp','Allo','Real',iGtmp,leng)
         Call dCopy_(leng,[Zero],0,Work(iHtmp),1)
         Call dCopy_(leng,[Zero],0,Work(iGtmp),1)
         Call Get_D1ao(ipD1ao,leng)
*
         NonEq=.False.
         First=.True.
         Dff=.False.
         Do_DFT=.True.
         Call Get_dScalar('PotNuc',PotNuc)
         Call DrvXV(Work(iHtmp),Work(iGtmp),Work(ipD1ao),
     &              PotNuc,leng,First,Dff,NonEq,lRF,
*
*------ Don't care about the last arguments: no (CAS-)DFT here I guess)
*
     &              'SCF',Zero,iCharge,iSpin,Work(ip_Dummy),
     &              Work(ip_Dummy),0,'1234',Do_DFT)
         Call Daxpy_(leng,One,Work(iHtmp),1,Work(itemp1),1)
*
*------ Hum, where the hell is FI (Fock Inactive) ???
*
*         Call Daxpy_(leng,One,Work(iGtmp),1,FI,1)
         Call GetMem('Gtmp','Free','Real',iGtmp,leng)
         Call GetMem('Htmp','Free','Real',iHtmp,leng)
         Call Free_Work(ipD1ao)
      End If
cnf
      ip=1
      ip2=0
      Do iS=1,nSym
        If (nBas(is).ne.0 .AND. nOrb(iS).ne.0) Then
           Call Square(work(itemp1+ip-1),
     &                   work(itemp2),
     &                   1,nBas(is),nBas(is))
           ip=ip+nBas(is)*(nBas(iS)+1)/2
           Call DGEMM_('T','N',
     &                 nOrb(iS),nBas(iS),nBas(iS),
     &                 1.0d0,Work(ipCMO+ip2),nBas(iS),
     &                 work(itemp2),nBas(iS),
     &                 0.0d0,Work(iTemp3),nOrb(iS))
           Call DGEMM_('N','N',
     &                 nOrb(is),nOrb(iS),nBas(iS),
     &                 1.0d0,Work(iTemp3),nOrb(iS),
     &                 Work(ipCMO+ip2),nBas(iS),
     &                 0.0d0,Work(kint1+ip2),nOrb(iS))
           ip2=ip2+nBas(is)**2
        End If
      End Do
      Call GetMem('HAM1','Free','REAL',itemp1,leng)
      Call GetMem('HAM2','Free','REAL',itemp2,leng)
      Call GetMem('HAM3','Free','REAL',itemp3,leng)

      Return
      End
