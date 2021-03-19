!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2011, Francesco Aquilante                              *
!***********************************************************************

      Subroutine Charge_GRID_IT(nSym,nBas,CMO,nCMO,OCCN,iDoIt,          &
     &                          long_prt)

!*********************************************************************
!
!  Author : F. Aquilante
!
!
!   Purpose: Compute Mulliken charges for each MO separately.
!            The analysis is performed ONLY for the occupied MOs
!            specified in GRID_IT input (this info is stored in iDoIt)
!
!   Note:  this functionality was requested by some Turbomole users
!          recently converted to MolCas. Its scientific value is not
!          too high in my opinion, therefore this subroutine is simply
!          a hack of the existing CHARGE_ and thus infinitely far from
!          efficient coding.
!
!                                            Toulouse, 28 Nov 2011
!
!*********************************************************************

      Implicit Real*8 (a-h,o-z)
      Integer nSym, nBas(nSym), nCMO, iDoIt(*)
      Real*8  CMO(nCMO), OCCN(*)
      Logical long_prt
#include "Molcas.fh"
#include "WrkSpc.fh"
      CHARACTER*(LENIN8) NAME(MxBas)


      MXTYP=0
      nTot1=0
      Do iSym = 1, nSym
         MxTyp=MxTyp+nBas(iSym)
         nTot1=nTot1+nBas(iSym)*(nBas(iSym)+1)/2
      End Do
      Call Get_cArray('Unique Basis Names',Name,(LENIN8)*MxTyp)
      Call GetMem('XOCC','ALLO','REAL',ipXocc,MXTYP)
      Call Get_iScalar('Unique atoms',nNUC)
      Call GetMem('QQ','ALLO','REAL',ipQQ,MXTYP*nNuc)
      Call GetMem('Ovrlp','Allo','Real',ipS,nTot1)
      iRc=-1
      iOpt=2
      iComp=1
      iSyLbl=1
      Call RdOne(iRc,iOpt,'Mltpl  0',iComp,Work(ipS),iSyLbl)
      If ( iRc.ne.0 ) then
         Write(6,*) 'charge_grid_it: iRc from Call RdOne not 0'
!         Write(6,*) 'Label = ',Label
         Write(6,*) 'iRc = ',iRc
         Call Abend
      Endif
      Write (6,*)
      Write (6,*)
      Write (6,*)
      Write (6,'(A)')       '         **************************'
      Call CollapseOutput(1,'       Charges per occupied MO ')
      Write (6,'(A)')       '         **************************'
      Write (6,*)
      Write (6,*)
      Write (6,*)

      Call FZero(Work(ipXocc),MXTYP)

      iCase=2
      jOcc=1
      Do iSym=1,nSym
         Do iOrb=1,nBas(iSym)

            If(IdoIt(jOcc).eq.1 .and. OCCN(jOcc).gt.0.0d0) Then

              Write (6,'(A,I4,A,I1,A,F6.4)')'          MO:',iOrb,       &
     &                                '      Symm.: ',iSym,             &
     &                                '      Occ. No.: ',OCCN(jOcc)

              lOcc=ipXocc+jOcc-1
              Work(lOcc)=OCCN(jOcc)

              Call FZero(Work(ipQQ),MxTYP*nNuc)
              Call One_CHARGE(NSYM,NBAS,Name,CMO,Work(ipXocc),Work(ipS),&
     &                        iCase,long_prt,                           &
     &                        MXTYP,Work(ipQQ),nNuc)
              Work(lOcc)=0.0d0
            EndIf

            jOcc=jOcc+1
         End Do
      End Do

      Call GetMem('XOCC','FRee','REAL',ipXocc,MXTYP)
      Call GetMem('Ovrlp','Free','Real',ipS,nTot1)
      Call GetMem('QQ','FREE','REAL',ipQQ,MXTYP*nNuc)
!
      Return
      End
