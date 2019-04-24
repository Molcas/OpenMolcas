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
* Copyright (C) 1990, Markus P. Fuelscher                              *
*               2013, Giovanni Li Manni                                *
************************************************************************
      Subroutine CASDFT_terms(CMO,F,FI,D1I,D1A,D1S)
*
*     This routine is a modification of SGFCIN, adapted to a CASDFT
*     implementation in which the CI step of a CASDFT calculation is
*         not corrected by DFT. DFT will play a role only in the Orbital
*         optimization step.
*     Purpose:
*     Generate the Fock-matrix for the frozen and inactive orbitals.
*     Compute also the core energy. Finally, transform the Fock-matrix
*     into the basis of the active orbitals.
*
*     M.P. Fuelscher, Lund, July 1990
*     GLM, Minneapolis,   May 2013
*
#ifdef _DMRG_
!     module dependencies
      use qcmaquis_interface_cfg
#endif
      Implicit Real*8 (A-H,O-Z)
*
#include "rasdim.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='CASDFT_Terms')
#include "rasscf.fh"
#include "WrkSpc.fh"
#include "rctfld.fh"
#include "pamint.fh"
#include "timers.fh"
#include "SysDef.fh"
*
      Dimension CMO(*) ,F(*) , FI(*) , D1I(*) , D1A(*), D1S(*)
      Character*8 Label
      Logical First, Dff, Do_DFT
      Parameter ( Zero=0.0d0 , One=1.0d0 )

      Call qEnter('CASDFT_Terms')

***********************************************************
C Local print level (if any)
***********************************************************
      IPRLEV=IPRLOC(3)
c      IPRLEV=100
      If ( IPRLEV.ge.DEBUG ) then
       Write(LF,*) 'Printing matrices in CASDFT_Terms'
       Write(LF,*)
       Write(LF,*) ' CMO in CASDFT_terms'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       ioff=1
       Do iSym = 1,nSym
         iBas = nBas(iSym)
         if(iBas.ne.0) then
           write(6,*) 'Sym =', iSym
           do i= 1,iBas
             write(6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
           end do
           iOff = iOff + (iBas*iBas)
         end if
       End Do

       Write(LF,*)
       Write(LF,*) ' D1I in AO basis in CASDFT_Terms'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       iOff=1
       Do iSym = 1,nSym
         iBas = nBas(iSym)
         call wrtmat(D1I(ioff),iBas,iBas, iBas, iBas)
         iOff = iOff + iBas*iBas
       End Do

       Write(LF,*)
       Write(LF,*) ' D1S in AO basis in CASDFT_Terms'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       iOff=1
       Do iSym = 1,nSym
         iBas = nBas(iSym)
         call wrtmat(D1S(ioff),iBas,iBas, iBas, iBas)
         iOff = iOff + iBas*iBas
       End Do


       Write(LF,*)
       Write(LF,*) ' D1A in AO basis in CASDFT_Terms'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       iOff=1
       Do iSym = 1,nSym
         iBas = nBas(iSym)
         call wrtmat(D1A(ioff),iBas,iBas, iBas, iBas)
         iOff = iOff + iBas*iBas
       End Do

      End If


*
***********************************************************
* Generate molecular charges
***********************************************************
      Call GetMem('Ovrlp','Allo','Real',iTmp0,nTot1+4)
      iRc=-1
      iOpt=2
      iComp=1
      iSyLbl=1
      Label='Mltpl  0'
      Call RdOne(iRc,iOpt,Label,iComp,Work(iTmp0),iSyLbl)
      Tot_Nuc_Charge=Work(iTmp0+nTot1+3)
      If ( iRc.ne.0 ) then
         Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
         Write(LF,*) 'Label = ',Label
         Write(LF,*) 'iRc = ',iRc
         Call QTrace
         Call Abend
      Endif
      Call GetMem('Ovrlp','Free','Real',iTmp0,nTot1+4)

      Tot_El_Charge=Zero
      Do iSym=1,nSym
         Tot_El_Charge=Tot_El_Charge
     &                -2.0D0*DBLE(nFro(iSym)+nIsh(iSym))
      End Do
      Tot_El_Charge=Tot_El_Charge-DBLE(nActEl)
      Tot_Charge=Tot_Nuc_Charge+Tot_El_Charge
*      If ( IPRLEV.ge.DEBUG ) then
*          write(6,*)
*       write(6,*) 'Total Charge :', Tot_Charge
*      end if
*
***********************************************************
* Load bare nuclei Hamiltonian
***********************************************************
      Call GetMem('Fcore','Allo','Real',iTmp1,nTot1)
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  6
      Label  = 'OneHam  '
      Call RdOne(iRc,iOpt,Label,iComp,Work(iTmp1),iSyLbl)
      If ( iRc.ne.0 ) then
         Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
         Write(LF,*) 'Label = ',Label
         Write(LF,*) 'iRc = ',iRc
         Call QTrace
         Call Abend
      Endif
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' OneHam in AO basis in CASDFT_Terms'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=0
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',Work(iTmp1+iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End IF
*
***********************************************************
* Load the nuclear repulsion energy
***********************************************************
      Call Get_dScalar('PotNuc',potNuc)
c      write(6,*)
c      write(6,*) 'PotNuc in casdft_terms.f:', PotNuc

*
***********************************************************
* Generate total density
***********************************************************
*
*      If ( IPRLEV.ge.DEBUG ) then
*       Call GetMem('D_MAT','Allo','Real',iTmp31,nBas*nBas)
*          call FZERO(Work(iTmp31),nBas*nBas)
*          Call Daxpy_(nBas*nBas,1.0D0,D1I,1,Work(iTmp31),1)
*          Call Daxpy_(nBas*nBas,1.0D0,D1A,1,Work(iTmp31),1)
*      Write(LF,*)
*       Write(LF,*) ' DMAT not folded in AO basis in CASDFT_Terms'
*       Write(LF,*) ' ---------------------'
*       Write(LF,*)
*          call wrtmat(Work(iTmp31),nBas,nBas, nBas, nBas)
*       Call GetMem('D_MAT','Free','Real',iTmp31,nBas*nBas)
*      End IF

      Call GetMem('DtmpI','Allo','Real',iTmp3,nTot1)
      Call GetMem('DtmpA','Allo','Real',iTmp4,nTot1)
      Call Fold(nSym,nBas,D1I,Work(iTmp3))
      Call Fold(nSym,nBas,D1A,Work(iTmp4))
      Call Daxpy_(nTot1,1.0D0,Work(iTmp4),1,Work(iTmp3),1)
      Call Put_D1ao(Work(iTmp3),nTot1)
      call xflush(6)
***********************************************************
* Generate spin-density
***********************************************************
        Call GetMem('DtmpS','Allo','Real',iTmp7,nTot1)
        Call Fold(nSym,nBas,D1S,Work(iTmp7))
        Call Put_D1Sao(Work(iTmp7),nTot1)
        Call GetMem('DtmpS','Free','Real',iTmp7,nTot1)
*
***********************************************************
* One- and two-electron type contributions
***********************************************************
*
      Call GetMem('htmp','Allo','Real',iTmp5,nTot1)
      Call GetMem('gtmp','Allo','Real',iTmp6,nTot1)
      Call dcopy_(nTot1,[0.0d0],0,Work(iTmp5),1)
      Call dcopy_(nTot1,[0.0d0],0,Work(iTmp6),1)
*
      First=.True.
      Dff=.False.
      Do_DFT=.True.

      Call Put_iArray('nFro',nFro,nSym)
      Call Put_iArray('nAsh',nAsh,nSym)
      Call Put_iArray('nIsh',nIsh,nSym)

      iCharge=Int(Tot_Charge)
c iTmp5 and iTmp6 are not updated in DrvXV...
      Call DrvXV(Work(iTmp5),Work(iTmp6),Work(iTmp3),
     &             PotNuc,nTot1,First,Dff,NonEq,lRF,
     &             KSDFT_TEMP,ExFac,iCharge,iSpin,D1I,D1A,
     &             nTot1,DFTFOCK,Do_DFT)

      Call Daxpy_(nTot1,1.0d0,Work(iTmp5),1,Work(iTmp1),1)
      Call Daxpy_(nTot1,1.0d0,Work(iTmp6),1,FI,1)

      Call GetMem('gtmp','Free','Real',iTmp6,nTot1)
      Call GetMem('htmp','Free','Real',iTmp5,nTot1)
      Call GetMem('DtmpA','Free','Real',iTmp4,nTot1)
      Call GetMem('DtmpI','Free','Real',iTmp3,nTot1)

***********************************************************
*     Compute energy contributions
***********************************************************
      Call GetMem('DoneI','Allo','Real',iTmp2,nTot1)

      Call Fold(nSym,nBas,D1I,Work(iTmp2))
*
      Eone = dDot_(nTot1,Work(iTmp2),1,Work(iTmp1),1)
      Call Get_dScalar('PotNuc',PotNuc_Ref)
      Eone = Eone + (PotNuc-PotNuc_Ref)
      Etwo = dDot_(nTot1,Work(iTmp2),1,FI,1)
      Call GetMem('DoneI','Free','Real',iTmp2,nTot1)
      EMY  = PotNuc_Ref+Eone+0.5d0*Etwo

      CASDFT_Funct = 0.0D0
      Call Get_dScalar('CASDFT energy',CASDFT_Funct)
      If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*) ' Nuclear repulsion energy :',PotNuc
*         Write(LF,*) ' Nuclear repulsion energy Ref :',PotNuc_Ref
         Write(LF,*) ' One-electron core energy :',Eone
         Write(LF,*) ' Two-electron core energy :',Etwo
         Write(LF,*) ' Total core energy        :',EMY
         Write(LF,*) ' CASDFT Energy            :',CASDFT_Funct
      End If

***********************************************************
* Printing matrices
***********************************************************
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FI matrix in CASDFT_Terms only 2-electron terms'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',FI(iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If

      Call DaXpY_(nTot1,One,Work(iTmp1),1,FI,1)

      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' Inactive Fock matrix in AO basis in CASDFT_terms'
        write(LF,*) '(it already contains OneHam and TwoEl contrib.)'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',FI(iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If

      Call GetMem('Fcore','Free','Real',iTmp1,nTot1)

***********************************************************
*     Transform FI to active orbital basis and move it over to F.
*     Remove also the symmetry blocking.
***********************************************************
***********************************************************
* Shall I add here the DFT contribution? Maybe not yet!
* I am commenting off... if needed we can always re-activate.
***********************************************************
      MXNB=0
      MXNA=0
      DO ISYM=1,NSYM
        MXNB=MAX(MXNB,NBAS(ISYM))
        MXNA=MAX(MXNA,NASH(ISYM))
      END DO
      CALL GETMEM('XXX0','ALLO','REAL',LX0,NTOT1)
      CALL GETMEM('XXX1','ALLO','REAL',LX1,NTOT1)
      CALL GETMEM('XXX2','ALLO','REAL',LX2,MXNB*MXNB)
      CALL GETMEM('XXX3','ALLO','REAL',LX3,MXNB*MXNA)
      CALL dcopy_(NTOT1,FI,1,WORK(LX1),1)
*      Call Get_dExcdRa(ipTmpFckI,nTmpFck)
*      CALL DaXpY_(NTOT1,1.0D0,Work(ipTmpFckI),1,WORK(LX1),1)
*     If ( IPRLEV.ge.DEBUG ) then
*         Write(LF,*)
*         Write(LF,*) ' Exchange Corr. in AO basis in CASDFT_Terms'
*         Write(LF,*) ' ---------------------'
*         Write(LF,*)
*         iOff=1
*         Do iSym = 1,nSym
*           iBas = nBas(iSym)
*           Call TriPrt(' ','(5G17.11)',Work(ipTmpFckI+ioff-1),iBas)
*           iOff = iOff + (iBas*iBas+iBas)/2
*         End Do
*     End If
*      Call Free_Work(ipTmpFckI)
*      If ( IPRLEV.ge.DEBUG ) then
*        Write(LF,*)
*        Write(LF,*) ' Modified FI in AO basis in CASDFT_Terms'
*        Write(LF,*) ' ---------------------'
*        Write(LF,*)
*        iOff=1
*        Do iSym = 1,nSym
*          iBas = nBas(iSym)
*          Call TriPrt(' ','(5G17.11)',Work(LX1+ioff-1),iBas)
*          iOff = iOff + (iBas*iBas+iBas)/2
*        End Do
*      End If
*
      CALL MOTRAC(CMO,WORK(LX1),WORK(LX2),WORK(LX3))
      CALL GETMEM('XXX3','FREE','REAL',LX3,MXNB*MXNA)
      CALL GETMEM('XXX2','FREE','REAL',LX2,MXNB*MXNB)
      CALL dcopy_(NACPAR,[ZERO],0,F,1)
      NTU=0
      ITU=0
      IADD=0

      IF (NACTEL.NE.0) THEN
         EMYN=EMY/DBLE(NACTEL)
      ELSE
         EMYN=0.0d0
      ENDIF
      DO NST=1,NSYM
        NAT=NASH(NST)
        IF(NAT.NE.0) THEN
          DO NT=1,NAT
            NTU=NTU+IADD
            DO NU=1,NT
              NTU=NTU+1
              ITU=ITU+1
              F(NTU)=WORK(LX1-1+ITU)
              IF(NT.EQ.NU) F(NTU)=F(NTU)+EMYN
              WORK(LX0-1+ITU) = F(NTU)
            END DO
          END DO
          IADD=IADD+NAT
        ENDIF
      END DO
#ifdef _DMRG_
      if(.not.doDMRG)then
        CALL CP_ONE_INT(WORK(LX0),ITU)
      end if
#endif
      CALL GETMEM('XXX1','FREE','REAL',LX1,NTOT1)
      CALL GETMEM('XXX0','FREE','REAL',LX0,NTOT1)

*     print h0
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*)
        Write(LF,*)
        Write(LF,*)
        Write(LF,*) ' Fock matrix in MO basis, h0, in CASDFT_TERMS'
        Write(LF,*) ' ------------'
        Call TriPrt(' ',' ',F,NAC)
      End If

      Call qExit('CASDFT_terms')

      Return
      End
