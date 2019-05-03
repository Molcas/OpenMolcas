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
* Copyright (C) Anders Ohrn                                            *
************************************************************************
*  Qfread
*
*> @brief
*>   Read in all data that comes from external Molcas routines and prepare various quantities, such as MME
*> @author A. Ohrn
*>
*> @details
*> This subroutine handles the interaction with the rest of Molcas.
*> Here orbitals and various matrices are stored and to some extent
*> modified to our purpose. We call on this subroutine even when
*> we are running classical stuff. The reason for this is that we
*> wish to collect some numbers, but really, we could easily have
*> constructed thing differently so that this subroutine only would
*> be called when quantum-classical stuff is running. For RASSI
*> implementation, we also read in and transform the transition density
*> matrix. At the end, we also call the routines that make the
*> multicenter multipole expansion.
*>
*> @note
*> Seward is mandatory for both SCF and RASSI. For SCF also,
*> Motra, Averd; for RASSI also, RASSCF and RASSI.
*>
*> @param[out] nAtomsCC Atoms in solvent
*> @param[out] Coord    Unique coordinates of the atoms in the molecule in the QM-region
*> @param[out] nBas     Number of basis functions in QM-region
*> @param[out] nBasCC   Like nBas but for a solvent molecule
*> @param[out] nCnC_C   Like nCnC, but for solvent
*> @param[out] nntyp    Number of basis-function types.
*> @param[out] nOcc     The total number of basis functions that belong to a certain basis-function type.
*> @param[out] natyp    Number of atoms of the i:th basis-function type
************************************************************************
*******JoseMEP the last three variables are included to the MEP calculation
      Subroutine Qfread(iQ_Atoms,nAtomsCC,Coord,nBas,nBasCC,nCnC_C
     &,nOcc,natyp,nntyp)
      Implicit Real*8 (a-h,o-z)

*-----------------------------------------------------------------------*
* Variables                                                             *
*-----------------------------------------------------------------------*
#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "qmcom.fh"
#include "qm1.fh"
#include "qm2.fh"
#include "integral.fh"
#include "numbers.fh"
#include "WrkSpc.fh"
#include "tratoc.fh"
#include "warnings.fh"

      Parameter (IndMax=nTraBuf) !nTraBuf definied in tratoc.fh

      Dimension Coord(MxAt*3),Chge(MxAt),CoordCC(3*3),ChgeCC(3)
      Dimension Cmo(MxBas**2),Cmo_S(MxBas**2),Occu(MxBas),Dummy(MxBas)

      Dimension nSh(MxAt),nfSh(MxAt,MxAngqNr),nCnC(MxBas),nCnC_C(MxBasC)
      Dimension nOcc(MxBas),natyp(MxAt),natypC(MxAt),iDumm(MxBas)
      Dimension nBas(MxSym),nBasCC(1),iCon(MxAt,MxPrCon)
      Dimension iC_icon(MxAt,MxPrCon)

      Character Line*120,BlLine*120,Title*100,OrbName*100,WhatGet*10
      Character StLine*120
      Dimension iDummy(1)

*-----------------------------------------------------------------------*
* Enter                                                                 *
*-----------------------------------------------------------------------*
      Call QEnter('QFREAD')

*------------------------------------------------------------------------*
* Print the joblabel. It is obtained in get_input.f                      *
*------------------------------------------------------------------------*
      If(ATitle) then
        lLine=Len(Line)
        Do 9990, i=1,lLine
          BlLine(i:i)=' '
          StLine(i:i)='*'
9990    Continue
        Write(6,*)
        Do 9991, i=1,6
          Line=BlLine
          If(i.eq.1.or.i.eq.6) Line=StLine
          If(i.eq.3) Line='Project:'
          If(i.eq.4) Write(Line,'(A72)')JobLab
          Call Center(Line)
          Write(6,*)'*'//Line//'*'
9991    Continue
        Write(6,*)
      Endif
      Write(6,*)'Auxiliary data being read and pre-processed.'
*----------------------------------------------------------------------*
* Collect some data from RUNFILE about the QM-region molecule.         *
*----------------------------------------------------------------------*
      Call Get_iScalar('nSym',nSym)
      If(nSym.ne.1) then !A minor restriction, no symmetry allowed,
                          !i.e. nSym=1.
        Write(6,*)
        Write(6,*)' QmStat does not run with symmetry!'
        Write(6,*)' The perturbation from the solvent breaks all symmet'
     &//'ry.'
        Call Quit(_RC_GENERAL_ERROR_)
      Endif
      Call Get_iScalar('Unique atoms',iQ_Atoms)
      If(iQ_Atoms.gt.MxAt) then
        Write(6,*)
        Write(6,*)'Maximum number of atoms exceeded. Increase MxAt in '
     &//'maxi.fh in QmStat source directory.'
        Call Quit(_RC_GENERAL_ERROR_)
      Endif
      Call Get_dArray('Nuclear charge',Chge,iQ_Atoms)
      Call Get_dArray('Unique Coordinates',Coord,3*iQ_Atoms)
      Call Get_dArray('Center of Mass',CT,3)
      Call Get_iArray('nBas',nBas,nSym)
      If(nBas(1).gt.MxBas) then
        Write(6,*)
        Write(6,*)'Maximum number of basis functions exceeded. Increase'
     &//' MxBas in maxi.fh in QmStat source directory.'
        Call Quit(_RC_GENERAL_ERROR_)
      Endif

*------------------------------------------------------------------------*
* Print elementary information about molecule.                           *
*------------------------------------------------------------------------*
      Write(6,*)
      Write(6,*)'     ------------------------------'
      Write(6,*)'     |       QM-region data       |'
      Write(6,*)'     ------------------------------'
      Write(6,*)
      Write(6,'(A,15X,I5)')'      Number of basis functions:'
     &                    ,(nBas(i),i=1,nSym)
      Write(6,'(A,3F10.6)')'      Centre of mass =',(CT(kk),kk=1,3)
      Call PrCoor

*----------------------------------------------------------------------*
* If the Qmtype is SCF we now want the orbitals. On the other hand if  *
* we are running RASSI, another route must be taken, hence here we     *
* inquire which QM-method that is used.                                *
*----------------------------------------------------------------------*
      If(QmType(1:3).eq.'SCF') then
*----------------------------------------------------------------------*
*      SSS    CCC  FFFFFF                                              *
*     SS    CC     FF                                                  *
*      SS  CC      FFFF                                                *
*       SS  CC     FF                                                  *
*     SSS    CCC   FF                                                  *
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
* Print information about orbitals and store coefficients in new       *
* variable.                                                            *
*----------------------------------------------------------------------*
        iLu=15
        iLu=IsFreeUnit(iLu)
        Write(OrbName,'(A)')'AVEORB'
        Write(WhatGet,'(A)')'CO'
        iWarn=1
        Call RdVec(OrbName,iLu,WhatGet,nSym,nBas,nBas,Cmo,Occu,Dummy
     &,iDumm,Title,iWarn,iErr)
        If(iErr.ne.0) then
          Write(6,*)
          Write(6,*)'Error when reading AVEORB'
          Write(6,*)
          Call Quit(_RC_IO_ERROR_READ_)
        Endif
        kaunter=0
        Call GetMem('OrbCoeffQ','Allo','Real',iV1,iOrb(1)*nBas(1))
        Do 102,j=1,iOrb(1)
          Do 103,k=1,nBas(1)
            Work(iV1+kaunter)=Cmo(k+(j-1)*nBas(1))
            kaunter=kaunter+1
103       Continue
102     Continue
*        Write(6,'(A,I4)')'      Number of Orbitals:',iOrb(1)
        Call Get_dScalar('PotNuc',PotNuc)

      Elseif(QmType(1:4).eq.'RASS') then
*----------------------------------------------------------------------*
*     RRRR      AA      SSS   SSS   II                                 *
*     RR  R    A  A    SS    SS     II                                 *
*     RR R    AAAAAA    SS    SS    II                                 *
*     RRR     AA  AA     SS    SS   II                                 *
*     RR R    AA  AA   SSS   SSS    II                                 *
*----------------------------------------------------------------------*
        Call TdmTrans(nBas)
      Endif
      Write(6,*)
      Write(6,*)
*----------------------------------------------------------------------*
* Compute various information about system. This we use for computing  *
* integrals later. GiveMeInfo collects stuff from seward, somtime with *
* some recomputations.                                                 *
*----------------------------------------------------------------------*
      Call GiveMeInfo(nBas(1),nntyp,natyp,BasOri,Icon,nPrimus,nBA_Q
     &,nCBoA_Q,nBonA_Q,ipE,ipC,nsh,nfsh,nSize,iPrint,MxAt,MxPrCon,MxBas
     &,MxAngqNr,ipACC,nACCSizeQ)
      iBas=0
      iAtom=0
      kold=1
      iold=1
      indold=0
      diff=0
      Do 149, i=1,nntyp
        nOcc(i)=0
149   Continue
      Do 150, i=1,nntyp
        na=natyp(i)
        Do 151, j=1,na
          ind=0
          jnd=0
          iAtom=iAtom+1
          ChaNuc(iAtom)=Chge(iAtom)
          info_atom(iAtom)=int(Chge(iAtom))
          nShj=nSh(i)
          Do 152, k=1,nShj
            nnaa=nfsh(i,k)
            Do 153, l=1,nnaa
              ibas=ibas+1
              indold=indold+1
              nOcc(i)=nOcc(i)+2*k-1
              nCnC(ibas)=nnaa
              ind=ind+1
              iQang(ibas)=k
              icont=Icon(i,ind)
              iCharOnBasQ(ibas)=int(Chge(iAtom))
              Do 1531, ix=1,2*k-1  !Here we construct an array of
                If(k.ne.kold) then !indeces which is used to put right
                  If(i.ne.iold) then !AO-overlap in right matrix pos.
                    Indold=Indold+nfsh(iold,kold)*(2*kold-2)
                    iold=i
                  Else
                    Indold=Indold+nfsh(i,kold)*(2*kold-2)
                  Endif
                  kold=k
                Endif
                iWoGehenQ(ibas,ix)=indold+nnaa*(ix-1)
1531          Continue
              Do 154, m=1,icont
                jnd=jnd+1
                alfa(ibas,m)=Work(ipE+i-1+MxAt*(jnd-1))
                cont(ibas,m)=Work(ipC+i-1+MxAt*(jnd-1))
154           Continue
153         Continue
152       Continue
151     Continue
150   Continue
      Kmax=ibas
      call dcopy_(nACCSizeQ,Work(ipACC),iONE,Trans,iONE)
      Call GetMem('AccTransa','Free','Real',ipACC,nACCSizeQ)
      Call GetMem('Exponents','Free','Real',ipE,nSize*MxAt) !Now we
      Call GetMem('ContrCoef','Free','Real',ipC,nSize*MxAt) !do not
                                             !need them, so deallocate.

*------------------------------------------------------------------------*
* Obtain and print information about solvent. This requires a renaming   *
* of the runfile.                                                        *
*------------------------------------------------------------------------*
      Call NameRun('WRUNFIL')
      Call Get_iScalar('nSym',nSymCC)
      If(nSymCC.ne.1) then
        Write(6,*)
        Write(6,*)' QmStat does not run with symmetry!'
        Call Quit(_RC_GENERAL_ERROR_)
      Endif
      Call Get_iScalar('Unique atoms',nAtomsCC)
      If(nAtomsCC.ne.3) then
        Write(6,*)
        Write(6,*)'Now now... what strange solvent molecule do you try'
     &//' to feed QmStat with?'
        Call Quit(_RC_GENERAL_ERROR_)
      Endif
      Call Get_dArray('Nuclear charge',ChgeCC,nAtomsCC)
      Call Get_dArray('Unique Coordinates',CoordCC,3*nAtomsCC)
      Call Get_iArray('nBas',nBasCC,nSymCC)
      If(nBasCC(1).gt.MxBasC) then
        Write(6,*)
        Write(6,*)'Number of solvent molecule basis functions exceeded.'
     &//' Increase MxBasC in maxi.fh in QmStat source directory.'
        Call Quit(_RC_GENERAL_ERROR_)
      Endif
      Write(6,*)
      Write(6,*)'     ------------------------------'
      Write(6,*)'     |   Solvent molecule data     |'
      Write(6,*)'     ------------------------------'
      Write(6,*)
      Write(6,'(A,15X,I5)')'      Number of basis functions:'
     &                    ,(nBasCC(i),i=1,nSymCC)
      Call PrCoor

*
*--- Collect information about the solvent orbitals.
*
      iLu=16
      iLu=IsFreeUnit(iLu)
      Write(OrbName,'(A)')'SOLORB'
      Write(WhatGet,'(A)')'CE'
      iWarn=1
      Call GetMem('OrbitalEnergy','Allo','Real',iOe,Sum(nBasCC))
      Call RdVec(OrbName,iLu,WhatGet,nSymCC,nBasCC,nBasCC,Cmo_S
     &,Dummy,Work(iOe),iDummy,Title,iWarn,iErr)
      Do 22, i=1,iOrb(2)
        c_orbene(i)=Work(iOe+i-1)
22    Continue
      Call GetMem('OrbitalEnergy','Free','Real',iOe,Sum(nBasCC))

*
*--- We should not need two solvent orbital vectors, so this should
*    be removed when the orbital rotation routine is fixed.
*
      Do 202,j=1,iOrb(2)
        Do 203,k=1,nBasCC(1)
          V3(k,j)=Cmo_S(k+(j-1)*nBasCC(1))
203     Continue
202   Continue
*      Write(6,'(A,I4)')'      Number of Orbitals:',iOrb(2)
      Write(6,*)
      Write(6,*)
*-----------------------------------------------------------------------*
* And now basis set information.                                        *
*-----------------------------------------------------------------------*
      Call GiveMeInfo(nBasCC(1),nntypC,natypC,SavOri,iC_Icon,mPrimus
     &,nBA_C,nCBoA_C,nBonA_C,ipE_C,ipC_C,nsh,nfsh,nSize,iPrint,MxAt
     &,MxPrCon,MxBas,MxAngqNr,ipACC,nACCSizeC)
      iBas=0
      iAtom=0
      mbash=0
      kold=1
      iold=1
      indold=0
      Do 250, i=1,nntypC  !Like the corresponding thing above for
        na=natypC(i)       !the QM-region.
        Do 251, j=1,na
          ind=0
          jnd=0
          nShj=nSh(i)
          iAtom=iAtom+1
          Do 252, k=1,nShj
            nnaa=nfsh(i,k)
            Do 253, l=1,nnaa
              ibas=ibas+1
              indold=indold+1
              nCnC_C(ibas)=nnaa
              ind=ind+1
              icont=iC_Icon(i,ind)
              iqn(ibas)=k
              iCharOnBasC(ibas)=int(ChgeCC(iAtom))
              Do 2531, ix=1,2*k-1  !Here we construct an array of
                If(k.ne.kold) then !indeces which is used to put right
                  If(i.ne.iold) then !AO-overlap in right matrix pos.
                    Indold=Indold+nfsh(iold,kold)*(2*kold-2)
                    iold=i
                  Else
                    Indold=Indold+nfsh(i,kold)*(2*kold-2)
                  Endif
                  kold=k
                Endif
                iWoGehenC(ibas,ix)=indold+nnaa*(ix-1)
2531          Continue
              Do 254, m=1,icont
                jnd=jnd+1
                beta(ibas,m)=Work(ipE_C+i-1+MxAt*(jnd-1))
                dont(ibas,m)=Work(ipC_C+i-1+MxAt*(jnd-1))
254           Continue
253         Continue
252       Continue
251     Continue
250   Continue
      Lmax=ibas
      If(nACCSizeC.gt.nACCSizeQ) then
        call dcopy_(nACCSizeC,Work(ipACC),iONE,Trans,iONE)
      Endif
      Call GetMem('AccTransa','Free','Real',ipACC,nACCSizeC)
      Call GetMem('Exponents','Free','Real',ipE_C,nSize*MxAt) !Now we
      Call GetMem('ContrCoef','Free','Real',ipC_C,nSize*MxAt) !do not
                                             !need them, so deallocate.
*
*----------------------------------------------------------------------*
* A small test to see if max-limits are violated.                      *
*----------------------------------------------------------------------*
      If(Kmax.gt.MxBB.or.Lmax.gt.MxBB) then
        Write(6,*)
        Write(6,*)'ERROR! MxBB too small!'
        Call Quit(_RC_INTERNAL_ERROR_)
      Endif
*----------------------------------------------------------------------*
* The multipoles and the Hamiltonian matrix are radically different    *
* between the QM-method alternatives, so once more an inquire.         *
*----------------------------------------------------------------------*
      Call NameRun('RUNFILE')
      If(QmType(1:3).eq.'SCF') then
        Call ScfHandM(Cmo,nBas,iQ_Atoms,nOcc,natyp,nntyp,Occu)
      Elseif(QmType(1:4).eq.'RASS') then
        Call RassiHandM(nBas,iQ_Atoms,nOcc,natyp,nntyp)
      Endif
*----------------------------------------------------------------------*
* Here is the end.                                                     *
*----------------------------------------------------------------------*
      Call QExit('QFREAD')
      Return
      End
