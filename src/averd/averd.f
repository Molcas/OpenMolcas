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
      Subroutine Averd(ireturn)
      Implicit Real*8 (a-h,o-z)
*
*
*-- Compute average density and corresponding natural orbitals. Two
*   possibilities exists, either construct average from input orbitals,
*   or from density matrices.
*
*   Author: Anders Ohrn.
*
*

*
*-- Include.
*
#include "mxdm.fh"
#include "real.fh"
#include "mxave.fh"
#include "WrkSpc.fh"

*
*-- Allocate.
*
      Dimension Wset(MxSets),nBas(MxSym)
      Logical PrOcc,PrEne,DensityBased
      Character Title*72, Fname*7,OLabel*10,Titorb*40,OrbFile*128
      Character BsLbl*4000,PLab*3
      Dimension Dummy(1),iDummy(7,8)


*
*-- Banner.
*
      ireturn=99

*
*-- Define defaults and initialize.
*
      Call Init_ave(Title,iPrint,Wset,Wsum,PrOcc,PrEne,DensityBased
     &         ,ThrOcc,Dummy(1),iDummy(1,1))

*
*-- Read input.
*
      Call Get_Averd_input(Title,Wset,iPrint,Nset,DensityBased,ThrOcc)

*
*-- Read some information from RUNFILE.
*
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      itBas=0
      Do 31, iSym=1,nSym
        itBas=itBas+nBas(isym)
31    Continue
      Call Get_cArray('Unique Basis Names',BsLbl,(LENIN8)*itBas)

*
*-- Some dimensions.
*
      lsmat=0
      ntot=0
      ntot2=0
      Do 70 i=1,nSym
        lsmat=lsmat+(nBas(i)*(nBas(i)+1))/2
        ntot=ntot+nBas(i)
        ntot2=ntot2+nBas(i)**2
70    Continue

*
*-- Read AO-basis overlap matrix.
*
      Call GetMem('Overlap','Allo','Real',iS,lsmat+4)
      OLabel='Mltpl  0'
      irc=0
      iopt=6
      icomp=1
      isyml=1
      Call RdOne(irc,iopt,OLabel,icomp,Work(iS),isyml)
      If(iprint.ge.99) then
        ind=0
        Do 72, iSym=1,nSym
          Call TriPrt('Overlap Matrix',' ',Work(iS+ind),nBas(iSym))
          ind=ind+nBas(iSym)*(nBas(iSym)+1)/2
72      Continue
      Endif

*
*-- Normalize weights.
*
      Do 80, iset=1,mxsets
        Wsum=Wsum+wset(iset)
80    Continue
      Do 81, iset=1,Nset
        Wset(iset)=Wset(iset)/Wsum
81    Continue

*
*-- Print some Bla Bla...
*
      If(iprint.ge.2) then
        Call Print_Input(Title,nSym,nBas,wSet,nSet)
      Endif

*
*-- Do the dirty work. Different paths for orbital- and density-based
*   averageing.
*

      Call GetMem('Density','Allo','Real',iDao,ntot2)
      call dcopy_(ntot2,[Zero],0,Work(iDao),1)
      If(.not.DensityBased) then
        Luinp=10
        Call GetMem('Orbitals','Allo','Real',iCMO,ntot2)
        Call GetMem('Occ','Allo','Real',iOcc,ntot)
        Do 90, iset=1,Nset
          Fname='NAT001'
          Write(Fname(4:6),'(I3.3)') iset
*------- Read orbital coefficients and occupation numbers.
          Call RdVec(Fname,Luinp,'CO',Nsym,nBas,nBas,Work(iCMO)
     &           ,Work(iOcc),Dummy,iDummy,Titorb,0,iErr)
          iC=0
          iO=0
          iD=0
*------- Up-date average density matrix.
          Do 93, isym=1,nSym
            kaunter=0
            Do 931, i=1,nBas(iSym)
              Do 932, j=1,nBas(iSym)
                Do 933, k=1,nBas(iSym)
                  Work(iDao+iD+kaunter)=Work(iDao+iD+kaunter)
     &  +Wset(iSet)*Work(iOcc+iO+k-1)*Work(iCMO+iC+i+(k-1)*nBas(iSym)-1)
     &                               *Work(iCMO+iC+j+(k-1)*nBas(iSym)-1)
933             Continue
                kaunter=kaunter+1
932           Continue
931         Continue
            iC=iC+nBas(isym)**2
            iD=iD+nBas(isym)**2
            iO=iO+nBas(isym)
93        Continue
*------- Print print print.
          If(iPrint.ge.5) then
            ThrO=1d-5
            Call Primo(Titorb,PrOcc,PrEne,ThrO,Dummy(1),nSym,nBas,nBas
     &                ,BsLbl,Dummy,Work(iOcc),Work(iCMO),-1)
          Endif
90      Continue
        Call GetMem('Orbitals','Free','Real',iCMO,ntot2)
        Call GetMem('Occ','Free','Real',iOcc,ntot)
      Else
        Call GetMem('DensityT','Allo','Real',iDtemp,lsmat)
        call dcopy_(lsmat,[Zero],0,Work(iDtemp),1)
        Do 95, iset=1,Nset
          Fname='RUN001'
          Write(Fname(4:6),'(I3.3)') iset
          Call NameRun(Fname)
*------- Collect density from runfile.
          Call Get_D1ao(ip_Dtmp,nDens)
          Call DaxPy_(lsmat,Wset(iset),Work(ip_Dtmp),1,Work(iDtemp),1)
          Call Free_Work(ip_Dtmp)
95      Continue
*----- Square the density matrix.
        iDt=0
        iDs=0
        Do 97, iSym=1,nSym
          nB=nBas(iSym)
          Call Dsq(Work(iDtemp+iDt),Work(iDao+iDs),1,nB,nB)
          iDt=iDt+nB*(nB+1)/2
          iDs=iDs+nB**2
97      Continue
        Call GetMem('DensityT','Free','Real',iDtemp,lsmat)
      Endif

*
*-- With the average density in store, lets orthogonalize (canonical),
*   then diagonalize to get natural orbitals.
*
      indT=0
      indS=0
      indB=0
      Call GetMem('NatOrbAcc','Allo','Real',iOrbs,ntot2)
      Call GetMem('NatOccAcc','Allo','Real',iOccs,ntot)
      Do 201, iSym=1,nSym
        nBT=nBas(iSym)*(nBas(iSym)+1)/2
        nBS=nBas(iSym)**2
        Call GetMem('EigV','Allo','Real',iVecs,nBS)
        Call GetMem('St','Allo','Real',iSt,nBT)
        Call GetMem('Si','Allo','Real',iSi,nBT)
        Call GetMem('Ss','Allo','Real',iSs,nBS)
        Call GetMem('Sp','Allo','Real',iSp,nBS)
        Call GetMem('AUX','Allo','Real',iAUX,nBS)
        Call GetMem('TransS','Allo','Real',iTrans,nBS)
        Call GetMem('TransSi','Allo','Real',iTrani,nBS)
        Call GetMem('OrthoDensS','Allo','Real',iOrtoD,nBS)
        Call GetMem('OrthoDensT','Allo','Real',iOrtoDt,nBT)
        Call GetMem('Occs','Allo','Real',iOccNat,nBas(iSym))
        call dcopy_(nBT,[Zero],0,Work(iSt),1)
        call dcopy_(nBT,[Zero],0,Work(iSi),1)
        kaunter=0
        Do 2011, iB1=1,nBas(iSym)
          Do 2012, iB2=1,nBas(iSym)
            Work(iVecs+kaunter)=Zero
            If(iB1.eq.iB2)Work(iVecs+kaunter)=One
            kaunter=kaunter+1
2012      Continue
2011    Continue
        Call Jacob(Work(iS+indT),Work(iVecs),nBas(iSym),nBas(iSym))
        Do 205, i=1,nBas(iSym)
          Sqroot=sqrt(Work(iS+indT+i*(i+1)/2-1))
          Work(iSt+i*(i+1)/2-1)=Sqroot
          Work(iSi+i*(i+1)/2-1)=One/Sqroot
205     Continue
        Call Square(Work(iSt),Work(iSs),1,nBas(iSym),nBas(iSym))
        Call Square(Work(iSi),Work(iSp),1,nBas(iSym),nBas(iSym))
        Call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One
     &            ,Work(iVecs),nBas(iSym),Work(iSs),nBas(iSym),Zero
     &            ,Work(iAUX),nBas(iSym))
        Call Dgemm_('N','T',nBas(iSym),nBas(iSym),nBas(iSym),One
     &            ,Work(iAUX),nBas(iSym),Work(iVecs),nBas(iSym),Zero
     &            ,Work(iTrans),nBas(iSym))
        Call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One
     &            ,Work(iVecs),nBas(iSym),Work(iSp),nBas(iSym),Zero
     &            ,Work(iAUX),nBas(iSym))
        Call Dgemm_('N','T',nBas(iSym),nBas(iSym),nBas(iSym),One
     &            ,Work(iAUX),nBas(iSym),Work(iVecs),nBas(iSym),Zero
     &            ,Work(iTrani),nBas(iSym))
        Call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One
     &            ,Work(iTrans),nBas(iSym),Work(iDao+indS),nBas(iSym)
     &            ,Zero,Work(iAUX),nBas(iSym))
        Call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One
     &            ,Work(iAUX),nBas(iSym),Work(iTrans),nBas(iSym),Zero
     &            ,Work(iOrtoD),nBas(iSym))
        kaunter=0
        Do 2051, iB1=1,nBas(iSym)
          Do 2052, iB2=1,nBas(iSym)
            Work(iVecs+kaunter)=Zero
            If(iB1.eq.iB2)Work(iVecs+kaunter)=One
            kaunter=kaunter+1
2052      Continue
2051    Continue
        kaunter=0
        Do 2053, i=1,nBas(iSym)
          Do 2054, j=1,i
            Work(iOrtoDt+kaunter)=Work(iOrtoD+i+(j-1)*nBas(iSym)-1)
            kaunter=kaunter+1
2054      Continue
2053    Continue
        Call Jacob(Work(iOrtoDt),Work(iVecs),nBas(iSym),nBas(iSym))
        Call Dgemm_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One
     &            ,Work(iTrani),nBas(iSym),Work(iVecs),nBas(iSym),Zero
     &            ,Work(iAUX),nBas(iSym))
        kaunt=0
        kaunter=0
        Do 208, i=1,nBas(iSym)
          Do 209, j=1,i
            If(i.eq.j) then
              Work(iOccNat+kaunt)=Work(iOrtoDt+kaunter)
              kaunt=kaunt+1
            Endif
            kaunter=kaunter+1
209       Continue
208     Continue
        Call Jacord3(Work(iOccNat),Work(iAUX),nBas(iSym),nBas(iSym))
        Call Add_Info('AVERD_OCC',Work(iOccNat),5,5)
        If(iPrint.ge.5) then
          Write(Titorb,'(A)')'All average orbitals in this irrep.'
          Thr=-1D0
          Call Primo(Titorb,.true.,.false.,Thr,Dum,1,nBas(iSym)
     &              ,nBas(iSym),BsLbl,Dummy,Work(iOccNat),Work(iAUX),-1)
        Endif
        call dcopy_(nBS,Work(iAUX),1,Work(iOrbs+indS),1)
        call dcopy_(nBas(iSym),Work(iOccNat),1,Work(iOccs+indB),1)
        Call GetMem('EigV','Free','Real',iVecs,nBS)
        Call GetMem('St','Free','Real',iSt,nBT)
        Call GetMem('Si','Free','Real',iSi,nBT)
        Call GetMem('Ss','Free','Real',iSs,nBS)
        Call GetMem('Sp','Free','Real',iSp,nBS)
        Call GetMem('AUX','Free','Real',iAUX,nBS)
        Call GetMem('TransS','Free','Real',iTrans,nBS)
        Call GetMem('TransSi','Free','Real',iTrani,nBS)
        Call GetMem('OrthoDensS','Free','Real',iOrtoD,nBS)
        Call GetMem('OrthoDensT','Free','Real',iOrtoDt,nBT)
        Call GetMem('Occs','Free','Real',iOccNat,nBas(iSym))
        indT=indT+nBT
        indS=indS+nBS
        indB=indB+nBas(iSym)
201   Continue
      Write(6,*)
      Write(6,*)
      Write(6,*)
      Write(6,*)'---Average natural orbital generation completed!---'
      Write(6,*)


*
*-- Write average orbitals to a file with the same format as
*   SCF-orbitals. To the outfile are orbital energies added just
*   to make the NEMO happy, the numbers are just bosh!
*
      Write(6,*)
      Write(6,*)
      Write(6,*)'Average orbitals put on AVEORB'
      Write(6,*)
      Write(6,*)'NB: Dummy orbital energies added to AVEORB for '
     &//'compatability reasons.'
      Write(6,*)'    They have no physical meaning.'
      Call GetMem('Zeros','Allo','Real',iZero,ntot)
      call dcopy_(ntot,[Zero],0,Work(iZero),1)
      LuOut=65
      LuOut=IsFreeUnit(LuOut)
      Title='Average Orbitals'
      OrbFile='AVEORB'
      Plab='COE'
      Call WrVec(OrbFile,LuOut,Plab,nSym,nBas,nBas,Work(iOrbs)
     &           ,Work(iOccs),Work(iZero),iDummy,Title)

*
*-- Say something about orbital occupation.
*
      Write(6,*)
      Write(6,*)
      Write(6,'(A)')' |  Average orbital occupation.'
      Write(6,'(A)')' |-----------------------------'
      Write(6,'(A,E18.8)')' |    Threshold: ',ThrOcc
      Write(6,*)
      nOrb=0
      iO=0
      Do 9991, iSym=1,nSym
        Do 9992, iB=1,nBas(iSym)
          If(Work(iOccs+iO+iB-1).lt.ThrOcc) GoTo 9992
          nOrb=nOrb+1
9992    Continue
        Write(6,9999)'      Symmetry:',iSym,'   Number of orbitals '
     &//'below threshold:',nOrb
        iO=iO+nBas(iSym)
9991  Continue
      Write(6,*)
9999  Format(A,I2,A,I4)
      Call GetMem('Zeros','Free','Real',iZero,ntot)
      Call GetMem('NatOrbAcc','Free','Real',iOrbs,ntot2)
      Call GetMem('NatOccAcc','Free','Real',iOccs,ntot)
      Call GetMem('Density','Free','Real',iDao,ntot2)
      Call GetMem('Overlap','Free','Real',iS,lsmat+4)

*
*-- Good Bye.
*
      ireturn=0
      Return
      End
