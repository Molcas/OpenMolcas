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
* Copyright (C) 1991,1993,1999, Roland Lindh                           *
*               1995, Martin Schuetz                                   *
************************************************************************
      SubRoutine Eval_Ints_New_(iiS,jjS,kkS,llS,TInt,nTInt,
     &                          iTOffs,nShi,nShj,nShk,nShl,
     &                          nShOffi,nShOffj,nShOffk,nShOffl,
     &                          Integ_Proc,
     &                          Dens,Fock,lDens,ExFac,nDens,
     &                          Ind,nInd,FckNoClmb,FckNoExch,
     &                          Thize,W2Disc,PreSch,Disc_Mx,Disc,
     &                          Quad_ijkl,DoIntegrals,DoFock)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals, parallel region          *
*          contains memory partitioning and loops over uncontracted    *
*          functions...                                                *
*                                                                      *
*  Input:                                                              *
*          iiS,jjS,kkS,llS     : shell indices                         *
*          TInt                : Computed Integrals                    *
*          nTInt               : dimension of TInt                     *
*          iTOffs              : iTOffs holds symmetry block offsets   *
*                                                                      *
*     nShi,nShj,          Dimensions used for blocks in Tint (input)   *
*     nshk,nshl:          Symmetry block isym,jsym,ksym,lsym for       *
*                         shells iS,jS,kS,lS starts at                 *
*                         iTOffs(ksym,jsym,isym)+1 and is dimensioned  *
*                         [nshl(lsym),nshk(ksym),nshj(jsym,nshi(isym)] *
*                         Note that l runs fastest! The dimensions     *
*                         must be larger or equal to the number of     *
*                         SAOs in the specified shells and symmetries, *
*                         otherwise chaos!!                            *
*                                                                      *
*     nShOffi,nShOffj,    Offsets of Integral symmetry blocks (input)  *
*     nShOffk,nShOffl:    An Integral (lso,kso|jso,iso) is placed at   *
*                         [lb,kb,jb,ib] where lb=lso-nShOffl(lsym),    *
*                         kb=kso-nShOffk(ksym) etc. Here lso,kso etc   *
*                         are the SAO labels within their symmetry.    *
*                         More explicitly, the Integral is stored in   *
*                         in Tint(ijkl), where                         *
*                                                                      *
*                         ijkl = iTOffs(ksym,jsym,isym)                *
*                           + (ib-1)*nshj(jsym)*nshk(ksym)*nshl(lsym)  *
*                           + (jb-1)*nshk(ksym)*nshl(lsym)             *
*                           + (kb-1)*nshl(lsym)                        *
*                           +  lb                                      *
*                                                                      *
*          Dens                : 1-particle density matrix             *
*          lDens               : length of density/Fock matrices       *
*          nDens               : # of density/Fock matrices            *
*          ExFac               : another scaling factor passed to      *
*                                Integ_Proc                            *
*          Ind,nInd            : auxiliary index list for Fock matrix  *
*                                construction (cf. d1ind_CpFck)        *
*          FckNoClmb           : no Coulomb contributions to Fock mat. *
*          FckNoExch           : no exchange contributions to Fock mat.*
*          Thize               : int threshold for disk write (SD)     *
*          W2Disc,PreSch       : booleans 1st iter / prescreening      *
*          Disc_Mx             : # ints to write on disk (semidirect)  *
*          iDisk               : act. position in file TMPINT          *
*                                                                      *
*  Auxiliary:                                                          *
*          ipMem1              : base pointer to Scratch space for     *
*                                Integral batch, which is further      *
*                                partitioned within this subroutine    *
*          MemMax              : amount of wrkspace for processing of  *
*                                integral batch                        *
*                                                                      *
*  Output: Fock                : 2el Hamiltonian                       *
*          Disc                : # ints written to disk (semidirect)   *
*                                                                      *
*                                                                      *
*  Local:                                                              *
*          Coor                : coordinates of four centers           *
*          iAngV               : angular momenta                       *
*          iCmpV               : # spherical components                *
*          iShelV,iShllV       : shell indices                         *
*          iAOV                : pointers to ??                        *
*          iStabs              : IDs of 4 unique centers, i.e. ptrs to *
*          Shijij,                           : swap booleans           *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter,QExit                                            *
*              GetMem                                                  *
*              Int_Setup                                               *
*              Dens_Info                                               *
*              MemRys                                                  *
*              PSOAO0                                                  *
*              Picky_                                                  *
*              TwoEl_NoSym                                             *
*              TwoEl_Sym                                               *
*              Integ_Proc                                              *
*                                                                      *
*             Roland Lindh / Martin Schuetz,                           *
*             Dept. of Theoretical Chemistry, University of Lund,      *
*             SWEDEN.                                                  *
*             Modified for k2 loop. August '91                         *
*             Modified for direct SCF. January '93                     *
*             Modified to minimize overhead for calculations with      *
*             small basis sets and large molecules. Sept. '93          *
*             parallel region split off in drvtwo.f, April '95         *
*             Total rehack May '99                                     *
************************************************************************
      use k2_setup
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
      External Integ_Proc
*     subroutine parameters
      Integer lDens
      Real*8  Coor(3,4),Thize,Fock(lDens,nDens),Dens(lDens,nDens),
     &        ExFac(nDens), Disc_Mx,Disc, TInt(nTInt), Tmax
      Integer iAngV(4),iCmpV(4),
     &        iShelV(4),iShllV(4),iAOV(4),iStabs(4),
     &        ipMem1,MemMax,
     &        iTOffs(8,8,8),Map4(4),
     &        nShi(0:7), nShj(0:7), nShk(0:7), nShl(0:7),
     &        nShOffi(0:7), nShOffj(0:7), nShOffk(0:7), nShOffl(0:7),
     &        Ind(nInd,nInd,2),kOp(4)
      Logical Shijij, W2Disc,PreSch,NoInts,FckNoClmb(nDens),
     &        FckNoExch(nDens), DoIntegrals,DoFock
*
#include "ndarray.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "k2.fh"
#include "setup.fh"
#include "status.fh"
*
      Common /ibas_ricd/ jbas_, lbas_
*     local variables to save
      Integer ipDDij,ipDDkl,ipDDik,ipDDil,ipDDjk,ipDDjl,
     &        iBsInc,jBsInc,kBsInc,lBsInc,iPrInc,jPrInc,kPrInc,lPrInc,
     &        ipMem2,
     &        Mem1,Mem2
      Save    ipDDij,ipDDkl,ipDDik,ipDDil,ipDDjk,ipDDjl,
     &        iBsInc,jBsInc,kBsInc,lBsInc,iPrInc,jPrInc,kPrInc,lPrInc,
     &        ipMem2,
     &        Mem1,Mem2
*     other local variables
      Integer iAOst(4), iPrimi,jPrimj,kPrimk,lPriml,
     &        iBasi,jBasj,kBask,lBasl,
     &  ipCffi,jpCffj,kpCffk,lpCffl,
     &  iBasn,jBasn,kBasn,lBasn,
     &  k2ij,nDCRR,k2kl,nDCRS, ipTmp,
     &  mDij,mDik,mDjk,mDkl,mDil,mDjl,
     &  mDCRij,mDCRik,mDCRjk,mDCRkl,mDCRil,mDCRjl,
     &  ipDij,ipDik,ipDjk,ipDkl,ipDil,ipDjl,
     &  ipZI,ipKab,ipP,nZeta,
     &  ipEta,ipEI,ipiEta,ipKcd,ipQ,nEta
      Integer   nSO,iBasAO,jBasAO,kBasAO,lBasAO,
     &  iS,jS,kS,lS,ijS,klS,ikS,ilS,jkS,jlS
      Logical IJeqKL
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      nElem(i)=(i+1)*(i+2)/2
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
      ipDDij=ip_Dummy
      mDCRij=1
      ipDDkl=ip_Dummy
      mDCRkl=1
      ipDDik=ip_Dummy
      ipDDil=ip_Dummy
      ipDDjk=ip_Dummy
      ipDDjl=ip_Dummy
*                                                                      *
************************************************************************
*                                                                      *
      If (ERI_Status.ne.Active) Then
         Call WarningMessage(2,
     &               'Eval_Ints_: Integral environment is not set up!')
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      iRout=9
      iPrint=nPrint(iRout)
*
      NoInts=.True.
      Tmax=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     If memory not allocated already at this point allocate!          *
*                                                                      *
      If (MemMax_int.eq.0) Then
        Call GetMem('MaxMem','Max','Real',iDum,MemMax)
        Call GetMem('MaxMem','Allo','Real',ipMem1,MemMax)
        MemMax_int=MemMax
        ipMem_int =ipMem1
      Else
        MemMax=MemMax_int
        ipMem1=ipMem_int
      End If
*
      Map4(1)=1
      Map4(2)=2
      Map4(3)=3
      Map4(4)=4
      iS_=Max(iiS,jjS)
      jS_=Min(iiS,jjS)
      kS_=Max(kkS,llS)
      lS_=Min(kkS,llS)
      If (iiS.ne.iS_) Then
          iTmp=Map4(1)
          Map4(1)=Map4(2)
          Map4(2)=iTmp
      End If
      If (kkS.ne.kS_) Then
          iTmp=Map4(3)
          Map4(3)=Map4(4)
          Map4(4)=iTmp
      End If
*     Write (*,*) ' -->',iS_,jS_,kS_,lS_,'<--'
*                                                                      *
************************************************************************
*                                                                      *
      Call Int_Setup(iSD,mSkal,iS_,jS_,kS_,lS_,
     &               Coor,Shijij,
     &               iAngV,iCmpV,iShelV,iShllV,iAOV,iStabs)
*                                                                      *
************************************************************************
*                                                                      *
      iPrimi   = nExp(iShllV(1))
      jPrimj   = nExp(iShllV(2))
      kPrimk   = nExp(iShllV(3))
      lPriml   = nExp(iShllV(4))
      iBasi    = nBasis(iShllV(1))
      jBasj    = nBasis(iShllV(2))
      kBask    = nBasis(iShllV(3))
      lBasl    = nBasis(iShllV(4))
      ipCffi   = ipCff(iShllV(1))
      jpCffj   = ipCff(iShllV(2))
      kpCffk   = ipCff(iShllV(3))
      lpCffl   = ipCff(iShllV(4))
      nZeta = iPrimi * jPrimj
      nEta = kPrimk * lPriml
      mDij=nZeta+1 ! Dummy initialize
      mDkl=nEta+1  ! Dummy initialize
*
      nHRRAB=iCmpV(1)*iCmpV(2)*
     &      (nabSz(iAngV(1)+iAngV(2)) -
     &       nabSz(Max(iAngV(1),iAngV(2))-1))
      nHRRCD=iCmpV(3)*iCmpV(4)*
     &      (nabSz(iAngV(3)+iAngV(4)) -
     &       nabSz(Max(iAngV(3),iAngV(4))-1))
      nHRRAB=nHRRAB*nIrrep
      nHRRCD=nHRRCD*nIrrep
      If (DoGrad_) Then
         ijCmp=nElem(iAngV(1))*nElem(iAngV(2))
         klCmp=nElem(iAngV(3))*nElem(iAngV(4))
      Else
         ijCmp=0
         klCmp=0
      End If
      mData1=nZeta*(nDArray+2*ijCmp)+nDScalar+nHRRAB
      mData2=nEta*(nDArray+2*klCmp)+nDScalar+nHRRCD
*                                                                      *
************************************************************************
*                                                                      *
*     partition memory for K2(ij)/K2(kl) temp spaces zeta,eta,kappa,P,Q
      ipZI  = ipZeta + nZeta
      ipKab = ipZI   + nZeta
      ipP   = ipKab  + nZeta
      ipEta = ipP    + nZeta*3
      ipEI  = ipEta  + nEta
      ipKcd = ipEI   + nEta
      ipQ   = ipKcd  + nEta
      ipiEta = ipiZet + nZeta
*                                                                      *
************************************************************************
*                                                                      *
*
*     No SO block in direct construction of the Fock matrix.
      nSO = MemSO2(iAngV(1),iAngV(2),iAngV(3),iAngV(4),
     &             iCmpV(1),iCmpV(2),iCmpV(3),iCmpV(4),
     &             iShelV(1),iShelV(2),iShelV(3),iShelV(4))
      If (nSO.eq.0) Then
        Return
*       Call GetMem('0Exit_DrvTwo0','Check','Real',iDum,iDum)
      End If
      If (.Not.DoIntegrals) nSO = 0
*
      iS = iShelV(1)
      jS = iShelV(2)
      kS = iShelV(3)
      lS = iShelV(4)
      ijS = iTri(iS,jS)
      klS = iTri(kS,lS)
      ikS = iTri(iS,kS)
      ilS = iTri(iS,lS)
      jkS = iTri(jS,kS)
      jlS = iTri(jS,lS)
*                                                                      *
************************************************************************
*                                                                      *
*-----Pick up pointers to k2 entities.
*
      k2ij  = IndK2(1,ijS)
      nDCRR = IndK2(2,ijS)
      k2kl  = IndK2(1,klS)
      nDCRS = IndK2(2,klS)
*                                                                      *
************************************************************************
*                                                                      *
*-----Pick up pointers to desymmetrized 1st order density
*     matrices. Observe that the desymmetrized 1st order
*     density matrices follows the contraction index.
*
      If (DoFock) Then
         ipTmp = ipDijs
         Nr_of_D=1
         Call Dens_Info(ijS,ipDij,ipDum,mDCRij,ipDDij,ipTmp,Nr_of_D)
         Call Dens_Info(klS,ipDkl,ipDum,mDCRkl,ipDDkl,ipTmp,Nr_of_D)
         Call Dens_Info(ikS,ipDik,ipDum,mDCRik,ipDDik,ipTmp,Nr_of_D)
         Call Dens_Info(ilS,ipDil,ipDum,mDCRil,ipDDil,ipTmp,Nr_of_D)
         Call Dens_Info(jkS,ipDjk,ipDum,mDCRjk,ipDDjk,ipTmp,Nr_of_D)
         Call Dens_Info(jlS,ipDjl,ipDum,mDCRjl,ipDDjl,ipTmp,Nr_of_D)
*
c        Write (*,*) ' Pointers to D=',
c    &                ipDij,ipDkl,ipDik,ipDil,ipDjk,ipDjl
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.99) Then
         Write (6,*) ' *** Centers ***'
         Write (6,'(3F7.3,6X,3F7.3)')
     &         ((Coor(i,j),i=1,3),j=1,2)
         Write (6,'(3F7.3,6X,3F7.3)')
     &         ((Coor(i,j),i=1,3),j=3,4)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Compute memory request for the primitives, i.e.
*     how much memory is needed up to the transfer
*     equation.
      Call MemRys(iAngV,MemPrm)
*                                                                      *
************************************************************************
*                                                                      *
*     Decide on the partioning of the shells based on the
*     available memory and the requested memory.
      Call PSOAO0(nSO,MemPrm,MemMax,
     &            iAngV,iCmpV,
     &            iBasi,iBsInc,jBasj,jBsInc,
     &            kBask,kBsInc,lBasl,lBsInc,
     &            iPrimi,iPrInc,jPrimj,jPrInc,
     &            kPrimk,kPrInc,lPriml,lPrInc,
     &            ipMem1,ipMem2,
     &            Mem1,Mem2,DoFock)
      If (iPrint.ge.59) Then
         Write (6,*) ' ************** Memory partioning **************'
         Write (6,*) ' ipMem1=',ipMem1
         Write (6,*) ' ipMem2=',ipMem2
         Write (6,*) ' Mem1=',Mem1
         Write (6,*) ' Mem2=',Mem2
         Write (6,*) ' iBasi,iBsInc=',iBasi,iBsInc
         Write (6,*) ' jBasj,jBsInc=',jBasj,jBsInc
         Write (6,*) ' kBasi,kBsInc=',kBask,kBsInc
         Write (6,*) ' lBasl,lBsInc=',lBasl,lBsInc
         Write (6,*) ' iPrimi,iPrInc=',iPrimi,iPrInc
         Write (6,*) ' jPrimj,jPrInc=',jPrimj,jPrInc
         Write (6,*) ' kPrimk,kPrInc=',kPrimk,kPrInc
         Write (6,*) ' lPriml,lPrInc=',lPriml,lPrInc
         Write (6,*) ' ***********************************************'
      End If
*                                                                      *
************************************************************************
*                                                                      *
      jbas_=jBasj
      lbas_=lBasl
*                                                                      *
************************************************************************
*                                                                      *
*     FacInt : scaling factor for the integrals passed down to         *
*              Integ_Proc                                              *
*
      FacInt=1.0D0
*                                                                      *
************************************************************************
*                                                                      *
*     These loops will partition the contraction loops if there is not
*     enough memory to store the whole SO/AO-block simultaneously. The
*     memory partitioning is determined by PSOAO0.
*
      Do iBasAO = 1, iBasi, iBsInc
         iBasn=Min(iBsInc,iBasi-iBasAO+1)
         iAOst(1) = iBasAO-1
*
         Do jBasAO = 1, jBasj, jBsInc
            jBasn=Min(jBsInc,jBasj-jBasAO+1)
            iAOst(2) = jBasAO-1
*
*---------- Move appropiate portions of the desymmetrized 1st
*           order density matrix.
*
            If (DoFock) Then
               Call Picky_(iBasi,iBsInc,iPrimi,iBasAO,iBasn,
     &                     jBasj,jBsInc,jPrimj,jBasAO,jBasn,
     &                     iCmpV(1),iCmpV(2),iShelV(1),iShelV(2),
     &                     mDCRij,ipDij,ipDDij,mDij,nIrrep)
            End If
*
            Do kBasAO = 1, kBask, kBsInc
               kBasn=Min(kBsInc,kBask-kBasAO+1)
               iAOst(3) = kBasAO-1
*
               If (DoFock) Then
                  Call Picky_(iBasi,iBsInc,iPrimi,iBasAO,iBasn,
     &                        kBask,kBsInc,kPrimk,kBasAO,kBasn,
     &                        iCmpV(1),iCmpV(3),iShelV(1),iShelV(3),
     &                        mDCRik,ipDik,ipDDik,mDik,nIrrep)
               End If
*
               If (DoFock) Then
                  Call Picky_(jBasj,jBsInc,jPrimj,jBasAO,jBasn,
     &                        kBask,kBsInc,kPrimk,kBasAO,kBasn,
     &                        iCmpV(2),iCmpV(3),iShelV(2),iShelV(3),
     &                        mDCRjk,ipDjk,ipDDjk,mDjk,nIrrep)
               End If
*
                Do lBasAO = 1, lBasl, lBsInc
                   lBasn=Min(lBsInc,lBasl-lBasAO+1)
                   iAOst(4) = lBasAO-1
*
                   If (DoFock) Then
                      Call Picky_(kBask,kBsInc,kPrimk,kBasAO,kBasn,
     &                            lBasl,lBsInc,lPriml,lBasAO,lBasn,
     &                            iCmpV(3),iCmpV(4),iShelV(3),iShelV(4),
     &                            mDCRkl,ipDkl,ipDDkl,mDkl,nIrrep)
                   End If
*
                   If (DoFock) Then
                      Call Picky_(iBasi,iBsInc,iPrimi,iBasAO,iBasn,
     &                            lBasl,lBsInc,lPriml,lBasAO,lBasn,
     &                            iCmpV(1),iCmpV(4),iShelV(1),iShelV(4),
     &                            mDCRil,ipDil,ipDDil,mDil,nIrrep)
                   End If
*
                   If (DoFock) Then
                      Call Picky_(jBasj,jBsInc,jPrimj,jBasAO,jBasn,
     &                            lBasl,lBsInc,lPriml,lBasAO,lBasn,
     &                            iCmpV(2),iCmpV(4),iShelV(2),iShelV(4),
     &                            mDCRjl,ipDjl,ipDDjl,mDjl,nIrrep)
                   End If
*                                                                      *
************************************************************************
*                                                                      *
*                 Compute SO/AO-integrals
*
                  If (nIrrep.eq.1) Then
*
                  Call TwoEl_NoSym_New(iS_,jS_,kS_,lS_,
     &                            Coor,
     &                            iAngV,iCmpV,iShelV,iShllV,
     &                            iAOV,iAOst,NoInts,
     &                            iStabs(1),iStabs(2),
     &                            iStabs(3),iStabs(4),
     &                            iPrimi,iPrInc,jPrimj,jPrInc,
     &                            kPrimk,kPrInc,lPriml,lPrInc,
     &                            Data_k2(k2ij),mData1,nDCRR,
     *                            Data_k2(k2kl),mData2,nDCRS,
     &                            IJeqKL,kOp,Disc_Mx,Disc,Thize,
     &                            Work(ipDDij),mDij,mDCRij,
     &                            Work(ipDDkl),mDkl,mDCRkl,
     & Work(ipDDik),mDik,mDCRik,Work(ipDDil),mDil,mDCRil,Work(ipDDjk),
     & mDjk,mDCRjk, Work(ipDDjl),mDjl,mDCRjl,Fock,Dens,lDens,
     &                            Work(ipCffi+(iBasAO-1)*iPrimi),iBasn,
     &                            Work(jpCffj+(jBasAO-1)*jPrimj),jBasn,
     &                            Work(kpCffk+(kBasAO-1)*kPrimk),kBasn,
     & Work(lpCffl+(lBasAO-1)*lPriml),lBasn, Work(ipFT),nFT,
     &Work(ipZeta),Work(ipZI),iWork(ipiZet),Work(ipKab),Work(ipP),nZeta,
     &Work(ipEta), Work(ipEI),iWork(ipiEta),Work(ipKcd),Work(ipQ),nEta,
     & Work(ipMem1),nSO,Work(ipMem2),Mem2,
     & Shijij,W2Disc,PreSch,Quad_ijkl,nHRRAB,nHRRCD,
     & DoIntegrals,DoFock,FckNoClmb(1),FckNoExch(1),Work(ipAux),nAux,
     & ExFac(1))
*
                  Else
*

                  Call TwoEl_Sym_New(iS_,jS_,kS_,lS_,
     &                            Coor,
     &                            iAngV,iCmpV,iShelV,iShllV,
     &                            iAOV,iAOst,NoInts,
     &                            iStabs(1),iStabs(2),
     &                            iStabs(3),iStabs(4),
     &                            iPrimi,iPrInc,jPrimj,jPrInc,
     &                            kPrimk,kPrInc,lPriml,lPrInc,
     &                            Data_k2(k2ij),mData1,nDCRR,
     &                            Data_k2(k2kl),mData2,nDCRS,
     &                            IJeqKL,kOp,Disc_Mx,Disc,Thize,
     &                            Work(ipDDij),mDij,mDCRij,
     &                            Work(ipDDkl),mDkl,mDCRkl,
     & Work(ipDDik),mDik,mDCRik,Work(ipDDil),mDil,mDCRil,Work(ipDDjk),
     & mDjk,mDCRjk, Work(ipDDjl),mDjl,mDCRjl,Fock,Dens,lDens,
     &                            Work(ipCffi+(iBasAO-1)*iPrimi),iBasn,
     &                            Work(jpCffj+(jBasAO-1)*jPrimj),jBasn,
     &                            Work(kpCffk+(kBasAO-1)*kPrimk),kBasn,
     & Work(lpCffl+(lBasAO-1)*lPriml),lBasn, Work(ipFT),nFT,
     &Work(ipZeta),Work(ipZI),iWork(ipiZet),Work(ipKab),Work(ipP),nZeta,
     &Work(ipEta), Work(ipEI),iWork(ipiEta),Work(ipKcd),Work(ipQ),nEta,
     & Work(ipMem1),nSO,Work(ipMem2),Mem2,
     & Shijij,W2Disc,PreSch,Quad_ijkl,nHRRAB,nHRRCD,
     & DoIntegrals,DoFock,FckNoClmb(1),FckNoExch(1),Work(ipAux),nAux,
     & ExFac(1))
*
                  End If
*                                                                      *
************************************************************************
*                                                                      *
*              Process SO/AO-integrals
*
                  nijkl=iBasn*jBasn*kBasn*lBasn
                  If (DoIntegrals.and..Not.NoInts) Then
*                    Get max AO/SO integrals
                     If (Petite) Then
                        n=nijkl*iCmpV(1)*iCmpV(2)*iCmpV(3)*iCmpV(4)
                        ip=ipMem2
                     Else
                        n=nijkl*nSO
                        ip=ipMem1
                     End If
                     Tmax=max(Tmax,
     &                        abs(Work(ip+iDAMax_(n,Work(ip),1)-1)))
                     If (Tmax.gt.CutInt) Then
                        Call Integ_Proc(iCmpV,iShelV,Map4,
     &                                  iBasn,jBasn,kBasn,lBasn,kOp,
     &                                  Shijij,IJeqKL,iAOV,iAOst,nijkl,
     &                                  Work(ipMem2),Work(ipMem1),nSO,
     &                                  iWork(ipiSOSym),mSkal,nSOs,
     &                                  TInt,nTInt,FacInt,
     &                                  iTOffs,nIrrep,
     &                                  nShi,nShj,nShk,nShl,
     &                                  nShOffi,nShOffj,nShOffk,nShOffl,
     &                                  Dens,Fock,lDens,ExFac,nDens,
     &                                  Ind,nInd,FckNoClmb,FckNoExch)
                     Else
                        Tmax=Zero
                     End If
                  End If
*
               End Do
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*
      Return
      End
