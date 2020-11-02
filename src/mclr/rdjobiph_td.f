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
      Subroutine RdJobIph_td(CIVec)
************************************************************************
*                                                                      *
*     Read the contents of the JOBIPH file.                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Use Arrays, only: CMO, G2t, G2sq, G1t
      Implicit Real*8 (a-h,o-z)
#include "Input.fh"
#include "Files_mclr.fh"
#include "glbbas_mclr.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
      Character*72 Line
      Real*8 rdum(1)
      Character(Len=1), Allocatable:: TempTxt(:)
      Real*8, Allocatable::  Tmp2(:), G2tts(:), G2tta(:)
      Real*8, Allocatable:: CIVec(:,:)

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

*----------------------------------------------------------------------*
*     Save the ROOT input parameter                                    *
*----------------------------------------------------------------------*
      kRoots=lRoots
*----------------------------------------------------------------------*
*     Read the table of disk adresses                                  *
*----------------------------------------------------------------------*
      Call DaName(LuJob,FnJob)
      iDisk=0
      Call iDaFile(LuJob,2,iToc,iTOCIPH,iDisk)
*----------------------------------------------------------------------*
*     Read the the system description                                  *
*----------------------------------------------------------------------*
      Call mma_allocate(TempTxt,LENIN8*MxOrb,Label='TempTxt')
      iDisk=iToc(1)
      Call WR_RASSCF_Info(LuJob,2,iDisk,
     &                    nActEl,iSpin,nSym,State_sym,nFro,
     &                    nIsh,nAsh,nDel,
     &                    nBas,MxSym,TempTxt,LENIN8*mxorb,
     &                    nConf,HeaderJP,144,
     &                    TitleJP,4*18*mxTit,PotNuc0,lRoots,
     &                    nRoots,iRoot,mxRoot,
     &                    nRs1,nRs2,nRs3,
     &                    nHole1,nElec3,iPt2,Weight)
      Call mma_deallocate(TempTxt)
*----------------------------------------------------------------------*
*     Overwrite the variable lroots if approriate                      *
*----------------------------------------------------------------------*
      If ( kRoots.ne.-1 ) then
         If ( iPt2.ne.0 ) then
            Write (6,*)
            Write (6,*) ' *** Error in subroutine RDJOBIPH_TD ***'
            Write (6,*) ' Pt2.ne.0'
            Write (6,*)
         Else if ( kRoots.gt.lRoots ) then
            Write (6,*)
            Write (6,*) ' *** Error in subroutine RDJOBIPH_TD ***'
            Write (6,*) ' kRoots.gt.lRoots'
            Write (6,*)
         End If
         lRoots=kRoots
         nRoots=1
      Else If ( nRoots.ne.1 ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine RDJOBIPH_TD ***'
         Write (6,*) ' nRoots.ne.1'
         Write (6,*)
      End If
*----------------------------------------------------------------------*
*     Precompute the total sum of variables and size of matrices       *
*----------------------------------------------------------------------*
      ntIsh=0
      ntItri=0
      ntIsqr=0
      ntAsh=0
      ntAtri=0
      ntAsqr=0
      ntBas=0
      ntBtri=0
      ntBsqr=0
      nna=0
      Length=0
      Do 10 iSym=1,nSym
         norb(isym)=nbas(isym)-ndel(isym)
         ntIsh=ntIsh+nIsh(iSym)
         ntItri=ntItri+nIsh(iSym)*(nIsh(iSym)+1)/2
         ntIsqr=ntIsqr+nIsh(iSym)*nIsh(iSym)
         ntAsh=ntAsh+nAsh(iSym)
         ntAtri=ntAtri+nAsh(iSym)*(nAsh(iSym)+1)/2
         ntAsqr=ntAsqr+nAsh(iSym)*nAsh(iSym)
         ntBas=ntBas+nBas(iSym)
         ntBtri=ntBtri+nBas(iSym)*(nBas(iSym)+1)/2
         ntBsqr=ntBsqr+nBas(iSym)*nBas(iSym)
         nA(iSym)=nna
         nnA=nnA+nAsh(isym)
         Length=Length+nbas(isym)*norb(isym)
10    Continue
*----------------------------------------------------------------------*
*     Load the orbitals used in the last macro iteration               *
*----------------------------------------------------------------------*
*
      Call mma_allocate(CMO,Length,Label='CMO')
      Call Get_CMO(CMO,Length)
*     iDisk=iToc(9)
*     IF(IPT2.EQ.0) iDisk=iToc(2)
*     Call dDaFile(LuJob,2,CMO,ntBsqr,iDisk)
      If ( .false. ) then
         jpCMO=1
         Do 15 iSym=1,nSym
            Write(Line,'(A,i2.2)') 'MO coefficients, iSym = ',iSym
            Call RecPrt(Line,' ',CMO(jpCMO),nBas(iSym),nBas(iSym))
            jpCMO=jpCMO+nBas(iSym)*nBas(iSym)
15       Continue
      End If
*----------------------------------------------------------------------*
*     Load the CI vector for the root lRoots                           *
*----------------------------------------------------------------------*
      Call mma_allocate(CIVec,nConf,1,Label='CIVec')
      iDisk=iToc(4)
      Do i=1,lroots-1
         Call dDaFile(LuJob,0,CIVec,nConf,iDisk)
      End Do
      Call dDaFile(LuJob,2,CIVec,nConf,iDisk)
      If (.false.) Call DVcPrt('CI coefficients',' ',CIVec,nConf)
*----------------------------------------------------------------------*
*     Load state energy                                                *
*----------------------------------------------------------------------*
      Call mma_allocate(Tmp2,mxRoot*mxIter,Label='Tmp2')
      iDisk=iToc(6)
      Call dDaFile(LuJob,2,Tmp2,mxRoot*mxIter,iDisk)
      ERASSCF(1)=0.0d0
      Do 20 iter=0,mxIter-1
         Temp=Tmp2(iter*mxRoot+lRoots)
         If ( Temp.ne.0.0D0 ) ERASSCF(1)=Temp
20    Continue
      Call mma_deallocate(Tmp2)
*     If ( debug ) Write(*,*) ' RASSCF energy =',ERASSCF(1)
*
      nAct  = 0
      nAct2 = 0
      nAct4 = 0
      Do iSym = 1, nSym
       nAct = nAct + nAsh(iSym)
       nAct2=nAct2+nAsh(iSym)**2
      End Do
      Do iS = 1, nSym
       Do jS = 1, nSym
        Do kS = 1, nSym
         lS=iEOr(iEOr(is-1,js-1),ks-1)+1
         nAct4=nAct4+nAsh(iS)*nAsh(jS)*nAsh(kS)*nAsh(lS)
        End Do
       End Do
      End Do
*-----------------------------------
* One electron dens - triang stor.
*-----------------------------------
      nG1 = nAct*(nAct+1)/2
      Call mma_allocate(G1t,nG1,Label='G1t')
      G1t(:)=0.0d0
      ipG1=ip_of_Work(G1t)
*
*--------------------------------------------------
* Triangular part of two electron dens,
* symmetric part
*--------------------------------------------------
      nG2=nG1*(nG1+1)/2

      Call mma_allocate(G2sq,nAct**4,Label='G2sq')
      Call mma_allocate(G2t,nG2,Label='G2t')
      Call mma_allocate(G2tts,nG2,Label='G2tts')
      Call mma_allocate(G2tta,nG2,Label='G2tta')
      iDisk=iToc(3)
      Do i=1,lroots-1
         Call dDaFile(LuJob,0,rdum,nG1,iDisk)
         Call dDaFile(LuJob,0,rdum,nG1,iDisk)
         Call dDaFile(LuJob,0,rdum,nG2,iDisk)
         Call dDaFile(LuJob,0,rdum,nG2,iDisk)
      End Do
      Call dDaFile(LuJob,2,G1t,nG1,iDisk)
      Call dDaFile(LuJob,0,rdum,nG1,iDisk)
      Call dDaFile(LuJob,2,G2tts,nG2,iDisk)
      Call dDaFile(LuJob,2,G2tta,nG2,iDisk)

      Do iB=1,nAct
         Do jB=1,iB
            iDij=iTri(ib,jB)
            Do kB=1,ib
               Do lB=1,kB
                  iDkl=iTri(kB,lB)
                  fact=1.0d00
                  if(iDij.ge.iDkl .and. kB.eq.lB) fact=2.0d00
                  if(iDij.lt.iDkl .and. iB.eq.jB) fact=2.0d00
                  iijkl=itri(iDij,iDkl)
                  G2t(iijkl)=Fact*G2tts(iijkl)
               End Do
            End Do
         End Do
      End Do
*
      ipg2   =ip_of_work(G2t)
      ipg2tmm=ip_of_work(G2t)
      ipg2tpp=ip_of_Work(G2t)

*---------------------------------------------------------------
* Rectangular part of the two el dens,
* symmetric and asymmetric contributions added
* ipG2sq = ipG2tts + ipG2tta
*-----------------------------------------------------------------
*
      Do iB=1,nAct
         Do jB=1,nact
            Factij=1.0d0
            If (ib.gt.jb) Factij=-1.0d0
            iDij=iTri(ib,jB)
            iDij2=ib+(jb-1)*NACT
            Do kB=1,nact
               Do lB=1,nact
                  Factkl=1.0d0
                  If (kb.gt.lb) Factkl=-1.0d0
                  iDkl=iTri(kB,lB)
                  iDkl2=kb+(lb-1)*NACT
                  fact=1.0d00
                  Fact2=Factij*Factkl
                  if(iDij.ge.iDkl .and. kB.eq.lB) fact=2.0d00
                  if(iDij.lt.iDkl .and. iB.eq.jB) fact=2.0d00
                  iijkl=itri(iDij,iDkl)
                  iijkl2=iDij2+nact**2*(iDkl2-1)
*
                  G2sq(iijkl2)=Fact*(G2tts(iijkl)+G2tta(iijkl)*Fact2)
               End Do
            End Do
         End Do
      End Do
      Call mma_deallocate(G2tts)
      Call mma_deallocate(G2tta)
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
