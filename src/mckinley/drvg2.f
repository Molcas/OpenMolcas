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
* Copyright (C) 1990, Roland Lindh                                     *
*               1995,1996, Anders Bernhardsson                         *
************************************************************************
      SubRoutine Drvg2(Hess,nhess,l_Grd,l_Hss)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals. The four outermost loops *
*          will controll the type of the two-electron integral, eg.    *
*          (ss|ss), (sd|pp), etc. The next four loops will generate    *
*          list of symmetry distinct centers that do have basis func-  *
*          tions of the requested type.                                *
*                                                                      *
* Called from: mckinley                                                *
*                                                                      *
* Input:                                                               *
*              nHess         : Size of gradient and hessian            *
*              l_Grd,l_Hss   : Boolean on/off for gradient/hessian     *
*                              generation                              *
*                                                                      *
* Calling    : QEnter                                                  *
*              StatP                                                   *
*              Drvk2                                                   *
*              DCopy   (ESSL)                                          *
*              Swap                                                    *
*              MemRg2 Calculate memory requirement for calc area       *
*              PSOAO1 Memory partioning                                *
*              PGet0                                                   *
*              TwoEl                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March 1990                                               *
*             Anders Bernhardsson 1995-1996                            *
************************************************************************
      use Real_Spherical
      use k2_setup
      use iSD_data
      use k2_arrays
      use pso_stuff
      use Basis_Info
      use Symmetry_Info, only: nIrrep, iOper
      use Sizes_of_Seward, only:S
      use Real_Info, only: CutInt
      Implicit Real*8 (A-H,O-Z)
      External Rsv_Tsk
#include "Molcas.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "disp2.fh"
#include "buffer.fh"
#include "etwas.fh"
#include "cputime.fh"
#include "nsd.fh"
#include "setup.fh"
*     Local arrays
      Real*8, Dimension(:), Allocatable :: DeDe2(:)
      Integer, Allocatable:: ipOffDA(:,:)
      Real*8  Coor(3,4), Hess(*)
      Integer iAngV(4), iCmpV(4), iShelV(4), iShllV(4),
     &        iAOV(4), iAOst(4), JndGrd(3,4,0:7), iFnc(4),
     &        JndHss(4,3,4,3,0:7)
      Logical Shik, Shjl, Shijij, JfGrd(3,4),lpick,
     &        JfHss(4,3,4,3), JfG(4),ltri,ldot, Rsv_Tsk,
     &        l_Hss,l_Grd,lGrad,n8,ldot2,new_fock,
     &        Post_Process
      Integer moip(0:7)
#ifdef _DEBUG_
      Character*40 format
#endif
      Real*8, Allocatable:: TMax(:,:)
      Integer, Allocatable:: Ind_ij(:,:)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
* - - - - - - P R O L O G
*
      Call StatusLine(' McKinley:',
     &                ' Computing 2-electron 2nd order derivatives')
*
      ipDij   = ip_Dummy
      ipDij2  = ip_Dummy
      ipDDij  = ip_Dummy
      ipDDij2 = ip_Dummy
      ipDkl   = ip_Dummy
      ipDkl2  = ip_Dummy
      ipDDkl  = ip_Dummy
      ipDDkl2 = ip_Dummy
      ipDik   = ip_Dummy
      ipDik2  = ip_Dummy
      ipDDik  = ip_Dummy
      ipDDik2 = ip_Dummy
      ipDil   = ip_Dummy
      ipDil2  = ip_Dummy
      ipDDil  = ip_Dummy
      ipDDil2 = ip_Dummy
      ipDjk   = ip_Dummy
      ipDjk2  = ip_Dummy
      ipDDjk  = ip_Dummy
      ipDDjk2 = ip_Dummy
      ipDjl   = ip_Dummy
      ipDjl2  = ip_Dummy
      ipDDjl  = ip_Dummy
      ipDDjl2 = ip_Dummy
      ipBuffer= ip_Dummy
      ipMOC   = ip_Dummy
      iFnc(1) = -99
      iFnc(2) = -99
      iFnc(3) = -99
      iFnc(4) = -99
      nDij=0
      nDkl=0
      nDik=0
      nDjl=0
      nDil=0
      nDjk=0
      mDCRij=0
      mDCRkl=0
      mDCRik=0
      mDCRjl=0
      mDCRil=0
      mDCRjk=0
      ipDijS =0
      ipDijS2=0
*
      Call CtrlMO(moip,nAco)
*
      iiii=0
      ipdex=0
      ndisp=0
      naco=0
      New_Fock=nirrep.eq.1
      Do iS=0,nIrrep-1
         nDisp=nDisp+ldisp(is)
         naco=naco+nAsh(is)
      End do
      n8=.true.
      Int_Direct=.true.
*
      call dcopy_(nHess,[Zero],0,Hess,1)
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*                                                                      *
************************************************************************
*                                                                      *
*-----Precompute k2 entities.
*
      lgrad=l_Grd
      lpick=lgrad.and.(.not.New_Fock)
      Pren = Zero
      Prem = Zero
      nIndK2 = S%nShlls*(S%nShlls+1)/2
      Call mma_allocate(IndK2,2,nIndk2)
      Call Drvk2_mck(ndede,new_Fock)
*
      Call StatP(0)
*                                                                      *
************************************************************************
*                                                                      *
*
*-----Allocate auxiliary array for symmetry transformation
*
      nAux = nIrrep**3
      If (nIrrep==1) nAux = 1
      Call mma_allocate(Aux,nAux,Label='Aux')
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate working area
*
      MxPrm = 0
      MxDij = 0
      MxBsC = 0
      Do iAng = 0, S%iAngMx
         MxPrm = Max(MxPrm,S%MaxPrm(iAng))
         Do 2900 iCnttp = 1,nCnttp
            iShll = dbsc(iCnttp)%iVal + iAng
            iPrim = Shells(iShll)%nExp
            If (iPrim.eq.0) Go To 2900
            If (Shells(iShll)%nBasis.eq.0) Go To 2900
            iBas  = Shells(iShll)%nBasis
            iCmp  = (iAng+1)*(iAng+2)/2
            MxBsC=Max(MxBsC,iBas*iCmp)
            MxDij= Max(MxDij,(iBas**2+1)*iCmp**2+iPrim**2+1)
 2900    Continue
      End Do
      MxDij = 6 * nIrrep * MxDij
      nZeta = MxPrm * MxPrm
      nEta  = MxPrm * MxPrm
      iii=nDens*10+10
      MemR=9*nZeta + 9*nEta +nZeta*nEta
      Call mma_allocate(Mem_INT,nZeta+nEta,Label='Mem_INT')
      ipIndZet=1
      ipIndEta=ipIndZet+nZeta
      Call mma_allocate(Mem_DBLE,MemR,Label='Mem_DBLE')
      ipZeta=1
*                                                                      *
************************************************************************
*                                                                      *
      If (lGrad) Then
*
*-----Calculate the size of memory needed for storing fock matrixes and
*     MO integrals and allocate it.
*
      nIndij=S%nShlls*(S%nShlls+1)/2
      nInt=0
      jDisp=0
      Do iIrrep=0,nIrrep-1
         Do iDisp=1,lDisp(iIrrep)
            jDisp=jDisp+1
            ipDisp(jDisp)=nInt+1
            Do jIrr=0,nIrrep-1
               kIrr=nrOpr(iEOr(iOper(iIrrep),iOper(jIrr)))
               If (jIrr.eq.kIrr) Then
                  nInt=nInt+nBas(jIrr)*(nBas(jIrr)+1)/2
               Else If (kIrr.lt.jIrr) Then
                  nInt=nInt+nBas(jIrr)*nBas(kIrr)
               End If
            End Do
*
            If (nMethod.eq.RASSCF) Then
               ipMO(jDisp,1)=nInt+1
               nInt=nInt+nMO(iIrrep)
               ipdisp2(jdisp)=nInt+1
               Do jIrr=0,nIrrep-1
                  kIrr=nrOpr(iEOr(iOper(iIrrep),iOper(jIrr)))
                  If (jIrr.eq.jIrr) Then
                     nInt=nInt+nBas(jIrr)*(nBas(jIrr)+1)/2
                  Else If (kIrr.lt.jIrr) Then
                     nInt=nInt+nBas(jIrr)*nBas(kIrr)
                  End If
               End Do
               ipMO(jDisp,2)=nInt+1-ipMO(jDisp,1)
            End If
*
         End Do
      End Do
      If (nMethod.eq.RASSCF) Then
         jDisp=0
         Do iIrrep=0,nIrrep-1
            Do iDisp=1,lDisp(iIrrep)
               jDisp=jDisp+1
               ipdisp3(jdisp)=nInt+1
               Do iS=0,nirrep-1
                  js=nrOpr(iEOr(iOper(is),iOper(iIrrep)))
                  nInt=nInt+nBas(iS)*nAsh(jS)
               End Do
            End Do
         End Do
      End If
      Call GetMem('Integrals','ALLO','REAL',ipInt,nInt)
      nTwo=0
      Do iIrrep=0,nIrrep-1
          nTwo=Max(nTwo,nFck(iIrrep))
      End Do
      If (Int_Direct) Then
         nTwo2=nInt
      Else
         nTwo2=nTwo
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Desymmetrize  densities.
*     Observe that the desymmetrized 1st order density matrices are
*     canonical, i.e. the relative order of the indices are canonically
*     ordered.
*
      ipDTemp=ip_Dummy
      ipDIN=ip_Dummy
      call dcopy_(nInt,[Zero],0,Work(ipInt),1)
      If (New_Fock) Then
         If (nmethod.ne.RASSCF) Then
            Call Get_D1ao_Var(ipDTemp,Length)
            Call DScal_(nDens,Half,Work(ipDTemp),1)
            ij=0
            Do i = 1, nBas(0)
               ij = ij + i
               Work(ipdtemp-1+ij)=Two*Work(ipDTemp-1+ij)
            End Do
         Else
            Call GetMem('DIN','Allo','Real',ipDIN,ndens)
            Call GetMem('DTemp','Allo','Real',ipDTemp,ndens)
            Call Din(Work(ipDIN))
            Call DScal_(nDens,Half,Work(ipDIN),1)
            ij=0
            Do i = 1, nBas(0)
               ij = ij + i
               Work(ipDIN-1+ij)=Two*Work(ipDIN-1+ij)
            End Do
            Call Dan(Work(ipDTemp))
            Call DScal_(nDens,Half,Work(ipDTemp),1)
            ij=0
            Do i = 1, nBas(0)
               ij = ij + i
               Work(ipDtemp-1+ij)=Two*Work(ipDTemp-1+ij)
            End Do
         End if
      Else
         mmdede=ndede
         Call mma_allocate(ipOffD,3,nIndij,label='ipOffD')
         Call mma_allocate(DeDe,mmDeDe+MxDij,label='DeDe')
         ipDijS = 1 + mmDeDe
         If (nMethod.ne.RASSCF) Then
            Call Get_D1ao_Var(ipDTemp,Length)
            Call DeDe_mck(Work(ipDTemp),nFck(0),ipOffD,nIndij,
     &                    Dede,mmDeDe,mDeDe,mIndij)
            Call GetMem('Dtemp','Free','Real',ipDTemp,ndens)
            ipDTemp=ip_Dummy
         Else
            Call mma_allocate(ipOffDA,3,nIndij,Label='ipOffDA')
            Call mma_allocate(DeDe2,mmDeDe+MxDij,label='DeDe2')
            ipDijS2 = 1 + mmDeDe
            Call GetMem('DIN','Allo','Real',ipDIN,ndens)
            Call GetMem('DTemp','Allo','Real',ipDTemp,ndens)
*
            Call Dan(Work(ipDTemp))
            Call DeDe_mck(Work(ipDTemp),nFck(0),ipOffD,nIndij,
     &                    DeDe,mmDeDe,mDeDe,mIndij)
            Call GetMem('Dtemp','Free','Real',ipDTemp,ndens)
            ipDTemp=ip_Dummy

            Call Din(Work(ipDIN))
            Call DeDe_mck(Work(ipDIN),nFck(0),ipOffDA,nIndij,
     &                    DeDe2,mmDeDe,mDeDe,mIndij)
            Call GetMem('DIN','Free','Real',ipDIN,ndens)
            ipDIN=ip_Dummy

            If (mDeDe.ne.nDeDe) Then
               Write (6,*) 'DrvG2: mDeDe.ne.nDeDe'
               Write (6,*) 'mDeDe,nDeDe=',mDeDe,nDeDe
               Call QTrace
               Call Abend
            End If
         End If
      End If
*
      nb=0
      Do is=0,nIrrep-1
         nb=nb+nBas(iS)
      End Do
*
      End If ! lGrad
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_iSD()
      Call Set_Basis_Mode('Valence')
      Call Nr_Shells(nSkal)
      Call Setup_iSD()
*
      nPairs=nSkal*(nSkal+1)/2
      nQuad=nPairs*(nPairs+1)/2
*                                                                      *
************************************************************************
*                                                                      *
*---  Compute entities for prescreening at shell level
*
      Call mma_allocate(TMax,nSkal,nSkal,Label='TMax')
      Call Shell_MxSchwz(nSkal,TMax)
      TMax_all=Zero
      Do iS = 1, nSkal
         Do jS = 1, iS
            TMax_all=Max(TMax_all,TMax(iS,jS))
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Create list of non-vanishing pairs
*
      Call mma_allocate(Ind_ij,2,nSkal*(nSkal+1)/2,Label='Ind_ij')
      nijS=0
      Do iS = 1, nSkal
         Do jS = 1, iS
            If (TMax_All*TMax(iS,jS).ge.CutInt) Then
               nijS = nijS + 1
               Ind_ij(1,nijS)=iS
               Ind_ij(2,nijS)=jS
            End If
         End Do
      End Do
      Call Init_Tsk(id_Tsk,nijS)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_MaxDBLE(MemMax)
      Call mma_allocate(Sew_Scr,MemMax-iii,Label='Sew_Scr')
      ipMem=1
      memmax=memmax-iii
*                                                                      *
************************************************************************
*                                                                      *
*     big loop over individual tasks, distributed over individual nodes
 10   Continue
*     make reservation of a task on global task list and get task range
*     in return. Function will be false if no more tasks to execute.
      If (.Not.Rsv_Tsk(id_Tsk,ijSh)) Go To 11
      iS = Ind_ij(1,ijSh)
      jS = Ind_ij(2,ijSh)
      Call CWTime(TCpu1,TWall1)
*                                                                      *
************************************************************************
*                                                                      *
*     Outer loops (ij) over angular momenta and centers
*
C     Do iS = 1, nSkal
         iShll  = iSD( 0,iS)
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iPrim  = iSD( 5,iS)
         iAO    = iSD( 7,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         iCnttp = iSD(13,iS)
         iCnt   = iSD(14,iS)
         Coor(1:3,1)=dbsc(iCnttp)%Coor(1:3,iCnt)
*
         iAngV(1) = iAng
         iShllV(1) = iShll
         iCmpV(1) = iCmp
         iShelV(1) = iShell
         iAOV(1) = iAO
*
C        Do jS = 1, iS
            jShll  = iSD( 0,jS)
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
            jPrim  = iSD( 5,jS)
            jAO    = iSD( 7,jS)
            mdcj   = iSD(10,jS)
            jShell = iSD(11,jS)
            jCnttp = iSD(13,jS)
            jCnt   = iSD(14,jS)
            Coor(1:3,2)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
            iAngV(2) = jAng
            iShllV(2) = jShll
            iCmpV(2) = jCmp
            iShelV(2) = jShell
            iAOV(2) = jAO
*
            ijAng = iAng + jAng
*
            nHrrab=0
            Do i=0,iAng+1
               Do j=0,jAng+1
                  If (i+j.le.iAng+jAng+1) Then
                     ijMax=Min(iAng,jAng)+1
                     nHrrab=nHrrab+ijMax*2+1
                  End If
               End Do
            End Do
            If (iShell.ge.jShell) Then
               ijShll = iShell*(iShell-1)/2+jShell
            Else
               ijShll = jShell*(jShell-1)/2+iShell
            End If
*                                                                      *
************************************************************************
*                                                                      *
*           Cltrls for MO transformation
*                                                                      *
************************************************************************
*                                                                      *
            If (nMethod.eq.RASSCF.and.l_Grd) Then
               iMemB=nACO**2*iCmp*iBas*jCmp*jBas*nDisp*nirrep
               If (iMemB.gt.MemMax) Then
                  Write (6,*) 'DrvG2: iMemB.gt.MemMax'
                  Write (6,*) 'iMemB=',iMemB
                  Write (6,*) 'MemMax=',MemMax
                  Write (6,*) 'Increase MOLCAS_MEM!'
                  Call QTrace()
                  Call Abend()
               End If
               Sew_Scr(1:iMemb)=Zero
            Else
               iMemb=0
            End If
*                                                                      *
************************************************************************
*                                                                      *
            Post_Process=.False.
            Do klSh = 1, nijS
               ks = Ind_ij(1,klSh)
               ls = Ind_ij(2,klSh)
*
               Aint=TMax(iS,jS)*TMax(kS,lS)
C              Write (*,*) 'is,js,ks,ls=',is,js,ks,ls
               If (AInt.lt.CutInt) Go To 400
*
C           Do kS = 1, nSkal
               kShll  = iSD( 0,kS)
               kAng   = iSD( 1,kS)
               kCmp   = iSD( 2,kS)
               kBas   = iSD( 3,kS)
               kPrim  = iSD( 5,kS)
               kAO    = iSD( 7,kS)
               mdck   = iSD(10,kS)
               kShell = iSD(11,kS)
               kCnttp = iSD(13,kS)
               kCnt   = iSD(14,kS)
               Coor(1:3,3)=dbsc(kCnttp)%Coor(1:3,kCnt)
*
               iAngV(3) = kAng
               iShllV(3) = kShll
               iCmpV(3) = kCmp
               iShelV(3) = kShell
               iAOV(3) = kAO
*
               Shik = iShell.eq.kShell
*
C              Do lS = 1, kS
                  lShll  = iSD( 0,lS)
                  lAng   = iSD( 1,lS)
                  lCmp   = iSD( 2,lS)
                  lBas   = iSD( 3,lS)
                  lPrim  = iSD( 5,lS)
                  lAO    = iSD( 7,lS)
                  mdcl   = iSD(10,lS)
                  lShell = iSD(11,lS)
                  lCnttp = iSD(13,lS)
                  lCnt   = iSD(14,lS)
                  Coor(1:3,4)=dbsc(lCnttp)%Coor(1:3,lCnt)
*
                  iAngV(4) = lAng
                  iShllV(4) = lShll
                  iCmpV(4) = lCmp
                  iShelV(4) = lShell
                  iAOV(4) = lAO
*
*
                  klAng = kAng + lAng
                  nHrrcd=0
                  Do i=0,kAng+1
                     Do j=0,lAng+1
                        If (i+j.le.kAng+lAng+1) Then
                           ijMax=Min(kAng,lAng)+1
                           nHrrcd=nHrrcd+ijMax*2+1
                        End If
                     End Do
                  End Do
*                                                                      *
************************************************************************
*                                                                      *
*-----------------Skip out if no symmetry and integral of
*                 one center type
*
                  If (kShell.ge.lShell) Then
                     klShll = kShell*(kShell-1)/2+lShell
                  Else
                     klShll = lShell*(lShell-1)/2+kShell
                  End If
*
*-----------------The code is working in such away that the MO needs
*                 upper and lower triangular parts of ij kl but hessian
*                 needs only lower, check if the integralbatch is lower
*                 or upper!!
*
                  lTri=iTri(iS,jS).ge.iTri(kS,lS)
                  If (.not.lTri.and.nMethod.ne.RASSCF) Goto 400
                  lDot=(lTri.and.l_Hss)
*
                  Shjl = jShell.eq.lShell
                  Shijij = Shik.and.Shjl
*                                                                      *
************************************************************************
*                                                                      *
              iCmpV(1)=icmp
              iCmpV(2)=jcmp
              iCmpV(3)=kcmp
              iCmpV(4)=lcmp
              iPrimi   = Shells(iShllV(1))%nExp
              jPrimj   = Shells(iShllV(2))%nExp
              kPrimk   = Shells(iShllV(3))%nExp
              lPriml   = Shells(iShllV(4))%nExp
              iBasi    = Shells(iShllV(1))%nBasis
              jBasj    = Shells(iShllV(2))%nBasis
              kBask    = Shells(iShllV(3))%nBasis
              lBasl    = Shells(iShllV(4))%nBasis
*                                                                      *
************************************************************************
*                                                                      *
*-------------Allocate memory for zeta, eta, kappa, P and Q.
*             Allocate also for Alpha, Beta , Gamma and Delta
*             in expanded form.
*
              nZeta = iPrimi * jPrimj
              nEta = kPrimk * lPriml
              MemR=9*nZeta + 9*nEta +nEta*nZeta
              ipZI  = ipZeta + nZeta
              ipKAB = ipZi   + nZeta
              ipP   = ipKAB  + nZeta
              ipxA  = ipP    + nZeta*3
              ipxB  = ipxA   + nZeta
              ipEta = ipxB   + nZeta
              ipEI  = ipEta  + nEta
              ipKCD = ipEI   + nEta
              ipQ   = ipKCD   + nEta
              ipxG  = ipQ    + nEta*3
              ipxD  = ipxG   + nEta
              ipxPre= ipxD   + nEta
*                                                                      *
************************************************************************
*                                                                      *
                  nab = nElem(iAng)*nElem(jAng)
                  ncd = nElem(kAng)*nElem(lAng)

                  ijS = iTri(iShell,jShell)
                  klS = iTri(kShell,lShell)
                  ikS = iTri(iShell,kShell)
                  ilS = iTri(iShell,lShell)
                  jkS = iTri(jShell,kShell)
                  jlS = iTri(jShell,lShell)
*                 If (.Not.l2DI) Then
*                    nab = 0
*                    ncd = 0
*                 End If
                  k2ij  = Indk2(1,ijS)
                  nDCRR = Indk2(2,ijS)
                  k2kl  = Indk2(1,klS)
                  nDCRS = Indk2(2,klS)
*
                  If (ltri) Then
*
*-------------------------------------------------------------------*
*
*                 Fix the 1st order density matrix
*
*-----------------Pick up pointers to desymmetrized 1st order density
*                 matrices. Observe that the desymmetrized 1st order
*                 density matrices follows the contraction index.
*
                  ipTmp =0
                  ipTmp2=0
                  If (lpick) Then
*
                  ipDij = ipOffD(1,ijS)
                  mDCRij= ipOffD(2,ijS)
                  nDij  = ipOffD(3,ijS)
*
                  ipTmp = ipDijs
                  If (nMethod.eq.RASSCF) Then
                     ipDij2 = ipOffDA(1,ijS)
                     ipTmp2= ipDijs2
                  End If
*
                  If (mDCRij.ne.0) Then
                     ipDDij = ipTmp
                     ipTmp = ipTmp + nDij*mDCRij
                     If (nMethod.eq.RASSCF) Then
                        ipDDij2=ipTmp2
                        ipTmp2= ipTmp2+ nDij*mDCRij
                     End If
                  Else
                     ipDDij = ip_Dummy
                  End If
*
                  ipDkl = ipOffD(1,klS)
                  If (nMethod.eq.RASSCF) ipDkl2 = ipOffDA(1,klS)
                  mDCRkl= ipOffD(2,klS)
                  nDkl  = ipOffD(3,klS)
                  If (mDCRkl.ne.0) Then
                     ipDDkl = ipTmp
                     ipTmp = ipTmp + nDkl*mDCRkl
                     If (nMethod.eq.RASSCF) Then
                       ipDDkl2=ipTmp2
                       ipTmp2= ipTmp2+ nDkl*mDCRkl
                     End If
                  Else
                     ipDDkl = ip_Dummy
                  End If
*
                  ipDik = ipOffD(1,ikS)
                  If (nMethod.eq.RASSCF) ipDik2 = ipOffDA(1,ikS)
                  mDCRik= ipOffD(2,ikS)
                  nDik  = ipOffD(3,ikS)
                  If (mDCRik.ne.0) Then
                     ipDDik = ipTmp
                     ipTmp = ipTmp + nDik*mDCRik
                    If (nMethod.eq.RASSCF) Then
                      ipDDik2=ipTmp2
                      ipTmp2= ipTmp2+ nDik*mDCRik
                    End If
                  Else
                     ipDDik = ip_Dummy
                  End If
*
                  ipDil = ipOffD(1,ilS)
                  If (nMethod.eq.RASSCF) ipDil2 = ipOffDA(1,ilS)
                  mDCRil= ipOffD(2,ilS)
                  nDil  = ipOffD(3,ilS)
                  If (mDCRil.ne.0) Then
                     ipDDil = ipTmp
                     ipTmp = ipTmp + nDil*mDCRil
                     If (nMethod.eq.RASSCF) Then
                      ipDDil2=ipTmp2
                      ipTmp2= ipTmp2+ nDil*mDCRil
                     End If
                  Else
                     ipDDil = ip_Dummy
                  End If
*
                  ipDjk = ipOffD(1,jkS)
                  If (nMethod.eq.RASSCF) ipDjk2 = ipOffDA(1,jkS)
                  mDCRjk= ipOffD(2,jkS)
                  nDjk  = ipOffD(3,jkS)
                  If (mDCRjk.ne.0) Then
                     ipDDjk = ipTmp
                     ipTmp = ipTmp + nDjk*mDCRjk
                     If (nMethod.eq.RASSCF) Then
                      ipDDjk2=ipTmp2
                      ipTmp2= ipTmp2 + nDjk*mDCRjk
                     End If
                  Else
                     ipDDjk = ip_Dummy
                  End If
*
                  ipDjl = ipOffD(1,jlS)
                  If (nMethod.eq.RASSCF) ipDjl2 = ipOffDA(1,jlS)
                  mDCRjl= ipOffD(2,jlS)
                  nDjl  = ipOffD(3,jlS)
                  If (mDCRjl.ne.0) Then
                     ipDDjl = ipTmp
                     ipTmp = ipTmp + nDjl*mDCRjl
                     If (nMethod.eq.RASSCF) Then
                      ipDDjl2=ipTmp2
                      ipTmp2= ipTmp2+ nDjl*mDCRjl
                     End If
                  Else
                     ipDDjl = ip_Dummy
                  End If
*
                  End If  ! If (lpick) Then
                  End If  ! If (ltri) Then
*                                                                      *
************************************************************************
*                                                                      *
*-----------------Compute total size of the second order density
*                 matrix in SO basis.
*
*----------------------------------------------------------------------*
                  nSO = MemSO2_P(iCmp,jCmp,kCmp,lCmp,
     &                           iAOV(1),iAOV(2),iAOV(3),iAOV(4))
                  ldot2=ldot
                  If (nSO.eq.0) ldot2=.false.
*
*-----------------Compute memory request for the primitives.
*
                  ider=2
                  if (.not.ldot2) iDer=1
                  Call MemRg2(iAngV,nRys,MemPrm,ider)
*
*----------------------------------------------------------------------*
*
*                 Calculate which derivatives that should be made.
*
*----------------------------------------------------------------------*
*
                  Call DerCtr(mdci,mdcj,mdck,mdcl,ldot2,JfGrd,
     &                        JndGrd,JfHss,JndHss,JfG,mBatch)
*

*----------------------------------------------------------------------*
*
*-----------------Decide on the partioning of the shells based on the
*                 available memory and the requested memory.
*
                  Call PSOAO2(nSO,MemPrm, MemMax,
     &                        iAngV, iCmpV, iAOV,iFnc,
     &                        iBasi,iBsInc, jBasj,jBsInc,
     &                        kBask,kBsInc, lBasl,lBsInc,
     &                        iPrimi,iPrInc,jPrimj,jPrInc,
     &                        kPrimk,kPrInc,lPriml,lPrInc,
     &                        nAco,
     &                        Mem1,Mem2,Mem3,Mem4,
     &                        MemX,MemPSO,
     &                        MemFck,nFT,memCMO2,MemFin,MemBuffer,
     &                        iMemB)

*
*----------------------------------------------------------------------*
*
*   Loop over basis function if we do not have enough of memory to
*   calculate them in one step.
*
*----------------------------------------------------------------------*
                  Do 500 iBasAO = 1, iBasi, iBsInc
                    iBasn=Min(iBsInc,iBasi-iBasAO+1)
                    iAOst(1) = iBasAO-1
*----------------------------------------------------------------------*
*
*----------------- Move appropiate portions of the desymmetrized 1st
*                  order density matrix.
*
*
*----------------------------------------------------------------------*
                  Do 510 jBasAO = 1, jBasj, jBsInc
                    jBasn=Min(jBsInc,jBasj-jBasAO+1)
                    iAOst(2) = jBasAO-1
                    If (lpick.and.nDij*mDCRij.ne.0) Then
                     Call Picky(DeDe(ipDij),iBasi,jBasj,
     &                         iPrimi*jPrimj,
     &                         iCmpV(1)*iCmpV(2),mDCRij,
     &                         iBasAO,iBasAO+iBasn-1,
     &                         jBasAO,jBasAO+jBasn-1,DeDe(ipDDij))
                    If (nMethod.eq.RASSCF)
     &              Call Picky(DeDe2(ipDij2),iBasi,jBasj,
     &                         iPrimi*jPrimj,
     &                         iCmpV(1)*iCmpV(2),mDCRij,
     &                         iBasAO,iBasAO+iBasn-1,
     &                         jBasAO,jBasAO+jBasn-1,DeDe2(ipDDij2))
                    End If
                    mDij = (iBasn*jBasn+1)*iCmpV(1)*iCmpV(2) +
     &                     iPrimi*jPrimj + 1
                    mDij = Min(nDij,mDij)
*
                  Do 520 kBasAO = 1, kBask, kBsInc
                    kBasn=Min(kBsInc,kBask-kBasAO+1)
                    iAOst(3) = kBasAO-1
                    If (lpick.and.nDik*mDCRik.ne.0) Then
                     Call Picky(DeDe(ipDik),iBasi,kBask,
     &                         iPrimi*kPrimk,
     &                         iCmpV(1)*iCmpV(3),mDCRik,
     &                         iBasAO,iBasAO+iBasn-1,
     &                         kBasAO,kBasAO+kBasn-1,DeDe(ipDDik))
                     If (nMethod.eq.RASSCF)
     &               Call Picky(DeDe2(ipDik2),iBasi,kBask,
     &                         iPrimi*kPrimk,
     &                         iCmpV(1)*iCmpV(3),mDCRik,
     &                         iBasAO,iBasAO+iBasn-1,
     &                         kBasAO,kBasAO+kBasn-1,DeDe2(ipDDik2))
                    End If
                    mDik = (iBasn*kBasn+1)*iCmpV(1)*iCmpV(3) +
     &                     iPrimi*kPrimk + 1
                    mDik = Min(nDik,mDik)
                    If (lpick.and.nDjk*mDCRjk.ne.0) Then
                     Call Picky(DeDe(ipDjk),jBasj,kBask,
     &                         jPrimj*kPrimk,
     &                         iCmpV(2)*iCmpV(3),mDCRjk,
     &                         jBasAO,jBasAO+jBasn-1,
     &                         kBasAO,kBasAO+kBasn-1,DeDe(ipDDjk))
                    If (nMethod.eq.RASSCF)
     &               Call Picky(DeDe2(ipDjk2),jBasj,kBask,
     &                         jPrimj*kPrimk,
     &                         iCmpV(2)*iCmpV(3),mDCRjk,
     &                         jBasAO,jBasAO+jBasn-1,
     &                         kBasAO,kBasAO+kBasn-1,DeDe2(ipDDjk2))
                    End If
                    mDjk = (jBasn*kBasn+1)*iCmpV(2)*iCmpV(3) +
     &                     jPrimj*kPrimk + 1
                    mDjk = Min(nDjk,mDjk)

                  Do 530 lBasAO = 1, lBasl, lBsInc
                    lBasn=Min(lBsInc,lBasl-lBasAO+1)
                    iAOst(4) = lBasAO-1
                    If (lpick.and.nDkl*mDCRkl.ne.0) Then
                    Call Picky(DeDe(ipDkl),kBask,lBasl,
     &                         kPrimk*lPriml,
     &                         iCmpV(3)*iCmpV(4),mDCRkl,
     &                         kBasAO,kBasAO+kBasn-1,
     &                         lBasAO,lBasAO+lBasn-1,DeDe(ipDDkl))
                    If (nMethod.eq.RASSCF)
     &              Call Picky(DeDe2(ipDkl2),kBask,lBasl,
     &                         kPrimk*lPriml,
     &                         iCmpV(3)*iCmpV(4),mDCRkl,
     &                         kBasAO,kBasAO+kBasn-1,
     &                         lBasAO,lBasAO+lBasn-1,DeDe2(ipDDkl2))
                    End If
                    mDkl = (kBasn*lBasn+1)*iCmpV(3)*iCmpV(4) +
     &                     kPrimk*lPriml + 1
                    mDkl = Min(nDkl,mDkl)
                    If (lpick.and.nDil*mDCRil.ne.0) Then
                    Call Picky(DeDe(ipDil),iBasi,lBasl,
     &                         iPrimi*lPriml,
     &                         iCmpV(1)*iCmpV(4),mDCRil,
     &                         iBasAO,iBasAO+iBasn-1,
     &                         lBasAO,lBasAO+lBasn-1,DeDe(ipDDil))
                    If (nMethod.eq.RASSCF)
     &              Call Picky(DeDe2(ipDil2),iBasi,lBasl,
     &                         iPrimi*lPriml,
     &                         iCmpV(1)*iCmpV(4),mDCRil,
     &                         iBasAO,iBasAO+iBasn-1,
     &                         lBasAO,lBasAO+lBasn-1,DeDe2(ipDDil2))
                    End If
                    mDil = (iBasn*lBasn+1)*iCmpV(1)*iCmpV(4) +
     &                     iPrimi*lPriml + 1
                    mDil = Min(nDil,mDil)
                    If (lpick.and.nDjl*mDCRjl.ne.0) Then
                    Call Picky(DeDe(ipDjl),jBasj,lBasl,
     &                         jPrimj*lPriml,
     &                         iCmpV(2)*iCmpV(4),mDCRjl,
     &                         jBasAO,jBasAO+jBasn-1,
     &                         lBasAO,lBasAO+lBasn-1,DeDe(ipDDjl))
                    If (nMethod.eq.RASSCF)
     &              Call Picky(DeDe2(ipDjl2),jBasj,lBasl,
     &                         jPrimj*lPriml,
     &                         iCmpV(2)*iCmpV(4),mDCRjl,
     &                         jBasAO,jBasAO+jBasn-1,
     &                         lBasAO,lBasAO+lBasn-1,DeDe2(ipDDjl2))
                    End If
                    mDjl = (jBasn*lBasn+1)*iCmpV(2)*iCmpV(4) +
     &                     jPrimj*lPriml + 1
                    mDjl = Min(nDjl,mDjl)
                    If (.not.lpick) Then
                     ipddjl1=ip_Dummy
                     ipddjl2=ip_Dummy
                     ipddil1=ip_Dummy
                     ipddil2=ip_Dummy
                     ipddkl1=ip_Dummy
                     ipddkl2=ip_Dummy
                     ipddij1=ip_Dummy
                     ipddij2=ip_Dummy
                     ipddik1=ip_Dummy
                     ipddik2=ip_Dummy
                     ipddjk1=ip_Dummy
                     ipddjk2=ip_Dummy
                    End If
*
*----------------------------------------------------------------------*
*
                    MEMCMO=nACO*(kCmp*kBasn+lCmp*lBasn)
*.................. MO tranformation buffer
                    ipBuffer = ipMem
                    ipMOC    = ipBuffer + MEMBUFFER
*.................. Area for the AO integrals
                    ipFin    = ipMOC    + MemCMO
*.................. Area for 2el density
                    ip_PP     = ipFin    + MemFin
                    ipMem2   = ip_PP     + Mem1         ! Work
                    ipMem3   = ipMem2   + Mem2          ! Work
                    ipMemX   = ipMem3   + Mem3          ! Work
*
*-------------------If MO transformation is performed in the standard way
*                   reserve memory for partial transfromed integrals
*
*
*-------------------Multilayer
*
                    ipMem4   = ipMem2   + Mem2 - Mem4

*
*----------------------------------------------------------------------*
*
*-------------------Get the 2nd order density matrix in SO basis.
*
*----------------------------------------------------------------------*
*

                    nijkl = iBasn*jBasn*kBasn*lBasn
                    Call Timing(dum,Time,Dum,Dum)
                    If (n8)
     &              Call PickMO(Sew_Scr(ipMOC),MemCMO,nAcO,iCmpV,
     &                          iBasAO,iBasn,jBasAO,jBasn,
     &                          kBasAO,kBasn,lBasAO,lBasn,iAOV)
                    If (ldot2)
     &               Call PGet0(iCmpV,
     &                         iBasn,jBasn,kBasn,lBasn,Shijij,
     &                         iAOV,iAOst,nijkl,Sew_Scr(ip_PP),nSO,
     &                         iFnc(1)*iBasn,iFnc(2)*jBasn,
     &                         iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,
     &                         Sew_Scr(ipMem2),Mem2,
     &                         iS,jS,kS,lS,nQuad,PMax)
                    Call Timing(dum,Time,Dum,Dum)
                    CPUStat(nTwoDens)=CPUStat(nTwoDens)+Time
*
*-------------------Compute gradients of shell quadruplet
*
                    ipD0=ip_of_Work(D0(1,1))
                    Call TwoEl_mck(Coor,iAngV,iCmpV,iShelV,iShllV,iAOV,
     &                   iAOst,mdci,mdcj,mdck,mdcl,nRys,
     &                   Data_k2(k2ij),nab,nDCRR,
     &                   Data_k2(k2kl),ncd,nDCRS,Pren,Prem,
     &            Shells(iShllV(1))%Exp,iPrimi,iPrInc,
     &            Shells(iShllV(2))%Exp,jPrimj,jPrInc,
     &            Shells(iShllV(3))%Exp,kPrimk,kPrInc,
     &            Shells(iShllV(4))%Exp,lPriml,lPrInc,
     &            Shells(iShllV(1))%pCff(1,iBasAO),iBasn,
     &            Shells(iShllV(2))%pCff(1,jBasAO),jBasn,
     &            Shells(iShllV(3))%pCff(1,kBasAO),kBasn,
     &            Shells(iShllV(4))%pCff(1,lBasAO),lBasn,
     &           Mem_DBLE(ipZeta),Mem_DBLE(ipZI),
     &           Mem_DBLE(ipP),Mem_DBLE(ipKab),nZeta,
     &           Mem_DBLE(ipEta), Mem_DBLE(ipEI),
     &           Mem_DBLE(ipQ),Mem_DBLE(ipKcd),nEta,
     &           Mem_DBLE(ipxA),Mem_DBLE(ipxB),
     &           Mem_DBLE(ipxG),Mem_DBLE(ipxD),
     &                   Mem_DBLE(ipxPre),
     &                   Hess, nhess,JfGrd,JndGrd,JfHss,JndHss,JfG,
     &                   Sew_Scr(ip_PP), nSO,Sew_Scr(ipMem2),Mem2,
     &                   Sew_Scr(ipMem3),Mem3,Sew_Scr(ipMem4),Mem4,
     &                   Aux,nAux,Sew_Scr(ipMemX),MemX,Shijij,
     &                   DeDe(ipDDij),DeDe2(ipDDij2),mDij,mDCRij,
     &                   DeDe(ipDDkl),DeDe2(ipDDkl2),mDkl,mDCRkl,
     &                   DeDe(ipDDik),DeDe2(ipDDik2),mDik,mDCRik,
     &                   DeDe(ipDDil),DeDe2(ipDDil2),mDil,mDCRil,
     &                   DeDe(ipDDjk),DeDe2(ipDDjk2),mDjk,mDCRjk,
     &                   DeDe(ipDDjl),DeDe2(ipDDjl2),mDjl,mDCRjl,
     &                   iCmpV,Sew_Scr(ipFin),MemFin,
     &                   Sew_Scr(ipMem2),Mem2+Mem3+MemX,nTwo2,nFT,
     &                   Mem_INT(ipIndEta),Mem_INT(ipIndZet),
     &                   Work(ipInt),ipd0,Sew_Scr(ipBuffer),MemBuffer,
     &                   lgrad,ldot2,n8,ltri,Work(ipDTemp),Work(ipDIN),
     &                   moip,nAco,Sew_Scr(ipMOC),MemCMO,new_fock)
                  Post_Process=.True.

*----------------------------------------------------------------------*
*
 530              Continue
 520              Continue
 510              Continue
 500              Continue
*
 400              Continue
C              End Do ! lS
C           End Do ! kS
            End Do ! klS
*
            If (nMethod.eq.RASSCF.and.Post_Process) Then
               ip1=ipMOC
               ip2=ip1+iCmp*iBas*naco
               ip3=ip2+nAco**2
               ip4=ip3+jcmp*jBas*naco
               ip5=ip4+iCmp*naco*iBas
               ip6=ip5+jcmp*jbas*naco
               Call CLR2(Sew_Scr(ipBuffer),Work(ipInt),
     &                   ibas,icmp,jbas,jcmp,iAOV(1),iAOV(2),
     &                   naco,ishelV,
     &                   Sew_Scr(ip1),Sew_Scr(ip2),Sew_Scr(ip3),
     &                   Sew_Scr(ip4),Sew_Scr(ip5),Sew_Scr(ip6))
            End If
*
C        End Do ! jS
C     End Do !  iS
*
      Call CWTime(TCpu2,TWall2)
      Call SavTim(4,TCpu2-TCpu1,TWall2-Twall1)
      Call SavStat(1,One,'+')
      Call SavStat(2,DBLE(nijs),'+')
      Go To 10
 11   Continue
*     End of big task loop
*                                                                      *
************************************************************************
*                                                                      *
* - - - - - - - E P I L O G
*                                                                      *
************************************************************************
*                                                                      *
      If (New_Fock) Then
         idd=0
         Do iS=0,nirrep-1
           Do iD=1,ldisp(is)
            idd=idd+1
            ip=ipInt-1+ipDisp(idd)
            Call DScal_(nDens,Half,work(ip),1)
            ij =ip-1
            Do i = 1, nBas(0)
             ij=ij+i
             Work(ij)=Two*Work(ij)
            End Do
           End Do
          End Do
         If (nmethod.eq.RASSCF) Then
         idd=0
         Do iS=0,nirrep-1
           Do iD=1,ldisp(is)
            idd=idd+1
            ip=ipInt-1+ipDisp2(idd)
            Call DScal_(nDens,Half,work(ip),1)
            ij =ip-1
            Do i = 1, nBas(0)
             ij=ij+i
             Work(ij)=Two*Work(ij)
            End Do
           End Do
          End Do

      End If
      End If
#ifdef _DEBUG_
      Call GADSum_SCAL(Pren)
      Call GADSum_SCAL(Prem)
      Write (Format,'(A,I2,A,I2,A)') '(A,F',
     &              3+Int(Log10(Pren)),
     &              '.0,A,F',
     &              3+Int(Log10(Prem)),
     &              '.0,A)'
      Write (6,Format)
     &   ' A total of', Pren,' entities were prescreened and',
     &                  Prem,' were kept.'
#endif
      Call mma_deallocate(Sew_Scr)
      Call Free_Tsk(id_Tsk)
*
*    YIPPIEEEE Finished OK fill it UP!!
*
      Call GADSum(Work(ipInt),nInt)
      jDisp=0
      Do iIrr=0,nIrrep-1
        Do iDisk=1,lDisp(iIrr)
         jDisp=jDisp+1
           Call WrDisk(Work(ipInt),nInt,jdisp,iIrr)
        End Do
      End Do
*
      Call mma_deallocate(Ind_ij)
      Call mma_deallocate(TMax)
      Call Free_iSD()
*
      If (.not.New_Fock) Then
         Call mma_deallocate(ipOffD)
         Call mma_deallocate(DeDe)
         If (nMethod.eq.RASSCF) Then
            Call mma_deallocate(DeDe2)
            Call mma_deallocate(ipOffDA)
         End If
      End If
*
      Call mma_deallocate(Mem_DBLE)
      Call mma_deallocate(Mem_INT)
*
      If (ipDIN.ne.ip_Dummy)
     &   Call GetMem('DIN','Free','Real',ipDIN,ndens)
      If (ipDTemp.ne.ip_Dummy)
     &   Call GetMem('DTemp','Free','Real',ipDTemp,ndens)
      Call GetMem('Integrals','Free','REAL',ipInt,nInt)
*
      Call mma_deallocate(Aux)
*
*-----Generate statistic of partioning
*
      Call mma_deallocate(IndK2)
      Call mma_deallocate(Data_k2)
*
      Return
      End
