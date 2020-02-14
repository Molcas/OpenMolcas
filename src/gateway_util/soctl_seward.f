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
      Subroutine SOCtl_Seward(Mamn,nMamn,DInf,nDInf,Info)
      Implicit Real*8 (a-h,o-z)
*
#include "itmax.fh"
#include "info.fh"
#include "rinfo.fh"
#include "real.fh"
#include "print.fh"
#include "stdalloc.fh"
*
      Character ChOper(0:7)*3,ChTemp*8,Mamn(nMamn)*(LENIN8)
      Character LP_Names(MxAtom)*(LENIN4)
      Character*60 Fmt
      Logical Type(0:7), lSkip, kECP, TstFnc, output, Get_BasisType
      Logical IsBasisAE
      Logical IsBasisANO
      Logical IsBasisUNK
      Integer Occ, Vir
      Parameter(Occ=1,Vir=0)
      Integer List(0:iTabMx), nFCore(0:7), nCore_Sh(0:iTabMx),
     &        List_AE(0:iTabMx)
      Integer jOffSO(0:7)
      Integer, Dimension(:), Allocatable :: Index, Index2, IndC, iCI,
     &                                      jCI, iOT, LPA, LPMM
      Real*8, Dimension(:), Allocatable :: LPQ
      Real*8, Dimension(:,:), Allocatable :: SM, LPC
      Real*8 DInf(nDInf)
      Character*(LENIN8) Clean_BName,ChTmp
      External Clean_BName

CSVC: the basis ids are tuples (c,n,l,m) with c the center index,
C     n the shell index, l the angmom value, and m the angmom component.
C     the angmom components of p are mapped (x,y,z) -> (1,-1,0)
C     examples: 3d1+ on atom 1: (1,3,2,1); 2py on atom 5: (5,2,1,-1)
CIFG: for Cartesian shells, l -> -l, m -> T(ly+lz)-(lx+ly), where T(n)=n*(n+1)/2
      integer :: llab,mlab
      integer, allocatable :: basis_ids(:,:), desym_basis_ids(:,:)
      integer, allocatable :: fermion_type(:)
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*     LVAL end MVAL dimensioned for L = MxAng
      dimension LVAL((MxAng+1)*(MxAng+1))
      dimension MVAL((MxAng+1)*(MxAng+1))
*                                                                      *
************************************************************************
*                                                                      *
      IsBasisAE=.false.
      IsBasisANO=.false.
      IsBasisUNK=.false.
      iRout=2
      iPrint = nPrint(iRout)
      Call qEnter('SOCtl')
cvv LP_NAMES was used later without initialization.
      do i=1,MxAtom
       LP_NAMES(i)(1:LENIN)='crap'
       LP_NAMES(i)(LENIN1:LENIN4)='crap'
      enddo

*                                                                      *
************************************************************************
*                                                                      *
*     Compute iBas, iBas_Aux, and iBas_Frag used for double checking
*     in SOCtl.
*     Compute cdMax, EtMax, IndS(iShell), nShlls, and
*     Ind_Shell(IndSOff(iCnttp,iCnt)).
*
      Call Misc_Seward(iBas,iBas_Aux,iBas_Frag,DInf,nDInf)
*                                                                      *
************************************************************************
*                                                                      *
*     initialize LVAL and MVAL
*     (note: this is wrong for Cartesian shells)
*
      k=0
      do i=0,MxAng
         do j=-i,i
            k=k+1
            lval(k)=i
            mval(k)=j
         enddo
      enddo
C     write(6,*) ' lval',k,(MxAng+1)**2
*     correct mval order for p-functions
      mval(2)=1
      mval(3)=-1
      mval(4)=0
      Call ICopy(MxAO,[-99],0,iCent,1)
      Call ICopy(MxAO,[-99],0,lnAng,1)
C     write(6,'(20i4)') (lval(i),i=1,k)
C     write(6,*) ' lval',k
C     write(6,'(20i4)') (mval(i),i=1,k)
*
      Call ICopy(1+iTabMx,[0],0,List   ,1)
      Call ICopy(1+iTabMx,[0],0,List_AE,1)
*
      isymunit=isfreeunit(58)
      call molcas_open(isymunit,'SYMINFO')
      rewind isymunit
      write(isymunit,'(A)') 'Symmetry information from seward'
      write(isymunit,'(A)')
     &'#of funct, unique centre, L, M , # of sym.ad.functions , Phases'
C     write(6,*) 'Symmetry info to file SYMINFO '
C     Show=.Not.Prprt
C     Show=Show.and..Not.Primitive_Pass
*                                                                      *
************************************************************************
*                                                                      *
*     Generate list of symmetry adapted or petite list basis functions
*                                                                      *
************************************************************************
*                                                                      *
      iSO = 0
      iSO_Aux=iBas+iBas_Frag
      iSO_Frag = 0
      iSO_Tot=0
      n2Tot = 0
      n2Max = 0
      nDim = 0
      m2Tot = 0
      iEMax = 0
      iAO=0
      lSkip=.False.
*
      Call ICopy(8,[0],0,nFCore,1)

      Call mma_Allocate(iCI,iBas,label='iCI')       ! Stuff for LoProp
      Call mma_Allocate(jCI,iBas,label='jCI')
!     Stuff for LocalDKH/X2C/BSS
      Call mma_Allocate(iOT,iBas,label='iOT')       ! Stuff for LoProp
      Call mma_Allocate(LPC,3,mCentr,label='LPC')
!     Stuff (not just) for LoProp
      Call mma_Allocate(LPQ,mCentr,label='LPQ')
!     Stuff (not just) for LoProp
      Call mma_Allocate(LPMM,mCentr,label='LPMM')
!     Stuff (not just) for LoProp
      Call mma_Allocate(LPA,mCentr,label='LPA')
      call mma_allocate(basis_ids,4,maxbfn+maxbfn_aux)
      call mma_allocate(desym_basis_ids,4,maxbfn+maxbfn_aux)
      call mma_allocate(fermion_type,maxbfn+maxbfn_aux)
*
      IsBasisAE  = Get_BasisType('AE_')
      IsBasisANO = Get_BasisType('ANO')
      IsBasisUNK = Get_BasisType('UNK')
      If (Show.and.iPrint.ge.6) then
         Write (6,*)
         Call CollapseOutput(1,'   SO/AO info:')
         Write (6,'(3X,A)')    '   -----------'
      End If
      If (Petite) Go To 199
*                                                                      *
************************************************************************
*                                                                      *
*---- Symmetry case.
*
      If (Show.and.iPrint.ge.6) then
         Write (6,*)
         Write (6,'(19x,a)')
     &            ' **************************************************'
         Write (6,'(19x,a)')
     &            ' ******** Symmetry adapted Basis Functions ********'
         Write (6,'(19x,a)')
     &            ' **************************************************'
         Write (6,*)
      End If
*
      Call mma_allocate(Index,5*iBas,label='Index')
      Call mma_allocate(Index2,5*iBas,label='Index2')
      Call ICopy(5*iBas,[0],0,Index,1)
      Call ICopy(5*iBas,[0],0,Index2,1)
      iCounter=0
      jCounter=0
      Call mma_Allocate(SM,iBas,iBas,label='SM')
      Call FZero(SM,iBas**2)
      Call mma_Allocate(IndC,2*mCentr)
      iAtoms=0
*
*     Loop over irreducible representations and symmetry operations,
*     respectively, for SO and Petite list, respectively.
*
      Do iIrrep = 0, 7
        jOffSO(iIrrep) = 0
      End Do
      Do i = 1, MxUnq
         IrrCmp(i) = 0
      End Do
      kIrrep=0
      Do 200 iIrrep = 0, nIrrep-1
         iOffSO(iIrrep) = iSO_Tot
         jOffSO(iIrrep) = iSO
         iAO = 0
         jSO = 0
         nBas(iIrrep) = 0
         nBas_Aux(iIrrep) = 0
         nBas_Frag(iIrrep) = 0
         nPrm(iIrrep) = 0
         Type(iIrrep)=.True.
*
*        Loop over distinct shell types
*
         mc  = 1
         iShell = 0
         If (iSkip(iIrrep).ne.0) Then
            Write (6,*)
            Write (6,*) ' All basis functions of Irrep', iIrrep+1,
     &                  ' are removed!'
            Write (6,*)
            lSkip=.True.
            Go To 2011
         End If
         iCnttp = 0
         Do 201 jCnttp = 1, nCnttp
*
*           Make sure that we process the dummy shell last
*
            If (jCnttp.eq.iCnttp_Dummy .and. jCnttp.ne.nCnttp) Then
               iCnttp = iCnttp + 2
            Else If (jCnttp.eq.nCnttp .and. iCnttp.eq.jCnttp) Then
               iCnttp = iCnttp_Dummy
            Else
               iCnttp = iCnttp + 1
            End If
*
            output = show .and. iPrint.ge.6
            If (AuxCnttp(iCnttp)) output=output .and. iPrint.ge.10
     &                            .and. iCnttp.ne.iCnttp_Dummy
            If (FragCnttp(iCnttp)) output=output .and. iPrint.ge.10
            kECP = ECP(iCnttp)
            lMax=nVal_Shells(iCnttp)-1
*
            Call OrbType(iAtmNr(iCnttp),List_AE,31)
            If (kECP) Then
*
*              ECP case
*
               Call ECP_Shells(iAtmNr(iCnttp),list)
*
*              No core to freeze!
*
               Call ICopy(lMax+1,[0],0,nCore_Sh,1)
            Else
*
*              Non-ECP case
*
*              Pick up the number of occupied orbitals in each shell type.
*
               Call ICopy(1+iTabMx,List_AE,1,List,1)
*
*              Pick up which orbitals should be frozen as default.
*
               If (Charge(iCnttp).ne.Zero) Then
                  Call Freeze_Default(iAtmNr(iCnttp),nCore_Sh,lMax)
               Else
*
*                 If there charge is zero we presume that these are
*                 ghost orbitals or something else. In any case
*                 we do not freeze any orbitals!
*
                  Call Freeze_Default(0             ,nCore_Sh,lMax)
               End If
            End If
*
*           Loop over distinct centers
*
            Do 202 iCnt = 1, nCntr(iCnttp)
               mdc = iCnt + mdciCnttp(iCnttp)
*
*              Loop over shells associated with this center
*              Start with s type shells
*
               kComp = 0
               kculf = 0
               iSh = ipVal(iCnttp) - 1
               If (nVal_Shells(iCnttp).lt.1) Then
                  Do iCo = 0, nIrrep/nStab(mdc)-1
                     iyy=Index_Center(mdc,iCo,IndC,iAtoms,mCentr)
                     iR=NrOpr(iCoSet(iCo,0,mdc),iOper,nIrrep)
                     ipxyz=(iCnt-1)*3+ipCntr(iCnttp)
                     XCoor=Dinf(ipxyz  )
                     If (iAnd(iOper(iR),1).ne.0) XCoor=-XCoor
                     YCoor=Dinf(ipxyz+1)
                     If (iAnd(iOper(iR),2).ne.0) YCoor=-YCoor
                     ZCoor=Dinf(ipxyz+2)
                     If (iAnd(iOper(iR),4).ne.0) ZCoor=-ZCoor
                     LPC(1,iyy)=XCoor
                     LPC(2,iyy)=YCoor
                     LPC(3,iyy)=ZCoor
                     LPQ(iyy)=Charge(iCnttp)
                     LPA(iyy)=iAtmnr(iCnttp)
                     LPMM(iyy)=IsMM(iCnttp)
                     LP_Names(iyy)=LblCnt(mdc)(1:LENIN)//':'
     &                       //ChOper(iOper(iR))
                  End Do
               End If
               Do 203 iAng = 0, nVal_Shells(iCnttp)-1
                  nCore=nCore_Sh(iAng)
                  iSh = iSh + 1
                  iShell = iShell + 1
                  If (nExp(iSh).eq.0) Go To 2033
                  If (nBasis(iSh).eq.0) Go To 2033
                  jComp = (iAng+1)*(iAng+2)/2
                  If(Prjct(iSh)) jComp = 2*iAng + 1
                  Do 204 iComp = 1, jComp
                     iAO = iAO + 1
                     If (iAO.gt.MxAO) Then
                        Call ErrTra
                        Write (6,*) ' Increase MxAO'
                        Call Abend
                     End If
                     lComp = kComp + iComp
                     lculf = kculf + icomp
*                    Get character of basis function
                     iChBs = iChBas(lComp)
                     If (Transf(iSh)) iChBs=iChBas(iSphCr(lComp))
*
*                    Skip if function not a basis of irreps.
*
                     If (.Not.TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                   nIrrep/nStab(mdc),iChTbl,iIrrep,iChBs,
     &                   nStab(mdc))) Go To 204
                     If(.not.FragShell(iSh) .and.
     &                  .not.AuxCnttp(iCnttp))
     &                 nFCore(iIrrep)=nFCore(iIrrep)+nCore
                     iEMax = Max(iEMax,IndS(iShell)+iComp)
                     If (IndS(iShell)+iComp.gt.MxUnq) Then
                        Call ErrTra
                        Write (6,*) ' Increase MxUnq'
                        Call Abend
                     End If
                     If (iSkip(iIrrep).eq.0) Then
                        IrrCmp(IndS(iShell)+iComp) =
     &                    iOr(IrrCmp(IndS(iShell)+iComp),2**iIrrep)
                     End If
                     If (output.and.Type(iIrrep)) Then
                        Write (6,*)
                        Write (6,'(10X,A,A)')
     &                    ' Irreducible representation : ',
     &                      lIrrep(iIrrep)
                        Write (6,'(10X,2A)')
     &                      ' Basis function(s) of irrep: ',
     &                       lBsFnc(iIrrep)
                        Write (6,*)
                        Write (6,'(A)')
     &                  ' Basis Label        Type   Center Phase'
                        Type(iIrrep)=.False.
                     End If
*
                     If (MaxBas(iAng).gt.0) iAOtSO(iAO,iIrrep) = jSO + 1
                     nPrm(iIrrep) = nPrm(iIrrep) + nExp(iSh)
                     m2Max = Max(m2Max,nExp(iSh)**2)
                     Do 205 iCntrc = 1, nBasis(iSh)
                        iSO_Tot = iSO_Tot + 1
                        If (AuxShell(iSh)) Then
                           iSO_Aux = iSO_Aux + 1
                           iSO_=iSO_Aux
                           nBas_Aux(iIrrep) = nBas_Aux(iIrrep) + 1
                        Else If (FragShell(iSh)) Then
                           iSO_Frag = iSO_Frag + 1
                           iSO_=iSO_Frag
                           nBas_Frag(iIrrep) = nBas_Frag(iIrrep) + 1
                        Else
                           iSO = iSO + 1
                           iSO_=iSO
                           nBas(iIrrep) = nBas(iIrrep) + 1
                        End If
                        If (iSO_.gt.nMamn) Then
                           Call qTrace
                           Write (6,*) ' iSO_.gt.nMamn'
                           Write (6,*) 'nMamn=',nMamn
                           Call Abend
                        End If
                        jSO = jSO + 1
*
                        ChTemp=LblCBs(lComp)
                        If (Transf(iSh)) ChTemp=LblSbs(lComp)
*
                        Call Name_to_lm(ChTemp,llab,mlab)
*
*                       Introduce a somewhat better labelling. Thnx LG!
*
                        If (IsBasisAE) Then
                           If (IsBasisANO) Then
                              Write (ChTemp(1:2),'(I2.2)') iAng+iCntrc
                           Else
                              If (nExp(iSh).eq.nBasis(iSh)) Then
                                 Write (ChTemp(1:1),'(A1)') '*'
                                 If (llab.ge.0)
     &                              Write(ChTemp(2:2),'(A1)') '0'
                              Else If (iCntrc.le.list(iAng)) Then
                                 Write (ChTemp(1:2),'(I2.2)')
     &                                 iAng+iCntrc
                              Else
                                 Write (ChTemp(1:1),'(A1)') '*'
                                 If (llab.ge.0)
     &                              Write(ChTemp(2:2),'(A1)') '0'
                              End If
                           End If
                        Else If (.Not.IsBasisUNK) Then
                           If (nExp(iSh).eq.nBasis(iSh)) Then
                              Write (ChTemp(1:1),'(A1)') '*'
                              If (llab.ge.0)
     &                           Write(ChTemp(2:2),'(A1)') '0'
                           Else If (iCntrc.le.list(iAng)) Then
                                 Write (ChTemp(1:2),'(I2.2)')
     &                                 iAng+iCntrc+
     &                                 (List_AE(iAng)-List(iAng))
                           Else
                              Write (ChTemp(1:1),'(A1)') '*'
                              If (llab.ge.0)
     &                           Write(ChTemp(2:2),'(A1)') '0'
                           End If
*
                        End If
                        ChTmp=Clean_BName(ChTemp,0)
*
                        If(output)
     &                  Write (6,'(I5,3X,A8,4X,A8,8(I3,4X,I2,4X))')
     &                        iSO_,LblCnt(mdc),ChTmp,
     &                        (mc+iCo,iPrmt(NrOpr(iCoSet(iCo,0,mdc),
     &                        iOper,nIrrep),iChbs)*
     &                        iChTbl(iIrrep,NrOpr(iCoSet(iCo,0,mdc),
     &                        iOper,nIrrep)),
     &                        iCo=0,nIrrep/nStab(mdc)-1 )
*
                        If (iSO_.gt.4*MxAO) Then
                           Write (6,*) 'iSO_.gt.2*MxAO'
                           Call Abend()
                        End If
                        iSOInf(1,iSO_)=iCnttp
                        iSOInf(2,iSO_)=iCnt
                        iSOInf(3,iSO_)=iAng
*
                        If (AuxShell(iSh).or.FragShell(iSh)) Go To 205
*
                        If (.Not.Primitive_Pass) Then
                           Write (isymunit,'(13(I4,4X))')
     &                        iSO_,mdc,LVAL(lculf),MVAL(lculf),
     &                        nIrrep/nStab(mdc),
     &                        (iPrmt(NrOpr(iCoSet(iCo,0,mdc),
     &                        iOper,nIrrep),iChbs)*
     &                        iChTbl(iIrrep,NrOpr(iCoSet(iCo,0,mdc),
     &                        iOper,nIrrep)),
     &                        iCo=0,nIrrep/nStab(mdc)-1 )
                        End If
*                                                                      *
************************************************************************
*                                                                      *
*---------------------------- Stuff (not just) for LoProp
*
                         Do iCo = 0, nIrrep/nStab(mdc)-1
                            ixxx = Index_NoSym(iCntrc,iComp,iAng,
     &                        mdc,iCo,Index,iCounter,iBas)
                            jxxx = Index_NoSym(iCntrc,iComp,iAng,
     &                        mdc,iirrep,Index2,jCounter,iBas)
                            fact =DBLE(iPrmt(NrOpr(iCoSet(iCo,0,mdc),
     &                        iOper,nIrrep),iChbs)*
     &                        iChTbl(iIrrep,NrOpr(iCoSet(iCo,0,mdc),
     &                        iOper,nIrrep)))
*
                            FacN = One/DBLE(nIrrep/nStab(mdc))
                            If (MolWgh.eq.1) Then
                               FacN= One
                            Else If (MolWgh.eq.2) Then
                               FacN= Sqrt(FacN)
                            End If
                            SM(ixxx,iSO)=Fact*FacN
                            iyy=Index_Center(mdc,iCo,IndC,iAtoms,mCentr)
*
                            iCI(ixxx)=iyy
                            jCI(jxxx)=icnt
*
                            If (iCntrc.le.list(iAng)) Then
                               iOT(ixxx)=Occ
                            Else
                               iOT(ixxx)=Vir
                            End If
*
                            iR=NrOpr(iCoSet(iCo,0,mdc),iOper,
     &                               nIrrep)
                            ipxyz=(iCnt-1)*3+ipCntr(iCnttp)
                            XCoor=Dinf(ipxyz  )
                            If (iAnd(iOper(iR),1).ne.0) XCoor=-XCoor
                            YCoor=Dinf(ipxyz+1)
                            If (iAnd(iOper(iR),2).ne.0) YCoor=-YCoor
                            ZCoor=Dinf(ipxyz+2)
                            If (iAnd(iOper(iR),4).ne.0) ZCoor=-ZCoor
                            LPC(1,iyy)=XCoor
                            LPC(2,iyy)=YCoor
                            LPC(3,iyy)=ZCoor
*
                            LPQ(iyy)=Charge(iCnttp)
                            LPMM(iyy)=IsMM(iCnttp)
                            LPA(iyy)=iAtmnr(iCnttp)
*
                            LP_Names(iyy)=LblCnt(mdc)(1:LENIN)//':'
     &                                    //ChOper(iOper(iR))
                            desym_basis_ids(1,ixxx) = iyy
                            desym_basis_ids(2,ixxx) = iCntrc
                            desym_basis_ids(3,ixxx) = llab
                            desym_basis_ids(4,ixxx) = mlab
                        End Do
*                                                                      *
************************************************************************
*                                                                      *
                        Mamn(iSO)=LblCnt(mdc)(1:LENIN)//ChTemp(1:8)
                        basis_ids(1,iSO) = mdc
                        basis_ids(2,iSO) = iCntrc
                        basis_ids(3,iSO) = llab
                        basis_ids(4,iSO) = mlab
                        fermion_type(iSO)=0
                        If (fMass(iCnttp).ne.1.0D0) fermion_type(iSO)=1
                        if (.Not.Primitive_Pass) then
                           kIrrep=kIrrep+1
                           icent(kIrrep)=mdc
                           lnang(kIrrep)=lval(lculf)
                           lmag(kIrrep)=mval(lculf)
                           lant(kIrrep)=nIrrep/nStab(mdc)
                        Endif
 205                 Continue
*
 204              Continue
 2033             kComp = kComp + (iAng+1)*(iAng+2)/2
                  kculf=kculf+ 2*iAng+1
 203           Continue
               mc = mc + nIrrep/nStab(mdc)
 202        Continue
*
 201     Continue
 2011    Continue
culf
         nrSym=nIrrep
         nrBas(iIrrep+1)=nBas(iIrrep)
*        write(6,*) ' nBas(iIrrep)', iIrrep, nBas(iIrrep)
         nDim = nDim + nBas(iIrrep)
         n2Tot = n2Tot + nBas(iIrrep)**2
         n2Max = Max(n2Max,nBas(iIrrep)**2)
         m2Tot = m2Tot + nPrm(iIrrep)**2
 200  Continue
*     If (lSkip) nDim = iBas
      If (iBas.ne.iSO .and.
     &    iBas_Aux.ne.iSO_Aux-iSO .and.
     &    .Not.lSkip) Then
         Write (6,*) 'iBas=',iBas
         Write (6,*) 'iBas_Aux=',iBas_Aux
         Write (6,*) 'iSO=',iSO
         Write (6,*) 'iSO_Aux=',iSO_Aux-iSO
         Write (6,*) 'iSO_Tot=',iSO_Tot
         Call ErrTra
         Call Abend
      End If
C redefine iOffSO array in case of Fragment AIEMP
      If (lFAIEMP) Then
        Do iIrrep = 0, nIrrep-1
          iOffSO(iIrrep) = jOffSO(iIrrep)
        End Do
      End IF
#ifdef _DEBUG_
      Call RecPrt('Symmetrization Matrix','(20F5.2)',SM,iBas,iBas)
#endif
      Call Put_dArray('SM',SM,iBas**2)
*
CSVC: basis IDs of both symmetric and non-symmetric case
      if (.not.Primitive_Pass) then
        Call Put_iArray('Fermion IDs',fermion_type,iSO)
        Call Put_iArray('Basis IDs',basis_ids,4*iSO)
        Call Put_iArray('Desym Basis IDs',desym_basis_ids,4*iBas)
      end if
*
      Call mma_deallocate(IndC)
      Call mma_deallocate(SM)
      Call mma_deallocate(Index)
      Call mma_deallocate(Index2)
      Go To 198
*                                                                      *
************************************************************************
*                                                                      *
*---- No symmetry case.
*
199   Continue
      If (Show.and.iPrint.ge.6) Then
         Write (6,*)
         Write (6,'(19x,a)')
     &            ' **************************************************'
         Write (6,'(19x,a)')
     &            ' ********** Petite list Basis Functions ***********'
         Write (6,'(19x,a)')
     &            ' **************************************************'
         Write (6,*)
      End If
*
      Do i = 1, MxUnq
         IrrCmp(i) = 1
      End Do
      kIrrep=0
      Do 300 iIrrep = 0, nIrrep-1
         iOffSO(iIrrep) = iSO_Tot
         iAO = 0
         jSO = 0
         nBas(iIrrep) = 0
         nBas_Aux(iIrrep) = 0
         nBas_Frag(iIrrep) = 0
         nPrm(iIrrep) = 0
         Type(iIrrep)=.True.
*
*        Loop over distinct shell types
*
         mc  = 1
         iShell = 0
         iCnttp = 0
         Do 301 jCnttp = 1, nCnttp
*
*           Make sure that we process the dummy shell last
*
            If (jCnttp.eq.iCnttp_Dummy .and. jCnttp.ne.nCnttp) Then
               iCnttp = iCnttp + 2
            Else If (jCnttp.eq.nCnttp .and. iCnttp.eq.jCnttp) Then
               iCnttp = iCnttp_Dummy
            Else
               iCnttp = iCnttp + 1
            End If
*
            output = show .and. iPrint.ge.6
            If (AuxCnttp(iCnttp).or.FragCnttp(iCnttp))
     &        output = output.and.iPrint.ge.10
     &                       .and.iCnttp.ne.iCnttp_Dummy
            kECP = ECP(iCnttp)
            lMax=nVal_Shells(iCnttp)-1
            Call OrbType(iAtmNr(iCnttp),List_AE,31)
            If (kECP) Then
               Call ECP_Shells(iAtmNr(iCnttp),list)
               Call ICopy(lmax+1,[0],0,nCore_Sh,1)
            Else
               Call ICopy(1+iTabMx,List_AE,1,List,1)
               If (Charge(iCnttp).ne.Zero) Then
                  Call Freeze_Default(iAtmNr(iCnttp),nCore_Sh,lMax)
               Else
                  Call Freeze_Default(0             ,nCore_Sh,lMax)
               End If
            End If
*
*           Loop over distinct centers
*
            Do 302 iCnt = 1, nCntr(iCnttp)
               mdc = iCnt + mdciCnttp(iCnttp)
*
*              Loop over shells associated with this center
*              Start with s type shells
*
               kComp = 0
               kculf = 0
               iSh = ipVal(iCnttp) - 1
               If (nVal_Shells(iCnttp).lt.1) Then
                  ipxyz=(iCnt-1)*3+ipCntr(iCnttp)
                  XCoor=Dinf(ipxyz  )
                  YCoor=Dinf(ipxyz+1)
                  ZCoor=Dinf(ipxyz+2)
                  LPC(1,mdc)=XCoor
                  LPC(2,mdc)=YCoor
                  LPC(3,mdc)=ZCoor
                  LPQ(mdc)=Charge(iCnttp)
                  LPMM(mdc)=IsMM(iCnttp)
                  LPA(mdc)=iAtmnr(iCnttp)
                  LP_Names(mdc)=LblCnt(mdc)(1:LENIN)//'    '
               End If
               Do 303 iAng = 0, nVal_Shells(iCnttp)-1
                  nCore=nCore_Sh(iAng)
                  iSh = iSh + 1
                  iShell = iShell + 1
                  If (nExp(iSh).eq.0) Go To 3033
                  If (nBasis(iSh).eq.0) Go To 3033
                  jComp = (iAng+1)*(iAng+2)/2
                  If(Prjct(iSh)) jComp = 2*iAng + 1
                  Do 304 iComp = 1, jComp
                     iAO = iAO + 1
                     If (iAO.gt.MxAO) Then
                        Call ErrTra
                        Write (6,*) ' Increase MxAO'
                        Call Abend
                     End If
                     lComp = kComp + iComp
                     lculf = kculf + iComp
*
*                    Skip if symmetry operator is not in the coset of
*                    this center.
*
                     Do 308 imc = 0, (nIrrep/nStab(mdc))-1
                        If (iCoSet(imc,0,mdc).eq.iOper(iIrrep))
     &                     Go To 307
 308                 Continue
                     Go To 304
 307                 Continue
                     If (IndS(iShell)+iComp.gt.MxUnq) Then
                        Call ErrTra
                        Write (6,*) ' Increase MxUnq'
                        Call Abend
                     End If
                     If (output.and.Type(iIrrep)) Then
                        Write (6,*)
                        Write (6,'(10X,2A)')
     &                      ' Basis functions generated by ',
     &                       ChOper(iIrrep)
                        Write (6,*)
                        Write (6,'(A)')
     &                  ' Basis Label        Type   Center'
                        Type(iIrrep)=.False.
                     End If
*
                     If (MaxBas(iAng).gt.0) iAOtSO(iAO,iIrrep) = jSO + 1
                     nPrm(iIrrep) = nPrm(iIrrep) + nExp(iSh)
                     m2Max = Max(m2Max,nExp(iSh)**2)
                     If(.not.FragShell(iSh) .and.
     &                  .not.AuxCnttp(iCnttp))
     &                 nFCore(0)=nFCore(0)+nCore
*
*                    Loop over contracted basis functions
*
                     Do 305 iCntrc = 1, nBasis(iSh)
                        iSO_Tot = iSO_Tot + 1
                        If (AuxShell(iSh)) Then
                           iSO_Aux = iSO_Aux + 1
                           iSO_=iSO_Aux
                           nBas_Aux(iIrrep) = nBas_Aux(iIrrep) + 1
                        Else If (FragShell(iSh)) Then
                           iSO_Frag = iSO_Frag + 1
                           iSO_=iSO_Frag
                           nBas_Frag(iIrrep) = nBas_Frag(iIrrep) + 1
                        Else
                           iSO = iSO + 1
                           iSO_=iSO
                           nBas(iIrrep) = nBas(iIrrep) + 1
                        End If
                        If (iSO_.gt.nMamn) Then
                           Call qTrace
                           Write (6,*) ' iSO_.gt.nMamn'
                           Write (6,*) 'nMamn=',nMamn
                           Call Abend
                        End If
                        jSO = jSO + 1
*
                        ChTemp=LblCBs(lComp)
                        If (Transf(iSh)) ChTemp=LblSbs(lComp)
*
                        Call Name_to_lm(ChTemp,llab,mlab)
*
*                       Introduce a somewhat better labelling. Thnx LG!
*
                        If (IsBasisAE) Then
                           If (IsBasisANO) Then
                              Write (ChTemp(1:2),'(I2.2)') iAng+iCntrc
                           Else
                              If (nExp(iSh).eq.nBasis(iSh)) Then
                                 Write (ChTemp(1:1),'(A1)') '*'
                                 If (llab.ge.0)
     &                              Write(ChTemp(2:2),'(A1)') '0'
                              Else If (iCntrc.le.list(iAng)) Then
                                 Write (ChTemp(1:2),'(I2.2)')
     &                                 iAng+iCntrc
                              Else
                                 Write (ChTemp(1:1),'(A1)') '*'
                                 If (llab.ge.0)
     &                              Write(ChTemp(2:2),'(A1)') '0'
                              End If
                           End If
                        Else If (.Not.IsBasisUNK) Then
                           If (nExp(iSh).eq.nBasis(iSh)) Then
                              Write (ChTemp(1:1),'(A1)') '*'
                              If (llab.ge.0)
     &                           Write(ChTemp(2:2),'(A1)') '0'
                           Else If (iCntrc.le.list(iAng)) Then
                                 Write (ChTemp(1:2),'(I2.2)')
     &                                 iAng+iCntrc+
     &                                 (List_AE(iAng)-List(iAng))
                           Else
                              Write (ChTemp(1:1),'(A1)') '*'
                              If (llab.ge.0)
     &                           Write(ChTemp(2:2),'(A1)') '0'
                           End If
                        End If
                        ChTmp=Clean_BName(ChTemp,0)
*
                        if(output) Write (6,'(I5,2X,A8,5X,A8,I3)')
     &                        iSO_,LblCnt(mdc),ChTmp,mc+imc
*
                        iSOInf(1,iSO_)=iCnttp
                        iSOInf(2,iSO_)=iCnt
                        iSOInf(3,iSO_)=iAng
*
                        If (AuxShell(iSh).or.FragShell(iSh)) Go To 305
                        Write (isymunit,'(13(I4,4X))')
     &                     iSO,mdc,LVAL(lculf),MVAL(lculf),
     &                     nIrrep/nStab(mdc),
     &                     (iPrmt(NrOpr(iCoSet(iCo,0,mdc),
     &                     iOper,nIrrep),iChbs)*
     &                     iChTbl(iIrrep,NrOpr(iCoSet(iCo,0,mdc),
     &                     iOper,nIrrep)),
     &                     iCo=0,nIrrep/nStab(mdc)-1 )
*                                                                      *
************************************************************************
*                                                                      *
*---------------------------- Stuff (not just) for LoProp
*
                        iCI(iSO)=mdc
                        jCI(iSO)=mdc
                        If (iCntrc.le.list(iAng)) Then
                           iOT(iSO)=Occ
                        Else
                           iOT(iSO)=Vir
                        End If
                        ipxyz=(iCnt-1)*3+ipCntr(iCnttp)
                        XCoor=Dinf(ipxyz  )
                        YCoor=Dinf(ipxyz+1)
                        ZCoor=Dinf(ipxyz+2)
                        LPC(1,mdc)=XCoor
                        LPC(2,mdc)=YCoor
                        LPC(3,mdc)=ZCoor
                        LPQ(mdc)=Charge(iCnttp)
                        LPMM(mdc)=IsMM(iCnttp)
                        LPA(mdc)=iAtmnr(iCnttp)
                        LP_Names(mdc)=LblCnt(mdc)(1:LENIN)//'    '
*                                                                      *
************************************************************************
*                                                                      *
                        Mamn(iSO)=LblCnt(mdc)(1:LENIN)//ChTemp(1:8)
                        basis_ids(1,iSO) = mdc
                        basis_ids(2,iSO) = iCntrc
                        basis_ids(3,iSO) = llab
                        basis_ids(4,iSO) = mlab
                        fermion_type(iSO)=0
                        If (fMass(iCnttp).ne.1.0D0) fermion_type(iSO)=1
                        If (.Not.Primitive_Pass) Then
                           kIrrep=kIrrep+1
                           icent(kIrrep)=mdc
                           lnang(kIrrep)=lval(lculf)
                           lmag(kIrrep)=mval(lculf)
                           lant(kIrrep)=nIrrep/nStab(mdc)
                        Endif
 305                 Continue
*
 304              Continue
 3033             kComp = kComp + (iAng+1)*(iAng+2)/2
                  kculf=kculf+ 2*iAng+1
 303           Continue
               mc = mc + nIrrep/nStab(mdc)
 302        Continue
*
 301     Continue
         nrSym=nIrrep
         nrBas(iIrrep+1)=nBas(iIrrep)
         nDim = nDim + nBas(iIrrep)
         n2Tot = n2Tot + nBas(iIrrep)**2
         n2Max = Max(n2Max,nBas(iIrrep)**2)
         m2Tot = m2Tot + nPrm(iIrrep)**2
 300  Continue
*
CSVC: basis IDs of non-symmetric case
      if (.not.Primitive_Pass) then
        Call Put_iArray('Fermion IDs',fermion_type,iBas)
        Call Put_iArray('Basis IDs',basis_ids,4*iBas)
      end if
*
      Do 310 jAO = 1, iAO
         Do 311 jIrrep = 0, nIrrep-1
            iAOtSO(jAO,jIrrep) = iAOtSO(jAO,jIrrep) + iOffSO(jIrrep)
 311     Continue
 310  Continue
*
*     Fix index list such that redundant operators will have the
*     same AO index.
*
      mc  = 1
      iShell = 0
      iAO = 0
      Do 401 iCnttp = 1, nCnttp
         kECP = ECP(iCnttp)
*
*        Loop over distinct centers
*
         Do 402 iCnt = 1, nCntr(iCnttp)
            mdc = iCnt + mdciCnttp(iCnttp)
            iChxyz=iChCnt(mdc)
*
*           Loop over shells associated with this center
*           Start with s type shells
*
            kComp = 0
            iSh = ipVal(iCnttp) - 1
            Do 403 iAng = 0, nVal_Shells(iCnttp)-1
               iSh = iSh + 1
               iShell = iShell + 1
               If (nExp(iSh).eq.0) Go To 4033
               If (nBasis(iSh).eq.0) Go To 4033
               jComp = (iAng+1)*(iAng+2)/2
               If(Prjct(iSh)) jComp = 2*iAng + 1
               Do 404 iComp = 1, jComp
                  iAO = iAO + 1
                  If (iAO.gt.MxAO) Then
                     Call ErrTra
                     Write (6,*) ' Increase MxAO'
                     Call Abend
                  End If
                  lComp = kComp + iComp
                  If (MaxBas(iAng).gt.0) Then
*
                     Do 408 imc = 0, (nIrrep/nStab(mdc))-1
                        itest1 = iAnd(iCoSet(imc,0,mdc),iChxyz)
                        Nr = NrOpr(iCoSet(imc,0,mdc),iOper,nIrrep)
                        Do 409 jIrrep = 0, nIrrep-1
                           itest2 = iAnd(iOper(jIrrep),iChxyz)
                           If (itest1.eq.itest2)
     &                        iAOtSO(iAO,jIrrep) = iAOtSO(iAO,Nr)
 409                    Continue
 408                 Continue
                  End If
 404           Continue
 4033          kComp = kComp + (iAng+1)*(iAng+2)/2
 403        Continue
            mc = mc + nIrrep/nStab(mdc)
 402     Continue
*
 401  Continue
*
 198  Continue
*                                                                      *
************************************************************************
*                                                                      *
      If (Show) Then
        If (iPrint.ge.6) Then ! std print has been moved to seward.f
*
*        Print out basis set information
*
         Fmt='(6X,A,T30,8I4)'
         Write(6,*)
         Write(6,'(6X,A)')'Basis set specifications :'
         Write(6,'(6X,A,T30,8(1X,A))')
     &         'Symmetry species',         (lIrrep(i),i=0,nIrrep-1)
         Write(6,Fmt)'Basis functions',          (nBas(i),i=0,nIrrep-1)
*
        End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Write info (not just) for LoProp
*
      If (.Not.Primitive_Pass) Then
         Call Put_cArray('LP_L',LP_Names(1),(LENIN4)*mCentr)
         Call Put_iArray('LP_A',LPA,mCentr)
         Call Put_dArray('LP_Q',LPQ,mCentr)
         Call Put_dArray('LP_Coor',LPC,3*mCentr)
         Call Put_iScalar('LP_nCenter',mCentr)
         Call Put_iArray('IsMM Atoms',LPMM,mCentr)
         Call Put_iArray('Center Index',iCI,iBas)
         Call Put_iArray('Orbital Type',iOT,iBas)
         Call Put_iArray('Non valence orbitals',nFCore,nIrrep)
      Else
         Call Put_iArray('Ctr Index Prim',jCI,iBas)
*
      End If
      call mma_deallocate(fermion_type)
      call mma_deallocate(desym_basis_ids)
      call mma_deallocate(basis_ids)
      Call mma_deallocate(LPA)
      Call mma_deallocate(LPMM)
      Call mma_deallocate(LPQ)
      Call mma_deallocate(LPC)
      Call mma_deallocate(iCI)
      Call mma_deallocate(jCI)
      Call mma_deallocate(iOT)
*
      If (Show.and.iPrint.ge.6) Then
         Call CollapseOutput(0,'   SO/AO info:')
         Write (6,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Mx_Unq=IEMax
      Mx_AO=iAO
*
#ifdef _DEBUG_
      Write (6,*) ' *** iAOtSO ***'
      Do 555 jAO = 1, iAO
         Write (6,*) (iAOtSO(jAO,jIrrep),jIrrep=0,nIrrep-1)
 555  Continue
      Write (6,*) ' *** IrrCmp ***'
      Do 556 iE = 1, iEMax
         Write (6,*) IrrCmp(iE)
 556  Continue
#endif
*
      write(isymunit,'(A3)') 'END'
      close(isymunit)
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit ('SOCtl')
      Return
      If (.False.) Call Unused_Integer(Info)
      End
