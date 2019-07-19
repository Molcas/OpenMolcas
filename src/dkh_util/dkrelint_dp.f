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
* Copyright (C) 2011, Daoling Peng                                     *
************************************************************************
      Subroutine DKRelint_DP
*
*     modified by D. Peng, ETH Zurich, October 2011
*
*     Interface/Driver routine for scalar relativistic
*           arbitrary-order DKH method,
*           exact decoupling X2C method &
*           exact decoupling BSS method.
*
      Implicit real*8(a-h,o-z)
#include "warnings.fh"
#include "itmax.fh"
#include "info.fh"
#include "rinfo.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "wldata.fh"
#include "oneswi.fh"
#include "RelLight.fh"
#include "relae.fh"
      integer ipaddr(3)
      Character*8 Label, pXpLbl
#ifdef MOLPRO
      character*(64) filename
#endif
      Integer nBas_prim(8), nBas_cont(8)
      Logical Debug
      Data Debug/.False./
      character*(3) paramtype
      integer relmethod,dkhorder,xorder,dkhparam
      logical delflag,DoFullLT
      integer stdout
      Dimension idum(1)

      If(IRFLAG1.eq.1)then
        Call DKRelint
        Return
      End IF
*                                                                      *
************************************************************************
*                                                                      *
      iRout=77
      iPrint=nPrint(iRout)
      Call QEnter('DKRelInt')
*
      If(Debug)Then
        idbg=6
      Else
        idbg=-1
      Endif
*
      stdout=6
*                                                                      *
************************************************************************
*                                                                      *
*     Save basis set info from contracted run
*
      if(iprint.ge.10) write(stdout,*) ' In DKRelInt', ncnttp
      kCof=0
      kAng=0
      kExp=0
      kC=0
*
*     Normalize coefficients
*
      do iCnttp=1,nCnttp
*
*       The none valence type shells comes at the end. When this block
*       is encountered stop the procedure.
*
        If(AuxCnttp(iCnttp) .or.
     &      FragCnttp(iCnttp) .or.
     &      nFragType(iCnttp).gt.0 ) Go To 999

*
        Do icnt = 1, nCntr(iCnttp)
        kC=kC+1
           Do iAngr=0,nAngr(kC)
               rI=DBLE(iAngr)+One+Half
              kAng=kAng+1
              Do iBas=1,nBasisr(kAng)
                 Sum=Zero
                 kExpi=kExp
                 kCofi=kCof
                 Do iExp=1,nPrimr(kAng)
                    kExpi=kExpi+1
                    kCofi=kCofi+1
                    rExpi=rExp(kExpi)
c                   write(stdout,'(a11,f20.8)') ' Exponents',rExpi
                    rCofi=rCof(kCofi)
                    kExpj=kExp
                    kCofj=kCof
                    Do jExp=1,nPrimr(kAng)
                       kExpj=kExpj+1
                       kCofj=kCofj+1
                       rExpj=rExp(kExpj)
                       rCofj=rCof(kCofj)
                       Sum=Sum+rCofi*rCofj*
     &                 (Two*sqrt(rExpi*rExpj)/(rExpi+rExpj))**rI
                    Enddo
                 Enddo
                 rNorm=One/sqrt(Sum)
                 if(iprint.ge.10) write(stdout,*) ' rNorm', kAng,rNorm
                   Do iExp=1,nPrimr(kAng)
                      rCof(kCof+iExp)=rCof(kCof+iExp)*rNorm
                     if(iprint.ge.10) then
                         write(stdout,'(a24,f20.6)')
     &                   ' normalized coefficients',
     &                   rCof(kCof+iExp)
                     endif
                   Enddo
                 kCof=kCof+nPrimr(kAng)
              Enddo
              kExp=kExp+nPrimr(kAng)
           Enddo
        Enddo
      Enddo
 999  Continue
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.10) Then
      i=0
      Do L=1,nrSym
         write(stdout,*) ' Irreducible representation', L
         Do ibas=1,nrBas(L)
            i=i+1
            Write (stdout,'(20i4)') i, icent(i),lnang(i),lmag(i)
         End Do
      End Do
      End If
*
      Call iCopy(8,nBas,1,nBas_Cont,1)
      nSym=nIrrep
#ifdef MOLPRO
      call icopy(8,nrbas_prim,1,nbas,1)
#else
*                                                                      *
************************************************************************
*                                                                      *
*     Close ONEINT
*
      iOpt=0
      Call ClsOne(iRC,iOpt)
*                                                                      *
************************************************************************
*                                                                      *
*     Open ONEREL
*
      iOpt = 0
      iRC = -1
      Lu_One=2
      Call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
      If (iRC.ne.0) Go To 9999
*
      Call OneBas('PRIM')
*                                                                      *
************************************************************************
*                                                                      *
      Call Get_iArray('nBas_Prim',nBas,nSym)
#endif
      Call iCopy(8,nBas,1,nBas_prim,1)
      If(iPrint.ge.10) then
         write(stdout,'(a,8i5)') ' Symmetries          ', nSym
         write(stdout,'(a,8i5)') ' Primitive basis fcns',
     &                          (nBas(i),i=0,nSym-1)
      Endif
*
*     Allocate memory for relativistic part
*
      VELIT=CLightAU
      iSizep=0
      iSizes=0
      iSizec=0
      iibas =0
      Do L=0,nSym-1
         iSizep=iSizep + nBas(L)*(nBas(L)+1)/2
         iSizes=iSizes + nBas(L)*nBas(L)
         iSizec=iSizec+nrBas(L+1)*(nrBas(L+1)+1)/2
         iibas =iibas + nbas(L)
      End Do
      If(iPrint.ge.10) write(stdout,*) ' iSizep', iSizep
*
      CALL GetMem('Kin     ','ALLO','REAL',iK,iSizep+4)
      CALL GetMem('SS      ','ALLO','REAL',iSS,iSizep+4)
      CALL GetMem('V       ','ALLO','REAL',iV,iSizep+4)
      CALL GetMem('pVp     ','ALLO','REAL',ipVp,iSizep+4)
*
      If (iprint.ge.20) write(stdout,*)
     &   '  indices', iss,ik,iv,ipvp
#ifdef MOLPRO
      call lesw(work(iss),iSizep,1,1101,0)
      call lesw(work(ik),iSizep,1,1401,0)
      Call lesw(Work(iv),iSizep,1,1411,0)
      call lesw(work(ipvp),iSizep,1,1412,0)
#else
      Label='Mltpl  0'
      iComp=1
      iOpt=0
      iRC = -1
      Call RdOne(iRC,iOpt,Label,1,Work(iSS),lOper)
      If (iRC.ne.0) then
         Write (stdout,*) 'DKRelInt: Error reading from ONEINT'
         Write (stdout,'(A,A)') 'Label=',Label
         Call Abend
      End If
      nComp=1
      ipaddr(1)=iSS
      If (iPrint.ge.20) Call PrMtrx(Label,[lOper],nComp,ipaddr,Work)
      Label='Attract '
      iRC = -1
      Call RdOne(iRC,iOpt,Label,1,Work(iV),lOper)
      If (iRC.ne.0) then
         Write (stdout,*) 'DKRelInt: Error reading from ONEINT'
         Write (stdout,'(A,A)') 'Label=',Label
         Call Abend
      End If
      Label='Kinetic '
      iRC = -1
      Call RdOne(iRC,iOpt,Label,1,Work(iK),lOper)
      If (iRC.ne.0) then
         Write (stdout,*) 'DKRelInt: Error reading from ONEINT'
         Write (stdout,'(A,A)') 'Label=',Label
         Call Abend
      End If
      Label='pVp     '
      iRC = -1
      Call RdOne(iRC,iOpt,Label,1,Work(ipVp),lOper)
      If (iRC.ne.0) then
         Write (stdout,*) 'DKRelInt: Error reading from ONEINT'
         Write (stdout,'(A,A)') 'Label=',Label
         Call Abend
      End If
*
      iOpt=0
      Call ClsOne(iRC,iOpt)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      If (IRELAE.ge.100) Then
         If (IRELAE.ge.1000) Then
*
*        IRELAE codes the method (DKH=1), the order (01 to 99), and the
*        parametrization (1 to 5), i.e. 3rd order DKH with sqrt
*        parametrization: 1-03-3, i.e. IRELAE=1033
*
            relmethod = 1
*
*           relmethod = 1  arbitrary order DKH
*                     = 2  exact decoupling X2C
*                     = 3  exact decoupling BSS
*
            xOrder = iRELAE/10000
            iTemp=iRELAE-1000 - xOrder*10000
            dkhorder = iTemp/10
            iTemp = iTemp - dkhorder*10
            dkhparam = iTemp
            If (xorder.gt.dkhorder) then
              Write(stdout,"(a,i3,a,i3)") " xorder was reduced from ",
     &                                     xorder," to ",dkhorder
              xorder=dkhorder
            End If
            If (iTemp.eq.1) Then
               paramtype='OPT'
            Else If (iTemp.eq.2) Then
               paramtype='EXP'
            Else If (iTemp.eq.3) Then
               paramtype='SQR'
            Else If (iTemp.eq.4) Then
               paramtype='MCW'
            Else If (iTemp.eq.5) Then
               paramtype='CAY'
            Else
               Write(stdout,*) 'dkrelint: Illegal parametrization!'
               Call Abend
            End If
c            write(stdout,"(A,A)") " Computing the",
c     &                             " arbitrary order DKH Hamiltonian"
         Else If (IRELAE.eq.101) Then
            relmethod = 2
            xorder = 1
c            write(stdout,"(A,A)") " Computing the",
c     &                            " exact decoupling X2C Hamiltonian"
         Else If (IRELAE.eq.102) Then
            relmethod = 3
            xorder = 1
c            write(stdout,"(A,A)") " Computing the",
c     &                            " exact decoupling BSS Hamiltonian"
         Else
            Write(stdout,*) 'dkrelint: Unknown method !'
            Call Abend
         End If
*                                                                      *
************************************************************************
*                                                                      *
         Call Allocate_Work(iK_Save,iSizep+4)
         Call Allocate_Work(iK_Done,iSizep+4)
         call dcopy_(iSizep+4,Work(iK),1,Work(iK_Save),1)
*
         Call Allocate_Work(iU_L,iSizes+4)
         Call Allocate_Work(iU_S,iSizes+4)
*
*        Read block information if do local transformation
*
         If (LDKroll) Then
           Call GetMem('Index  ','ALLO','INTE',indx,iibas+4)
           Call xdr_indx(iibas,iWork(indx))
CDP           write(6,*) "radild : ",radild
           DoFullLT=.true.
           if(radiLD.eq.0.d0) DoFullLT=.false.
           if(DoFullLT)then
             if (relmethod.eq.1.and.xorder.eq.0) then
               xorder = dkhorder
             end if
             Write(6,"(A)") "   DLU Local Transformation"
           else
             Write(6,"(A)") "   DLH Local Transformation"
           end if
         End If
*
*        Do the Hamiltonian seperately
*
         k=0
         ks=0
         kz=0
*
         Do L = 0, nSym-1
            If (L.eq.nSym-1) delflag=.TRUE.
            n=nBas(L)
            iSize=n*(n+1)/2
            If (iSize.eq.0) Go To 911
*                                                                      *
************************************************************************
*                                                                      *
            If (LDKroll) Then
              Call GetMem('InfoLoc','ALLO','INTE',iLoc,n+4)
              Call GetMem('MapLoc ','ALLO','INTE',iMap,n+4)
              call xdr_info_local(n,iWork(indx+kz),nbl,iWork(iLoc),
     &                            iWork(iMap) )
CDP              write(6,"(a,i1,i5,a,99i4)") '   Sym: ',L+1,n,
CDP     &                    '  = Local ',(iWork(iLoc+i),i=0,nbl-1)
              Call XDR_Local_Ham(n,isize,n*n,relmethod,dkhparam,
     &                     dkhorder,xorder,Work(iSS+k),Work(iK+k),
     &                     Work(iV+k),Work(ipVp+k),Work(iU_L+ks),
     &                     Work(iU_S+ks),iWork(indx+kz),nbl,
     &                     iWork(iLoc),iWork(iMap),DoFullLT,clightau)
              Call GetMem('InfoLoc','FREE','INTE',iLoc,n+4)
              Call GetMem('MapLoc ','FREE','INTE',iMap,n+4)
            Else
              Call XDR_Ham(n,isize,n*n,relmethod,dkhparam,dkhorder,
     &                     xorder,Work(iSS+k),Work(iK+k),Work(iV+k),
     &                     Work(ipVp+k),Work(iU_L+ks),Work(iU_S+ks),
     &                     clightau)
            End If
*                                                                      *
************************************************************************
*                                                                      *
            ks=ks+n*n
            kz=kz+n
 911        k=k+isize
         End Do
*                                                                      *
************************************************************************
*                                                                      *
         call dcopy_(iSizep+4,Work(iK),1,Work(iK_Done),1)
*                                                                      *
************************************************************************
*                                                                      *
         If (xOrder.le.0) Go To 912
*                                                                      *
************************************************************************
*                                                                      *
*        Pick up the number of property integrals to process.
*
#ifdef MOLPRO
#else
         numb_props=nProp_Int(.False.,iWork(ip_iDummy),0)
         Call Allocate_iWork(ipInd,4*numb_props)
         numb_props=nProp_Int(.True.,iWork(ipInd),numb_props)
         Do iProps=1,numb_props
            ipOp   =ipInd + (iProps-1)*4
            ip_MEF =ipInd + (iProps-1)*4 + 1
            ipiComp=ipInd + (iProps-1)*4 + 2
            ipjCent=ipInd + (iProps-1)*4 + 3
*
            call dcopy_(iSizep+4,Work(iK_Save),1,Work(iK),1)
*                                                                      *
************************************************************************
*                                                                      *
*           Read the property integrals
*
*           Open ONEREL
*
            iOpt = 0
            iRC = -1
            Lu_One=2
            Call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
            If (iRC.ne.0) Go To 9999
*
            Call OneBas('PRIM')
*
            iMEF = iWork(ip_MEF)
            If (iWork(ipOp).eq.1) Then
               Write (Label,'(a,i2)') 'MLTPL ',iMEF
            Else If (iWork(ipOp).eq.2) Then
               jCent = iWork(ipjCent)
               Write (Label,'(a,i1,i5)') 'EF',iMEF,jCent
            Else If (iWork(ipOp).eq.3) Then
               jCent = iWork(ipjCent)
               Write (Label,'(a,i5)') 'Cnt',jCent
            Else
               Write (6,*) 'DKRelInt: illegal property!'
               Call Abend
            End If
            iComp = iWork(ipiComp)
C           Write (6,*)
C           Write (6,*) 'Label=',Label
C           Write (6,*) 'iComp=',iComp
C           Write (6,*)
*
            iOpt=1
            iRC = -1
            lOper=-1
            Call iRdOne(iRC,iOpt,Label,iComp,idum,lOper)
            If (iRC.eq.0) nInt=idum(1)
C           Write (6,*) 'lOper=',lOper
            CALL GetMem('X       ','ALLO','REAL',iX,nInt+4)
            iRC = -1
            iOpt=0
            Call RdOne(iRC,iOpt,Label,iComp,Work(iX),lOper)
            If (iRC.ne.0) then
               Write (stdout,*) 'DKRelInt: Error reading from ONEREL'
               Write (stdout,'(A,A)') 'Label=',Label
               Write (stdout,'(A,A)') 'iRC=',iRC
               Call Abend
            End If
            Call CmpInt(Work(iX),nInt,nBas_Prim,nSym,lOper)
            If (nInt.eq.0) Then
               iOpt=0
               Call ClsOne(iRC,iOpt)
               Go To 666
            End If
*
            If (iWork(ipOp).eq.1) Then
               Write (pXpLbl,'(A,I2)') 'pMp   ', iMEF
            Else If (iWork(ipOp).eq.2) Then
               Write (pXpLbl,'(A,I1,I5)') 'PP',iMEF,jCent
            Else If (iWork(ipOp).eq.3) Then
               Write (pXpLbl,'(A,I2)') 'pCp   ',jCent
            End If
            iOpt=1
            iRC = -1
            Call iRdOne(iRC,iOpt,pXpLbl,iComp,idum,lOper)
            If (iRC.eq.0) nInt=idum(1)
            CALL GetMem('pXp     ','ALLO','REAL',ipXp,nInt+4)
            iOpt=0
            iRC = -1
            Call RdOne(iRC,iOpt,pXpLbl,iComp,Work(ipXp),lOper)
            If (iRC.ne.0) then
               Write (stdout,*) 'DKRelInt: Error reading from ONEREL'
               Write (stdout,'(A,A)') 'pXpLbl=',pXpLbl
               Write (stdout,'(A,A)') 'iRC=',iRC
               Call Abend
            End If
            Call CmpInt(Work(ipXp),nInt,nBas_Prim,nSym,lOper)
*
            iOpt=0
            Call ClsOne(iRC,iOpt)
*
            Call GetMem('Core','Max','Real',iDum(1),Mem_Available)
C           Write (6,*) 'Mem_Available=',Mem_Available
            delflag=.FALSE.
            k=0
            ks=0
            kz=0
            Do L = 0, nSym-1
               If (L.eq.nSym-1 .and.
     &             iProps.eq.numb_props) delflag=.TRUE.
               n=nBas(L)
               iSize=n*(n+1)/2
               If (iSize.eq.0) Go To 91
*
*              Skip if the propetry operator does not have a total
*              symmetric component!
*
               If(iAnd(1,lOper).eq.0) Go To 91
*                                                                      *
************************************************************************
*                                                                      *
               Call XDR_Prop(n,isize,n*n,relmethod,dkhparam,dkhorder,
     &                       xorder,Work(iSS+k),Work(iK+k),Work(iV+k),
     &                       Work(ipVp+k),Work(iX+k),Work(ipXp+k),
     &                       Work(iU_L+ks),Work(iU_S+ks),clightau)
               ks=ks+n*n
               kz=kz+n
 91            k=k+isize
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*           Put the picture change corrected integral back to the
*           ONEINT file. Primitives in Work(iX).
*
*           First contract the result, store in Work(ip_Prop)
*
            Call Allocate_Work(ip_Prop,iSizec+4)
            Call repmat(idbg,Work(iX),Work(ip_Prop))
*
*                                                                      *
************************************************************************
*                                                                      *
*           Read the contracted property integrals from OneInt
*
            iOpt = 0
            iRC = -1
            Lu_One=2
            Call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
            If (iRC.ne.0) Go To 9999
*
            iOpt=1
            iRC = -1
            lOper=-1
            Call iRdOne(iRC,iOpt,Label,iComp,idum,lOper)
            If (iRC.eq.0) nInt=idum(1)
            CALL GetMem('Y       ','ALLO','REAL',iY,nInt+4)
            iRC = -1
            iOpt=0
            Call RdOne(iRC,iOpt,Label,iComp,Work(iY),lOper)
C           Write (6,*) 'Y1=',DDot_(nInt,Work(iY),1,
C    &                                  1.0D0,0)
            If (iRC.ne.0) then
               Write (stdout,*) 'DKRelInt: Error reading from ONEINT'
               Write (stdout,'(A,A)') 'Label=',Label
               Call Abend
            End If
*
*           Put the picture change corrected blocks in. Note that this
*           is just the diagonal symmetry blocks.
*
            Call Cp_Prop_Int(Work(iY),nInt,Work(ip_Prop),iSizec,
     &                       nrBas,nIrrep,lOper)
*
*           Now write it back to disc
*
            iOpt=0
            Call WrOne(iRC,iOpt,Label,iComp,Work(iY),lOper)
*
            iOpt=0
            Call ClsOne(iRC,iOpt)
*
            Call Free_Work(iY)
            Call Free_Work(ip_Prop)
*
            Call GetMem('pXp     ','FREE','REAL',ipXp,iSizep+4)
 666        Continue
            Call GetMem('X       ','FREE','REAL',iX,iSizep+4)
*                                                                      *
************************************************************************
*                                                                      *
         End Do ! iProps
*
         Call Free_iWork(ipInd)
#endif
 912     Continue
*                                                                      *
************************************************************************
*                                                                      *
         call dcopy_(iSizep+4,Work(iK_Done),1,Work(iK),1)
         Call Free_Work(iK_Done)
         Call Free_Work(iK_Save)
         Call Free_Work(iU_L)
         Call Free_Work(iU_S)
         If (LDKroll) Then
           Call GetMem('Index  ','FREE','INTE',indx,iibas+4)
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*        The Old code, no option for property integrals.
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*        Loop over the symmetry blocks
*
         epsilon=1.d-10
         delflag=.FALSE.
         k=0
         Do L = 0, nSym-1
            If (L.eq.nSym-1) delflag=.TRUE.
            n=nBas(L)
            iSize=n*(n+1)/2
            If (iSize.eq.0) goto 9
*           Allocate
*
            CALL GetMem('P       ','ALLO','REAL',iP,isize+4)
            CALL GetMem('G       ','ALLO','REAL',iG,isize+4)
            CALL GetMem('Ev2     ','ALLO','REAL',iEv2,n*n+4)
            CALL GetMem('Eig     ','ALLO','REAL',iEig,n*n+4)
            CALL GetMem('Sinv    ','ALLO','REAL',iSinv,n*n+4)
            CALL GetMem('Ew      ','ALLO','REAL',iEw,n+4)
            CALL GetMem('E       ','ALLO','REAL',iE,n+4)
            CALL GetMem('Aa      ','ALLO','REAL',iAa,n+4)
            CALL GetMem('Rr      ','ALLO','REAL',iRr,n+4)
            CALL GetMem('Tt      ','ALLO','REAL',iTt,n+4)
            CALL GetMem('Re1r    ','ALLO','REAL',iRe1r,n*n+4)
            CALL GetMem('Auxi    ','ALLO','REAL',iAuxi,n*n+4)
            CALL GetMem('Twrk4   ','ALLO','REAL',iTwrk4,n*200+4)
            Length    = N*N+4
            Length2   = iSize+4
            iDim      = N
            If (IRELAE .EQ. 0) Then
               Length    = 1
               Length2   = 1
               iDim      = 1
            End If
            CALL GetMem('Even1   ','ALLO','REAL',iEven1,Length)
            CALL GetMem('Pvpt    ','ALLO','REAL',iPvpt,Length2)
            CALL GetMem('Bu      ','ALLO','REAL',iBu,Length2)
*
*           call to package relsew
*
            CALL SCFCLI(idbg,epsilon,Work(iSS+k),Work(iK+k),Work(iV+k),
     &                  Work(ipVp+k),n,iSize,VELIT,Work(iBu),Work(iP),
     &                  Work(iG),Work(iEv2),Work(iEig),Work(iSinv),
     &                  Work(iEw),
     &                  Work(iE),Work(iAa),Work(iRr),Work(iTt),
     &                  Work(iPvpt),Work(iEven1),
     &                  Work(iRe1r),Work(iAuxi),
     &                  Work(iTwrk4),iDim)


            CALL GetMem('Bu      ','FREE','REAL',iBu,Length2)
            CALL GetMem('P       ','FREE','REAL',iP,isize+4)
            CALL GetMem('G       ','FREE','REAL',iG,isize+4)
            CALL GetMem('Ev2     ','FREE','REAL',iEv2,n*n+4)
            CALL GetMem('Eig     ','FREE','REAL',iEig,n*n+4)
            CALL GetMem('Sinv    ','FREE','REAL',iSinv,n*n+4)
            CALL GetMem('Ew      ','FREE','REAL',iEw,n+4)
            CALL GetMem('E       ','FREE','REAL',iE,n+4)
            CALL GetMem('Aa      ','FREE','REAL',iAa,n+4)
            CALL GetMem('Rr      ','FREE','REAL',iRr,n+4)
            CALL GetMem('Tt      ','FREE','REAL',iTt,n+4)
            CALL GetMem('Pvpt    ','FREE','REAL',iPvpt,Length2)
            CALL GetMem('Even1   ','FREE','REAL',iEven1,Length)
            CALL GetMem('Re1r    ','FREE','REAL',iRe1r,n*n+4)
            CALL GetMem('Auxi    ','FREE','REAL',iAuxi,n*n+4)
            CALL GetMem('Twrk4   ','FREE','REAL',iTwrk4,n*200+4)
 9          k=k+isize
         End Do
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate arrays in contracted basis
*
*
      Call GetMem('H       ','ALLO','REAL',iH,iSizec+4)
      Call FZero(Work(iH),iSizec+4)
      Call GetMem('H_o     ','ALLO','REAL',iH_nr,iSizec+4)
      Call FZero(Work(iH_nr),iSizec+4)
      Call GetMem('H_temp  ','ALLO','REAL',iH_temp,iSizec+4)
      Call FZero(Work(iH_temp),iSizec+4)
c
#ifdef MOLPRO
c... store relativistic H
      Call repmat(idbg,Work(ik),Work(iH_temp))
      call fperm(Work(iH_temp),Work(ih))
      Call writem(Work(iH),isizec+2,1,1200,0,'H0')
      Call writem(Work(iH),isizec+2,1,1210,0,'H01')
c... store V=H-T
      Call lesw(Work(iss),iSizec,1,1400,0)
      Call daxpy_(iSizec,-one,Work(iss),1,Work(iH),1)
      Call writem(Work(iH),isizec+2,1,1410,0,'POT')
c... reset contracted basis size
      Call iCopy(8,nBas_Cont,1,nBas,1)
      CALL GetMem('V       ','FREE','REAL',iV,iSizep+4)
      CALL GetMem('SS      ','FREE','REAL',iSS,iSizep+4)
      CALL GetMem('Kin     ','FREE','REAL',iK,iSizep+4)
#else
*
*     Note: in combination with ECPs V is only based on the effective
*           charges of the atoms. In the primitive basis, however, we
*           have temporarily introduced the actual atomic charges. We
*           have to fix this now. Hence the somewhat strange way in
*           which the DKH corrected Hamiltonian is computed.
*
*     Compute stripped non-relativistic H (iH_Temp)
*                                                                      *
************************************************************************
*                                                                      *
*     Open ONEREL
*
      iOpt = 0
      iRC = -1
      Lu_One=2
      Call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
      If (iRC.ne.0) Go To 9999
*
      Call OneBas('PRIM')
*
      Label='Kinetic '
      Call RdOne(iRC,iOpt,Label,1,Work(iSS),lOper)
      If (iRC.ne.0) then
         Write (stdout,*) 'DKRelInt: Error reading from ONEINT'
         Write (stdout,'(A,A)') 'Label=',Label
         Call Abend
      End If
      Label='Attract '
      Call RdOne(iRC,iOpt,Label,1,Work(iV),lOper)
      If (iRC.ne.0) then
         Write (stdout,*) 'DKRelInt: Error reading from ONEINT'
         Write (stdout,'(A,A)') 'Label=',Label
         Call Abend
      End If
*     Add Kinetic and Attraction term
      Call DaXpY_(iSizep+4,One,Work(iSS),1,Work(iV),1)
      If (iPrint.ge.20) Then
         Call iSwap(8,nBas,1,nBas_Prim,1)
         Call PrMtrx('Attract+Kinetic (prim)',[lOper],nComp,[iV],Work)
         Call PrMtrx('Kinetic (prim)',[lOper],nComp,[iSS],Work)
         Call iSwap(8,nBas,1,nBas_Prim,1)
      End If
      call dcopy_(4,[Zero],0,Work(iH_Temp+iSizec),1)
*     Contract and store in iH_temp
      Call repmat(idbg,Work(iV),Work(iH_temp))
*
      CALL GetMem('V       ','FREE','REAL',iV,iSizep+4)
      CALL GetMem('SS      ','FREE','REAL',iSS,iSizep+4)
*
*     Close ONEREL and re-open ONEINT
*
      iOpt = 0
      iRC = -1
      Call ClsOne(iRC,iOpt)
      If (iRC.ne.0) Go To 9999
      iOpt = 0
      iRC = -1
      Lu_One=2
      Call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
      If (iRC.ne.0) Go To 9999
*
      If (iPrint.ge.20) Then
         Call iSwap(8,nBas,1,nBas_Cont,1)
         Call PrMtrx('iH_temp (cont)',[lOper],nComp,[iH_temp],Work)
         Call iSwap(8,nBas,1,nBas_Cont,1)
      End If
*
*     Transform DKH Hamiltonian to contracted basis (iH)
*
      If (iPrint.ge.20) Then
         Call iSwap(8,nBas,1,nBas_Prim,1)
         Call PrMtrx('iK (prim)',[lOper],nComp,[iK],Work)
         Call iSwap(8,nBas,1,nBas_Prim,1)
      End If
      Call repmat(idbg,Work(iK),Work(iH))
      If (iPrint.ge.20) Then
         Call iSwap(8,nBas,1,nBas_Cont,1)
         Call PrMtrx('iH (cont)',[lOper],nComp,[iH],Work)
         Call iSwap(8,nBas,1,nBas_Cont,1)
      End If
*
      CALL GetMem('Kin     ','FREE','REAL',iK,iSizep+4)
*
      iOpt = 0
      iRC = -1
      Label='OneHam 0'
      Call RdOne(iRC,iOpt,Label,1,Work(iH_nr),lOper)
      If (iRC.ne.0) then
         Write (stdout,*) 'DKRelInt: Error reading from ONEINT'
         Write (stdout,'(A,A)') 'Label=',Label
         Call Abend
      End If
      iOpt = 0
      iRC = -1
*
*     final Hamiltonian computed as H(nrel) + ( Hrel(s) - Hnrel(s))
*     where (s) is stripped and with full charge
*
      Call DaXpY_(iSizec+4,-One,Work(iH_temp),1,Work(iH),1)
      Call DaXpY_(iSizec+4,One,Work(iH_nr),1,Work(iH),1)
*
      Call Get_iArray('nBas',nBas,nSym)
      If(iPrint.ge.10) then
         write(stdout,'(a11,10i5)') ' Symmetries', nSym
         write(stdout,'(a11,10i5)') ' Contracted',(nBas(i),i=0,nSym-1)
      Endif
      Label='OneHam 0'
      lOper=1
      nComp=1
      ipaddr(1)=iH
      If (iPrint.ge.20) Call PrMtrx(Label,[lOper],nComp,ipaddr,Work)
*
*     Replace 1-el Hamiltonian on ONEINT
*
      iRC = -1
      Call WrOne(iRC,iOpt,Label,1,Work(iH),lOper)
      Label='OneHam  '
      Call WrOne(iRC,iOpt,Label,1,Work(iH),lOper)
      If (iRC.ne.0) then
         Write (stdout,*) 'DKRelInt: Error reading from ONEINT'
         Write (stdout,'(A,A)') 'Label=',Label
         Call Abend
      End If
*
      If (IRELAE .eq. 23) Then   ! IORA
*
*     Replace overlap on ONEINT
*
         Call repmat(idbg,Work(ipVp),Work(iH))
         iRC = -1
         Label='Mltpl  0'
         Call WrOne(iRC,iOpt,Label,1,Work(iH),lOper)
         If (iRC.ne.0) then
            Write (stdout,*) 'DKInt: Error reading from ONEINT'
            Write (stdout,'(A,A)') 'Label=',Label
            Call Abend
         End If
      End If
#endif
*
      CALL GetMem('OneHam  ','FREE','REAL',iH,iSizec+4)
      CALL GetMem('H_o     ','FREE','REAL',iH_nr,iSizec+4)
      CALL GetMem('H_temp  ','FREE','REAL',iH_temp,iSizec+4)
      CALL GetMem('pVp     ','FREE','REAL',ipVp,iSizep+4)
*
      Call QExit('DKRelInt')
      Return
*
 9999 Continue
      Write (stdout,*) ' *** Error in subroutine DKRelInt ***'
      Write (stdout,*) '     Abend in subroutine OpnOne'
      Call Abend
      End
