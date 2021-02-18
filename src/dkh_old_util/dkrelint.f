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
* Copyright (C) 2005, Alexander Wolf                                   *
*               2005, Markus Reiher                                    *
************************************************************************
      Subroutine DKRelint
      use Basis_Info
      use Temporary_Parameters, only: force_out_of_core
      use DKH_Info
      use Symmetry_Info, only: nIrrep
c
c     modified by A. Wolf and M. Reiher, Uni Bonn, Feb. and March 2005
c       (extended for use of generalized arbitrary-order DKH)
c       NB: If the standard 2nd order DKH is wanted,
c           the original routines by Hess are called!
c
      Implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "warnings.fh"
#include "rinfo.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "wldata.fh"
#include "oneswi.fh"
#include "RelLight.fh"
#include "dkhparameters.fh"
      Integer DKHMemMax, DKHMemCheck, xOrder_save
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
      integer dkhorder,xorder,xorder_dum
      logical dkhscfflg
      Integer snumber,tnumber,unumber,nbasp,nbaso
      Integer Get_sNumber, Get_tNumber, Get_uNumber
      logical no_hamil,no_prop,no_s,no_u,LDKpert
      dimension nInt(1)
*                                                                      *
************************************************************************
*                                                                      *
      iRout=77
      iPrint=nPrint(iRout)
*
      If(Debug)Then
        idbg=6
      Else
        idbg=-1
      Endif
*
      stdout=6
      dbgunit=32
      dkhunit1=11
      dkhunit2=12
      dkhunit3=13
      dkhunit4=14
      dkhunit5=15
      outunit1=21
      outunit2=22
      eigunit=29
#ifdef MOLPRO
      call tmpnm(dkhunit1,filename)
      call molcas_open (dkhunit1, filename)
      call tmpnm(dkhunit2,filename)
      call molcas_open (dkhunit2, filename)
      call tmpnm(dkhunit3,filename)
      call molcas_open (dkhunit3, filename)
      call tmpnm(dkhunit4,filename)
      call molcas_open (dkhunit4, filename)
      call tmpnm(dkhunit5,filename)
      call molcas_open (dkhunit5, filename)
      call tmpnm(dbgunit,filename)
      call molcas_open (dbgunit, filename)
      call tmpnm(outunit1,filename)
      call molcas_open (outunit1, filename)
      call tmpnm(outunit2,filename)
      call molcas_open (outunit2, filename)
      call tmpnm(eigunit,filename)
      call molcas_open (eigunit, filename)
#endif
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
        If(dbsc(iCnttp)%Aux .or.
     &     dbsc(iCnttp)%Frag .or.
     &     dbsc(iCnttp)%nFragType.gt.0 ) Go To 999

*
        Do icnt = 1, dbsc(iCnttp)%nCntr
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
      iSizec=0
      Do L=0,nSym-1
         iSizep=iSizep + nBas(L)*(nBas(L)+1)/2
         iSizec=iSizec+nrBas(L+1)*(nrBas(L+1)+1)/2
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
      If (IRELAE.ge.1000) Then
          DKH_Verbose=iPrint.ge.6
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*        Use the new code of M. Reiher et al.
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*        memory management:
*        Out_Of_Core=.False.  (a lot of mem available -> store all
*                              square matrices in RAM)
*        Out_Of_Core=.True.   (export all matrices to disk: causes
*                              quite some I/O)
*                                                                      *
************************************************************************
*                                                                      *
         Out_Of_Core=.false.
*
*
*        IRELAE codes the method (DKH=1), the order (02 to 14), and the
*        parametrization (1 to 5), i.e. 3rd order DKH with sqrt
*        parametrization: 1-03-3, i.e. IRELAE=1033
*
         xOrder = iRELAE/10000
         iTemp=iRELAE-1000 - xOrder*10000
         dkhorder = iTemp/10
         iTemp = iTemp - dkhorder*10

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
*
*        Choose/adjust  "maxoperators" and "maxuops" for parser1
*        appropriately. Through common block!
*
         xorder_dum=xorder
         xOrder=Min(xOrder,dkhOrder) ! xorder never larger than dkhorder
         if (xorder_dum.ne.xorder) then
           write(stdout,*) "xorder was reduced from ",xorder_dum,
     *                     " to ",xorder
         end if
CMR  xOrder=0 may be used in order to assess the PCE, however, the DKH
CMR  code cannot process this case yet and currently we just skip DKH if
CMR  xorder=0
CMR  In any case, the following line must be skipped
CMR         xOrder=Max(2,xOrder)        ! oxoder never smaller than 2
         Call adjust_param (dkhorder,xorder)
*
*        Call symbolic parser procedure
*
*        dkhscfflg can be used as an input flag; however, it is set
*        FALSE so that the property operator is treated
*        perturbatively as in 4-comp.
*        1st-order response, i.e. the U_i do not depend on the property
*        (see part III).
*
         dkhscfflg=.False.
*
*                                                                      *
************************************************************************
*                                                                      *
         Call Allocate_Work(iK_Save,iSizep+4)
         Call Allocate_Work(iK_Done,iSizep+4)
         call dcopy_(iSizep+4,Work(iK),1,Work(iK_Save),1)
*
*        Do the Hamiltonian seperately
*
         xOrder_Save = xOrder
         xOrder=0
         No_Hamil=.False.
         No_Prop=.True.
CMR THE FOLLOWING HAS TO BE OTPIMIZED (mem VS. cpu); AND THEY MAY BE
CMR TREATED SEPARATELY ESPECIALLY IN COMBINATION WITH Out_Of_Core
CMR - IF THAT IS .TRUE. WE SHOULD WRITE XORDER.EQ.0
         if (dkhorder.le.4) then
           no_s=.true.
           no_u=.true.
         else
           no_s=.false.
           no_u=.false.
         end if
*
         Call DKHparser_symbolic(dkhorder,xorder,paramtype,dkhscfflg,
     &                           no_prop,no_s,no_u)
         snumber=Get_sNumber(dkhunit3)
         tnumber=Get_tNumber(dkhunit4)
         unumber=Get_uNumber(dkhunit5)
*
         no_hamil=.false.
         no_prop=.true.
         k=0
*
         Do L = 0, nSym-1
            n=nBas(L)
            iSize=n*(n+1)/2
            If (iSize.eq.0) Go To 911
*                                                                      *
************************************************************************
*                                                                      *
            Call Allocate_Work(ip1,iSizep+4)
            Call Allocate_Work(ip2,iSizep+4)
            Call GetMem('Core','Max','Real',iDum,Mem_Available)
*                                                                      *
************************************************************************
*                                                                      *
*           Determine how much memory is needed for DKHinf
*
            nbasp=1
            nbaso=n
#ifdef MOLPRO
#else
            Call get_natoms_all(nAtom)
#endif
            LDKpert=.false.
            maxsiz=n
            nblock=1
            If (LDKroll) Then
               Call GetMem('Index2 ','ALLO','INTE',indx2,4*nAtom)
               Call GetMem('Index  ','ALLO','INTE',indx,n)
               Call GetMem('Coord  ','ALLO','REAL',icoord,nAtom*3)
               Call calc_indx(iWork(indx2),iWork(indx),Work(icoord),
     &                        n,nAtom,maxsiz,nblock)
               Call GetMem('Coord  ','FREE','REAL',icoord,nAtom*3)
               Call GetMem('Index  ','FREE','INTE',indx,n)
*
               If ((nCtrLD.ne.0).or.(radiLD.ne.0.0d0)) LDKpert=.true.
            Else
               indx2=ip_iDummy
            End If
            If (LDKroll.and.(.not.LDKpert)) Then
               DKHMemCheck=6+maxsiz*(6+3*n)+n+3*isize+
     &                     maxsiz*maxsiz*(13+snumber+tnumber+unumber)+
     &                     maxsiz*(maxsiz+1)*2
               DKHMemCheck=DKHMemCheck+4*maxsiz*maxsiz+4
            Else
               DKHMemCheck=6+(9+snumber+tnumber+unumber)*n*n+6*n
c                       ^= 6*nbasp*nbasp
               DKHMemCheck=DKHMemCheck+isize*3+n*n*3
               DKHMemCheck=DKHMemCheck+4*n*n+4
            End If
c           Write (6,*) 'DKHMemCheck=',DKHMemCheck
c           Write (6,*) 'Mem_Available=',Mem_Available

            Out_Of_Core = DKHMemCheck.gt.Mem_Available
c           Write (6,*) 'Out_of_Core=',Out_of_Core
            Out_Of_Core = Out_Of_Core .or. Force_Out_Of_Core
c           Write (6,*) 'Out_of_Core=',Out_of_Core
*
            If (Out_Of_Core) Then
***            DKH matrices are written to disk.
               nbaso=1
               dkhmemmax=6+3*n*n+6*n+6+snumber+tnumber+unumber
c          6*nbasp*nbasp=^           ^= 6*nbaso*nbaso
               dkhmemmax=dkhmemmax+isize*3+n*n*3
               If ((dkhmemmax.gt.Mem_Available).or.(LDKroll)) Then
                  Write (6,*) 'DKRelInt: Insufficient memory(1)'
                  Write (6,*) 'dkhmemmax.gt.Mem_Available'
                  Write (6,*) 'dkhmemmax=',dkhmemmax
                  Write (6,*) 'Mem_Available=',Mem_Available
                  Call Quit(_RC_MEMORY_ERROR_)
               End If
            Else
               If (LDKroll.and.(.not.LDKpert)) Then
                  dkhmemmax=DKHMemCheck-4*maxsiz*maxsiz-4
               Else
                  dkhmemmax=DKHMemCheck-4*n*n-4
               End If
            End if
C           Write (6,*) 'dkhmemmax=',dkhmemmax
*
            CALL GetMem('DKHmem  ','ALLO','REAL',iDKHmem,dkhmemmax+4)
*                                                                      *
************************************************************************
*                                                                      *
            Call dkhparser_driver(n,isize,dkhscfflg,dkhorder,xorder,
     &                            Work(iSS+k), Work(iK+k),Work(iV+k),
*    &                            Work(ipVp+k),Work(ip_Dummy),
*    &                            Work(ip_Dummy),clightau,paramtype,
     &                            Work(ipVp+k),Work(ip1),Work(ip2),
     &                            clightau,paramtype,dkhmemmax,
     &                            Work(iDKHmem),no_hamil,no_prop,
     &                            nbasp,nbaso,LDKroll,iWork(indx2),
     &                            nAtom,maxsiz,nblock,LDKpert)
*                                                                      *
************************************************************************
*                                                                      *
            CALL GetMem('DKHmem  ','FREE','REAL',iDKHmem,dkhmemmax+4)
            Call Free_Work(ip2)
            Call Free_Work(ip1)
 911        k=k+isize
         End Do
*                                                                      *
************************************************************************
*                                                                      *
         call dcopy_(iSizep+4,Work(iK),1,Work(iK_Done),1)
*                                                                      *
************************************************************************
*                                                                      *
         xOrder=xOrder_Save
         If (xOrder.le.0) Go To 912
*                                                                      *
************************************************************************
*                                                                      *
*        Loop over all property matrices, do not process the Hamiltonian
*        here.
*
         No_Hamil=.True.
         No_Hamil=.False.
CMR THE FOLLOWING HAS TO BE OTPIMIZED (mem VS. cpu); AND S AND U MAY BE
CMR TREATED SEPARATELY ESPECIALLY IN COMBINATION WITH Out_Of_Core
CMR - IF THAT IS .TRUE. WE SHOULD WRITE XORDER.EQ.0
         if (xorder.le.1) then
           no_s=.true.
           no_u=.true.
         else
           no_s=.false.
           no_u=.false.
         end if
         Call DKHparser_symbolic(dkhorder,xorder,paramtype,dkhscfflg,
     &                           no_prop,no_s,no_u)
         snumber=Get_sNumber(dkhunit3)
         tnumber=Get_tNumber(dkhunit4)
         unumber=Get_uNumber(dkhunit5)
*
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
            Call iRdOne(iRC,iOpt,Label,iComp,nInt,lOper)
C           Write (6,*) 'lOper=',lOper
            CALL GetMem('X       ','ALLO','REAL',iX,nInt(1)+4)
            iRC = -1
            iOpt=0
            Call RdOne(iRC,iOpt,Label,iComp,Work(iX),lOper)
            If (iRC.ne.0) then
               Write (stdout,*) 'DKRelInt: Error reading from ONEREL'
               Write (stdout,'(A,A)') 'Label=',Label
               Write (stdout,'(A,A)') 'iRC=',iRC
               Call Abend
            End If
            Call CmpInt(Work(iX),nInt(1),nBas_Prim,nSym,lOper)
            If (nInt(1).eq.0) Then
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
            Call iRdOne(iRC,iOpt,pXpLbl,iComp,nInt,lOper)
            CALL GetMem('pXp     ','ALLO','REAL',ipXp,nInt(1)+4)
            iOpt=0
            iRC = -1
            Call RdOne(iRC,iOpt,pXpLbl,iComp,Work(ipXp),lOper)
            If (iRC.ne.0) then
               Write (stdout,*) 'DKRelInt: Error reading from ONEREL'
               Write (stdout,'(A,A)') 'pXpLbl=',pXpLbl
               Write (stdout,'(A,A)') 'iRC=',iRC
               Call Abend
            End If
            Call CmpInt(Work(ipXp),nInt(1),nBas_Prim,nSym,lOper)
*
            iOpt=0
            Call ClsOne(iRC,iOpt)
*
            Call GetMem('Core','Max','Real',iDum,Mem_Available)
C           Write (6,*) 'Mem_Available=',Mem_Available
            no_hamil=.true.
            no_prop=.false.
            k=0
            Do L = 0, nSym-1
               n=nBas(L)
               iSize=n*(n+1)/2
               If (iSize.eq.0) Go To 91
*
*              Skip if the propetry operator does not have a total
*              symmetric component!
*
               If(iAnd(1,lOper).eq.0) Go To 91
*
*              Determine how much memory is needed for DKHinf
*
               nbasp=n
               nbaso=n
               If (LDKroll) Then
                  DKHMemCheck=maxsiz*(6+3*n)+n+3*isize+
     &    maxsiz*maxsiz*(19+snumber+tnumber+unumber)+
     &    maxsiz*(maxsiz+1)*3
                  DKHMemCheck=DKHMemCheck+8*maxsiz*maxsiz+4
               Else
                  DKHMemCheck=(15+snumber+tnumber+unumber)*n*n+6*n
                  DKHMemCheck=DKHMemCheck+isize*5+n*n*3
                  DKHMemCheck=DKHMemCheck+8*n*n+4
               End If
C              Write (6,*) 'DKHMemCheck=',DKHMemCheck
               Out_Of_Core = DKHMemCheck.gt.Mem_Available
C              Write (6,*) 'Out_of_core=',Out_of_core
               Out_Of_Core = Out_Of_Core .or. Force_Out_Of_Core
C              Write (6,*) 'Out_of_core=',Out_of_core
*
               If (Out_Of_Core) Then
***               DKH matrices are written to disk.
                  nbasp=1
                  nbaso=1
                  dkhmemmax=6+3*n*n+6*n+6+snumber+tnumber+unumber
c             6*nbasp*nbasp=^           ^= 6*nbaso*nbaso
                  dkhmemmax=dkhmemmax+isize*5+n*n*3
                  If ((dkhmemmax.gt.Mem_Available).or.(LDKroll)) Then
                     Write (6,*) 'DKRelInt: Insufficient memory(2)'
                     Write (6,*) 'dkhmemmax.gt.Mem_Available'
                     Write (6,*) 'dkhmemmax=',dkhmemmax
                     Write (6,*) 'Mem_Available=',Mem_Available
                     Call Quit(_RC_MEMORY_ERROR_)
                  End If
               Else ! In-core
***               DKH matrices are kept in RAM.
                  If (LDKroll) Then
                      DKHMemCheck=DKHMemCheck-8*maxsiz*maxsiz-4
                   Else
                      dkhmemmax=DKHMemCheck-8*n*n-4
                   End If
               End If
C              Write (6,*) 'dkhmemmax=',dkhmemmax
               CALL GetMem('DKHmem  ','ALLO','REAL',iDKHmem,dkhmemmax+4)
*                                                                      *
************************************************************************
*                                                                      *
               Call dkhparser_driver(n,isize,dkhscfflg,dkhorder,xorder,
     &                               Work(iSS+k), Work(iK+k),Work(iV+k),
     &                               Work(ipVp+k),Work(iX+k),
     &                               Work(ipXp+k),clightau,paramtype,
     &                               dkhmemmax,Work(iDKHmem),no_hamil,
     &                               no_prop,nbasp,nbaso,LDKroll,
     &                               iWork(indx2),nAtom,maxsiz,nblock,
     &                               LDKpert)
*                                                                      *
************************************************************************
*                                                                      *
               CALL GetMem('DKHmem  ','FREE','REAL',iDKHmem,dkhmemmax+4)
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
            Call iRdOne(iRC,iOpt,Label,iComp,nInt,lOper)
            CALL GetMem('Y       ','ALLO','REAL',iY,nInt(1)+4)
            iRC = -1
            iOpt=0
            Call RdOne(iRC,iOpt,Label,iComp,Work(iY),lOper)
C           Write (6,*) 'Y1=',DDot_(nInt(1),Work(iY),1,
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
            Call Cp_Prop_Int(Work(iY),nInt(1),Work(ip_Prop),iSizec,
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
         If (LDKroll) Call GetMem('Index2 ','FREE','INTE',indx2,4*nAtom)
*                                                                      *
************************************************************************
*                                                                      *
         call dcopy_(iSizep+4,Work(iK_Done),1,Work(iK),1)
         Call Free_Work(iK_Done)
         Call Free_Work(iK_Save)
         If (DKH_Verbose) Then
            Write (stdout,7766)
 7766       Format (43X,'...  done.',//2X)
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
         k=0
         Do L = 0, nSym-1
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

*
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
      close(dkhunit1, status='DELETE')
      close(dkhunit2, status='DELETE')
      close(dkhunit3, status='DELETE')
      close(dkhunit4, status='DELETE')
      close(dkhunit5, status='DELETE')
      close(dbgunit, status='DELETE')
      close(outunit1, status='DELETE')
      close(outunit2, status='DELETE')
      close(eigunit, status='DELETE')
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
      Return
*
 9999 Continue
      Write (stdout,*) ' *** Error in subroutine DKRelInt ***'
      Write (stdout,*) '     Abend in subroutine OpnOne'
      Call Abend
      End
