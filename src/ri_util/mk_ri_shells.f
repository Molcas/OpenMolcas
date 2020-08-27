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
* Copyright (C) Roland Lindh                                           *
************************************************************************
      Subroutine Mk_RI_Shells(LuRd)
************************************************************************
*                                                                      *
*    Objective: To expand the data for the auxiliary functions         *
*                                                                      *
* Called from: RdCtl                                                   *
*                                                                      *
* Calling    : GetBS                                                   *
*                                                                      *
*     Author: Roland Lindh                                             *
*                                                                      *
************************************************************************
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
#include "real.fh"
#include "print.fh"
      Logical Hit, IfTest
      Character*13 DefNm
      Character*80 Ref(2), BSLbl, BSLB*180
      Character*80 atom,type,author,basis,CGTO, Aux
      Character*80 atomb
      Character *256 Basis_lib, Fname
      Character*180, Allocatable :: STDINP(:) !CGGn
      Character*180 Line, Get_Ln
      External Get_Ln
      Integer StrnLn
      External StrnLn
      Logical Quit_On_Error
      common /getlnQOE/ Quit_On_Error
      Character*180 Get_Ln_Quit

      Integer BasisTypes(4), nDel(MxAng)
      Data DefNm/'basis_library'/
*                                                                      *
************************************************************************
*                                                                      *
#include "getbs_interface.fh"
*                                                                      *
************************************************************************
*                                                                      *
      Call qEnter('Mk_RI_Shells')
      iRout = 2
      iPrint = nPrint(iRout)
*
      Call mma_allocate(STDINP,mxAtom*2,label='STDINP')
      IfTest=.False.
*     IfTest=.True.
*
*     Add the auxiliary basis set
*
      BasisTypes(1)=0
      BasisTypes(2)=0
      BasisTypes(3)=0
      BasisTypes(4)=0
      iShll = Mx_Shll - 1
      lSTDINP=0
      mCnttp=nCnttp
*
*     Branch to special loop for reading external RICD basis sets
*     since these have a different infrastructure.
*
      If (iRI_Type.eq.5) Go To 1000

      Do iCnttp = 1, mCnttp
         If (dbsc(iCnttp)%Frag.or.dbsc(iCnttp)%nVal.eq.0) cycle
         mdc = dbsc(iCnttp)%mdci
         nCnttp=nCnttp+1
*
         If (nCnttp.gt.Mxdbsc) Then
            Call WarningMessage(2,'Error in Mk_RI_Shells')
            Write (6,*) 'Mk_RI_Shells: Increase Mxdbsc'
            Call Abend()
         End If
*
*        Resolve the name of the valence basis and find the name of
*        the appropriate auxiliary basis set.
*
         Bsl(nCnttp)=Bsl_Old(iCnttp)
*
         Hit=.True.
         Call Decode(Bsl(nCnttp),atom,1,Hit)
         Hit=.True.
         Call Decode(Bsl(nCnttp),type,2,Hit)
         Hit=.True.
         Call Decode(Bsl(nCnttp),author,3,Hit)
         Hit=.True.
         Call Decode(Bsl(nCnttp),basis,4,Hit)
         Hit=.True.
         Call Decode(Bsl(nCnttp),CGTO,5,Hit)
         Hit=.False.
         Call Decode(Bsl(nCnttp),Aux,6,Hit)
         If (.Not.Hit) Aux = ' '
*
         n=Index(Atom,' ')-1
         Bsl(nCnttp)(1:n+1)=atom(1:n)//'.'
         nn = n + 1
*
         n=Index(Type,' ')-1
         Bsl(nCnttp)(nn+1:nn+n+5)=Type(1:n)//'.....'
*
*        Modify basis set library correctly
*
         Indx=Index(Bsl(nCnttp),' ')
         BSLbl=' '
         BSLbl(1:Indx-1)=Bsl(nCnttp)(1:Indx-1)
         Call WhichMolcas(Basis_lib)
         If (Basis_lib(1:1).ne.' ') Then
            ib=index(Basis_lib,' ')-1
            If(ib.lt.1)
     &       Call SysAbendMsg('rdCtl','Too long PATH to MOLCAS',' ')
            If (iRI_Type.eq.1) Then
               Fname=Basis_lib(1:ib)//'/basis_library/j_Basis'
            Else If (iRI_Type.eq.2) Then
               Fname=Basis_lib(1:ib)//'/basis_library/jk_Basis'
            Else If (iRI_Type.eq.3) Then
               Fname=Basis_lib(1:ib)//'/basis_library/c_Basis'
            Else
               Call WarningMessage(2,'Error in Mk_RI_Shells')
               Write (6,*) 'Wrong iRI_Type!'
               Write (6,*) 'iRI_Type=',iRI_Type
               Call Abend()
            End If
         Else
            If (iRI_Type.eq.1) Then
               Fname=DefNm//'/j_Basis'
            Else If (iRI_Type.eq.2) Then
               Fname=DefNm//'/jk_Basis'
            Else If (iRI_Type.eq.3) Then
               Fname=DefNm//'/c_Basis'
            Else
               Call WarningMessage(2,'Error in Mk_RI_Shells')
               Write (6,*) 'Wrong iRI_Type!'
               Write (6,*) 'iRI_Type=',iRI_Type
               Call Abend()
            End If
         End If
*
         If (Show.and.iPrint.ge.6) Then
            Write (6,*)
            Write (6,*)
            Write(6,'(1X,A,I5,A,A)')
     &              'Basis Set ',nCnttp,' Label: ', BSLbl(1:Indx-1)
            Write(6,'(1X,A,A)') 'Basis set is read from library:',
     &                          Fname(1:index(Fname,' '))
         End if
*
         jShll = iShll
         SODK(nCnttp)=.False.
         Bsl_Old(nCnttp)=Bsl(nCnttp)
         Call GetBS(Fname,Bsl(nCnttp),Indx-1,lAng,iShll,
     &              MxAng,Charge(nCnttp),
     &              iAtmNr(nCnttp),BLine,Ref,PAM2(nCnttp),
     &              NoPairL(nCnttp),SODK(nCnttp),
     &              CrRep(nCnttp),nProj,nAIMP,UnNorm,nDel,
     &               nVal,   nPrj,   nSRO,   nSOC,  nPP,
     &              ipVal_, ipPrj_, ipSRO_, ipSOC_,ipPP_,
     &              LuRd,BasisTypes,
     &              STDINP,lSTDINP,.False.,.true.,' ',nCnttp)
         dbsc(nCnttp)%Aux=.True.
*
         Charge(nCnttp)=Zero
*
         If (Show.and.iPrint.ge.6 .and.
     &      Ref(1).ne.BLine .and. Ref(2).ne.Bline) Then
            Write (6,'(1x,a)')  'Basis Set Reference(s):'
            If (Ref(1).ne.BLine) Write (6,'(5x,a)') Ref(1)
            If (Ref(2).ne.BLine) Write (6,'(5x,a)') Ref(2)
            Write (6,*)
            Write (6,*)
         End If
         lPAM2 = lPAM2 .or. PAM2(nCnttp)
         dbsc(nCnttp)%ECP=(nPrj+nSRO+nSOC+
     &                dbsc(nCnttp)%nM1+dbsc(nCnttp)%nM2).ne.0
         lPP=lPP .or. nPP.ne.0
         lECP = lECP .or. dbsc(nCnttp)%ECP
         lNoPair = lNoPair .or. NoPairL(nCnttp)
*
         iAngMx=Max(iAngMx,lAng)
*        No transformation needed for s and p shells
         Shells(jShll+1)%Transf=.False.
         Shells(jShll+1)%Prjct =.False.
         Shells(jShll+2)%Transf=.False.
         Shells(jShll+2)%Prjct =.False.
         pChrg(nCnttp)=pChrg(iCnttp)
         Fixed(nCnttp)=Fixed(iCnttp)
         dbsc(nCnttp)%Parent_iCnttp=iCnttp
C        pChrg(nCnttp)=.False.
C        Fixed(nCnttp)=.False.
         dbsc(nCnttp)%iVal = ipVal_
         dbsc(nCnttp)%iPrj = ipPrj_
         ipSRO(nCnttp) = ipSRO_
         ipSOC(nCnttp) = ipSOC_
         ipPP(nCnttp)  = ipPP_
         dbsc(nCnttp)%nVal = nVal
         dbsc(nCnttp)%nPrj = nPrj
         nSRO_Shells(nCnttp) = nSRO
         nSOC_Shells(nCnttp) = nSOC
         nPP_Shells(nCnttp)  = nPP
         nTot_Shells(nCnttp) = nVal+nPrj+nSRO+nSOC+nPP
         lAux = lAux .or. dbsc(nCnttp)%Aux
         Do iSh = jShll+1, iShll
            Shells(iSh)%nBasis=Shells(iSh)%nBasis_c
            Call mma_deallocate(Shells(iShll)%pCff)
            Call mma_allocate(Shells(iShll)%pCff,
     &                        Shells(iSh)%nExp,Shells(iSh)%nBasis,
     &                        Label='pCff')
            Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
            Shells(iSh)%Aux=.True.
         End Do
*
         nCnt = dbsc(iCnttp)%nCntr
         dbsc(nCnttp)%nCntr=nCnt
         dbsc(nCnttp)%mdci =mdc
         Call mma_allocate(dbsc(nCnttp)%Coor,3,nCnt,Label='dbsc:C')
         dbsc(nCnttp)%Coor(:,:)=dbsc(iCnttp)%Coor(:,:)
*
         Mx_Shll=iShll+1
         Max_Shells=Mx_Shll
         Mx_mdc=mdc
*
      End Do
      Go To 1100
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Specially designed loop to read a RICD auxiliary basis set from
*     an external library.
*
 1000 Continue
*
      Lu_lib=17
      Lu_lib=IsFreeUnit(Lu_lib)
      call molcas_open(Lu_lib,'RICDLIB')
*
      Do iCnttp = 1, mCnttp
         If (dbsc(iCnttp)%Frag.or.dbsc(iCnttp)%nVal.eq.0) cycle
         mdc = dbsc(iCnttp)%mdci
*
         Hit=.True.
         Call Decode(Bsl_Old(iCnttp),atom,1,Hit)
         Type=' '
         Author=' '
         basis=' '
         CGTO=' '
         Aux=' '
         If (IfTest) Then
            Write (6,*) 'Bsl_Old=',Bsl_Old(iCnttp)
            Write (6,*) 'Atom=',Atom
         End If
*
         Indx=Index(Bsl_Old(iCnttp),' ')
         BSLbl=' '
         BSLbl(1:Indx-1)=Bsl_Old(iCnttp)(1:Indx-1)
*
*        Find the basis set
*
         ReWind(Lu_lib)
*
*        Loop over the basis set library to find the correct label
*
         If (IfTest) Write (6,*) ' Locate basis set label in library'
   10    BSLB = Get_Ln_Quit(Lu_lib,0)
         If (Quit_On_Error) Then
            iLast3=StrnLn(BsLbl)
            Call WarningMessage(2,
     &          'The requested basis set label: '//
     &          BsLbl(:iLast3)//';'//
     &          'was not found in basis library: '//'RICDLIB')
            Call Abend()
         End If
*
         Call UpCase(BSLB)
         If (BSLB(1:1).ne.'/') Go To 10
         If (IfTest) Write(6,*) 'BSLB=',BSLB
         n=Index(BSLB,' ')
         Do i=n,80
            BSLB(i:i)='.'
         End Do
         Hit=.True.
         Call Decode(BSLB(2:80),atomb,1,Hit)
         If (atomb.ne.atom) Go To 10
         If (IfTest) Write(6,*) 'atomb=',atomb
*
*        Now we should have found the correct basis set label!
*
         nSet=-1

         Do While (nSet.ne.0)
            Line=Get_Ln(Lu_lib)
            If (IfTest) Then
               Write(6,*) 'nSet=',nSet
               Write(6,*) 'Line=',Line
            End If
            Call Get_I1(2,lAng)
            If (nSet.eq.-1) Call Get_I1(3,nSet)
            If (IfTest) Write(6,*) 'lAng,nSet=',lAng,nSet
*
            Line=Get_Ln(Lu_lib)
            Line=Get_Ln(Lu_lib)
*
            nCnttp=nCnttp+1
            If (nCnttp.gt.Mxdbsc) Then
               Call WarningMessage(2,'Error in Mk_RI_Shells')
               Write (6,*) 'Mk_RI_Shells: Increase Mxdbsc'
               Call Abend()
            End If
            If (Show.and.iPrint.ge.6) Then
               Write (6,*)
               Write (6,*)
               Write(6,'(1X,A,I5,A,A)')
     &                 'Basis Set ',nCnttp,' Label: ', BSLb
               Write(6,'(1X,A)') 'Basis set is read from the workdir.'
            End if
*
            SODK(nCnttp)=.False.
            Bsl(nCnttp)=BSLB(2:80)
            Bsl_Old(nCnttp)=Bsl(nCnttp)
*
*           Loop over the angular shells
*
            jShll = iShll
            Do iAng = 0, lAng
               iShll=iShll+1
               Line=Get_Ln(Lu_lib)
               Call Get_I1(1,nPrim)
               Call Get_I1(2,nCntrc)
               Call Get_I1(3,iSph)
               If (IfTest) Then
                  Write (6,*) 'iAng=',iAng
                  Write (6,*) 'nPrim=',nPrim
                  Write (6,*) 'nCntrc=',nCntrc
                  Write (6,*) 'iSph=',iSph
               End If
*
*              Read Gaussian exponents
*
               Call mma_Allocate(Shells(iShll)%Exp,nPrim,Label='ExpRI')
               Shells(iShll)%nExp=nPrim
               Shells(iShll)%nBasis_C = nCntrc
               iEnd = iStrt - 1
               If (nPrim.gt.0) then
                  If (IfTest) Write(6,*) ' Read gaussian exponents'
                  Call Read_v(Lu_lib,Shells(iShll)%Exp,1,nPrim,1,Ierr)
                  If (Ierr.ne.0) Then
                     Call WarningMessage(2,
     &                     'GetBS: Error while reading the exponents')
                     Call Quit_OnUserError()
                  End If
                  If (IfTest) Write(6,*) ' Done with exponents'
               If (iPrint.ge.99.or.IfTest)
     &           Call RecPrt(' Exponents',' ',Shells(iShll)%Exp,nPrim,1)
               End If
               iStrt = iEnd + 1
*
*              Read contraction coefficients. Storage of coefficients
*              for both contracted and uncontracted case.
*
               Call mma_allocate(Shells(iShll)%Cff_c,nPrim,nCntrc,2,
     &                           Label='Cff_c')
               Call mma_allocate(Shells(iShll)%pCff,nPrim,nCntrc,
     &                           Label='pCff')
               Shells(iShll)%nBasis=nCntrc
               Call mma_allocate(Shells(iShll)%Cff_p,nPrim,nPrim,2,
     &                           Label='Cff_p')
               iEnds= iEnd
               iEnds= iEnds
               iEndc = iStrt - 1
               iEnd  = iStrt - 1
*              Read contraction coefficients
*              Observe that the matrix will have nPrim rows and
*              nCntrc columns
               If (IfTest) Write (6,*) ' Read coefficients'
*
*              Read in coeffs. in GC format, as the standard case
*
               If (IfTest) Write (6,*) ' Standard case'
               If (nPrim*nCntrc.gt.0) Then
                  Shells(iShll)%Cff_c(:,:,:)=Zero
*
*              Note that we now change the order!!!
               Read (Lu_lib,*) ((Shells(iShll)%Cff_c(i,j,2),
     &                           j=1,nCntrc),i=1,nPrim)
               If (Ierr.ne.0) Then
                  Call WarningMessage(2,
     &                   'GetBS: Error reading coeffs in GC format')
                  Call Quit_OnUserError()
               End If
               If (IfTest)
     &            Call RecPrt(' Coeffs',' ',Shells(iShll)%Cff_c(1,1,2),
     &                        nPrim,nCntrc)
               Shells(iShll)%Cff_c(:,:,1) = Shells(iShll)%Cff_c(:,:,2)
               If (IfTest)
     &            Call RecPrt(' Coeffs',' ',Shells(iShll)%Cff_c(1,1,1),
     &                        nPrim,nCntrc)
*
*              Put in unit matrix of uncontracted set
*
               Shells(iShll)%Cff_p(:,:,1)=Zero
               Do i=1,nPrim
                  Shells(iShll)%Cff_p(i,i,1)=One
               End Do
*
               Shells(iShll)%Cff_p(:,:,2)=Shells(iShll)%Cff_p(:,:,1)
               Call Nrmlz(Shells(iShll)%Exp,nPrim,
     &                    Shells(iShll)%Cff_p(1,1,1),nPrim ,iAng)

               Shells(iShll)%Cff_c(:,:,2)=Shells(iShll)%Cff_c(:,:,1)
               Call Fix_Coeff(nPrim,nCntrc,Shells(iShll)%Cff_c(1,1,2),
     &                                     Shells(iShll)%Cff_p(1,1,1),
     &                                     'F')
               Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
               End If
*
               iEnd =iEnds
               If (iSph.eq.0) Then
                  Shells(iShll)%Transf=.False.
                  Shells(iShll)%Prjct =.False.
               Else If (iSph.eq.1) Then
                  Shells(iShll)%Transf=.False.
                  Shells(iShll)%Prjct =.True.
               Else If (iSph.eq.2) Then
                  Shells(iShll)%Transf=.True.
                  Shells(iShll)%Prjct =.False.
               Else
                  Shells(iShll)%Transf=.True.
                  Shells(iShll)%Prjct =.True.
               End If

               Shells(iShll)%nBasis=Shells(iShll)%nBasis_C
               Shells(iShll)%Aux=.True.
*
            End Do ! iAng
*
            dbsc(nCnttp)%Aux=.True.
            Charge(nCnttp)=Zero
            PAM2(nCnttp)=.False.
            lPAM2 = lPAM2 .or. PAM2(nCnttp)
            nVal=lAng+1
            nPrj=0
            nSRO=0
            nSOC=0
            nPP=0
            dbsc(nCnttp)%ECP=.False.
            lECP = lECP .or. dbsc(nCnttp)%ECP
            lPP=lPP .or. nPP.ne.0
            NoPairL(nCnttp)=.False.
            lNoPair = lNoPair .or. NoPairL(nCnttp)
            iAngMx=Max(iAngMx,lAng)
*
            dbsc(nCnttp)%nOpt = 0
            dbsc(nCnttp)%iVal = jShll + 1
            dbsc(nCnttp)%iPrj = -1
            ipSRO(nCnttp) = -1
            ipSOC(nCnttp) = -1
            ipPP(nCnttp)  = -1
*
            dbsc(nCnttp)%nVal = nVal
            dbsc(nCnttp)%nPrj = nPrj
            nSRO_Shells(nCnttp) = nSRO
            nSOC_Shells(nCnttp) = nSOC
            nPP_Shells(nCnttp)  = nPP
            nTot_Shells(nCnttp) = nVal+nPrj+nSRO+nSOC+nPP
            lAux = lAux .or. dbsc(nCnttp)%Aux
*
            nCnt = dbsc(iCnttp)%nCntr
            dbsc(nCnttp)%nCntr=nCnt
            dbsc(nCnttp)%mdci =mdc
            Call mma_allocate(dbsc(nCnttp)%Coor,3,nCnt,Label='dbsc:C')
            dbsc(nCnttp)%Coor(:,:)=dbsc(iCnttp)%Coor(:,:)
*
            Mx_Shll=iShll+1
            Max_Shells=Mx_Shll
            Mx_mdc=mdc
*
            nSet=nSet-1
            If (nSet.ne.0) Line=Get_Ln(Lu_lib)
         End Do ! Do While (nSet.ne.0)
*
      End Do ! iCnttp
*
      Close(Lu_lib)
*                                                                      *
************************************************************************
*                                                                      *
*     Add the final DUMMY SHELL!
*
 1100 Continue
      Call Mk_Dummy_Shell()
      Call mma_deallocate(STDINP)
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('Mk_RI_Shells')
      Return
      End
