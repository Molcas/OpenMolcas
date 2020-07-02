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
* Copyright (C) 1990,2020,  Roland Lindh                               *
*               1990, IBM                                              *
*               1993, Per Boussard                                     *
************************************************************************
      SubRoutine GetECP(lUnit,ipExp,ipCff,nExp,nBasis,MxShll,iShll,
     &                  BLine,ipBk,CrRep,nProj,ipAkl,ip_Occ,
     &                  ipPP,nPP,UnNorm,DInf,nDInf,nCnttp)
************************************************************************
*                                                                      *
*    Objective: To read ECP information, excluding the valence basis-  *
*               set. This means that we read (from input stream) the   *
*               effective charge, the model potential terms (M1 and    *
*               M2), projection parameters and projection orbital      *
*               basis-set (exponents and coefficients)                 *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : RecPrt, Rdbsl                                           *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*                                                                      *
*     Modified: Per Boussard -93.                                      *
************************************************************************
      Use Basis_Info, only: dbsc
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
*for ip_Dummy:
      Character*180 Line, Get_Ln
      Character*(*) BLine
*     External Get_Ln
      Real*8, Dimension(:), Allocatable :: Scrt1, Scrt2
      Integer ipExp(MxShll), ipCff(MxShll), ipBk(MxShll),
     &        nExp(MxShll), nBasis(MxShll),
     &        ipAkl(MxShll), ip_Occ(MxShll), mPP(2)
      Logical UnNorm
      Real*8 DInf(nDInf)
      Integer nCnttp
*                                                                      *
************************************************************************
*                                                                      *
      iRout=6
      iPrint = nPrint(iRout)
      Call QEnter('GetECP')
*
      ip_Dummy=-1
      iStrt = ipExp(iShll+1)
*                                                                      *
************************************************************************
*                                                                      *
*---- Process Pseudo Potentials
*
      Line=Get_Ln(lUnit)
*
      If (Index(Line,'PP').ne.0) Then
*
         ipPP=iShll+1
*
C        Write (6,*) 'Line=',Line
         If (Index(Line,'PPSO').ne.0) Then
           Call Get_i(4,mPP,2)
         Else
           Call Get_i(4,mPP,1)
           mPP(2)=0
         Endif
C        Write (6,*) 'mPP=',mPP
         nPP = 1 + mPP(1) + mPP(2)
*
         Do iPP = 0, mPP(1)
            iShll = iShll + 1
C           Write (6,*) 'iPP,nPP=',iPP,nPP
            If (iShll.gt.MxShll) Then
               Call WarningMessage(2,
     &                     'Abend in GetECP: Increase MxShll')
               Call Abend()
            End If
*
*---------- Pick up the number of terms in the shell
            Line=Get_Ln(lUnit)
C           Write (6,*) 'Line=',Line
            Call Get_i1(1,kcr)
C           Write (6,*) 'kcr,iShll=',kcr,iShll
            nExp(iShll) = kcr
            ipExp(iShll) = iStrt
*
            Do jcr = 1, kcr
               Line=Get_Ln(lUnit)
*
               Call Get_I1(1,ncr)
C              Write (6,*) 'ncr=',ncr
               DInf(iStrt)=DBLE(ncr)
               iStrt=iStrt+1
               Call Get_F1(2,zcr)
C              Write (6,*) 'zcr=',zcr
               DInf(iStrt)=zcr
               iStrt=iStrt+1
               Call Get_F1(3,ccr)
C              Write (6,*) 'ccr=',ccr
               DInf(iStrt)=ccr
               iStrt=iStrt+1
            End Do
*
         End Do
         Do iPP = 1, mPP(2)
            iShll = iShll + 1
C           Write (6,*) 'iPP,nPP=',iPP,nPP
            If (iShll.gt.MxShll) Then
               Call ErrTra
               Write (6,*) 'Abend in GetECP: Increase MxShll'
               Call Abend
            End If
*
*---------- Pick up the number of terms in the shell
            Line=Get_Ln(lUnit)
C           Write (6,*) 'Line=',Line
            Call Get_i1(1,kcr)
C           Write (6,*) 'kcr,iShll=',kcr,iShll
            nExp(iShll) = kcr
            ipExp(iShll) = iStrt
*
            Do jcr = 1, kcr
               Line=Get_Ln(lUnit)
*
               Call Get_I1(1,ncr)
C              Write (6,*) 'ncr=',ncr
               ncr=ncr+1000
               DInf(iStrt)=DBLE(ncr)
               iStrt=iStrt+1
               Call Get_F1(2,zcr)
C              Write (6,*) 'zcr=',zcr
               DInf(iStrt)=zcr
               iStrt=iStrt+1
               Call Get_F1(3,ccr)
C              Write (6,*) 'ccr=',ccr
               DInf(iStrt)=ccr
               iStrt=iStrt+1
            End Do
*
         End Do
         iEnd = iStrt - 1
         ipExp(iShll+1) = iEnd + 1
C        Write (6,*) 'Done'
*
         Go To 999
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Write (6,*) ' Reading M1'
      If (Index(Line,'M1').eq.0) Then
         Call WarningMessage(2,
     &              'ERROR: Keyword M1 expected, offending line : '
     &            //Line)
         Call Quit_OnUserError()
      Endif
*     Read(Line,*)nM1
      Line=Get_Ln(lUnit)
      Call Get_i1(1,nM1)
*     Write (*,*) ' nM1=',nM1
      dbsc(nCnttp)%nM1=nM1
      Call mma_allocate(dbsc(nCttp)%M1xp,nM1,Label='dbsc:M1xp')
      Call mma_allocate(dbsc(nCttp)%M1cf,nM1,Label='dbsc:M1cf')
* Note: all broadcasting of DInf will be done at end of getbs!!
*     If (nM1.gt.0) Read(lUnit,*) (dbsc(nCnttp)%M1xp(i),i=1,nM1)
      If (nM1.gt.0) Call Read_v(lUnit,dbsc(nCntp)%M1xp,1,nM1,1,ierr)
*     If (nM1.gt.0) Read(lUnit,*) (dbsc(nCnttp)%M1cf(i),i=1,nM1)
      If (nM1.gt.0) Call Read_v(lUnit,dbsc(nCnttp)%M1cf,1,nM1,1,ierr)
*                                                                      *
************************************************************************
*                                                                      *
*     Write (*,*) ' Reading M2'
      Line=Get_Ln(lUnit)
      If (Index(Line,'M2').eq.0) Then
         Call WarningMessage(2,
     &              'ERROR: Keyword M2 expected, offending line : '
     &            //Line)
         Call Quit_OnUserError()
      Endif
*     Read(Line,*)nM2
      Line=Get_Ln(lUnit)
      Call Get_i1(1,nM2)
      dbsc(nCnttp)%nM2=nM2
      Call mma_allocate(dbsc(nCttp)%M2xp,nM2,Label='dbsc:M2xp')
      Call mma_allocate(dbsc(nCttp)%M2cf,nM2,Label='dbsc:M2cf')
*     Write (*,*) ' nM2=',nM2
*     If (nM2.gt.0) Read(lUnit,*) (dbsc(nCnttp)%M2xp(i),i=1,nM2)
      If (nM2.gt.0) Call Read_v(lUnit,dbsc(nCnttp)%M2xp,1,nM2,1,ierr)
*     If (nM2.gt.0) Read(lUnit,*) (dbsc(nCnttp)%M2cf(i),i=1,nM2)
      If (nM2.gt.0) Call Read_v(lUnit,dbsc(nCnttp)%M2cf,1,nM2,1,ierr)
      ipExp(iShll+1) = iEnd + 1
*                                                                      *
************************************************************************
*                                                                      *
*     Read core-repulsion parameter
*
*     Write (6,*) ' Reading Core-repulsion'
      Line=Get_Ln(lUnit)
      If (Index(Line,'COREREP').eq.0) Then
         Call WarningMessage(2,
     &              'ERROR: Keyword COREREP expected, offending line : '
     &            //Line)
         Call Quit_OnUserError()
      Endif
*     Read(Line,*)CrRep
      Line=Get_Ln(lUnit)
      Call Get_F1(1,CrRep)
*                                                                      *
************************************************************************
*                                                                      *
*     Now read Projection parameters
*
*     Write (6,*) ' Reading projection operator'
      Line=Get_Ln(lUnit)
      If (Index(Line,'PROJOP').eq.0) Then
         Call WarningMessage(2,
     &              'ERROR: Keyword PROJOP expected, offending line : '
     &            //Line)
         Call Quit_OnUserError()
      Endif
*
*     Now read projection basis set
*     Loop over each shell type (s,p,d,etc....)
*     The ECP is assumed to contain s,p,d and f functions,
*     Read(Line,*,Err=993) nProj
      Line = Get_Ln(lUnit)
      Call Get_I1(1,nProj)
*
      If (nProj.lt.0) Go To 999
      Do 10 iAng = 0, nProj
*        Write (6,*) ' iAng=',iAng
         n_Elec=2*(2*iAng+1)
         iShll = iShll + 1
         If (iShll.gt.MxShll) Then
            Call WarningMessage(2,
     &                  'Abend in GetECP: Increase MxShll')
            Call Quit_OnUserError()
         End If
*        Read (Line,*,Err=993) nPrim, nCntrc
         Line = Get_Ln(lUnit)
         Call Get_I1(1,nPrim)
         Call Get_i1(2,nCntrc)
*
         nExp(iShll) = nPrim
         nBasis(iShll) = nCntrc
*
*------- Check if occupation number is included on the line
*
         ip_Occ(iShll)=iStrt
         iSS=1
         Call NxtWrd(Line,iSS,iEE)
         iSS=iEE+1
         Call NxtWrd(Line,iSS,iEE)
         iSS=iEE+1
         Call NxtWrd(Line,iSS,iEE)
         If (iEE.gt.0) Then
            Do i = 1, nCntrc
               Call Get_i1(2+i,n_Occ)
               DInf(iStrt)=DBLE(n_Occ)/DBLE(n_Elec)
               iStrt=iStrt+1
*              Write (*,*) 'n_Occ=',n_Occ
            End Do
         Else
            Do i = 1, nCntrc
*              Write (*,*) 'n_Occ=',n_Elec
               DInf(iStrt)=One
               iStrt=iStrt+1
            End Do
         End If
*
*        Read "orbital energies"
*
*        Write (6,*) ' Reading Bk'
         ipBk(iShll) = iStrt
         ipAkl(iShll) = ip_Dummy
         iEnd = iStrt + nCntrc - 1
*        If (nCntrc.gt.0) Read (lUnit,*,Err=992) (DInf(i),i=iStrt,iEnd)
         If (nCntrc.gt.0) Call Read_v(lUnit,DInf,iStrt,iEnd,1,ierr)
         If (ierr.ne.0) goto 992
         iStrt = iEnd + 1
*
*        Read gaussian EXPonents
*
*        Write (6,*) ' Reading Exponents'
         ipExp(iShll) = iStrt
         nExp(iShll) = nPrim
         iEnd = iStrt + nPrim - 1
*        If (nPrim.gt.0) Read (lUnit,*,Err=992) (DInf(i),i=iStrt,iEnd)
         If (nPrim.gt.0) Call Read_v(lUnit,DInf,iStrt,iEnd,1,ierr)
         If (ierr.ne.0) goto 992
         iStrt = iEnd + 1
*        Call RecPrt(' Exponents',DInf(iStrt),nPrim,1)
         ipCff(iShll) = iStrt
         nBasis(iShll) = nCntrc
         iEnd = iStrt + nPrim*nCntrc - 1
*
*        Read contraction coefficients
*        Observe that the matrix will have nPrim rows and
*        nCntrc columns
*
*        Write (6,*) ' Reading coefficients'
         Do 20 iPrim = 0, nPrim-1
*           Read (lUnit,*,Err=991) (DInf(i),i=iStrt+iPrim,iEnd,nPrim)
            Call Read_v(lUnit,DInf,iStrt+iPrim,iEnd,nPrim,ierr)
            If(ierr.ne.0) goto 991
 20      Continue
*
*--------Renormalize
*
         call mma_allocate(Scrt1,nPrim**2)
         call mma_allocate(Scrt2,nPrim*nCntrc)
         Call Nrmlx(DInf(ipExp(iShll)),nPrim,DInf(ipCff(iShll)),
     &                   nCntrc, Scrt1,nPrim**2,
     &                           Scrt2, nPrim*nCntrc,iAng)
         call mma_deallocate(Scrt1)
         call mma_deallocate(Scrt2)
*
*        Duplicate!
*
         iStrt = iEnd + 1
         iEnd = iStrt + nPrim*nCntrc - 1
         Call DCopy_(nPrim*nCntrc,DInf(ipCff(iShll)),1,DInf(iStrt),1)
*
         If (.Not.UnNorm) Then
            If (nExp(iShll)*nBasis(iShll).ge.1) Then
               Call Nrmlz(DInf(ipExp(iShll)),nExp(iShll),
     &                    DInf(ipCff(iShll)),nBasis(iShll),iAng)
            End If
         End If
*
         iStrt = iEnd + 1
         If (iShll.lt.MxShll) ipExp(iShll+1) = iEnd + 1
*
 10   Continue
 999  Continue
*
      Call QExit('GetECP')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_character(BLine)
*
*                                                                      *
************************************************************************
*                                                                      *
*998  Continue
*        Call WarningMessage(2,
*    &      'Abend in GetBS: Basis set label not found in library')
*        Call Quit_OnUserError()
*997  Continue
*        Call WarningMessage(2,
*    &      'Abend in GetBS: Error on opening the basis set library')
*        Call Quit_OnUserError()
*996  Continue
*        Call WarningMessage(2,
*    &      'Abend in GetBS: No library file available with this name')
*        Call Quit_OnUserError()
*995  Continue
*        Call WarningMessage(2,
*    &      'Abend in GetBS: Error while reading from basis library')
*        Call Quit_OnUserError()
*993  Continue
*        Call WarningMessage(2,
*    &      'Abend in GetBS: Error while reading the primitive line'
*    &                //' Line reads:'//Line)
*        Call Quit_OnUserError()
 992  Continue
         Call WarningMessage(2,
     &      'Abend in GetBS: Error while reading the exponents')
         Call Quit_OnUserError()
 991  Continue
         Call WarningMessage(2,
     &      'Abend in GetBS: Error while reading the coefficients')
         Call Quit_OnUserError()
      End
