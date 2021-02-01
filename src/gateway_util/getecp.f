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
#include "compiler_features.h"
#ifdef _IN_MODULE_
      SubRoutine GetECP(lUnit,iShll,nProj,UnNorm)
      Use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
      Integer lUnit
      Integer iShll
      Integer nProj
      Logical UnNorm
*     Local variables
      Character(LEN=180) Line, Get_Ln
*     External Get_Ln
      Integer mPP(2)
      Real*8, Dimension(:), Allocatable :: Scrt1, Scrt2
*                                                                      *
************************************************************************
*                                                                      *
      Interface
         Subroutine Read_v(lunit,Work,istrt,iend,inc,ierr)
         Integer lUnit, iStrt, Inc, iErr
         Real*8 Work(iend)
         End Subroutine Read_v
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
*---- Process Pseudo Potentials
*
      Line=Get_Ln(lUnit)
*
      If (Index(Line,'PP').ne.0) Then
*
         dbsc(nCnttp)%iPP=iShll+1
*
C        Write (6,*) 'Line=',Line
         If (Index(Line,'PPSO').ne.0) Then
           Call Get_i(4,mPP,2)
         Else
           Call Get_i(4,mPP,1)
           mPP(2)=0
         Endif
C        Write (6,*) 'mPP=',mPP
         dbsc(nCnttp)%nPP = 1 + mPP(1) + mPP(2)
*
         Do iPP = 0, mPP(1)
            iShll = iShll + 1
C           Write (6,*) 'iPP,dbsc(nCnttp)%nPP=',iPP,dbsc(nCnttp)%nPP
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
            Shells(iShll)%nExp = 3*kcr
            Call mma_allocate(Shells(iShll)%Exp,3*kcr,Label='Exp')
*
            iStrt=1
            Do jcr = 1, kcr
               Line=Get_Ln(lUnit)
*
               Call Get_I1(1,ncr)
C              Write (6,*) 'ncr=',ncr
               Shells(iShll)%Exp(iStrt)=DBLE(ncr)
               iStrt=iStrt+1
               Call Get_F1(2,zcr)
C              Write (6,*) 'zcr=',zcr
               Shells(iShll)%Exp(iStrt)=zcr
               iStrt=iStrt+1
               Call Get_F1(3,ccr)
C              Write (6,*) 'ccr=',ccr
               Shells(iShll)%Exp(iStrt)=ccr
               iStrt=iStrt+1
            End Do
*
         End Do
         Do iPP = 1, mPP(2)
            iShll = iShll + 1
C           Write (6,*) 'iPP,dbsc(nCnttp)%nPP=',iPP,dbsc(nCnttp)%nPP
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
            Shells(iShll)%nExp=3*kcr
            Call mma_allocate(Shells(iShll)%Exp,3*kcr,Label='Exp')
*
            iStrt=1
            Do jcr = 1, kcr
               Line=Get_Ln(lUnit)
*
               Call Get_I1(1,ncr)
C              Write (6,*) 'ncr=',ncr
               ncr=ncr+1000
               Shells(iShll)%Exp(iStrt)=DBLE(ncr)
               iStrt=iStrt+1
               Call Get_F1(2,zcr)
C              Write (6,*) 'zcr=',zcr
               Shells(iShll)%Exp(iStrt)=zcr
               iStrt=iStrt+1
               Call Get_F1(3,ccr)
C              Write (6,*) 'ccr=',ccr
               Shells(iShll)%Exp(iStrt)=ccr
               iStrt=iStrt+1
            End Do
*
         End Do
C        Write (6,*) 'Done'
*
         Return
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     M1 section                                                       *
*                                                                      *
*     Write (6,*) ' Reading M1'
      If (Index(Line,'M1').eq.0) Then
         Call WarningMessage(2,
     &              'ERROR: Keyword M1 expected, offending line : '
     &            //Line)
         Call Quit_OnUserError()
      Endif
      Line=Get_Ln(lUnit)
      Call Get_i1(1,nM1)
      dbsc(nCnttp)%nM1=nM1
      If (nM1.gt.0) Then
         Call mma_allocate(dbsc(nCnttp)%M1xp,nM1,Label='dbsc:M1xp')
         Call mma_allocate(dbsc(nCnttp)%M1cf,nM1,Label='dbsc:M1cf')
         Call Read_v(lUnit,dbsc(nCnttp)%M1xp(1),1,nM1,1,ierr)
         Call Read_v(lUnit,dbsc(nCnttp)%M1cf(1),1,nM1,1,ierr)
       End If
*                                                                      *
************************************************************************
*                                                                      *
*     M2 section                                                       *
*                                                                      *
*     Write (*,*) ' Reading M2'
      Line=Get_Ln(lUnit)
      If (Index(Line,'M2').eq.0) Then
         Call WarningMessage(2,
     &              'ERROR: Keyword M2 expected, offending line : '
     &            //Line)
         Call Quit_OnUserError()
      Endif
      Line=Get_Ln(lUnit)
      Call Get_i1(1,nM2)
      dbsc(nCnttp)%nM2=nM2
      If (nM2.gt.0) Then
         Call mma_allocate(dbsc(nCnttp)%M2xp,nM2,Label='dbsc:M2xp')
         Call mma_allocate(dbsc(nCnttp)%M2cf,nM2,Label='dbsc:M2cf')
         Call Read_v(lUnit,dbsc(nCnttp)%M2xp,1,nM2,1,ierr)
         Call Read_v(lUnit,dbsc(nCnttp)%M2cf,1,nM2,1,ierr)
      End If
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
      Line=Get_Ln(lUnit)
      Call Get_F1(1,dbsc(nCnttp)%CrRep)
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
      If (nProj.lt.0) Return
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
         Shells(iShll)%nExp=nPrim
         Shells(iShll)%nBasis = nCntrc
*
*------- Check if occupation number is included on the line
*
         iSS=1
         Call NxtWrd(Line,iSS,iEE)
         iSS=iEE+1
         Call NxtWrd(Line,iSS,iEE)
         iSS=iEE+1
         Call NxtWrd(Line,iSS,iEE)
         Call mma_allocate(Shells(iShll)%Occ,nCntrc,Label='Occ')
         If (iEE.gt.0) Then
            Do i = 1, nCntrc
               Call Get_i1(2+i,n_Occ)
               Shells(iShll)%Occ(i)=DBLE(n_Occ)/DBLE(n_Elec)
*              Write (*,*) 'n_Occ=',n_Occ
            End Do
         Else
            Shells(iShll)%Occ(:)=One
         End If
*
*        Read "orbital energies"
*
*        Write (6,*) ' Reading Bk'
         Call mma_allocate(Shells(iShll)%Bk,nCntrc,Label='Bk')
         Shells(iShll)%nBk=nCntrc
         If (nCntrc.gt.0) Call Read_v(lUnit,Shells(iShll)%Bk,
     &                                1,nCntrc,1,ierr)
         If (ierr.ne.0) goto 992
*
*        Read gaussian EXPonents
*
*        Write (6,*) ' Reading Exponents'
         Call mma_allocate(Shells(iShll)%Exp,nPrim,Label='Exp')
         Shells(iShll)%nExp=nPrim
         If (nPrim.gt.0) Call Read_v(lUnit,Shells(iShll)%Exp,1,nPrim,1,
     &                               ierr)
         If (ierr.ne.0) goto 992
*        Call RecPrt(' Exponents',Shells(iShll)%nExp,nPrim,1)
         Call mma_allocate(Shells(iShll)%Cff_c,nPrim,nCntrc,2,
     &                     Label='Cff_c')
         Call mma_allocate(Shells(iShll)%pCff,nPrim,nCntrc,
     &                     Label='pCff')
         Call mma_allocate(Shells(iShll)%Cff_p,nPrim,nPrim,2,
     &                     Label='Cff_p')
         Shells(iShll)%Cff_p(:,:,:)=Zero  ! dummy assign
         Shells(iShll)%nBasis=nCntrc
*
*        Read contraction coefficients
*        Observe that the matrix will have nPrim rows and
*        nCntrc columns
*
*        Write (6,*) ' Reading coefficients'
         Do iPrim = 1, nPrim
            Call Read_v(lUnit,Shells(iShll)%Cff_c(1,1,1),
     &                  iPrim,nPrim*nCntrc,nPrim,ierr)
            If(ierr.ne.0) goto 991
         End Do
*
*--------Renormalize
*
         call mma_allocate(Scrt1,nPrim**2)
         call mma_allocate(Scrt2,nPrim*nCntrc)
         Call Nrmlx(Shells(iShll)%Exp,nPrim,
     &              Shells(iShll)%Cff_c(1,1,1),nCntrc,
     &                      Scrt1,nPrim**2,
     &                      Scrt2, nPrim*nCntrc,iAng)
         call mma_deallocate(Scrt1)
         call mma_deallocate(Scrt2)
*
*        Duplicate!
*
         Shells(iShll)%Cff_c(:,:,2)=Shells(iShll)%Cff_c(:,:,1)
*
         If (.Not.UnNorm) Then
            nExpi=Shells(iShll)%nExp
            If (nExpi*Shells(iShll)%nBasis.ge.1) Then
               Call Nrmlz(Shells(iShll)%Exp,nExpi,
     &                    Shells(iShll)%Cff_c(1,1,1),
     &                    Shells(iShll)%nBasis,iAng)
            End If
         End If
         Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
*
 10   Continue
*
      Return
*                                                                      *
************************************************************************
*                                                                      *
 992  Continue
         Call WarningMessage(2,
     &      'Abend in GetBS: Error while reading the exponents')
         Call Quit_OnUserError()
 991  Continue
         Call WarningMessage(2,
     &      'Abend in GetBS: Error while reading the coefficients')
         Call Quit_OnUserError()
      End Subroutine GetECP

#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
      dummy_empty_procedure(GetECP)

#endif
