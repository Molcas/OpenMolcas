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
* Copyright (C) 1995, Martin Schuetz                                   *
************************************************************************
      SubRoutine ExpKap(kapOV,U,mynOcc)
************************************************************************
*                                                                      *
*     purpose: U = exp(kappa)                                          *
*                                                                      *
*     method:  Either via ETaylr for larger orbital rotations, or via  *
*              ESchur for smaller orbital rotations.                   *
*                                                                      *
*     output:                                                          *
*       U       : unitary matrix to transform old CMOs                 *
*                                                                      *
*     called from: RotMOS                                              *
*                                                                      *
*     calls to: ?????                                                  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1995                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
*
*     declaration subroutine parameters
      Real*8 kapOV(nOV),U(nOFS)
      Integer mynOcc(8)
*
      Parameter (Thresh = 5.0D-3)
c     Parameter (Thresh = 0.0D0)
c     Parameter (Thresh = 99999.99D99)
*
      Real*8 Cpu1,Tim1,Tim2,Tim3
*
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
*
      ddd=abs(kapOV(IDAMAX_(nOV,kapOV,1)))
      If (ddd.lt.Thresh) Then
         Call ESchur(kapOV,U,mynOcc)
      Else
         Call ETaylr(kapOV,U,mynOcc)
      End If
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call NrmClc(kapOV,nOV,'ExpKap','kapOV')
#endif
*
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld(10) = TimFld(10) + (Cpu2 - Cpu1)
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      SubRoutine ETaylr(kapOV,U,mynOcc)
************************************************************************
*                                                                      *
*     purpose: U = exp(kappa)                                          *
*                                                                      *
*     method:  The series in x(t)x (where x is the occ/virt part of    *
*              kappa) are interpreted as parts of Taylor expansions of *
*              cosines and sines.                                      *
*              cost: diagonalization of one (nocc)x(nocc) matrix       *
*              exakt method, no adding up of matrix products necessary *
*                                                                      *
*     output:                                                          *
*       U       : unitary matrix to transform old CMOs                 *
*                                                                      *
*     called from: ExpKap                                              *
*                                                                      *
*     calls to: ?????                                                  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1995                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*
*     declaration subroutine parameters
      Real*8 kapOV(nOV),U(nOFS)
      Integer mynOcc(8)
*
*     declaration local variables
      Integer iSym,lAi,nOrbmF,nVrt,nOccmF,ii,ioffs,iAof,iUof,iptr
      Real*8 quot,dlambd,dlamb2,Uii,A2ii,A3ii
      Real*8 Cpu1,Tim1,Tim2,Tim3
      Real*8, Dimension(:), Allocatable:: UNew, A1, A2, A3, Usm, Dg,
     &                                   Scratch
*
      Parameter (THLDEP = 1.0D-14)
*
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
*
* Dummy initialize (prevent stupid compiler complaints):
      A3ii=0.0D0
      A2ii=0.0D0
      Uii=0.0D0
      ioffs=1
      ivoffs=1
      Do iSym=1,nSym
* initialize some pointers and sizes of intermediate matrices
        nOrbmF=nOrb(iSym)-nFro(iSym)
        If (nOrbmF.lt.1) GoTo 200
        nVrt=nOrb(iSym)-mynOcc(iSym)
        nOccmF=mynOcc(iSym)-nFro(iSym)
        If ((nVrt.lt.1).OR.(nOccmF.lt.1)) GoTo 200
        lAi=nOccmF*nOccmF
* allocate and init memory for U=exp(kappa) matrix
        Call mma_allocate(UNew,nOrbmF*nOrbmF,Label='UNew')
        Call FZero(UNew,nOrbmF**2)
* allocate space for intermediate matrices A,A2,A3,Usm
        Call mma_allocate(A2,lAi,Label='A2')
        Call mma_allocate(A3,lAi,Label='A3')
        Call mma_allocate(Usm,lAi,Label='Usm')
        Call mma_allocate(A1,lAi,Label='A1')
        Call mma_allocate(Dg,nOccmF,Label='Dg')
        call dcopy_(lAi,[Zero],0,A2,1)
        call dcopy_(lAi,[Zero],0,A3,1)
        call dcopy_(lAi,[Zero],0,Usm,1)
        call dcopy_(lAi,[Zero],0,A1,1)
        call dcopy_(nOccmF,[One],0,Usm,nOccmF+1)
*       compute A=x(T)x and diagonalize it
        Call MxMt(kapOV(ivoffs),nVrt,1,kapOV(ivoffs),1,nVrt,
     &            A1,nOccmF,nVrt)
*
        Call mma_allocate(Scratch,nOccmF**2,Label='Scratch')
*
        Dummy=0.0D0
        iDum=0
        Call Diag_Driver('V','A','L',nOccmF,A1,Scratch,
     &                   nOccmF,Dummy,Dummy,iDum,iDum,Dg,
     &                   Usm,nOccmF,1,0,'J',nFound,iErr)
*
        Call mma_deallocate(Scratch)
*-- compute matrices A1,A2,A3 using eigenvalues of A and
*-- backtransformation with eigenvectors U thereafter
        Do ii=0,nOccmF-1
          iAof=ii*(nOccmF+1)
          iUof=ii*(nOrbmF+1)
          dlamb2=Dg(1+ii)
C---Replaced code sequence, PAM Apr 2002
C          If (abs(dlamb2).lt.THLDEP) Then
C*           apply Cauchy criterium for limiting case...
C*           A1 goes directly into U=exp(kappa) matrix
C            UNew(1+iUof)=One
C            A2(1+iAof)=One
C            A3(1+iAof)=-Half
C          Else If (dlamb2.gt.Zero) Then
C            dlambd=sqrt(dlamb2)
C            quot=DCOS(dlambd)
C*           A1 goes directly into U=exp(kappa) matrix
C            UNew(1+iUof)=quot
C            A2(1+iAof)=DSIN(dlambd)/dlambd
C            quot=quot/dlamb2
C            quot=quot-(One-Half*dlamb2)/dlamb2
C            A3(1+iAof)=quot
C          Else
C*           large, negative eigenvalue for x(T)x,
C*           cannot proceede
C            Write (*,*) 'ExpKap: large, negative eigenvalue for x(T)x'
C            Write (*,*) 'dlamb2=',dlamb2
C            Call Abend()
C          End If
C---Replaced by:-----------

CPAM05: Deactivate the following line if small imaginary rotations arising
C from numerical noise should be allowed (avoiding discontinuous derivatives.)
          If(dlamb2.lt.0.0D0) dlamb2=0.0d00

          If(dlamb2.lt.-THLDEP) then
            Write (6,*) 'ExpKap: negative eigenvalue for x(T)x'
            Write (6,*) 'dlamb2=',dlamb2
            Call Abend()
          Else If(dlamb2.lt.0.16d0) then
            A3ii=dlamb2*(151200-dlamb2*(5040-dlamb2*(90-dlamb2)))/
     &                                                        3628800
            A2ii=1-(dlamb2*(60480-dlamb2*(3024-
     &                                    dlamb2*(72-dlamb2))))/362880
            Uii=1.0D0+dlamb2*(A3ii-0.5D0)
          Else
            dlambd=SQRT(dlamb2)
            quot=COS(dlambd)
            Uii=quot
            A2ii=SIN(dlambd)/dlambd
            A3ii=(quot-(One-Half*dlamb2))/dlamb2
          End If
          UNew(1+iUof)=Uii
          A2(1+iAof)=A2ii
          A3(1+iAof)=A3ii
C---End of replacement, PAM Apr 2002
        End Do
* Diagonal is no longer used...
        Call mma_deallocate(Dg)
* Here comes backtransformation. Eigenvectors in Usm.
        Call DGEMM_('N','N',
     &              nOccmF,nOccmF,nOccmF,
     &              1.0d0,Usm,nOccmF,
     &                    UNew,nOrbmF,
     &              0.0d0,A1,nOccmF)
        Call DGEMM_('N','T',
     &              nOccmF,nOccmF,nOccmF,
     &              1.0d0,A1,nOccmF,
     &                    Usm,nOccmF,
     &              0.0d0,UNew,nOrbmF)
        Call DGEMM_('N','N',
     &              nOccmF,nOccmF,nOccmF,
     &              1.0d0,Usm,nOccmF,
     &                    A2,nOccmF,
     &              0.0d0,A1,nOccmF)
        Call DGEMM_('N','T',
     &              nOccmF,nOccmF,nOccmF,
     &              1.0d0,A1,nOccmF,
     &                    Usm,nOccmF,
     &              0.0d0,A2,nOccmF)
        Call DGEMM_('N','N',
     &              nOccmF,nOccmF,nOccmF,
     &              1.0d0,Usm,nOccmF,
     &                    A3,nOccmF,
     &              0.0d0,A1,nOccmF)
        Call DGEMM_('N','T',
     &              nOccmF,nOccmF,nOccmF,
     &              1.0d0,A1,nOccmF,
     &                    Usm,nOccmF,
     &              0.0d0,A3,nOccmF)
* dispose this stuff...
        Call mma_deallocate(A1)
        Call mma_deallocate(Usm)
*-- assemble rest of U=exp(kappa) matrix:
*--1) occ/virt (OV) Block as xA2 (nVrt x nOccmF):
        iptr=1+nOccmF
* iptr points to appropriate position within U=exp(kappa)
        Call DGEMM_('N','N',
     &              nVrt,nOccmF,nOccmF,
     &              1.0d0,kapOV(ivoffs),nVrt,
     &                    A2,nOccmF,
     &              0.0d0,UNew(iptr),nOrbmF)
* dispose A2
        Call mma_deallocate(A2)
*--2) virt/virt (VV) Block as 1-0.5xx(T)+xA3x(T)
*       initiate VV block in U=exp(kappa) as identity matrix
        iptr=1+nOrbmF*nOccmF+nOccmF
        call dcopy_(nVrt,[One],0,UNew(iptr),nOrbmF+1)
* compute VV block in U=exp(kappa) as 1-0.5xx(T)+xA3x(T)
        Call DGEMM_('N','T',nVrt,nVrt,nOccmF,
     &             -Half,kapOV(ivoffs),nVrt,
     &                   kapOV(ivoffs),nVrt,
     &               One,UNew(iptr),nOrbmF)
        Call DGEMM_('N','T',
     &              nOccmF,nVrt,nOccmF,
     &              1.0d0,A3,nOccmF,
     &                    kapOV(ivoffs),nVrt,
     &              0.0d0,UNew(iptr-nOccmF),nOrbmF)
        Call DGEMM_('N','N',nVrt,nVrt,nOccmF,
     &              One,kapOV(ivoffs),nVrt,
     &                  UNew(iptr-nOccmF),nOrbmF,
     &              One,UNew(iptr),nOrbmF)
* dispose stuff...
        Call mma_deallocate(A3)
*--3) virt/occ (VO) Block as -A2x(T)
* compute VO block as -A2x(T) from OV block by
* just transpose xA2 and scale with -1.0
        iptr=iptr-nOccmF
        Call DGETMO(UNew(1+nOccmF),nOrbmF,nVrt,nOccmF,
     &              UNew(iptr),nOrbmF)
cVV: GCC 3.0.4 had a problem with next loop..
c        Do ii=0,nOccmF-1
         ii=0
         if(nOccmF.lt.1) goto 610
600      continue
          Call DSCAL_(nVrt,-One,UNew(iptr+ii),nOrbmF)
          ii=ii+1
          if(ii.le.nOccmF-1) goto 600
c        End Do
610      Continue
*       copy symblock to U
        call dcopy_(nOrbmF*nOrbmF,UNew,1,U(ioffs),1)
        ioffs=ioffs+nOrbmF*nOrbmF
        ivoffs=ivoffs+nOccmF*nVrt
*       free memory
        Call mma_deallocate(UNew)
  200   Continue
      End Do
*
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld(11) = TimFld(11) + (Cpu2 - Cpu1)
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      SubRoutine ESchur(kapOV,U,mynOcc)
************************************************************************
*                                                                      *
*     purpose: U = exp(kappa)                                          *
*                                                                      *
*     method:  Reduction of kappa matrix to (2nocc)x(2nocc) size by    *
*              Transformation to Schur basis                           *
*                                                                      *
*     output:                                                          *
*       U       : unitary matrix to transform old CMOs                 *
*                                                                      *
*     called from: ExpKap                                              *
*                                                                      *
*     calls to: ?????                                                  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1995                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*
*     declaration subroutine parameters
      Real*8 kapOV(nOV),U(nOFS)
      Integer mynOcc(8)
*
*     declaration local variables
      Integer iSym,lkap,lkap2,lrkap,lrkap2,lrdmat,nOccmF,nOrbmF,nVrt,
     &        irmofs,iptr,iptr2,iptr3,iptr4,ivoffs,ii,ij,nzero
      Real*8 vprod,Fact
      Real*8 Cpu1,Tim1,Tim2,Tim3
      Real*8, Dimension(:), Allocatable:: rKap, rMat, rPro, rSum, rAux,
     &                                    HTra, TSum
*     declarations of functions
      real*8 ddot_
*
      Parameter (THLDEP = 1.0D-18)
*
* PAM 2008 Hardcoded EPS:
      EPS=1.0D-15
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
*
      ioffs=1
      ivoffs=1
      Do iSym=1,nSym
*       initialize some pointers and sizes of intermediate matrices
        nOrbmF=nOrb(iSym)-nFro(iSym)
        nVrt=nOrb(iSym)-mynOcc(iSym)
        nOccmF=mynOcc(iSym)-nFro(iSym)
        If ((nVrt.lt.1).OR.(nOccmF.lt.1)) GoTo 200
*       size of reduced kappa matrix kap_red
        lrkap=2*nOccmF
        lrkap2=lrkap*lrkap
        lkap=nOrbmF
        lkap2=lkap*lkap
*       size of transformation matrix V with kap = V kap_red V(T)
        lrdmat=lrkap*nOrbmF
*       offset in V
        irmofs=nOccmF*nOrbmF
*       allocate space for intermediate matrices V & kap_red
        Call mma_allocate(rKap,lrkap2,Label='rKap')
        Call mma_allocate(rMat,lrdmat,Label='rMat')
*       clear kap_red and V
        call dcopy_(lrkap2,[Zero],0,rKap,1)
        call dcopy_(lrdmat,[Zero],0,rMat,1)
*       set occ-occ part of V to identity matrix
        call dcopy_(nOccmF,[One],0,rMat,nOrbmF+1)
*       copy kap vectors from kapOV to appropriate positions in V
*       and perform Gram-Schmidt orthonormalization
        nzero=0
        iptr=1+irmofs+nOccmF
        Do ii=1,nOccmF
          call dcopy_(nVrt,kapOV(ivoffs),1,rMat(iptr),1)
          iptr2=1+irmofs+nOccmF
          Do ij=1,ii-1
            iptr3=(ii-1)*lrkap+nOccmF+ij
            iptr4=lrkap*nOccmF+(ij-1)*lrkap+ii
            vprod=DDOT_(nVrt,rMat(iptr),1,rMat(iptr2),1)
            rKap(iptr3)= vprod
            rKap(iptr4)=-vprod
            call daxpy_(nVrt,-vprod,rMat(iptr2),1,rMat(iptr),1)
            iptr2=iptr2+nOrbmF
          End Do
*         normalize, or delete in case of linear dependence
          vprod=DNORM2(nVrt,rMat(iptr),1)
          iptr3=(ii-1)*lrkap+nOccmF+ii
          iptr4=lrkap*nOccmF+(ii-1)*lrkap+ii
          if (vprod.gt.THLDEP) Then
             Call DSCAL_(nVrt,One/vprod,rMat(iptr),1)
             rKap(iptr3)= vprod
             rKap(iptr4)=-vprod
          Else
             call dcopy_(nVrt,[Zero],0,rMat(iptr),1)
             rKap(iptr3)=Zero
             rKap(iptr4)=Zero
             nzero=nzero+1
          End If
          ivoffs=ivoffs+nVrt
          iptr=iptr+nOrbmF
        End Do
*       compute exp(kappa_red)-1
        Call mma_allocate(rSum,lrkap2,Label='rSum')
        Call mma_allocate(rPro,lrkap2,Label='rPro')
        Call mma_allocate(rAux,lrkap2,Label='rAux')
*       clear kappa_red_sum
        call dcopy_(lrkap2,[Zero],0,rSum,1)
*       now sum over exp expansion (until machine precision)
        call dcopy_(lrkap2,rKap,1,rPro,1)
        Call DGEADD(rSum,lrkap,'N',rPro,lrkap,'N',
     &              rSum,lrkap,lrkap,lrkap)
        Do kexp=2,MxKp2U
          Fact=One/DBLE(kexp)
          call dcopy_(lrkap2,rPro,1,rAux,1)
          call dcopy_(lrkap2,[Zero],0,rPro,1)
          Call DGEMM_('N','N',lrkap,lrkap,lrkap,
     &                 Fact,rAux,lrkap,
     &                      rKap,lrkap,
     &                 Zero,rPro,lrkap)
          Call DGEADD(rSum,lrkap,'N',rPro,lrkap,'N',
     &                rSum,lrkap,lrkap,lrkap)
*         if addup smaller than machine precision, exit loop
          If (abs(rPro(IDAMAX_(lrkap2,rPro,1))).le.eps)
     &      GoTo 100
        End Do
  100   Continue
        If (kexp.ge.MxKp2U) Then
           Write (6,*) 'ExpKap: kexp.ge.MxKp2U'
           Write (6,*) 'kexp=',kexp
           Write (6,*) 'MxKp2U=',MxKp2U
           Call Abend()
        End If
*       we don't use these anymore...
        Call mma_deallocate(rAux)
        Call mma_deallocate(rPro)
        Call mma_deallocate(rKap)
*       but this...
        Call mma_allocate(HTra,lrdmat,Label='HTra')
*       transform back: U = exp(kap) = V exp(kap_red) V(T)
        Call DGEMM_('N','N',
     &              nOrbmF,lrkap,lrkap,
     &              1.0d0,rMat,nOrbmF,
     &                    rSum,lrkap,
     &              0.0d0,HTra,nOrbmF)
*       we don't use this anymore...
        Call mma_deallocate(rSum)
*       but this...
        Call mma_allocate(TSum,lkap2,Label='TSum')
        Call DGEMM_('N','T',
     &              nOrbmF,nOrbmF,lrkap,
     &              1.0d0,HTra,nOrbmF,
     &                    rMat,nOrbmF,
     &              0.0d0,TSum,nOrbmF)
*       we don't use this anymore...
        Call mma_deallocate(rMat)
*       Add identity matrix to exp(kappa)-1
        Do ii=1,lkap
          iptr=(ii-1)*lkap+ii
          TSum(iptr)= TSum(iptr)+One
        End Do
*       copy symblock to U
        call dcopy_(lkap2,TSum,1,U(ioffs),1)
        ioffs=ioffs+lkap2
*       free memory
        Call mma_deallocate(TSum)
        Call mma_deallocate(HTra)
  200   Continue
      End Do
*
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld(12) = TimFld(12) + (Cpu2 - Cpu1)
      Return
      End
