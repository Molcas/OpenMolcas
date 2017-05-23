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
************************************************************************
*                                                                      *
* This routine computes the expectation value <s^2> for a UHF          *
* wavefunction.                                                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Parameters:                                                          *
* Ca     -  Alpha orbitals, input.                                     *
* Cb     -  Beta orbitals, input.                                      *
* S      -  Overlap matrix, input.                                     *
* nAlpha -  Number of alpha orbitals occupied per irrep, input.        *
* nBeta  -  Number of beta orbitals occupied per irrep, input.         *
* nBas   -  Number of basis functions per irrep, input.                *
* nOrb   -  Number of orbitals per irrep, input.                       *
* nSym   -  Number of irreps, input.                                   *
* s2     -  <s^2>, output.                                             *
*                                                                      *
************************************************************************
      Subroutine s2calc(Ca,Cb,S,nAlpha,nBeta,nBas,nOrb,nSym,s2)
      Implicit None
#include "stdalloc.fh"
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Real*8  Ca(*)
      Real*8  Cb(*)
      Real*8  S(*)
      Integer nAlpha(*)
      Integer nBeta(*)
      Integer nBas(*)
      Integer nOrb(*)
      Integer nSym
      Real*8  s2
*----------------------------------------------------------------------*
* Local variables.                                                     *
*----------------------------------------------------------------------*
      Real*8  sz
      Real*8  sb
      real*8  so
      Integer iSym
      Integer npSmat
      Integer npHalf
      Integer npTfrm
      Integer idxCMO
      Integer idxOvl
      Integer i
      Real*8, Dimension(:), Allocatable :: Smat, Half, Tfrm
*----------------------------------------------------------------------*
* Debug printing stuff                                                 *
*----------------------------------------------------------------------*
*define _DEBUG_
#ifdef _DEBUG_
      Write(6,'(a,8i5)') 'nSym  . . . . . .',nSym
      Write(6,'(a,8i5)') 'nBas  . . . . . .',(nBas(iSym),iSym=1,nSym)
      Write(6,'(a,8i5)') 'nOrb  . . . . . .',(nOrb(iSym),iSym=1,nSym)
      Write(6,'(a,8i5)')
     &   'nAlpha  . . . . .',(nAlpha(iSym),iSym=1,nSym)
      Write(6,'(a,8i5)')
     &   'nBeta . . . . . .',(nBeta(iSym),iSym=1,nSym)
      i = 1
      Do iSym = 1, nSym
         Call RecPrt('Ca',' ',Ca(i),nBas(iSym),nBas(iSym))
         Call RecPrt('Cb',' ',Cb(i),nBas(iSym),nBas(iSym))
         Call RecPrt('S ',' ',S (i),nBas(iSym),nBas(iSym))
         i = i + nBas(iSym)**2
      End Do
#endif
*----------------------------------------------------------------------*
* Initialize                                                           *
*----------------------------------------------------------------------*
      s2=0.0d0
*----------------------------------------------------------------------*
* Compute <sz> and N_beta                                              *
*----------------------------------------------------------------------*
      sz=0.0d0
      sb=0.0d0
      Do iSym=1,nSym
         sz=sz+0.5d0*(nAlpha(iSym)-nBeta(iSym))
         sb=sb+1.0d0*nBeta(iSym)
      End Do
      s2=s2 + sz*(sz+1.0d0) + sb
#ifdef _DEBUG_
      Write(6,'(a,f12.6)') 'sz  . . . . . . .',sz
      Write(6,'(a,f12.6)') 'sb  . . . . . . .',sb
      Write(6,'(a,f12.6)') 's2  . . . . . . .',s2
#endif
*----------------------------------------------------------------------*
* Compute size of scratch matrices, the allocate.                      *
*----------------------------------------------------------------------*
      npSmat=0
      npHalf=0
      npTfrm=0
      Do iSym=1,nSym
         npSmat=Max(npSmat,nBas(iSym)*nBas(iSym))
         npHalf=Max(npHalf,nAlpha(iSym)*nBas(iSym))
         npTfrm=Max(npTfrm,nAlpha(iSym)*nBeta(iSym))
      End Do
#ifdef _DEBUG_
      Write(6,'(a,i5)') 'npSmat  . . . . .',npSmat
      Write(6,'(a,i5)') 'npHalf  . . . . .',npHalf
      Write(6,'(a,i5)') 'npTfrm  . . . . .',npTfrm
#endif
      If(npTfrm.eq.0) Return
      Call mma_allocate(Smat,npSmat)
      Call mma_allocate(Half,npHalf)
      Call mma_allocate(Tfrm,npTfrm)
*----------------------------------------------------------------------*
* Compute sum_a sum_b [ C_alpha* S C_beta ]^2                          *
*----------------------------------------------------------------------*
*     Call xxdGemul(a,lda,'N/T', b,ldb,'N/T', c,ldc, l,m,n)
*     C=A*B, m is contraction index
      so=0.0d0
      idxCMO=1
      idxOvl=1
      Do iSym=1,nSym
*        Write(6,'(a,i1)') 'Irrep ',iSym
         If(nAlpha(iSym)*nBeta(iSym).gt.0) Then
            Call Square(S(idxOvl),Smat,1,nBas(iSym),nBas(iSym))
            Call DGEMM_('T','N',
     &                  nAlpha(iSym),nBas(iSym),nBas(iSym),
     &                  1.0d0,Ca(idxCMO),nBas(iSym),
     &                        Smat,nBas(iSym),
     &                  0.0d0,Half,nAlpha(iSYm))
*           Call RecPrt('Half transform','(12f12.6)',
*    &                  Half,nAlpha(iSym),nBas(iSym))
            Call DGEMM_('N','N',
     &                  nAlpha(iSym),nBeta(iSym),nBas(iSym),
     &                  1.0d0,Half,nAlpha(iSym),
     &                        Cb(idxCMO),nBas(iSym),
     &                  0.0d0,Tfrm,nAlpha(iSym))
*           Call RecPrt('Transform','(12f12.6)',
*    &                  Tfrm,nAlpha(iSym),nBeta(iSym))
            Do i=1,nAlpha(iSym)*nBeta(iSym)
               so=so+Tfrm(i)*Tfrm(i)
            End Do
         End If
         idxCMO=idxCMO+nBas(iSym)*nOrb(iSym)
         idxOvl=idxOvl+nBas(iSym)*(nBas(iSym)+1)/2
      End Do
      s2=s2 - so
*     Write(6,'(a,f12.6)') 'so  . . . . . . .',so
*     Write(6,'(a,f12.6)') 's2  . . . . . . .',s2
*----------------------------------------------------------------------*
* Deallocate scratch matrix.                                           *
*----------------------------------------------------------------------*
      Call mma_deallocate(Tfrm)
      Call mma_deallocate(Half)
      Call mma_deallocate(Smat)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
