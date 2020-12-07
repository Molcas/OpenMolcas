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
      Subroutine GF(nX,mInter,nInter,Tmp1,Tmp2,EVec,EVal,RedM,
     &              iNeg,dDipM,mTR,Smmtrc,nAtom,DipM)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Logical Smmtrc(3,nAtom)
      Real*8 dDipM(3,nInter+mTR), DipM(3), Tmp1(nX**2), Tmp2(nX**2),
     &       EVec(2*mInter,mInter), EVal(2*mInter), RedM(mInter)
      Real*8, Allocatable:: G(:), GInv(:), F(:)
*                                                                      *
************************************************************************
*                                                                      *
      iRout=138
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt('GF: dDipM',' ',dDipM,3,nInter)
      Call RecPrt('GF: DipM',' ',DipM,3,1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Note that all calculations will be done in the cartesian basis!  *
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute harmonic frequencies
*                                                                      *
************************************************************************
*                                                                      *
*     Generate the G matrix (mass tensor in cartesians)
*
      Call mma_allocate(G,nX**2,Label='G')
      Call mma_allocate(GInv,nX**2,Label='GInv')
      Call Mk_G(G,GInv,nX)
*
*
*     Get the force constant matrix in cartesians
*
      Call mma_allocate(F,nX**2,Label='F')
      Call Get_H(F,nX)
*                                                                      *
************************************************************************
*                                                                      *
*     Form the GF-matrix (actually G^(1/2)FG^(1/2))
*
      Call GF_Mult(G,F,Tmp2,mInter)  ! Result in Tmp2
      Call mma_deallocate(F)
*
*     Compute the frequencies and harmonic eigenfunctions in
*     Cartesians.
*
      Call GF_Harmonic_Frequencies(G,GInv,Tmp1,Tmp2,
     &                             EVec,EVal,RedM,iNeg,nX,mInter)
*
      Call mma_deallocate(G)
      Call mma_deallocate(GInv)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Compute the dipole moment derivative in Cartesians.
*
      Call Get_dDipM(dDipM,DipM,mInter,nInter)
*                                                                      *
************************************************************************
*                                                                      *
*    Transform from cartesian to normal coordinates
*
      Do iNC = 1, mInter
         call dcopy_(mInter,EVec(1,iNC),2,Tmp2,1)
         ix = (iNC-1)*3 + 1
         iy = (iNC-1)*3 + 2
         iz = (iNC-1)*3 + 3
         Tmp1(ix) = Zero
         Tmp1(iy) = Zero
         Tmp1(iz) = Zero
         i = 0
         Do iAtom = 1, nAtom
            Do ixyz = 1, 3
               If (Smmtrc(ixyz,iAtom)) Then
                  i = i + 1
                  Tmp1(ix) = Tmp1(ix) + dDipM(1,i)*Tmp2(i)
                  Tmp1(iy) = Tmp1(iy) + dDipM(2,i)*Tmp2(i)
                  Tmp1(iz) = Tmp1(iz) + dDipM(3,i)*Tmp2(i)
               End If
            End Do
         End Do
      End Do
      call dcopy_(3*mInter,Tmp1,1,dDipM,1)
#ifdef _DEBUGPRINT_
      Call RecPrt('dDipM(normal coord.)',' ',dDipM,3,mInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
