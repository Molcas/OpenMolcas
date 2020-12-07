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
      Subroutine Get_dDipM(dDipM,DipM,mInter,nInter)
************************************************************************
*                                                                      *
*     Objective: to compute the dipole moment derivative in Cartesians *
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      use Slapaf_Info, only: Cx
      Implicit Real*8 (a-h,o-z)
#include "info_slapaf.fh"
#include "stdalloc.fh"
      Real*8 dDipM(3,mInter), DipM(3)
      Logical Found
      Real*8, Allocatable:: Tmp2(:), BOld(:), TROld(:)
*
      nX=3*nsAtom
*
      Call mma_allocate(Tmp2,nX**2,Label='Tmp2')
      Call mma_allocate(BOld,nX*nInter,Label='BOld')
      Call Qpg_dArray('BMxOld',Found,nBMx)
      If (Found.and.(nBMx.eq.nX*nInter)) Then
         Call Get_dArray('BMxOld',BOld,nX*nInter)
      Else
         Call Get_dArray('BMtrx',BOld,nX*nInter)
      End If
      If (mTROld.gt.0) Then
         Call mma_allocate(TROld,nX*mTROld,Label='TROld')
         Call Qpg_dArray('TROld',Found,nTR)
         If (Found.and.(nTR.eq.nX*mTROld)) Then
            Call Get_dArray('TROld',TROld,nX*mTROld)
         Else
            Call Get_dArray('TR',TROld,nX*mTROld)
         End If
      Else
         Call mma_allocate(TROld,1,Label='TROld')
      End If
*
      Call Get_dDipM_(nX,BOld,TROld,mInter,nInter,Degen,
     &                Tmp2,dDipM,mTROld,Cx,Smmtrc,nsAtom,DipM)
*
      Call mma_deallocate(TROld)
      Call mma_deallocate(BOld)
      Call mma_deallocate(Tmp2)
*
      Return
      End
      Subroutine Get_dDipM_(nX,BMtrx,TRVec,mInter,nInter,Degen,
     &                     Tmp2,dDipM,mTR,Coor,Smmtrc,nAtom,DipM)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Logical Smmtrc(3,nAtom)
      Real*8 TRVec(nX,mTR), Degen(3,nAtom), BMtrx(nX,nInter),
     &       Tmp2(nX**2), Coor(3,nAtom), dDipM(3,nInter+mTR), DipM(3)
*
      Real*8 CM(3)
      Parameter ( thr = 1.0D-12 )
#ifdef _DEBUGPRINT_
      Call RecPrt('TRVec',' ',TRVec,nX,mTR)
      Call RecPrt('BMtrx',' ',BMtrx,nX,nInter)
      Call RecPrt('dDipM',' ',dDipM,3,nInter+mTR)
      Call RecPrt('DipM',' ',DipM,3,1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Analysis of the translational and rotational modes.
*     Compute the dipole moment derivative with respect to the
*     translational and rotational modes.
*
      Do i = 1, 3
         CM(i)=Zero
         rNorm=Zero
         Do iAtom = 1, nAtom
            rNorm=rNorm+Degen(i,iAtom)
            If (Smmtrc(i,iAtom)) Then
               CM(i) = CM(i) + Degen(i,iAtom)*Coor(i,iAtom)
            End If
         End Do
         CM(i)=CM(i)/rNorm
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Loop over the rotations and translations
*
      iTR=nInter+mTR
      Do iX = mTR, 1, -1
         Tx=Zero
         Ty=Zero
         Tz=Zero
         Rx=Zero
         Ry=Zero
         Rz=Zero
         Do i = 1, nAtom
            Tx = Tx + TRVec((i-1)*3+1,iX)*Degen(1,i)
            Ty = Ty + TRVec((i-1)*3+2,iX)*Degen(2,i)
            Tz = Tz + TRVec((i-1)*3+3,iX)*Degen(3,i)
            Rx = Rx +(TRVec((i-1)*3+2,iX)*(Coor(3,i)-CM(3)) -
     &                TRVec((i-1)*3+3,iX)*(Coor(2,i)-CM(2)))
     &              * Degen(1,i)
            Ry = Ry +(TRVec((i-1)*3+3,iX)*(Coor(1,i)-CM(1)) -
     &                TRVec((i-1)*3+1,iX)*(Coor(3,i)-CM(3)))
     &              * Degen(2,i)
            Rz = Rz +(TRVec((i-1)*3+1,iX)*(Coor(2,i)-CM(2)) -
     &                TRVec((i-1)*3+2,iX)*(Coor(1,i)-CM(1)))
     &              * Degen(3,i)
         End Do
#ifdef _DEBUGPRINT_
         Write (6,*) 'Tx,Ty,Tz=',Tx,Ty,Tz
         Write (6,*) 'Rx,Ry,Rz=',Rx,Ry,Rz
#endif
*
         If (Rx**2+Ry**2+Rz**2.lt.thr .and.
     &       Tx**2+Ty**2+Tz**2.gt.thr) Then
*
*-----------Translation, dipole moment invariant to translation
*
*
#ifdef _DEBUGPRINT_
            Write (6,*) 'Translation'
#endif
            call dcopy_(3,[Zero],0,dDipM(1,iTR),1)
*
         Else If (Tx**2+Ty**2+Tz**2.lt.thr .and.
     &            Rx**2+Ry**2+Rz**2.gt.thr) Then
            rNorm=(Rx**2+Ry**2+Rz**2)
*
*           Rotation, dipole moment variant to rotation
*
            If (rNorm.gt.thr) Then
*
*              General axis
*
               dDipM(1,iTR)= (DipM(2)*Rz - DipM(3)*Ry)/rNorm
               dDipM(2,iTR)= (DipM(3)*Rx - DipM(1)*Rz)/rNorm
               dDipM(3,iTR)= (DipM(1)*Ry - DipM(2)*Rx)/rNorm
            Else
               Call WarningMessage(2,' GF: too small rNorm!')
               Call Abend()
*
            End If
         End If
         iTR = iTR - 1
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('dDipM(Original)',' ',dDipM,3,nInter+mTR)
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
*     Now backtransform to cartesian coordinates.
*
*     dmu/dx = dq/dx  dmu/dq
*
      Do ix = 1, 3
*
         jj = 0
         jjj = 0
         Do jAtom = 1, nAtom
            Do jx = 1, 3
               jjj = jjj + 1
               If (Smmtrc(jx,jAtom)) Then
                  jj = jj + 1
                  ij=(jj-1)*3+ix
*
                  tmp_ij=0.0D0
                  Do k = 1, nInter
                     tmp_ij = tmp_ij + dDipM(ix,k) * BMtrx(jjj,k)
                  End Do
                  Do k = 1, mTR
                     tmp_ij = tmp_ij + dDipM(ix,nInter+k) * TRVec(jjj,k)
                  End Do
                  Tmp2(ij)=tmp_ij
*
               End If
            End Do
         End Do
*
      End Do
      call dcopy_(3*mInter,Tmp2,1,dDipM,1)
#ifdef _DEBUGPRINT_
      Call RecPrt('dDipM(cartesian)',' ',dDipM,3,mInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
