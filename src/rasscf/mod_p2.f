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
      Subroutine Mod_P2(P2mo,nP2Act,D1mo,nD1mo,DS1mo,ExFac,nDet)
      Implicit Real*8 (A-H,O-Z)
#include "nq_info.fh"
#include "real.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='MOD_P2  ')
      Real*8 P2mo(nP2Act),D1mo(nD1mo), DS1mo(nD1mo)
      Integer iOff_Ash(0:7)
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      iOff_=0
      Do iIrrep = 0, mIrrep-1
         iOff_Ash(iIrrep)=iOff_
         iOff_=iOff_+nAsh(iIrrep)
      End Do
*
************************************************************************
*
*---   active space Ptvxy
*
      IF (nDet.eq.1) THEN
*
      P2Act=DBLE(nP2Act)
      Call Put_Temp('nP2Act  ',P2Act,1)
      Call Put_Temp('P2_RAW  ',P2mo,nP2Act)
*
      Do iIrrep = 0, mIrrep-1
         Do jIrrep = 0, mIrrep-1
         ijIrrep=iEor(iIrrep,jIrrep)
            Do kIrrep = 0, mIrrep-1
            ijkIrrep=iEor(ijIrrep,kIrrep)
*
               Do k_ = 1, nASh(kIrrep)
                  k=iOff_Ash(kIrrep)+k_
                  Do l_ = 1, nASh(ijkIrrep)
                     l=iOff_Ash(ijkIrrep)+l_
                     If (l.gt.k) Go To 100
                     kl   = iTri(k,l)
                     Do i_ = 1, nASh(iIrrep)
                        i=iOff_Ash(iIrrep)+i_
                        il   = iTri(i,l)
                        ik   = iTri(i,k)
                        Do j_ = 1, nASh(jIrrep)
                           j=iOff_Ash(jIrrep)+j_
                           If (j.gt.i) Go To 200
                           ij   = iTri(i,j)
                           If (kl.gt.ij) Go To 200
                           ijkl = iTri(ij,kl)
                           jk   = iTri(j,k)
                           jl   = iTri(j,l)
*
                           Fact=One
                           If (iIrrep.eq.jIrrep) Then
                              If (k.eq.l) Fact=Two
                           End If
                           P2mo(ijkl) = Fact * P2mo(ijkl)
*
                           If (iIrrep.eq.ijkIrrep) Then
                              P2mo(ijkl) = P2mo(ijkl)
     &                                   + (1.0d0-ExFac)*
     &                                   (0.25d0* D1mo(jk)* D1mo(il)+
     &                                    0.25d0*DS1mo(jk)*DS1mo(il))
                           End If
*
                           If (iIrrep.eq.kIrrep) Then
                              P2mo(ijkl) = P2mo(ijkl)
     &                                   + (1.0d0-ExFac)*
     &                                   (0.25d0* D1mo(jl)* D1mo(ik)+
     &                                    0.25d0*DS1mo(jl)*DS1mo(ik))
                           End If
*
                           P2mo(ijkl) = P2mo(ijkl)/Fact
*
 200                       Continue
                        End Do
                     End Do
 100                 Continue
                  End Do
               End Do
*
            End Do
         End Do
      End Do
      Call Put_Temp('P2_KS   ',P2mo,nP2Act)
      ELSE
         Write(LF,*) ' Not implemented yet!!! nDet=',nDet
         Call Abend()
      END IF
*
      RETURN
      END
