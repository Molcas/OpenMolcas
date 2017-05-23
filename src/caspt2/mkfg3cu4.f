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
* Copyright (C) 2014, Naoki Nakatani                                   *
************************************************************************
#ifdef _ENABLE_BLOCK_DMRG_
      Subroutine MKFG3CU4(IFF,G1,F1,G2,F2,G3,F3,idxG3,W3)
*
* Load 1-el, 2-el, and 3-el density matrices to resp. G1, G2, and G3
* and compute 1-el to 4-el contractions of Fock operator F1, F2, and F3.
* For F3, use cumulant reconstruction from 1-el, 2-el, and 3-el density matrices.
*
* Written by N. Nakatani, Oct. 2014
*
      IMPLICIT NONE
*
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"
*
      INTEGER, INTENT(IN) :: IFF
      REAL*8, INTENT(OUT) :: G1(NLEV,NLEV),G2(NLEV,NLEV,NLEV,NLEV)
      REAL*8, INTENT(OUT) :: F1(NLEV,NLEV),F2(NLEV,NLEV,NLEV,NLEV)
      REAL*8, INTENT(OUT) :: G3(*), F3(*)
      INTEGER*1, INTENT(IN) :: idxG3(6,*)
      REAL*8, INTENT(IN) :: W3(NLEV,NLEV,NLEV,NLEV)
*
      REAL*8  G1SUM
      INTEGER IT,IU,IV,IX,IY,IZ,IW
      INTEGER JT,JU,JV,JX,JY,JZ
      INTEGER IZSYM,IYZSYM,IXYZSYM,IVXYZSYM
      INTEGER IG3

      REAL*8, EXTERNAL :: CU4F3H
*
*
      If(NACTEL.GT.1) Then
* load 2-el density matrix
        Call block_load2pdm(nlev,G2,jstate,jstate)
* compute 1-el density matrix from 2-el density matrix
        Do iu=1,nlev
          Do it=1,nlev
            G1sum=0.0D0
            If(ism(it).EQ.ism(iu)) Then
              Do iw=1,nlev
                G1sum=G1sum+G2(iw,iw,it,iu)
              End Do
              G1(it,iu)=G1sum/(NACTEL-1)
            End If
          End Do
        End Do
      Else
* special case for NACTEL = 1
        Call block_load1pdm(nlev,G1,jstate,jstate)
      End If
*
      Do iz=1,nlev
        izSym=ism(iz)
        Do iy=1,nlev
          iyzSym=Mul(ism(iy),izSym)
          If(IFF.NE.0.AND.iyzSym.EQ.1) Then
            Do iw=1,nlev
              F1(iy,iz)=F1(iy,iz)+G2(iw,iw,iy,iz)*EPSA(iw)
            End Do
          End If
        End Do
      End Do

* skip 3RDM part if NACTEL <= 2
      If(NACTEL.LE.2) GoTo 999

      Do iz=1,nlev
        izSym=ism(iz)
        Do iy=1,nlev
          iyzSym=Mul(ism(iy),izSym)
* load 3PDM of which is G3(:,:,:,:,iy,iz)
          Call block_load3pdm2f(nlev,W3,jstate,jstate,iy,iz)
          If(IFF.NE.0) Then
            Do ix=1,nlev
              ixyzSym=Mul(ism(ix),iyzSym)
              Do iv=1,nlev
                ivxyzSym=Mul(ism(iv),ixyzSym)
                If(ivxyzSym.EQ.1) Then
                  Do iw=1,nlev
                    F2(iv,ix,iy,iz)=F2(iv,ix,iy,iz)
     &                             +W3(iw,iw,iv,ix)*EPSA(iw)
                  End Do
                End If
              End Do
            End Do
          End If

          Do iG3=1,NG3
            jt=idxG3(1,iG3)
            ju=idxG3(2,iG3)
            jv=idxG3(3,iG3)
            jx=idxG3(4,iG3)
            jy=idxG3(5,iG3)
            jz=idxG3(6,iG3)
            If(iy.EQ.jy.AND.iz.EQ.jz) Then
              G3(iG3)=W3(jt,ju,jv,jx)
              If(IFF.NE.0) Then
* CU4F3 Contrib. :: + G1(lT,lT)*G3(iP,iQ,jP,jQ,kP,kQ)
                F3(iG3)=F3(iG3)+EASUM*G3(iG3)
                Do iw=1,nlev
                  F3(iG3)=F3(iG3)
* CU4F3 Contrib. :: - 0.5D0*G1(iP,lT)*G3(lT,iQ,jP,jQ,kP,kQ)
     &                   -0.5D0*G1(jt,iw)*W3(iw,ju,jv,jx)*EPSA(iw)
* CU4F3 Contrib. :: - 0.5D0*G1(lT,iQ)*G3(iP,lT,jP,jQ,kP,kQ)
     &                   -0.5D0*G1(iw,ju)*W3(jt,iw,jv,jx)*EPSA(iw)
                End Do
              End If
            End If

            If(IFF.NE.0.AND.iy.EQ.jy.AND.iz.EQ.jz) Then
              Do iw=1,nlev
                F3(iG3)=F3(iG3)
* CU4F3 Contrib. :: - 0.5D0*G1(jP,lT)*G3(lT,jQ,iP,iQ,kP,kQ)
     &                 -0.5D0*G1(jv,iw)*W3(iw,jx,jt,ju)*EPSA(iw)
* CU4F3 Contrib. :: - 0.5D0*G1(lT,jQ)*G3(jP,lT,iP,iQ,kP,kQ)
     &                 -0.5D0*G1(iw,jx)*W3(jv,iw,jt,ju)*EPSA(iw)
              End Do
            End If

            If(IFF.NE.0.AND.iy.EQ.jv.AND.iz.EQ.jx) Then
              Do iw=1,nlev
                F3(iG3)=F3(iG3)
* CU4F3 Contrib. :: - 0.5D0*G1(kP,lT)*G3(lT,kQ,iP,iQ,jP,jQ)
     &                 -0.5D0*G1(jy,iw)*W3(iw,jz,jt,ju)*EPSA(iw)
* CU4F3 Contrib. :: - 0.5D0*G1(lT,kQ)*G3(kP,lT,iP,iQ,jP,jQ)
     &                 -0.5D0*G1(iw,jz)*W3(jy,iw,jt,ju)*EPSA(iw)
              End Do
            End If
          End Do
        End Do
      End Do
*
      If(IFF.NE.0) Then
        Do iG3=1,NG3
          it=idxG3(1,iG3)
          iu=idxG3(2,iG3)
          iv=idxG3(3,iG3)
          ix=idxG3(4,iG3)
          iy=idxG3(5,iG3)
          iz=idxG3(6,iG3)
          F3(iG3)=F3(iG3)+CU4F3H(nlev,EPSA,EASUM,
     &                    G1,G2,F1,F2,it,iu,iv,ix,iy,iz)
        End Do
      End If

 999  Return
      End
#elif defined (NAGFOR)
c Some compilers do not like empty files
      Subroutine empty_MKFG3CU4()
      End
#endif
