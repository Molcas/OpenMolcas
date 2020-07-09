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
* Copyright (C) 1990, Roland Lindh                                     *
*               1995, Martin Schuetz                                   *
************************************************************************
      Subroutine PLFInd_Clmb_dtraf
     &                 (AOInt,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,
     &                  iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp,
     &                  iadSO)
************************************************************************
*  object: to sift and index the SO integrals. The SO integrals are    *
*          moved to the proper positions in symmetry blocks, defined   *
*          by a shell quadruplet for 1st quarter direct transformation *
*                                                                      *
*          the indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
*          april '90                                                   *
*          Martin G. Schuetz, Dep. of Theoretical Chemistry            *
*          University of Lund, Sweden                                  *
*          august '95                                                  *
************************************************************************
      use index_arrays, only: iShOff, nShBF
      Implicit Real*8 (A-H,O-Z)
*
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
*
      Real*8 AOInt(ijkl,iCmp,jCmp,kCmp,lCmp)
      Integer iShell(4), iAO(4), kOp(4), iAOst(4), iadSO
      Logical Shijij
      Integer indSO(4), indSOb(4), indSOo(4), nSOiSh(4),
     &        IntOff, IntPtr
      Logical usShij
*     statement functions
c     ITgOff(i,j,k,l) =
c    &             ((indSOo(j)*(indSOo(j)+1)/2+indSOo(i))*nSOiSh(k)
c    &            +indSOo(k))*nSOiSh(l)
c     ISqOff(i,j,k,l) =
c    &             ((indSOo(j)*nSOiSh(i)+indSOo(i))*nSOiSh(k)
c    &            +indSOo(k))*nSOiSh(l)
*
      Call qEnter('PLFInd')
      iRout = 109
      iPrint = nPrint(iRout)
*     and set up corresponding logicals...
      usShij = iShell(1).eq.iShell(2)
*     compute 1st index of components and # functions of
*     the shells ...unscrambled...
      Do i = 1, 4
        nSOiSh(i) = nShBF(0,iShell(i))
        indSOb(i) = iShOff(0,iShell(i))
      End Do
*
      If (usShij) Then
*       slow indices ij triangular...
        Do i1 = 1, iCmp
          Do i2 = 1, jCmp
            Do i3 = 1, kCmp
              Do i4 = 1, lCmp
*
*               Unfold the way the eight indices have been reordered.
                iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
                jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
                kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
                lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
                nijkl = 0
                Do lAOl = 0, lBas-1
                  lSOl = lSO + lAOl
                  Do kAOk = 0, kBas-1
                    kSOk = kSO + kAOk
                    Do jAOj = 0, jBas-1
                      jSOj = jSO + jAOj
                      Do iAOi = 0, iBas-1
                        nijkl = nijkl + 1
                        iSOi = iSO + iAOi
*
*                       unscramble SO function indices ...
                        indSO(4)=lSOl
                        indSO(3)=kSOk
                        indSO(2)=jSOj
                        indSO(1)=iSOi

*                       ij triangular here ...
                        idum=indSO(1)
                        indSO(1)=Min(idum,indSO(2))
                        indSO(2)=Max(idum,indSO(2))
*                       compute SO function offset relative to 1st
*                       index of component of shell in irrep ...
                        indSOo(1)=indSO(1)-indSOb(1)
                        indSOo(2)=indSO(2)-indSOb(2)
                        indSOo(3)=indSO(3)-indSOb(3)
                        indSOo(4)=indSO(4)-indSOb(4)
*
*                       compute offset and position of integral
c                       IntOff=ITgOff(1,2,3,4)
                        IntOff=
     &                    ((indSOo(2)*(indSOo(2)+1)/2+indSOo(1))
     &                     *nSOiSh(3)
     &                     +indSOo(3))*nSOiSh(4)
                        IntPtr=IntOff+indSOo(4)+1
                        Work(iadSO+IntPtr-1)=AOInt(nijkl,i1,i2,i3,i4)
*
                      End Do
                    End Do
                  End Do
                End Do
*
              End Do
            End Do
          End Do
        End Do
      Else
        Do i1 = 1, iCmp
          Do i2 = 1, jCmp
            Do i3 = 1, kCmp
              Do i4 = 1, lCmp
*
*               Unfold the way the eight indices have been reordered.
                iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
                jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
                kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
                lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
*
                nijkl = 0
                Do lAOl = 0, lBas-1
                  lSOl = lSO + lAOl
*                 unscramble SO function indices ...
*                 and compute SO function offset relative to 1st
*                 index of component of shell in irrep ...
                  indSO(4)=lSOl
                  indSOo(4)=indSO(4)-indSOb(4)
                  Do kAOk = 0, kBas-1
                    kSOk = kSO + kAOk
                    indSO(3)=kSOk
                    indSOo(3)=indSO(3)-indSOb(3)
                    Do jAOj = 0, jBas-1
                      jSOj = jSO + jAOj
                      indSO(2)=jSOj
                      indSOo(2)=indSO(2)-indSOb(2)
                      Do iAOi = 0, iBas-1
                        nijkl = nijkl + 1
                        iSOi = iSO + iAOi
                        indSO(1)=iSOi
                        indSOo(1)=indSO(1)-indSOb(1)
*
*                       compute offset and position of integral
c                       IntOff=ISqOff(1,2,3,4)
                        IntOff=
     &                    ((indSOo(2)*nSOiSh(1)+indSOo(1))
     &                     *nSOiSh(3)
     &                     +indSOo(3))*nSOiSh(4)
                        IntPtr=IntOff+indSOo(4)+1
                        Work(iadSO+IntPtr-1)=AOInt(nijkl,i1,i2,i3,i4)
*
                      End Do
                    End Do
                  End Do
                End Do
*
              End Do
            End Do
          End Do
        End Do
      End If
*
*     Call GetMem(' Exit PLFIdS','CHECK','REAL',iDum,iDum)
      Call qExit('PLFInd')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(Shijij)
      End
