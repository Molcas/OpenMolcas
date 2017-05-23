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

      SUBROUTINE CHO_eval_twxy(irc,ipScr,ipChoV,ipInt,nAorb,
     &                         JSYM,NUMV,DoReord)

      Implicit Real*8 (a-h,o-z)
      Integer irc,ipChoV(*),ipInt,nAorb(*),JSYM,NUMV,iAorb(8)
      Integer ipScr(8,8)
      Logical DoReord

      parameter (zero = 0.0D0, one = 1.0D0)

#include "WrkSpc.fh"
#include "cholesky.fh"
#include "choorb.fh"

C ************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
C ************************************************
      iTri(i,j) = Max(i,j)*(Max(i,j)-3)/2 + i + j
C ************************************************


      If (NumV .lt. 1) Return

C --- Computing the integrals (TT|TT),(TW,TW) and (TW|XY)
C ---------------------------------------------------------
C --- (tw|xy)  <-  (tw|xy)  +  sum_J  L(tw,#J) * L(xy,#J)
C==========================================================

      Do iSymy=1,nSym

         iSymx=MulD2h(iSymy,JSYM)

         Nxy  = nAorb(iSymx)*nAorb(iSymy)
     &        + Min(0,JSYM-2)*nAorb(iSymx)*(nAorb(iSymy)-1)/2

         If (iSymx.ge.iSymy.and.Nxy.gt.0) then

            Do iSymw=iSymy,nSym   ! iSymw.ge.iSymy

               iSymt=MulD2h(iSymw,JSYM)

               Ntw  = nAorb(iSymt)*nAorb(iSymw)
     &              + Min(0,JSYM-2)*nAorb(iSymt)*(nAorb(iSymw)-1)/2

               If (iSymt.ge.iSymw.and.Ntw.gt.0) then

                  CALL DGEMM_('N','T',Ntw,Nxy,NumV,
     &                       ONE,Work(ipChoV(iSymw)),Ntw,
     &                       WORK(ipChoV(iSymy)),Nxy,ONE,
     &                       Work(ipScr(iSymw,iSymy)),Ntw)


               End If

            End Do

         End If

      End Do


C --- Reorder to the LT-storage as required by the CI-routines
C ---
C --- LT-storage defined by  iTri(iTri(t,w),iTri(x,y))
C --- but done through explicit loops
C ------------------------------------------------------------
      IF (DoReord) THEN

         iAorb(1)= 0
         Do iSym = 2,nSym
            iAorb(iSym) = iAorb(iSym-1) + nAorb(iSym-1)
         End Do

         If (JSYM.eq.1) Then

            Do iSymy=1,nSym   ! iSymx=iSymy

               Nxy = nAorb(iSymy)*(nAorb(iSymy)+1)/2

               If (Nxy.gt.0) then

                  Do iSymw=iSymy,nSym  ! iSymt=iSymw

                     Ntw = nAorb(iSymw)*(nAorb(iSymw)+1)/2

                     If (Ntw.gt.0) then

                        Do iy=1,nAorb(iSymy)
                         iyG = iAorb(iSymy) + iy  !global index
                         Do ix=iy,nAorb(iSymy)
                          ixG  = iAorb(iSymy) + ix
                          ixy  = ix*(ix-1)/2 + iy
                          ixyG = ixG*(ixG-1)/2 + iyG ! global index
                          Do iw=1,nAorb(iSymw)
                           iwG  = iAorb(iSymw) + iw
                           Do it=iw,nAorb(iSymw)

                              itG  = iAorb(iSymw) + it
                              itw  = it*(it-1)/2 + iw
                              itwG = itG*(itG-1)/2 + iwG  !global index

                              iScr = ipScr(iSymw,iSymy) + Ntw*(ixy-1)
     &                             + itw - 1

                              iRes = ipInt + iTri(itwG,ixyG) - 1

                              Work(iRes) = Work(iScr)

                           End Do
                          End Do
                         End Do
                        End Do

                     End If

                  End Do

               End If

            End Do

         Else  ! Jsym.ne.1

            Do iSymy=1,nSym

               iSymx=MulD2h(iSymy,JSYM)
               Nxy = nAorb(iSymx)*nAorb(iSymy)

               If (iSymx.gt.iSymy.and.Nxy.gt.0) then

                  Do iSymw=iSymy,nSym  ! iSymw.ge.iSymy

                     iSymt=MulD2h(iSymw,JSYM)
                     Ntw = nAorb(iSymt)*nAorb(iSymw)

                     If (iSymt.gt.iSymw.and.Ntw.gt.0) then

                        Do iy=1,nAorb(iSymy)
                         iyG = iAorb(iSymy) + iy  !global index
                         Do ix=1,nAorb(iSymx)
                          ixG  = iAorb(iSymx) + ix
                          ixy  = nAorb(iSymx)*(iy-1) + ix
                          ixyG = ixG*(ixG-1)/2 + iyG ! global index
                          Do iw=1,nAorb(iSymw)
                           iwG  = iAorb(iSymw) + iw
                           Do it=1,nAorb(iSymt)

                              itG  = iAorb(iSymt) + it
                              itw  = nAorb(iSymt)*(iw-1) + it
                              itwG = itG*(itG-1)/2 + iwG  !global index

                              iScr = ipScr(iSymw,iSymy) + Ntw*(ixy-1)
     &                             + itw - 1

                              iRes = ipInt + iTri(itwG,ixyG) - 1

                              Work(iRes) = Work(iScr)

                           End Do
                          End Do
                         End Do
                        End Do

                     End If

                  End Do

               End If

            End Do


         EndIf

      ENDIF



      irc=0

      Return
      END

**************************************************************
