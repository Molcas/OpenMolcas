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

      SUBROUTINE CHO_rassi_twxy(irc,Scr,ChoV,ipInt,nAorb,JSYM,NUMV,
     &                          DoReord)

      use Data_Structures, only: Laq_type, twxy_type
      Implicit Real*8 (a-h,o-z)
      Integer irc,ipInt,nAorb(*),JSYM,NUMV,iAorb(8)
      Type (Laq_Type) ChoV
      Type (twxy_type) Scr
      Logical DoReord

#include "real.fh"
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

         If (Nxy.gt.0) then

            Do iSymw=iSymy,nSym   ! iSymw.ge.iSymy (particle symmetry)

               iSymt=MulD2h(iSymw,JSYM)

               Ntw  = nAorb(iSymt)*nAorb(iSymw)

               If (Ntw.gt.0) then

                  CALL DGEMM_('N','T',Ntw,Nxy,NumV,
     &                       ONE,ChoV%pA(iSymt)%A,Ntw,
     &                           ChoV%pA(iSymx)%A,Nxy,
     &                       One,Scr%pA(iSymw,iSymy)%A,Ntw)


               End If

            End Do

         End If

      End Do


C --- Reorder to the storage required by the RASSI program
C ---
C --- There is no permutational symmetry but only particle
C --- symmetry in the (tw|xy) integrals
C ------------------------------------------------------------
      IF (DoReord) THEN

         iAorb(1)= 0
         Do iSym = 2,nSym
            iAorb(iSym) = iAorb(iSym-1) + nAorb(iSym-1)
         End Do

         nTA = iAorb(nSym) + nAorb(nSym) ! total # active orbitals


         Do iSymy=1,nSym

            iSymx=MulD2h(iSymy,JSYM)

            Nxy = nAorb(iSymx)*nAorb(iSymy)

            If (Nxy<=0) Cycle

            Do iSymw=iSymy,nSym

               iSymt=MulD2h(iSymw,JSYM)

               Ntw = nAorb(iSymt)*nAorb(iSymw)

               If (Ntw<=0) Cycle

               Do iy=1,nAorb(iSymy)

                   iyG = iAorb(iSymy) + iy  !global index

                   Do ix=1,nAorb(iSymx)

                    ixG  = iAorb(iSymx) + ix

                    ixy  = ix + nAorb(iSymx)*(iy-1)

                    ixyG = nTA*(iyG-1) + ixG ! global index

                    Do iw=1,nAorb(iSymw)

                     iwG  = iAorb(iSymw) + iw

                     Do it=1,nAorb(iSymt)

                        itG  = iAorb(iSymt) + it

                        itw  = it + nAorb(iSymt)*(iw-1)

                        itwG = nTA*(iwG-1) + itG ! global index

                        iRes = ipInt + iTri(itwG,ixyG) - 1

                        Work(iRes) = Scr%pA(iSymw,iSymy)%A(itw,ixy)

                     End Do

                    End Do

                   End Do

               End Do

            End Do

         End Do


      ENDIF

      irc=0

      Return
      END

**************************************************************
