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

      SUBROUTINE REORD_Pmat(ipDA2,ipPmat,ipDSA2)

      Implicit Real*8 (a-h,o-z)

      Integer  ipDA2,ipPmat
      Integer  nnA(8,8),ipDSA2(8,8,8)
      Integer  iAorb(8),cho_irange
      External cho_irange
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"

C ************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
C ************************************************
      iTri(i,j) = Max(i,j)*(Max(i,j)-3)/2 + i + j
C ************************************************


      Call set_nnA(nSym,nAsh,nnA)

C --------------------------------------------------
C --- P[tw],xy :  stored as squared in [tw] and
C ---             lower triangular (packed) in (xy)
C ---             and blocked by symmetry
C --------------------------------------------------
      nPmat=0
      Do JSYM=1,nSym
         Do iSymy=1,nSym
            iSymx=MulD2h(iSymy,JSYM)
            if (iSymx.le.iSymy) then
            Do iSymw=1,nSym
               iSymt=MulD2h(iSymw,JSYM)
               ipDSA2(iSymw,iSymy,JSYM) = ipPmat + nPmat
               nPmat = nPmat
     &               + nAsh(iSymt)*nAsh(iSymw)*nnA(iSymx,iSymy)
            End Do
            endif
         End do
      End Do

      iAorb(1)= 0
      Do iSym = 2,nSym
         iAorb(iSym) = iAorb(iSym-1) + nAsh(iSym-1)
      End Do

      nAcOrb = iAorb(nSym) + nAsh(nSym)

C --- Copy out elements in the reorder 2-el density matrix
C --- From Canonical storage ( t>=w, x>=y, (tw)>=(xy) )
C --- to the cholesky storage P[tw],xy
C --- Diagonal (xy) elements are already scaled by a factor 1/2
C --------------------------------------------------------------

      Do iyG=1,nAcOrb    ! global index in total # of active orbitals

         iSymy = CHO_IRANGE(iyG,iAorb,nSym,.FALSE.)
         iy = iyG - iAorb(iSymy)

         Do ixG=iyG,nAcOrb

            ixyG = ixG*(ixG-1)/2 + iyG
            iSymx = CHO_IRANGE(ixG,iAorb,nSym,.FALSE.)
            ix = ixG - iAorb(iSymx)
            if (iSymx.eq.iSymy) then
               ixy = iTri(ix,iy)
            else
               ixy  = nAsh(iSymx)*(iy-1) + ix
            endif

            Do iwG=1,nAcOrb

               iSymw = CHO_IRANGE(iwG,iAorb,nSym,.FALSE.)
               iw = iwG - iAorb(iSymw)

               Do itG=iwG,nAcOrb

                  iSymt = CHO_IRANGE(itG,iAorb,nSym,.FALSE.)
                  it = itG - iAorb(iSymt)

                  iSymTW = MulD2h(iSymt,iSymw)
                  iSymXY = MulD2h(iSymx,iSymy)

                  itwG = itG*(itG-1)/2 + iwG

                  If (iSymTW.eq.iSymXY .and. itwG.ge.ixyG) Then

                     kfrom = ipDA2 + iTri(itwG,ixyG) - 1

                     Ntw = nAsh(iSymt)*nAsh(iSymw)
                     itw1 = nAsh(iSymt)*(iw-1) + it

                     kto1 = ipDSA2(iSymw,iSymy,iSymXY) - 1
     &                    + Ntw*(ixy-1) + itw1

                     Work(kto1) = Work(kfrom)

                     itw2 = nAsh(iSymw)*(it-1) + iw

                     kto2 = ipDSA2(iSymt,iSymy,iSymXY) - 1
     &                    + Ntw*(ixy-1) + itw2

                     Work(kto2) = Work(kfrom)

                     if (iSymt.eq.iSymw) then
                        itw = iTri(it,iw)
                     else
                        itw  = nAsh(iSymt)*(iw-1) + it
                     endif
                     Nxy = nAsh(iSymx)*nAsh(iSymy)
                     ixy1 = nAsh(iSymx)*(iy-1) + ix

                     kto3 = ipDSA2(iSymy,iSymw,iSymXY) - 1
     &                    + Nxy*(itw-1) + ixy1

                     Work(kto3) = Work(kfrom)

                     ixy2 = nAsh(iSymy)*(ix-1) + iy

                     kto4 = ipDSA2(iSymx,iSymw,iSymXY) - 1
     &                    + Nxy*(itw-1) + ixy2

                     Work(kto4) = Work(kfrom)

                  End If

               End Do

            End Do

         End Do

      End Do


      Return
      END
