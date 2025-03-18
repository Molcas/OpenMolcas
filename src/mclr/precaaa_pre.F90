!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996, Anders Bernhardsson                              *
!***********************************************************************
      Subroutine Precaaa_Pre(ActInt,A_J,Scr)
!
      use MCLR_Data, only: nA
      use input_mclr, only: ntAsh,nSym,nAsh,nIsh,nBas
      Implicit None
      Real*8 ActInt(ntAsh,ntAsh,ntAsh,ntAsh),A_J(*),Scr(*)
!
      Integer iSym, iA, iAabs, iAtot, iB, iBabs, iBtot, jSym, iC,       &
     &        iCabs, iCtot, iD, iDabs, iDtot
      Real*8 Val

      Do iSym = 1, nSym
        Do iA = 1, nAsh(iSym)
          iAabs = iA + nA(iSym)
          iAtot = iA + nIsh(iSym)
          Do iB = 1, nAsh(iSym) ! iA
            iBabs = iB + nA(iSym)
            iBtot = iB + nIsh(iSym)
            Do jSym = 1, iSym
              Call Coul(iSym,iSym,jSym,jSym,iAtot,iBtot,A_J,Scr)
              Do iC = 1, nAsh(jSym)
                iCabs = iC + nA(jSym)
                iCtot = iC + nIsh(jSym)
                Do iD = 1, nAsh(jSym) ! iC
                  iDabs = iD + nA(jSym)
                  iDtot = iD + nIsh(jSym)
                  Val = A_J(iCtot+nBas(jSym)*(iDtot-1))
                  ActInt(iAabs,iBabs,iCabs,iDabs) = Val
!                 ActInt(iBabs,iAabs,iCabs,iDabs) = Val
!                 ActInt(iAabs,iBabs,iDabs,iCabs) = Val
!                 ActInt(iBabs,iAabs,iDabs,iCabs) = Val
!                 ActInt(iCabs,iDabs,iAabs,iBabs) = Val
!                 ActInt(iCabs,iDabs,iBabs,iAabs) = Val
!                 ActInt(iDabs,iCabs,iAabs,iBabs) = Val
!                 ActInt(iDabs,iCabs,iBabs,iAabs) = Val
                End Do
              End Do
            End Do
          End Do
        End Do
      End Do
!
      End Subroutine Precaaa_Pre
