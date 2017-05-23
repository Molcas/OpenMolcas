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
* Copyright (C) 1996, Markus P. Fuelscher                              *
************************************************************************
      Subroutine Upd_FA(PUVX,F,D,ExFac)
************************************************************************
*                                                                      *
*     compute FIA, FAA, and FAS from the integral set (pu!vx)          *
*                                                                      *
*     calling arguments:                                               *
*     PUVX    : input, array of real                                   *
*               ERIs with indices (pu!vx)                              *
*               (three active, one general)                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

      Dimension PUVX(*), F(*), D(*)

#include "rasdim.fh"
#include "general.fh"

      Integer case, state_symmetry
      Integer   off_PUVX, off_Dmat, off_Fmat
      Dimension off_PUVX(mxSym), off_Dmat(mxSym), off_Fmat(mxSym)

      iTri(i)=(i*i-i)/2

*     nasty, but necessary
      state_symmetry=lSym

*     generate offsets

      iStack = 0
      Do iSym = 1,nSym
         off_Dmat(iSym) = iStack
         iAsh = nAsh(iSym)
         iStack = iStack+ (iAsh*iAsh+iAsh)/2
      End Do

      iStack = 0
      Do iSym = 1,nSym
         off_Fmat(iSym) = iStack
         iOrb = nOrb(iSym)
         iStack = iStack+ (iOrb*iOrb+iOrb)/2
      End Do

      iStack = 0
      Do iSym = 1,nSym
        off_PUVX(iSym) = iStack
        iOrb = nOrb(iSym)
        Do jSym = 1,nSym
          jAsh = nAsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym
            kAsh = nAsh(kSym)
            Do lSym = 1,kSym
              lAsh = nAsh(lSym)
              klSym = 1 + ieor(kSym-1,lSym-1)
              If ( ijSym.eq.klSym) then
                kl_Orb_pairs = kAsh*lAsh
                If ( kSym.eq.lSym ) kl_Orb_pairs = (kAsh*kAsh+kAsh)/2
                iStack = iStack + iOrb*jAsh*kl_Orb_pairs
              End If
            End Do
          End Do
        End Do
      End Do

*     clear the subbocks FIA, FAA and FSA

      Do iSym = 1,nSym
        iOrb  = nOrb(iSym)
        iAsh  = nAsh(iSym)
        iIsh  = nIsh(iSym)
        iFoff = off_Fmat(iSym)
        Do iU = iIsh+1,iIsh+iAsh
          Do iP = 1,iU
            iPU = iP + iTri(iU)
            F(iFoff+iPU) = 0.0d0
          End Do
        End Do
        Do iU = iIsh+iAsh+1,iOrb
          Do iP = iIsh+1,iIsh+iAsh
            iPU = iP + iTri(iU)
            F(iFoff+iPU) = 0.0d0
          End Do
        End Do
      End Do

*     generate the subblocks FIA, FAA and FSA

      Do iSym = 1,nSym
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iPUVX = off_PUVX(iSym)
        Do jSym = 1,iSym
          jOrb = nOrb(jSym)
          jAsh = nAsh(jSym)
          jIsh = nIsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym
            kAsh = nAsh(kSym)
            Do lSym = 1,kSym
              lAsh = nAsh(lSym)
              klSym = 1 + ieor(kSym-1,lSym-1)

*             find cases
              case = 4
              If ( iSym.eq.jSym ) case = case-2
              If ( iSym.eq.kSym ) case = case-1

              If ( ijSym.eq.klSym .and.
     &             iAsh*jAsh*kAsh*lAsh.ne.0 ) then

                Goto (100,200,300,400) case

*               symmetry case (II!II)
100             Continue
                iFoff = off_Fmat(iSym)
                iDoff = off_Dmat(iSym)
                Do iV = 1,kAsh
                  Do iX = 1,iV
                    iVX = iTri(iV) + iX
                    DVX = 2.0D0*D(iDoff+iVX)
                    If ( iX.eq.iV ) DVX = D(iDoff+iVX)
                    Do iU = 1,jAsh
                      iUV = iTri(iU) + iV
                      If ( iV.gt.iU ) iUV  = iTri(iV) + iU
                      DUV = ExFac*0.5D0*D(iDoff+iUV)
                      iUX = iTri(iU) + iX
                      If ( iX.gt.iU ) iUX  = iTri(iX) + iU
                      DUX = ExFac*0.5D0*D(iDoff+iUX)
                      If ( iX.eq.iV ) then
                        DUV = ExFac*0.5D0*DUV
                        DUX = ExFac*0.5D0*DUX
                      End If
                      iPUVX = off_PUVX(iSym)
*                     inactive/active block
                      Do iP = 1,iIsh
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPU  = iTri(iIsh+iU) + iP
                        F(iFoff+iPU) = F(iFoff+iPU) + DVX*Temp
                        iPV  = iTri(iIsh+iV) + iP
                        F(iFoff+iPV) = F(iFoff+iPV) - DUX*Temp
                        iPX  = iTri(iIsh+iX) + iP
                        F(iFoff+iPX) = F(iFoff+iPX) - DUV*Temp
                      End Do
*                     active/active block and iP<=(iIsh+iU)
                      Do iP = iIsh+1,iIsh+iU
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPU  = iTri(iIsh+iU) + iP
                        F(iFoff+iPU) = F(iFoff+iPU) + DVX*Temp
                        iPV  = iTri(iIsh+iV) + iP
                        If ( iP.gt.(iIsh+iV) )
     &                  iPV  = iTri(iP) + iIsh + iV
                        F(iFoff+iPV) = F(iFoff+iPV)-ExFac*0.5D0*DUX*Temp
                        If ( iP.eq.(iIsh+iV) )
     &                  F(iFoff+iPV) = F(iFoff+iPV)-ExFac*0.5D0*DUX*Temp
                        iPX  = iTri(iIsh+iX) + iP
                        If ( iP.gt.(iIsh+iX) )
     &                  iPX  = iTri(iP) + iIsh + iX
                        F(iFoff+iPX) = F(iFoff+iPX)-ExFac*0.5D0*DUV*Temp
                        If ( iP.eq.(iIsh+iX) )
     &                  F(iFoff+iPX) = F(iFoff+iPX)-ExFac*0.5D0*DUV*Temp
                      End Do
*                     active/active block and iP>(iIsh+iU)
                      Do iP = iIsh+iU+1,iIsh+iAsh
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPV  = iTri(iIsh+iV) + iP
                        If ( iP.gt.(iIsh+iV) )
     &                  iPV  = iTri(iP) + iIsh + iV
                        F(iFoff+iPV) = F(iFoff+iPV)-ExFac*0.5D0*DUX*Temp
                        If ( iP.eq.(iIsh+iV) )
     &                  F(iFoff+iPV) = F(iFoff+iPV)-ExFac*0.5D0*DUX*Temp
                        iPX  = iTri(iIsh+iX) + iP
                        If ( iP.gt.(iIsh+iX) )
     &                  iPX  = iTri(iP) + iIsh + iX
                        F(iFoff+iPX) = F(iFoff+iPX)-ExFac*0.5D0*DUV*Temp
                        If ( iP.eq.(iIsh+iX) )
     &                  F(iFoff+iPX) = F(iFoff+iPX)-ExFac*0.5D0*DUV*Temp
                      end Do
*                     active/secondary block
                      Do iP = iIsh+iAsh+1,iOrb
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPU  = iTri(iP) + iIsh + iU
                        F(iFoff+iPU) = F(iFoff+iPU) + DVX*Temp
                        iPV  = iTri(iP) + iIsh + iV
                        F(iFoff+iPV) = F(iFoff+iPV) - DUX*Temp
                        iPX  = iTri(iP) + iIsh + iX
                        F(iFoff+iPX) = F(iFoff+iPX) - DUV*Temp
                      End Do
                      off_PUVX(iSym) = off_PUVX(iSym) + iOrb
                    End Do
                  End Do
                End Do
                Goto 500

*               symmetry case (II!KK)
200             Continue
                iFoff = off_Fmat(iSym)
                kDoff = off_Dmat(kSym)
                Do iV = 1,kAsh
                  Do iX = 1,iV
                    iVX = iTri(iV) + iX
                    DVX = 2.0D0*D(kDoff+iVX)
                    If ( iX.eq.iV ) DVX = D(kDoff+iVX)
                    Do iU = 1,jAsh
                      iPUVX = off_PUVX(iSym)
*                     inactive/active block
                      Do iP = 1,iIsh
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPU  = iTri(iIsh+iU) + iP
                        F(iFoff+iPU) = F(iFoff+iPU) + DVX*Temp
                      End Do
*                     active/active block and iP<=(iIsh+iU)
                      Do iP = iIsh+1,iIsh+iU
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPU  = iTri(iIsh+iU) + iP
                        F(iFoff+iPU) = F(iFoff+iPU) + DVX*Temp
                      End Do
                      iPUVX = iPUVX + iAsh - iU
*                     active/secondary block
                      Do iP = iIsh+iAsh+1,iOrb
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPU  = iTri(iP) + iIsh + iU
                        F(iFoff+iPU) = F(iFoff+iPU) + DVX*Temp
                      End Do
                      off_PUVX(iSym) = off_PUVX(iSym) + iOrb
                    End Do
                  End Do
                End Do
                Goto 500

*               symmetry case (IJ!IJ)
300             Continue
                iFoff = off_Fmat(iSym)
                jFoff = off_Fmat(jSym)
                iDoff = off_Dmat(iSym)
                jDoff = off_Dmat(jSym)
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                    Do iU= 1,jAsh
                      iUX = iTri(iU) + iX
                      If ( iX.gt.iU ) iUX  = iTri(iX) + iU
                      DUX = ExFac*0.5D0*D(jDoff+iUX)
                      iPUVX = off_PUVX(iSym)
*                     inactive/active block
                      Do iP = 1,iIsh
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPV  = iTri(iIsh+iV) + iP
                        F(iFoff+iPV) = F(iFoff+iPV) - DUX*Temp
                      End Do
*                     active/active block and iP<=(iIsh+iV)
                      Do iP = iIsh+1,iIsh+iV
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPV  = iTri(iIsh+iV) + iP
                        F(iFoff+iPV) = F(iFoff+iPV) - DUX*Temp
                      End Do
                      iPUVX = iPUVX + iAsh - iV
*                     active/secondary block
                      Do iP = iIsh+iAsh+1,iOrb
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPV  = iTri(iP) + iIsh + iV
                        F(iFoff+iPV) = F(iFoff+iPV) - DUX*Temp
                      End Do
                      off_PUVX(iSym) = off_PUVX(iSym) + iOrb
                    End Do
                    Do iU= 1,iAsh
                      iUV = iTri(iU) + iV
                      If ( iV.gt.iU ) iUV  = iTri(iV) + iU
                      DUV = ExFac*0.5D0*D(iDoff+iUV)
                      iPUVX = off_PUVX(jSym)
*                     inactive/active block
                      Do iP = 1,jIsh
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPX  = iTri(jIsh+iX) + iP
                        F(jFoff+iPX) = F(jFoff+iPX) - DUV*Temp
                      End Do
*                     active/active block and iP<=(jIsh+iX)
                      Do iP = jIsh+1,jIsh+iX
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPX  = iTri(jIsh+iX) + iP
                        F(jFoff+iPX) = F(jFoff+iPX) - DUV*Temp
                      End Do
                      iPUVX = iPUVX + jAsh - iX
*                     active/secondary block
                      Do iP = jIsh+jAsh+1,jOrb
                        iPUVX = iPUVX + 1
                        Temp = PUVX(iPUVX)
                        iPX  = iTri(iP) + jIsh + iX
                        F(jFoff+iPX) = F(jFoff+iPX) - DUV*Temp
                      End Do
                      off_PUVX(jSym) = off_PUVX(jSym) + jOrb
                    End Do
                  End Do
                End Do
                Goto 500

*               symmetry case (IJ!KL)
400             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                    off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                    off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do

500             Continue
              End If

            End Do
          End Do
        End Do
      End Do

*     nasty, but necessary
      lSym=state_symmetry

      Return
      End
