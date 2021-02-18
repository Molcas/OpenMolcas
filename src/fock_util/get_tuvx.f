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
      Subroutine Get_TUVX(PUVX,TUVX)
************************************************************************
*                                                                      *
*     extract the subset of ERIs (tu!vx) from the set (pu!vx)          *
*                                                                      *
*     calling arguments:                                               *
*     PUVX    : input, array of real                                   *
*               ERIs with indices (pu!vx)                              *
*               (three active, one general)                            *
*     TUVX    : output, array of real                                  *
*               ERIs with indices (tu!vx)                              *
*               (four active)                                          *
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

      Dimension PUVX(*), TUVX(*)

#include "rasdim.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='XXXXXXXX')

      Integer off_Ash(mxSym), off_PUVX(mxSym)

      iTri(i) = (i*i-i)/2

*     generate offsets
      iStack = 0
      Do iSym = 1,nSym
        off_Ash(iSym) = iStack
        iStack = iStack + nAsh(iSym)
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

*     select integrals with all 4 indices active

      Do iSym = 1,nSym
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iPUVX = off_PUVX(iSym)
        Do jSym = 1,nSym
          jAsh = nAsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym
            kAsh = nAsh(kSym)
            lSym = 1 + ieor(ijSym-1,kSym-1)
            lAsh = nAsh(lSym)

            If ( lSym.le.kSym .and.
     &           iAsh*jAsh*kAsh*lAsh.ne.0 ) then


c           Write(LF,*)'sym(p,w,x,y),offset= ',isym,jsym,ksym,lsym,iPUVX
c           call recprt('(pw|xy)','(1P,5D16.8)',PUVX(iPUVX+1),iorb*jAsh,
c     &                    kAsh*lAsh+Min(ijSym-2,0)*kAsh*(lAsh-1)/2)

              Do iV = 1,kAsh
                lMax = lAsh
                If ( kSym.eq.lSym ) lMax = iV
                Do iX = 1,lMax
                  Do iU = 1,jAsh
                    Do iP = 1,iOrb
                      iT = iP - iIsh
                      iPUVX=iPUVX+1
                      If ( iT.gt.0 .and. iT.le.iAsh ) then
                        iiT = iT + off_Ash(iSym)
                        iiU = iU + off_Ash(jSym)
                        If ( iiU.gt.iiT ) then
                          iiT = iU + off_Ash(jSym)
                          iiU = iT + off_Ash(iSym)
                        End If
*
                        iTU = iiU + iTri(iiT)
                        iiV = iV + off_Ash(kSym)
                        iiX = iX + off_Ash(lSym)
                        If ( iiX.gt.iiV ) then
                          iiV = iX + off_Ash(lSym)
                          iiX = iV + off_Ash(kSym)
                        End If
                        iVX = iiX + iTri(iiV)
                        If ( iVX.gt.iTU ) then
                          iTemp = iTU
                          iTU = iVX
                          iVX = iTemp
                        End If
                        iTUVX = iVX + iTri(iTU)
                        TUVX(iTUVX) = PUVX(iPUVX)
                      End If

                    End Do
                  End Do
                End Do
              End Do
            End If

          End Do
        End Do
      End Do

      Return
      End
