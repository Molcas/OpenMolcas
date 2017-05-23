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
      Subroutine Calc_PUVX2(PUVX,nPUVX,TabMO,mAO,nCoor,nTabMOs,
     &                     dF_dRho,ndF_dRho,nD,Weights)
      Implicit Real*8 (A-H,O-Z)
      Dimension PUVX(nPUVX)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_PUVX(mxSym),off_Bas(mxSym)
      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),
     &       dF_dRho(ndF_dRho,nCoor)
*
      lsym_tmp=lsym
*
*      Check dimensions: This is inconsistent! RL
*
      If(ndF_dRho.eq.3.or.ndF_dRho.eq.5) Then
      Else
         Call WarningMessage(2,'Calc_PUVX2: Dim. error!!!')
         Write(6,*) 'ndF_Rho:',ndF_dRho
         Call Abend()
      End If
*     generate offsets
      iStack  = 0
      iStack1 = 0
      Do iSym = 1,nSym
        off_Ash(iSym)    = iStack
        off_Bas(iSym)    = iStack1
        off_BasAsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)
        iStack1 = iStack1 + nBas(iSym)
        iStack  = iStack  + nAsh(iSym)
      End Do
*
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
      NrInt=iStack
*
      If (nPUVX.ne.NrInt) Then
         Call WarningMessage(2,
     &              ' Wrong number of two electron DFT int.!!!')
         Call Abend()
      End If
*

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
              Do iV = 1,kAsh
                jV = iV + off_BasAsh(kSym)
                lMax = lAsh
                If ( kSym.eq.lSym ) lMax = iV
                Do iX = 1,lMax
                  jX = iX + off_BasAsh(lSym)
                  Do iU = 1,jAsh
                  jU = iU + off_BasAsh(jSym)
                    Do iP = 1,iOrb
                      iT = iP - iIsh
                      iPUVX=iPUVX+1
                      jP = iP     +   off_Bas(iSym)
      If(ndF_dRho/nD.eq.4) Then
************************************************************************
                      Do iGrid=1,nCoor
                      PUVX(iPUVX) = PUVX(iPUVX) +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*dF_dRho(2,iGrid) +
*
     &                  (TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)+

     &                   TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)+

     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)+

     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX))*
     &                   Weights(iGrid)*dF_dRho(4,iGrid) +
*
     &                  (TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)+

     &                   TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)+

     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)+

     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX))*
     &                   Weights(iGrid)*dF_dRho(6,iGrid) +
*
     &                  (TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)+

     &                   TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)+

     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)+

     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX))*
     &                   Weights(iGrid)*dF_dRho(8,iGrid)
                      End Do
********************************************************************
      Else
********************************************************************
                      Do iGrid=1,nCoor
                      PUVX(iPUVX) = PUVX(iPUVX) +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*dF_dRho(2,iGrid)
                      End Do
*
********************************************************************
      End If
                    End Do
                  End Do
                End Do
              End Do
            End If

          End Do
        End Do
      End Do
*
      lsym=lsym_tmp
*
      Return
      End
