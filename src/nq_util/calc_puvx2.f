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



      Subroutine Calc_OTPUVX(PUVX,TabMO,mAO,nCoor,nTabMOs
     &                       ,P2_ontop,nP2_ontop,Rho,nRho,
     &                      dF_dRho,ndF_dRho,RhoI,RhoA,mRho,
     &                      Weights,D1MO,nD1MO,nIrrep)
      Implicit Real*8 (A-H,O-Z)
      Dimension PUVX(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_PUVX(mxSym),off_Bas(mxSym),
     &        off_ish(mxSym),off_BasIsh(mxSym),
     &        off_BasVsh(mxSym)
      Integer   off_Dmat, off_Fmat
      Dimension off_Dmat(mxSym), off_Fmat(mxSym)

      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),P2_ontop(nP2_ontop,nCoor),
     &       dF_dRho(ndF_dRho,nCoor),Rho(nRho,nCoor),
     &       RhoI(mRho,nCoor),RhoA(mRho,nCoor)
      Integer nIrrep
      Real*8 D1MO(nD1MO)
      Integer count_tmp
      Real*8 DVX
      real*8 Fact
      Integer case
      iTrii(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      iTri(i)=(i*i-i)/2
*
      thrsrho=1.0d-15
      thrsrho2=1.0d-15
      thrspi=1.0d-30
*
      lsym_tmp=lsym

      Call Unused_real_array(RhoI)
      Call Unused_real_array(RhoA)

      iStack  = 0
      iStack1 = 0
      iStack2 = 0
      off_ish(:) = 0
      off_Ash(:) = 0
      off_Bas(:) = 0
      off_BasAsh(:) = 0
      off_BasIsh(:) = 0
      off_BasVsh(:) = 0
      ntot1 = 0

      Do iSym = 1,nSym
        off_ish(isym)    = iStack2
        off_Ash(iSym)    = iStack
        off_Bas(iSym)    = iStack1
        off_BasVsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)+nAsh(iSym)
        off_BasAsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)
        off_BasIsh(iSym) = iStack1+nFro(iSym)
        ntot1 = iTrii(nBas(iSym),nBas(iSym)) + ntot1
        iStack2 = iStack2 + nIsh(iSym)
        iStack1 = iStack1 + nBas(iSym)
        iStack  = iStack  + nAsh(iSym)
      End Do



!count the number of tmp_pot:
      count_tmp = 0
      do isym=1,nsym
        do jsym=1,nsym
        count_tmp = count_tmp + nIsh(isym)*(nIsh(jsym)+nAsh(jsym))**2
        end do
      end do

*
!Calculate PUVX offsets
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
*

************************************************************************
* BUILD PUVX (needed to construct FA)
************************************************************************
      Fact = 1.00d0
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
                    Do iGrid=1,nCoor
                       dTot=Rho(1,iGrid)+Rho(2,iGrid)
                     ratio = 0.0d0
                  if(dTot.ge.thrsrho.and.
     &               P2_ontop(1,iGrid).ge.thrspi) then
              ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                      if((1.0d0-ratio).gt.thrsrho2) then
                      Zeta  = sqrt(1.0d0-ratio)
                      PUVX(iPUVX) = PUVX(iPUVX) +
     &                   Fact*(1.0D0/(Zeta*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*dble(nIrrep)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      else
                     PUVX(iPUVX) = PUVX(iPUVX) +0.0D0
                     end if
                   end if
                      End Do
                    End Do
                  End Do
                End Do
              End Do
            End If
          End Do
        End Do
      End Do

!Construction of Fock matrix pieces
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
!I think I can generate the contributions from the potentials to the new
!FI and FA terms at this level.  Then there will be no need to store a
!large number of V_TE potentials.
!      CALL GETMEM('FI_V','ALLO','REAL',ifiv,ntot1)
!      CALL GETMEM('FA_V','ALLO','REAL',ifav,ntot1)
!      Call Get_dArray('FI_V',Work(ifiv),ntot1)
!      Call Get_dArray('FA_V',Work(ifav),ntot1)

      Do iSym = 1,nSym !sym for p
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iPUVX = off_PUVX(iSym)
        Do jSym = 1,nSym !sym for u
          jOrb = nOrb(jSym)
          jAsh = nAsh(jSym)
          jIsh = nIsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym !Sym for v
            korb = nOrb(kSym)
            kAsh = nAsh(kSym)
            kIsh = nIsh(kSym)
            Do lSym = 1,kSym !sym for x
              lOrb = nOrb(lSym)
              lAsh = nAsh(lSym)
              lIsh = nIsh(lSym)
              klSym = 1 + ieor(kSym-1,lSym-1)

*             find cases
              case = 4
              If ( iSym.eq.jSym ) case = case-2
              If ( iSym.eq.kSym ) case = case-1

              If ( ijSym.eq.klSym .and.
!     &             iAsh*jAsh*kAsh*lAsh.ne.0 ) then
     &             iOrb*jOrb*kOrb*lOrb.ne.0 ) then

                Goto (100,200,300,400) case

!Symmetry case (II|II)
100             Continue
                iFoff = off_Fmat(iSym)
                iDoff = off_Dmat(iSym)
                Do iV = 1,kIsh
                  jV = iV + off_basIsh(kSym)
                  iX=iV
                  jX=jV
                  DVX=2.0D0
                  If ( iX.eq.iV ) DVX = DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !V
!ACTIVE CONTRIBUTIONS
                Do iV = 1,kAsh
                  jV = iV + off_basAsh(ksym)
                do iX=1,iV
                    jX = iX + off_basAsh(ksym)
                    !iVX = iTri(iV) + iX
                    !jVX = iTri(jV) + jX
                    iVX = iTri(iV+off_Ash(ksym)) + iX + off_Ash(ksym)
                    DVX=D1MO(iVX)
                    !DVX=D1MO(iDoff+iVX)!Why *2?
                  If ( iX.eq.iV ) DVX = DVX*0.5
          !          DVX=2.0D0*D1MO(iDoff+iVX)!Why *2?
          !          If (iX.eq.iV) DVX=D1MO(iDoff+iVX)
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iU + off_basAsh(jSym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(iSym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !X
                end do !V
                Goto 500
!Symmetry case (II|KK)
200             Continue
                iFoff = off_Fmat(iSym)
                kDoff = off_Dmat(kSym)
                Do iV = 1,kIsh
                  jV = iV + off_basIsh(ksym)
                  iX=iV
                  jX=jV
                  DVX=2.0D0
                  If ( iX.eq.iV ) DVX = DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+nAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !V
!ACTIVE CONTRIBUTIONS
                Do iV = 1,kAsh
                  jV = iV + off_basAsh(ksym)
                  do iX=1,iV
                    jX = iX + off_basAsh(ksym)
                    iVX = iTri(iV+off_Ash(ksym)) + iX + off_Ash(ksym)
                    DVX=D1MO(iVX)
                    !DVX=D1MO(kDoff+iVX)
                  If ( iX.eq.iV ) DVX = DVX*0.5
            !        DVX=2.0D0*D1MO(iDoff+iVX)
            !        If (iX.eq.iV) DVX=D1MO(iDoff+iVX)
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+nAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !X
                end do !V
                Goto 500
*               symmetry case (IJ!IJ)
300             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                   ! off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                   ! off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do
                Goto 500
*               symmetry case (IJ!KL)
400             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                  !  off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                  !  off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do
                Goto 500

500             Continue
              End If

            End Do
          End Do
        End Do
      End Do

!      Call Put_dArray('FI_V',Work(ifiv),ntot1)
!      Call Put_dArray('FA_V',Work(ifav),ntot1)
!      CALL GETMEM('FI_V','FREE','REAL',ifiv,ntot1)
!      CALL GETMEM('FA_V','FREE','REAL',ifav,ntot1)
      lsym=lsym_tmp
*
      Return
      End


      Subroutine Calc_OTOE(OE,TabMO,mAO,nCoor,nTabMOs
     &                       ,P2_ontop,nP2_ontop,Rho,nRho,
     &                      dF_dRho,ndF_dRho,RhoI,RhoA,mRho,Weights,
     &                      nIrrep)
      Implicit Real*8 (A-H,O-Z)
      Dimension OE(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_Bas(mxSym)!,off_Ish(mxSym)!,off_Fmat(mxSym),
!     &        off_Dmat(mxSym)
      Integer off_basIsh(mxSym),off_Fmat(mxSym)
      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),P2_ontop(nP2_ontop,nCoor),
     &       dF_dRho(ndF_dRho,nCoor),Rho(nRho,nCoor),
     &       RhoI(mRho,nCoor),RhoA(mRho,nCoor)
      Integer VX,ix,jx,iv,jv
      Integer nIrrep
*
      !iTri(i) = (i*i-i)/2
       iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      thrsrho=1.0d-15
      thrsrho2=1.0d-15
*
      Call Unused_real_array(RhoI)
      Call Unused_real_array(RhoA)

      lsym_tmp=lsym
      iStack  = 0
      iStack1 = 0
      Do iSym = 1,nSym
        off_Ash(iSym)    = iStack
        off_Bas(iSym)    = iStack1
        off_BasAsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)
        off_BasIsh(iSym) = iStack1+nFro(iSym)
        iStack1 = iStack1 + nBas(iSym)
        iStack  = iStack  + nAsh(iSym)
      End Do

      iStack = 0
      do iSym=1,nSym
        off_Fmat(iSym) = iStack
        iorb = nOrb(iSym)
        iStack = iStack + (iOrb*iOrb + iOrb)/2
      end do
*

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!OE pieces - both orbs must be in the same irrep, eh?
      Do iSym = 1,nSym
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        Do iV = 1,iOrb!iIsh+iAsh+iVsh
          jV = iV + off_BasIsh(iSym)
          do iX = 1,iV
            jX = iX + off_BasIsh(iSym)
            VX = off_Fmat(iSym) + iTri(iV,iX)
            fact=1.0d0
            Do iGrid = 1, nCoor
            dTot=Rho(1,iGrid)+Rho(2,iGrid)
            ratio = 0.0d0
                ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                if (dTot.ge.thrsrho) then
                 if((1.0d0-ratio).gt.thrsrho2) then
!                if (.false.) then
                  Zeta  = sqrt(1.0d0-ratio)
                  OE(VX)=OE(VX)+TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                    *(dF_dRho(1,iGrid)*(0.5D0*(1+Zeta)
     &                    +2.0D0*P2_ontop(1,iGrid)/(Zeta*dTot**2))
     &                    +dF_dRho(2,iGrid)*(0.5D0*(1-Zeta)
     &                    -2.0D0*P2_ontop(1,iGrid)/(Zeta*dTot**2)))
     &                    *Weights(iGrid)*Dble(nIrrep)

                  else
                    OE(VX)=OE(VX)+TabMO(1,iGrid,jV)*
     &                     TabMO(1,iGrid,jX)*(dF_dRho(1,iGrid)
     &                     +dF_dRho(2,iGrid))*Weights(iGrid)*0.5
     *                     *dble(nIrrep)
                  end if
                 end if
            End Do     ! iGrid
          End do
        End Do
      End Do
           lsym=lsym_tmp
      Return
      End

      Subroutine Calc_OTOEf(OE,TabMO,mAO,nCoor,nTabMOs
     &                       ,P2_ontop,nP2_ontop,Rho,nRho,
     &                      dF_dRho,ndF_dRho,RhoI,RhoA,mRho,Weights,
     &                      nIrrep)
      Implicit Real*8 (A-H,O-Z)
      Dimension OE(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_Bas(mxSym)
!     &        off_Dmat(mxSym)
      Integer off_basIsh(mxSym),off_Fmat(mxSym)
      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),P2_ontop(nP2_ontop,nCoor),
     &       dF_dRho(ndF_dRho,nCoor),Rho(nRho,nCoor),
     &       RhoI(mRho,nCoor),RhoA(mRho,nCoor)
!      Integer count_tmp
      Integer VX,ix,jx,iv,jv
      Integer nIrrep
*
      !iTri(i) = (i*i-i)/2
       iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      thrsrho=1.0d-15
      thrsrho2=1.0d-15
      thrsrho3=0.9000000000d0
      thrsrho4=1.1500000000d0
      Ab1=-4.756065601d+2
      Bb1=-3.794733192d+2
      Cb1=-8.538149682d+1
*
      Call Unused_real_array(RhoI)
      Call Unused_real_array(RhoA)

      lsym_tmp=lsym
      iStack  = 0
      iStack1 = 0
      Do iSym = 1,nSym
        off_Ash(iSym)    = iStack
        off_Bas(iSym)    = iStack1
        off_BasAsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)
        off_BasIsh(iSym) = iStack1+nFro(iSym)
        iStack1 = iStack1 + nBas(iSym)
        iStack  = iStack  + nAsh(iSym)
      End Do

      iStack = 0
      do iSym=1,nSym
        off_Fmat(iSym) = iStack
        iorb = nOrb(iSym)
        iStack = iStack + (iOrb*iOrb + iOrb)/2
      end do
*

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!OE pieces - both orbs must be in the same irrep, eh?
      Do iSym = 1,nSym
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        Do iV = 1,iOrb!iIsh+iAsh
          jV = iV + off_BasIsh(iSym)
          do iX = 1,iV
            jX = iX + off_BasIsh(iSym)
            VX = off_Fmat(iSym) + iTri(iV,iX)
            fact=1.0d0
            Do iGrid = 1, nCoor
            dTot=Rho(1,iGrid)+Rho(2,iGrid)
            ratio = 0.0d0
                ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                if (dTot.ge.thrsrho) then
                 if(((1.0d0-ratio).gt.thrsrho2)
     &               .and.(ratio.lt.thrsrho3)) then
                  Zeta  = sqrt(1.0d0-ratio)
                  OE(VX)=OE(VX)+TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                    *(dF_dRho(1,iGrid)*(0.5D0*(1+Zeta)
     &                    +2.0D0*P2_ontop(1,iGrid)/(Zeta*dTot**2))
     &                    +dF_dRho(2,iGrid)*(0.5D0*(1-Zeta)
     &                    -2.0D0*P2_ontop(1,iGrid)/(Zeta*dTot**2)))
     &                    *Weights(iGrid)*nIrrep
                  else if((ratio.ge.thrsrho3)
     &              .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (10.0D0*Ab1*(ratio-1.15d0)**4.0d0) +
     &                        (8.0D0*Bb1*(ratio-1.15d0)**3.0d0) +
     &                        (6.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                     OE(VX)=OE(VX)+TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                    *(dF_dRho(1,iGrid)*(0.5D0*(1.0D0+Zeta)
     &                    -2.0D0*P2_ontop(1,iGrid)/(dTot**2)*Deriv)
     &                    +dF_dRho(2,iGrid)*(0.5D0*(1.0D0-Zeta)
     &                    +2.0D0*P2_ontop(1,iGrid)/(dTot**2)*Deriv))
     &                    *Weights(iGrid)*nIrrep

                  else
                    OE(VX)=OE(VX)+TabMO(1,iGrid,jV)*
     &                     TabMO(1,iGrid,jX)*(dF_dRho(1,iGrid)
     &                     +dF_dRho(2,iGrid))*Weights(iGrid)*0.5*nIrrep
                  end if
                 end if
            End Do     ! iGrid
          End do
        End Do
      End Do

           lsym=lsym_tmp

      Return
      End

      Function Delta(x,y)
      Integer Delta
      Integer x,y

      if (x.eq.y) then
        Delta = 1
      else
        Delta = 0
      end if
      end Function


      Subroutine Calc_OTPUVX_ft(PUVX,TabMO,mAO,nCoor,nTabMOs
     &                       ,P2_ontop,nP2_ontop,Rho,nRho,
     &                      dF_dRho,ndF_dRho,RhoI,RhoA,mRho,
     &                      Weights,D1MO,nD1MO,nIrrep)
      Implicit Real*8 (A-H,O-Z)
      Dimension PUVX(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer nIrrep
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_PUVX(mxSym),off_Bas(mxSym),
     &        off_ish(mxSym),off_BasIsh(mxSym),
     &        off_BasVsh(mxSym)
      Integer   off_Dmat, off_Fmat
      Dimension off_Dmat(mxSym), off_Fmat(mxSym)

      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),P2_ontop(nP2_ontop,nCoor),
     &       dF_dRho(ndF_dRho,nCoor),Rho(nRho,nCoor),
     &       RhoI(mRho,nCoor),RhoA(mRho,nCoor)
      Real*8 D1MO(nD1MO)
              Integer count_tmp
      Integer p,q
      Real*8 DVX
      real*8 Fact
      Integer case
      iTrii(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      iTri(i)=(i*i-i)/2
*
      thrsrho=1.0d-15
      thrsrho2=1.0d-15
      thrsrho3=0.9000000000d0
      thrsrho4=1.1500000000d0
      Ab1=-4.756065601d+2
      Bb1=-3.794733192d+2
      Cb1=-8.538149682d+1
*
      Call Unused_real_array(RhoI)
      Call Unused_real_array(RhoA)
      lsym_tmp=lsym

      iStack  = 0
      iStack1 = 0
      iStack2 = 0
      off_ish(:) = 0
      off_Ash(:) = 0
      off_Bas(:) = 0
      off_BasAsh(:) = 0
      off_BasIsh(:) = 0
      off_BasVsh(:) = 0

      Do iSym = 1,nSym
        off_ish(isym)    = iStack2
        off_Ash(iSym)    = iStack
        off_Bas(iSym)    = iStack1
        off_BasVsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)+nAsh(iSym)
        off_BasAsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)
        off_BasIsh(iSym) = iStack1+nFro(iSym)
        ntot1 = iTrii(nBas(iSym),nBas(iSym))
        iStack2 = iStack2 + nIsh(iSym)
        iStack1 = iStack1 + nBas(iSym)
        iStack  = iStack  + nAsh(iSym)
      End Do



!count the number of tmp_pot:
      count_tmp = 0
      do isym=1,nsym
        do jsym=1,nsym
        count_tmp = count_tmp + nIsh(isym)*(nIsh(jsym)+nAsh(jsym))**2
        end do
      end do

*
!Calculate PUVX offsets
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
*

************************************************************************
* BUILD PUVX (needed to construct FA)
************************************************************************
      Fact = 1.00d0
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
                    Do iGrid=1,nCoor
                       dTot=Rho(1,iGrid)+Rho(2,iGrid)
                     ratio = 0.0d0
                  if(dTot.ge.thrsrho) then
              ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                    if((1.0d0-ratio).gt.thrsrho2.and.
     &                 (ratio.lt.thrsrho3)) then
                      Zeta  = sqrt(1.0d0-ratio)
                      PUVX(iPUVX) = PUVX(iPUVX) +
     &                   Fact*(1.0D0/(Zeta*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*nIrrep*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      else if((ratio.ge.thrsrho3).and.
     &                       (ratio.le.thrsrho4)) then
                        Ze_int = ratio - 1.15d0
                        PUVX(iPUVX) = PUVX(iPUVX) +
     &                  2.0d0/dTot*(5.0d0*Ab1*Ze_int**4.0d0
     &                + 4.0d0*Bb1*Ze_int**3.0d0
     &                + 3.0d0*Cb1*Ze_int**2.0d0)
     &                * (dF_dRho(1,iGrid)- dF_dRho(2,iGrid))
     &                * TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                  TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
     &                 * Weights(iGrid)*nIrrep

                      else
                        PUVX(iPUVX) = PUVX(iPUVX) + 0.0D0
                      end if
                      else
                     PUVX(iPUVX) = PUVX(iPUVX) + 0.0D0
                     end if
                      End Do
                    End Do
                  End Do
                End Do
              End Do
            End If

          End Do
        End Do
      End Do
*

************************************************************************
* BUILD PUVX_TMP (needed to construct FI)
************************************************************************
!We need integrals (kk|pq) and (pk|qk) where k is inact, pq are
!occupied.  p and q must be from the same irrep.
      CALL GETMEM('PUVX_TMP','ALLO','REAL',iTMPP,count_tmp)
      Call Get_dArray('TEP_I',Work(iTMPP),count_tmp)
!AMS - I don't actually need this for the gradients.  It may be useful
!for CI-MC-PDFT calculations, though.

!read the intermediate values from the runfile:


      iPUVX = -1
      do isym=1,nsym
        do jsym=1,nsym
          do k=1+off_bas(isym),nish(isym)+off_bas(isym)
            !Case 1: p and q are inactive
            do p=1+off_bas(jsym),nIsh(jsym)+off_bas(jsym)
              do q=1+off_bas(jsym),nIsh(jsym)+off_bas(jsym)
              iPUVX = iPUVX + 1
                Do iGrid=1,nCoor
                   dTot=Rho(1,iGrid)+Rho(2,iGrid)
                   ratio = 0.0d0
                  if(P2_ontop(1,iGrid).eq.0.0D0)then
                  Work(iTMPP+iPUVX) = Work(iTMPP+iPUVX) + 0.0D0
                  else if(dTot.ge.thrsrho.and.
     &              P2_ontop(1,iGrid).ge.thrsrho) then
                    ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                    if((1.0d0-ratio).gt.thrsrho2) then
                      Zeta  = sqrt(1.0d0-ratio)
                       Work(iTMPP+iPUVX) = Work(iTMPP+iPUVX) +
     &                       (1.0D0/(Zeta*dTot))
     &                       *TabMO(1,iGrid,k)*TabMO(1,iGrid,k)*
     &                       TabMO(1,iGrid,p)*TabMO(1,iGrid,q)*
     &                       Weights(iGrid)*
     &                       (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                    else
                       Work(iTMPP+iPUVX) = Work(iTMPP+iPUVX) +0.0D0
                    end if
                  end if
                end do
              end do
            end do
            !Case 2: p is inactive, q is active
            do p=1+off_bas(jsym),nIsh(jsym)+off_bas(jsym)
              do q=1+off_basash(jsym),nAsh(jsym)+off_basash(jsym)
              iPUVX = iPUVX + 1
                Do iGrid=1,nCoor
                   dTot=Rho(1,iGrid)+Rho(2,iGrid)
                   ratio = 0.0d0
                  if(P2_ontop(1,iGrid).eq.0.0D0)then
                  Work(iTMPP+iPUVX) = Work(iTMPP+iPUVX) + 0.0D0
                  else if(dTot.ge.thrsrho.and.
     &              P2_ontop(1,iGrid).ge.thrsrho) then
                    ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                    if((1.0d0-ratio).gt.thrsrho2) then
                      Zeta  = sqrt(1.0d0-ratio)
                       Work(iTMPP+iPUVX) = Work(iTMPP+iPUVX) +
     &                       (1.0D0/(Zeta*dTot))
     &                       *TabMO(1,iGrid,k)*TabMO(1,iGrid,k)*
     &                       TabMO(1,iGrid,p)*TabMO(1,iGrid,q)*
     &                       Weights(iGrid)*
     &                       (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                    else
                       Work(iTMPP+iPUVX) = Work(iTMPP+iPUVX) +0.0D0
                    end if
                  end if
                end do
              end do
            end do
            !Case 3: p and q are active
            do p=1+off_basash(jsym),nAsh(jsym)+off_basash(jsym)
              do q=1+off_basash(jsym),nAsh(jsym)+off_basash(jsym)
              iPUVX = iPUVX + 1
                Do iGrid=1,nCoor
                   dTot=Rho(1,iGrid)+Rho(2,iGrid)
                   ratio = 0.0d0
                  if(P2_ontop(1,iGrid).eq.0.0D0)then
                  Work(iTMPP+iPUVX) = Work(iTMPP+iPUVX) + 0.0D0
                  else if(dTot.ge.thrsrho.and.
     &              P2_ontop(1,iGrid).ge.thrsrho) then
                    ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                    if((1.0d0-ratio).gt.thrsrho2) then
                      Zeta  = sqrt(1.0d0-ratio)
                       Work(iTMPP+iPUVX) = Work(iTMPP+iPUVX) +
     &                       (1.0D0/(Zeta*dTot))
     &                       *TabMO(1,iGrid,k)*TabMO(1,iGrid,k)*
     &                       TabMO(1,iGrid,p)*TabMO(1,iGrid,q)*
     &                       Weights(iGrid)*
     &                       (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                    else
                       Work(iTMPP+iPUVX) = Work(iTMPP+iPUVX) +0.0D0
                    end if
                  end if
                end do
              end do
            end do
          end do
        end do
      end do

!      Call Put_dArray('TEP_I',Work(iTMPP),count_tmp)
      CALL GETMEM('PUVX_TMP','Free','REAL',iTMPP,count_tmp)
      lsym=lsym_tmp

!OK - potential modifications - for the iTMPP terms calculated above,
!maybe I really only want the third case - the active-active case.
!These are the only terms that need to be considered for the CI fock
!part.

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
!I think I can generate the contributions from the potentials to the new
!FI and FA terms at this level.  Then there will be no need to store a
!large number of V_TE potentials.
!      CALL GETMEM('FI_V','ALLO','REAL',ifiv,ntot1)
!      CALL GETMEM('FA_V','ALLO','REAL',ifav,ntot1)
!      Call Get_dArray('FI_V',Work(ifiv),ntot1)
!      Call Get_dArray('FA_V',Work(ifav),ntot1)

      Do iSym = 1,nSym !sym for p
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iPUVX = off_PUVX(iSym)
        Do jSym = 1,nSym !sym for u
          jOrb = nOrb(jSym)
          jAsh = nAsh(jSym)
          jIsh = nIsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym !Sym for v
            kAsh = nAsh(kSym)
            kIsh = nIsh(kSym)
            Do lSym = 1,kSym !sym for x
              lAsh = nAsh(lSym)
              lIsh = nIsh(lSym)
              klSym = 1 + ieor(kSym-1,lSym-1)

*             find cases
              case = 4
              If ( iSym.eq.jSym ) case = case-2
              If ( iSym.eq.kSym ) case = case-1

              If ( ijSym.eq.klSym .and.
     &             iAsh*jAsh*kAsh*lAsh.ne.0 ) then

                Goto (100,200,300,400) case

!Symmetry case (II|II)
100             Continue
                iFoff = off_Fmat(iSym)
                iDoff = off_Dmat(iSym)
                Do iV = 1,kIsh
                  jV = iV + off_basIsh(kSym)
                  iX=iV
                  jX=jV
                  DVX=2.0D0
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      !iPU  = iTri(jU) + jP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      !iPU  = iTri(jU) + jP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      !iPU  = iTri(jU) + jP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      !iPU  = iTri(jU) + jP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho) then
!                        if(dTot.ge.thrsrho.and.
!     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      !iPU  = iTri(jU) + jP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                end do !V
!ACTIVE CONTRIBUTIONS
                Do iV = 1,kAsh
                  jV = iV + off_basAsh(ksym)
                do iX=1,iV
                    jX = iX + off_basAsh(ksym)
                    iVX = iTri(iV) + iX
                    DVX=2.0D0*D1MO(iDoff+iVX)
                    If (iX.eq.iV) DVX=D1MO(iDoff+iVX)
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iU + off_basAsh(jSym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(iSym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho) then
!                        if(dTot.ge.thrsrho.and.
!     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                end do !X
                end do !V
                Goto 500
!Symmetry case (II|KK)
200             Continue
                iFoff = off_Fmat(iSym)
                iDoff = off_Dmat(iSym)
                Do iV = 1,kIsh
                  jV = iV + off_basIsh(ksym)
                  iX=iV
                  jX=jV
                  DVX=2.0D0
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+nAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                end do !V
!ACTIVE CONTRIBUTIONS
                Do iV = 1,kAsh
                  jV = iV + off_basAsh(ksym)
                  do iX=1,iV
                    jX = iX + off_basAsh(ksym)
                    iVX = iTri(iV) + iX
                    DVX=2.0D0*D1MO(iDoff+iVX)
                    If (iX.eq.iV) DVX=D1MO(iDoff+iVX)
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+nAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*nIrrep*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) + DVX*V_PUVX
                    end do
                  end do
                end do !X
                end do !V
                Goto 500
*               symmetry case (IJ!IJ)
300             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                   ! off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                   ! off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do
                Goto 500
*               symmetry case (IJ!KL)
400             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                  !  off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                  !  off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do
                Goto 500

500             Continue
              End If

            End Do
          End Do
        End Do
      End Do
!      Call Put_dArray('FI_V',Work(ifiv),ntot1)
!      Call Put_dArray('FA_V',Work(ifav),ntot1)
!      CALL GETMEM('FI_V','FREE','REAL',ifiv,ntot1)
!      CALL GETMEM('FA_V','FREE','REAL',ifav,ntot1)
      lsym=lsym_tmp
      Return
      End

      Subroutine Calc_OTPUVX_FTLSDA2(PUVX,TabMO,mAO,nCoor,nTabMOs
     &                       ,P2_ontop,nP2_ontop,Rho,nRho,
     &                    dF_dRho,ndF_dRho,RhoI,RhoA,mRho,Weights,
     &                    D1MO,nD1MO,nIrrep)
      Implicit Real*8 (A-H,O-Z)
      Dimension PUVX(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_PUVX(mxSym),off_Bas(mxSym),
     &        off_ish(mxSym),off_BasIsh(mxSym),
     &        off_BasVsh(mxSym)
      Integer   off_Dmat, off_Fmat
      Dimension off_Dmat(mxSym), off_Fmat(mxSym)


      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),P2_ontop(nP2_ontop,nCoor),
     &       dF_dRho(ndF_dRho,nCoor),Rho(nRho,nCoor),
     &       RhoI(mRho,nCoor),RhoA(mRho,nCoor)
      Integer nIrrep
      Real*8 D1MO(nD1MO)
      Real*8 DVX
      real*8 Fact
      Integer case
      iTrii(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      iTri(i)=(i*i-i)/2
*
      thrsrho=1.0d-15
      thrsrho2=1.0d-15
      thrspi=1.0d-30
      thrsCEH=1.0d-24
      thrsrho3=0.9000000000d0
      thrsrho4=1.1500000000d0
      Ab1=-4.756065601d+2
      Bb1=-3.794733192d+2
      Cb1=-8.538149682d+1
*
      Call Unused_real_array(RhoI)
      Call Unused_real_array(RhoA)
      lsym_tmp=lsym

      iStack  = 0
      iStack1 = 0
      iStack2 = 0
      off_ish(:) = 0
      off_Ash(:) = 0
      off_Bas(:) = 0
      off_BasAsh(:) = 0
      off_BasIsh(:) = 0
      off_BasVsh(:) = 0
      ntot1 = 0

      Do iSym = 1,nSym
        off_ish(isym)    = iStack2
        off_Ash(iSym)    = iStack
        off_Bas(iSym)    = iStack1
        off_BasVsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)+nAsh(iSym)
        off_BasAsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)
        off_BasIsh(iSym) = iStack1+nFro(iSym)
        ntot1 = iTrii(nBas(iSym),nBas(iSym)) + ntot1
        iStack2 = iStack2 + nIsh(iSym)
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
              lAsh = nAsh(lsym)
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
*
      Fact = 1.000d0
        Do iSym = 1,nSym
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iPUVX = off_PUVX(iSym)
        Do jSym = 1,nSym
          jAsh = nAsh(jSym)
          jIsh = nIsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym
            kAsh = nAsh(kSym)
!            kIsh = nIsh(kSym)
            lSym = 1 + ieor(ijSym-1,kSym-1)
            lAsh = nAsh(lSym)
!            lIsh = nIsh(lSym)

            If ( lSym.le.kSym .and.
     &           iAsh*jAsh*kAsh*lAsh.ne.0 ) then
!              Do iV = kIsh+1,kIsh+kAsh
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

                    Do iGrid=1,nCoor
                       dTot=Rho(1,iGrid)+Rho(2,iGrid)
                     ratio = 0.0d0
                  if(dTot.ge.thrsrho.and.
     &               P2_ontop(1,iGrid).ge.thrspi) then
                    ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                    if(((1.0d0-ratio).gt.thrsrho2)
     &                 .and.(ratio.lt.thrsrho3)) then
                      Zeta  = sqrt(1.0d0-ratio)
                      PUVX(iPUVX) = PUVX(iPUVX) +
     &                    (1.0D0/(Zeta*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*dble(nIrrep)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                     elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      PUVX(iPUVX) = PUVX(iPUVX) +
     &                    (1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*dble(nIrrep)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                     else
                      PUVX(iPUVX) = PUVX(iPUVX) +0.0D0
                     end if
                     !end standard translation
                   end if

                      End Do
                    End Do
                  End Do
                End Do
              End Do
            End If

          End Do
        End Do
      End Do

!Construction of Fock matrix pieces
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
!I think I can generate the contributions from the potentials to the new
!FI and FA terms at this level.  Then there will be no need to store a
!large number of V_TE potentials.
!      CALL GETMEM('FI_V','ALLO','REAL',ifiv,ntot1)
!      CALL GETMEM('FA_V','ALLO','REAL',ifav,ntot1)
!      Call Get_dArray('FI_V',Work(ifiv),ntot1)
!      Call Get_dArray('FA_V',Work(ifav),ntot1)

      Do iSym = 1,nSym !sym for p
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iPUVX = off_PUVX(iSym)
        Do jSym = 1,nSym !sym for u
          jOrb = nOrb(jSym)
          jAsh = nAsh(jSym)
          jIsh = nIsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym !Sym for v
            korb = nOrb(kSym)
            kAsh = nAsh(kSym)
            kIsh = nIsh(kSym)
            Do lSym = 1,kSym !sym for x
              lOrb = nOrb(lSym)
              lAsh = nAsh(lSym)
              lIsh = nIsh(lSym)
              klSym = 1 + ieor(kSym-1,lSym-1)

*             find cases
              case = 4
              If ( iSym.eq.jSym ) case = case-2
              If ( iSym.eq.kSym ) case = case-1

              If ( ijSym.eq.klSym .and.
!     &             iAsh*jAsh*kAsh*lAsh.ne.0 ) then
     &             iOrb*jOrb*kOrb*lOrb.ne.0 ) then

                Goto (100,200,300,400) case

!Symmetry case (II|II)
100             Continue
                iFoff = off_Fmat(iSym)
                iDoff = off_Dmat(iSym)
                Do iV = 1,kIsh
                  jV = iV + off_basIsh(kSym)
                  iX=iV
                  jX=jV
                  DVX=2.0D0
                  If ( iX.eq.iV ) DVX = DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                     P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                       .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                        elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                        else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                       .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                        V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                     P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                       .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                   Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                   Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !V
!ACTIVE CONTRIBUTIONS
                Do iV = 1,kAsh
                  jV = iV + off_basAsh(ksym)
                do iX=1,iV
                    jX = iX + off_basAsh(ksym)
                    !iVX = iTri(iV) + iX
                    !jVX = iTri(jV) + jX
                    iVX = iTri(iV+off_Ash(ksym)) + iX + off_Ash(ksym)
                    DVX=D1MO(iVX)
                    !DVX=D1MO(iDoff+iVX)!Why *2?
                  If ( iX.eq.iV ) DVX = DVX*0.5
          !          DVX=2.0D0*D1MO(iDoff+iVX)!Why *2?
          !          If (iX.eq.iV) DVX=D1MO(iDoff+iVX)
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iU + off_basAsh(jSym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(iSym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                       .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                   !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !X
                end do !V
                Goto 500
!Symmetry case (II|KK)
200             Continue
                iFoff = off_Fmat(iSym)
                kDoff = off_Dmat(kSym)
                Do iV = 1,kIsh
                  jV = iV + off_basIsh(ksym)
                  iX=iV
                  jX=jV
                  DVX=2.0D0
                  If ( iX.eq.iV ) DVX = DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+nAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                       .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                     !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                   end do
                end do !V
!ACTIVE CONTRIBUTIONS
                Do iV = 1,kAsh
                  jV = iV + off_basAsh(ksym)
                  do iX=1,iV
                    jX = iX + off_basAsh(ksym)
                    iVX = iTri(iV+off_Ash(ksym)) + iX + off_Ash(ksym)
                    DVX=D1MO(iVX)
                    !DVX=D1MO(kDoff+iVX)
                  If ( iX.eq.iV ) DVX = DVX*0.5
            !        DVX=2.0D0*D1MO(iDoff+iVX)
            !        If (iX.eq.iV) DVX=D1MO(iDoff+iVX)
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                       .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+nAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                        end do
                      end do
                    !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrspi) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if(((1.0d0-ratio).gt.thrsrho2)
     &                        .and.(ratio.lt.thrsrho3)) then
                            Zeta  = sqrt(1.0d0-ratio)
                            V_PUVX = V_PUVX +
     &                      Fact*(1.0D0/(Zeta*dTot))
     &                      *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                      TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                      Weights(iGrid)*
     &                      (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                      V_PUVX = V_PUVX +
     &                    Fact*(1.0D0/(2.0D0*dTot))
     &                   *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   (20.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                    +16.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                    +12.0D0*Cb1*(ratio-1.15d0)**2.0d0)*
     &                   Weights(iGrid)*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !X
                end do !V
                Goto 500
*               symmetry case (IJ!IJ)
300             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                   ! off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                   ! off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do
                Goto 500
*               symmetry case (IJ!KL)
400             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                  !  off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                  !  off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do
                Goto 500

500             Continue
              End If

            End Do
          End Do
        End Do
      End Do

!      Call Put_dArray('FI_V',Work(ifiv),ntot1)
!      Call Put_dArray('FA_V',Work(ifav),ntot1)
!      CALL GETMEM('FI_V','FREE','REAL',ifiv,ntot1)
!      CALL GETMEM('FA_V','FREE','REAL',ifav,ntot1)
*
      lsym=lsym_tmp
*
      Return
      End

      Subroutine Calc_OTOE_FTLSDA(OE,TabMO,mAO,nCoor,nTabMOs
     &                       ,P2_ontop,nP2_ontop,Rho,nRho,
     &                      dF_dRho,ndF_dRho,RhoI,RhoA,mRho,Weights,
     &                      nIrrep)
      Implicit Real*8 (A-H,O-Z)
      Dimension OE(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_PUVX(mxSym),off_Bas(mxSym)!,off_Fmat(mxSym),
!     &        off_Dmat(mxSym)
      Integer nIrrep
      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),P2_ontop(nP2_ontop,nCoor),
     &       dF_dRho(ndF_dRho,nCoor),Rho(nRho,nCoor),
     &       RhoI(mRho,nCoor),RhoA(mRho,nCoor)
*
      !iTri(i) = (i*i-i)/2
       iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      thrsrho=1.0d-15
      thrsrho2=1.0d-15
      thrsrho3=0.9000000000d0
      thrsrho4=1.1500000000d0
      thrspi=1.0d-30
      thrsCEH=1.0d-24
      Ab1=-4.756065601d+2
      Bb1=-3.794733192d+2
      Cb1=-8.538149682d+1
*
      Call Unused_real_array(RhoI)
      Call Unused_real_array(RhoA)
      lsym_tmp=lsym
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

      Do iSym = 1,nSym
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iPUVX = off_PUVX(iSym)
        Do jSym = 1,nSym
          jAsh = nAsh(jSym)
          jIsh = nIsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym
            kAsh = nAsh(kSym)
            kIsh = nIsh(kSym)
            lSym = 1 + ieor(ijSym-1,kSym-1)
            lAsh = nAsh(lSym)
            lIsh = nIsh(lSym)

            If ( lSym.le.kSym .and.
     &           iAsh*jAsh*kAsh*lAsh.ne.0 ) then
              Do iV = 1,kIsh+kAsh
                jV = iV + off_BasAsh(kSym)
                lMax = lIsh+lAsh
                If ( kSym.eq.lSym ) lMax = iV
                Do iX = 1,lMax
                  jX = iX + off_BasAsh(lSym)

                  iVX=iTri(iV + off_BasAsh(kSym) ,
     *                    iX + off_BasAsh(lSym) )

                  Do iGrid = 1, nCoor
*CEH -OE addition
                    dTot=Rho(1,iGrid)+Rho(2,iGrid)
                    ratio = 0.0d0
                     if(dTot.ge.thrsrho.and.
     &                 P2_ontop(1,iGrid).ge.thrspi) then
                      ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                    if(((1.0d0-ratio).gt.thrsrho2)
     &                 .and.(ratio.lt.thrsrho3)) then
                       Zeta  = sqrt(1.0d0-ratio)
                     OE(iVX)=OE(iVX)+TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                    *(dF_dRho(1,iGrid)*(0.5D0*(1.0D0+Zeta)
     &                    +2.0D0*P2_ontop(1,iGrid)/(Zeta*dTot**2))
     &                    +dF_dRho(2,iGrid)*(0.5D0*(1.0D0-Zeta)
     &                    -2.0D0*P2_ontop(1,iGrid)/(Zeta*dTot**2)))
     &                    *Weights(iGrid)
                     elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (10.0D0*Ab1*(ratio-1.15d0)**4.0d0) +
     &                        (8.0D0*Bb1*(ratio-1.15d0)**3.0d0) +
     &                        (6.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                     OE(iVX)=OE(iVX)+TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                    *(dF_dRho(1,iGrid)*(0.5D0*(1.0D0+Zeta)
     &                    -2.0D0*P2_ontop(1,iGrid)/(dTot**2)*Deriv)
     &                    +dF_dRho(2,iGrid)*(0.5D0*(1.0D0-Zeta)
     &                    +2.0D0*P2_ontop(1,iGrid)/(dTot**2)*Deriv))
     &                    *Weights(iGrid)*nIrrep
                       else
                         OE(iVX)=OE(iVX)+TabMO(1,iGrid,jV)*
     &                     TabMO(1,iGrid,jX)*(dF_dRho(1,iGrid)
     &                     +dF_dRho(2,iGrid))*Weights(iGrid)*0.5D0
     &                     *nIrrep
                       end if
                      end if
*End CEH -OE addition
                  End Do     ! iGrid


                  Do iU = 1,jIsh+jAsh
                  jU = iU + off_BasAsh(jSym)
                    Do iP = 1,iOrb
                      iT = iP - iIsh
                      iPUVX=iPUVX+1
                      jP = iP     +   off_Bas(iSym)

                    End Do
                  End Do
                End Do
              End Do
            End If

          End Do
        End Do
      End Do
           lsym=lsym_tmp

      Return
      End
      Subroutine Calc_OTOE_ft(OE,TabMO,mAO,nCoor,nTabMOs
     &                       ,P2_ontop,nP2_ontop,Rho,nRho,
     &                      dF_dRho,ndF_dRho,RhoI,RhoA,mRho,Weights,
     &                      nIrrep)
      Implicit Real*8 (A-H,O-Z)
      Dimension OE(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_Bas(mxSym)!,off_Fmat(mxSym),
!     &        off_Dmat(mxSym)
      Integer off_basIsh(mxSym),off_Fmat(mxSym)
      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),P2_ontop(nP2_ontop,nCoor),
     &       dF_dRho(ndF_dRho,nCoor),Rho(nRho,nCoor),
     &       RhoI(mRho,nCoor),RhoA(mRho,nCoor)
      Integer nIrrep
      Integer VX,ix,jx,iv,jv
      real*8 zeta,zeta_in
*
      !iTri(i) = (i*i-i)/2
       iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)

      thrsrho=1.0d-15
      thrsrho2=1.0d-15
      thrsrho3=0.9000000000d0
      thrsrho4=1.1500000000d0
      thrspi=1.0d-30
      Ab1=-4.756065601d+2
      Bb1=-3.794733192d+2
      Cb1=-8.538149682d+1
*
      Call Unused_real_array(RhoI)
      Call Unused_real_array(RhoA)
      lsym_tmp=lsym
      iStack  = 0
      iStack1 = 0
      Do iSym = 1,nSym
        off_Ash(iSym)    = iStack
        off_Bas(iSym)    = iStack1
        off_BasAsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)
        off_BasIsh(iSym) = iStack1+nFro(iSym)
        iStack1 = iStack1 + nBas(iSym)
        iStack  = iStack  + nAsh(iSym)
      End Do

      iStack = 0
      do iSym=1,nSym
        off_Fmat(iSym) = iStack
        iorb = nOrb(iSym)
        iStack = iStack + (iOrb*iOrb + iOrb)/2
      end do
*

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!OE pieces - both orbs must be in the same irrep, eh?
      Do iSym = 1,nSym
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        Do iV = 1,iOrb!iIsh+iAsh
          jV = iV + off_BasIsh(iSym)
          do iX = 1,iV
            jX = iX + off_BasIsh(iSym)
            VX = off_Fmat(iSym) + iTri(iV,iX)
            Do iGrid = 1, nCoor
*CEH -OE addition
            dTot=Rho(1,iGrid)+Rho(2,iGrid)
            ratio = 0.0d0
                if (dTot.ge.thrsrho.and.
     &               P2_ontop(1,iGrid).ge.thrspi) then
                ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  if(((1.0d0-ratio).gt.thrsrho2).and.
     &            (ratio.lt.thrsrho3)) then
                  Zeta  = sqrt(1.0d0-ratio)
                  OE(VX)=OE(VX)+TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                    *(dF_dRho(1,iGrid)*(0.5D0*(1.0d0+Zeta)
     &                    +2.0D0*P2_ontop(1,iGrid)/(Zeta*dTot**2))
     &                    +dF_dRho(2,iGrid)*(0.5D0*(1.0d0-Zeta)
     &                    -2.0D0*P2_ontop(1,iGrid)/(Zeta*dTot**2)))
     &                    *Weights(iGrid)*nIrrep
                   else if((ratio.ge.thrsrho3).and.
     &                     (ratio.le.thrsrho4)) then
                     Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &         (Bb1*(ratio-1.15d0)**4.0d0) + (Cb1*(ratio-1.15d0)**3.0d0)
                     Zeta_in = 5.0d0*Ab1*(ratio - 1.15d0)**4.0d0
     &               + 4.0d0*Bb1*(ratio - 1.15d0)**3.0d0
     &               + 3.0d0*Cb1*(ratio - 1.15d0)**2.0d0
                     OE(VX) = OE(VX)
     &               + TabMO(1,iGrid,jV) * TabMO(1,iGrid,jX)
     &               *(dF_dRho(1,iGrid) * ((1+Zeta)/2 -Ratio*(Zeta_in))
     &               + dF_dRho(2,iGrid) * ((1-Zeta)/2 -Ratio*(Zeta_in)))
     &               * Weights(iGrid)*nIrrep
                   else
                    OE(VX)=OE(VX)+TabMO(1,iGrid,jV)*
     &                     TabMO(1,iGrid,jX)*(dF_dRho(1,iGrid)
     &                     +dF_dRho(2,iGrid))*Weights(iGrid)*0.5d0

                   end if
                end if

                !end if
!              end if
*End CEH -OE addition
            End Do     ! iGrid
          End do
        End Do
      End Do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           lsym=lsym_tmp

      Return
      End



      Subroutine Calc_OTPUVXGGA(PUVX,TabMO,mAO,nCoor,nTabMOs
     &                       ,P2_ontop,nP2_ontop,Rho,nRho,
     &                      dF_dRho,ndF_dRho,RhoI,RhoA,mRho,
     &                      Weights,D1MO,nD1MO,nIrrep)
      Implicit Real*8 (A-H,O-Z)
      Dimension PUVX(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer nIrrep
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_PUVX(mxSym),off_Bas(mxSym),
     &        off_ish(mxSym),off_BasIsh(mxSym),
     &        off_BasVsh(mxSym)
      Integer   off_Dmat, off_Fmat
      Dimension off_Dmat(mxSym), off_Fmat(mxSym)

      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),P2_ontop(nP2_ontop,nCoor),
     &       dF_dRho(ndF_dRho,nCoor),Rho(nRho,nCoor),
     &       RhoI(mRho,nCoor),RhoA(mRho,nCoor)
      Real*8 D1MO(nD1MO)
      Integer count_tmp
      Real*8 DVX
      real*8 Fact
      Real*8 time1,time2
      Integer case
      iTrii(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      iTri(i)=(i*i-i)/2
*
      thrsrho=1.0d-15
      thrsrho2=1.0d-15
      thrspi=1.0d-30
*
      lsym_tmp=lsym

      Call Unused_real_array(RhoI)
      Call Unused_real_array(RhoA)

      Call CPU_Time(time1)

      iStack  = 0
      iStack1 = 0
      iStack2 = 0
      off_ish(:) = 0
      off_Ash(:) = 0
      off_Bas(:) = 0
      off_BasAsh(:) = 0
      off_BasIsh(:) = 0
      off_BasVsh(:) = 0
      ntot1 = 0

      Do iSym = 1,nSym
        off_ish(isym)    = iStack2
        off_Ash(iSym)    = iStack
        off_Bas(iSym)    = iStack1
        off_BasVsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)+nAsh(iSym)
        off_BasAsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)
        off_BasIsh(iSym) = iStack1+nFro(iSym)
        ntot1 = iTrii(nBas(iSym),nBas(iSym)) + ntot1
        iStack2 = iStack2 + nIsh(iSym)
        iStack1 = iStack1 + nBas(iSym)
        iStack  = iStack  + nAsh(iSym)
      End Do



!count the number of tmp_pot:
      count_tmp = 0
      do isym=1,nsym
        do jsym=1,nsym
        count_tmp = count_tmp + nIsh(isym)*(nIsh(jsym)+nAsh(jsym))**2
        end do
      end do

*
!Calculate PUVX offsets
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
*

************************************************************************
* BUILD PUVX potentials
!Note - I really only need the TUVX potentials.  Area for speedup.
************************************************************************
      Fact = 1.00d0
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
!                        MO1=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
                  Do iU = 1,jAsh
                  jU = iU + off_BasAsh(jSym)
                    Do iP = 1,iOrb
                      iT = iP - iIsh
                      iPUVX=iPUVX+1
                      jP = iP     +   off_Bas(iSym)
                    Do iGrid=1,nCoor
                       dTot=Rho(1,iGrid)+Rho(2,iGrid)
                     ratio = 0.0d0
                  if(dTot.ge.thrsrho.and.
     &               P2_ontop(1,iGrid).ge.thrspi) then
              ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                      if((1.0d0-ratio).gt.thrsrho2) then
                      Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
!                      interm1 = 1/(zeta*dTot)
!                      interm2 = interm1/dTot
                      PUVX(iPUVX) = PUVX(iPUVX) +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*dble(nIrrep)*(
     &                   1.0D0/(Zeta*dTot)*
!     &                   interm1*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
!     &                    +RHOPx/interm2
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
!     &                    +RHOPx/interm2
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
!     &                    +RHOPy/interm2
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
!     &                    +RHOPy/interm2
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
!     &                    +RHOPz/interm2
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
!     &                    +RHOPz/interm2
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      else
                     PUVX(iPUVX) = PUVX(iPUVX) +0.0D0
                     end if
                   end if
                      End Do
                    End Do
                  End Do
                End Do
              End Do
            End If
          End Do
        End Do
      End Do

      Call CPU_Time(time2)
      PUVX_Time = PUVX_time + (time2-time1)
!Construction of Fock matrix pieces
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
!I think I can generate the contributions from the potentials to the new
!FI and FA terms at this level.  Then there will be no need to store a
!large number of V_TE potentials.
!      CALL GETMEM('FI_V','ALLO','REAL',ifiv,ntot1)
!      CALL GETMEM('FA_V','ALLO','REAL',ifav,ntot1)
!      Call Get_dArray('FI_V',Work(ifiv),ntot1)
!      Call Get_dArray('FA_V',Work(ifav),ntot1)

      Do iSym = 1,nSym !sym for p
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iPUVX = off_PUVX(iSym)
        Do jSym = 1,nSym !sym for u
          jOrb = nOrb(jSym)
          jAsh = nAsh(jSym)
          jIsh = nIsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym !Sym for v
            kOrb = nOrb(kSym)
            kAsh = nAsh(kSym)
            kIsh = nIsh(kSym)
            Do lSym = 1,kSym !sym for x
              lOrb = nOrb(lSym)
              lAsh = nAsh(lSym)
              lIsh = nIsh(lSym)
              klSym = 1 + ieor(kSym-1,lSym-1)

*             find cases
              case = 4
              If ( iSym.eq.jSym ) case = case-2
              If ( iSym.eq.kSym ) case = case-1

              If ( ijSym.eq.klSym .and.
!     &             iAsh*jAsh*kAsh*lAsh.ne.0 ) then
     &             iOrb*jOrb*kOrb*lOrb.ne.0 ) then

                Goto (100,200,300,400) case

!Symmetry case (II|II)
100             Continue
                iFoff = off_Fmat(iSym)
                iDoff = off_Dmat(iSym)
                Do iV = 1,kIsh
                  jV = iV + off_basIsh(kSym)
                  iX=iV
                  jX=jV
                  DVX=2.0D0
                  If ( iX.eq.iV ) DVX = DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !V

      Call CPU_Time(time1)
      FI_time = FI_Time + (time1-time2)
!ACTIVE CONTRIBUTIONS
                Do iV = 1,kAsh
                  jV = iV + off_basAsh(ksym)
                do iX=1,iV
                    jX = iX + off_basAsh(ksym)
                    iVX = iTri(iV+off_Ash(ksym)) + iX + off_Ash(ksym)
                    DVX=D1MO(iVX)
                    !iVX = iTri(iV) + iX
                    !DVX=2.0D0*D1MO(iDoff+iVX)!Why *2?
                    If (iX.eq.iV) DVX=DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iU + off_basAsh(jSym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(iSym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !X
                end do !V
      Call CPU_Time(time2)
      FA_time = FA_Time + (time2-time1)
                Goto 500
!Symmetry case (II|KK)
200             Continue
                iFoff = off_Fmat(iSym)
                iDoff = off_Dmat(kSym)
                Do iV = 1,kIsh
                  jV = iV + off_basIsh(ksym)
                  iX=iV
                  jX=jV
                  DVX=2.0D0
                  If ( iX.eq.iV ) DVX = DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+nAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !V
      Call CPU_Time(time1)
      FI_time = FI_Time + (time1-time2)
!ACTIVE CONTRIBUTIONS
                Do iV = 1,kAsh
                  jV = iV + off_basAsh(ksym)
                  do iX=1,iV
                    jX = iX + off_basAsh(ksym)
                    iVX = iTri(iV+off_Ash(ksym)) + iX + off_Ash(ksym)
                    DVX=D1MO(iVX)
                    !iVX = iTri(iV) + iX
                    !DVX=2.0D0*D1MO(iDoff+iVX)
                    If (iX.eq.iV) DVX=DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+nAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP + iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !X
                end do !V
      Call CPU_Time(time2)
      FA_time = FA_Time + (time2-time1)
                Goto 500
*               symmetry case (IJ!IJ)
300             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                   ! off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                   ! off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do
                Goto 500
*               symmetry case (IJ!KL)
400             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                  !  off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                  !  off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do
                Goto 500

500             Continue
              End If

            End Do
          End Do
        End Do
      End Do

!      Call Put_dArray('FI_V',Work(ifiv),ntot1)
!      Call Put_dArray('FA_V',Work(ifav),ntot1)
!      CALL GETMEM('FI_V','FREE','REAL',ifiv,ntot1)
!      CALL GETMEM('FA_V','FREE','REAL',ifav,ntot1)
      lsym=lsym_tmp
*
      Return
      End

      Subroutine Calc_OTPUVXGGA_ft(PUVX,TabMO,mAO,nCoor,nTabMOs
     &                       ,P2_ontop,nP2_ontop,Rho,nRho,
     &                      dF_dRho,ndF_dRho,RhoI,RhoA,mRho,
     &                      Weights,D1MO,nD1MO,nIrrep)
      Implicit Real*8 (A-H,O-Z)
      Dimension PUVX(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer nIrrep
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_PUVX(mxSym),off_Bas(mxSym),
     &        off_ish(mxSym),off_BasIsh(mxSym),
     &        off_BasVsh(mxSym)
      Integer   off_Dmat, off_Fmat
      Dimension off_Dmat(mxSym), off_Fmat(mxSym)

      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),P2_ontop(nP2_ontop,nCoor),
     &       dF_dRho(ndF_dRho,nCoor),Rho(nRho,nCoor),
     &       RhoI(mRho,nCoor),RhoA(mRho,nCoor)
      Real*8 D1MO(nD1MO)
      real*8 V_puvx
      Integer count_tmp
      Real*8 DVX
      real*8 Fact
      Integer case
      iTrii(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      iTri(i)=(i*i-i)/2
*
      thrsrho=1.0d-15
      thrsrho2=1.0d-15
      thrspi=1.0d-30
      thrsCEH=1.0d-24
      thrsrho3=0.9000000000d0
      thrsrho4=1.1500000000d0
      Ab1=-4.756065601d+2
      Bb1=-3.794733192d+2
      Cb1=-8.538149682d+1
*
      lsym_tmp=lsym

      Call Unused_real_array(RhoI)
      Call Unused_real_array(RhoA)

      iStack  = 0
      iStack1 = 0
      iStack2 = 0
      off_ish(:) = 0
      off_Ash(:) = 0
      off_Bas(:) = 0
      off_BasAsh(:) = 0
      off_BasIsh(:) = 0
      off_BasVsh(:) = 0
      ntot1 = 0

      Do iSym = 1,nSym
        off_ish(isym)    = iStack2
        off_Ash(iSym)    = iStack
        off_Bas(iSym)    = iStack1
        off_BasVsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)+nAsh(iSym)
        off_BasAsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)
        off_BasIsh(iSym) = iStack1+nFro(iSym)
        ntot1 = iTrii(nBas(iSym),nBas(iSym)) + ntot1
        iStack2 = iStack2 + nIsh(iSym)
        iStack1 = iStack1 + nBas(iSym)
        iStack  = iStack  + nAsh(iSym)
      End Do



!count the number of tmp_pot:
      count_tmp = 0
      do isym=1,nsym
        do jsym=1,nsym
        count_tmp = count_tmp + nIsh(isym)*(nIsh(jsym)+nAsh(jsym))**2
        end do
      end do

*
!Calculate PUVX offsets
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
*

************************************************************************
* BUILD PUVX potentials
!Note - I really only need the TUVX potentials.  Area for speedup.
************************************************************************
      Fact = 1.00d0
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
                    Do iGrid=1,nCoor
                       dTot=Rho(1,iGrid)+Rho(2,iGrid)
                     ratio = 0.0d0
                  if(dTot.ge.thrsrho.and.
     &               P2_ontop(1,iGrid).ge.thrspi) then
              ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                      d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      PUVX(iPUVX) = PUVX(iPUVX) +
     &                   Weights(iGrid)*dble(nIrrep)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                       d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                       d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                    PUVX(iPUVX) = PUVX(iPUVX) +
     &                   Weights(iGrid)*dble(nIrrep)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                     PUVX(iPUVX) = PUVX(iPUVX) +0.0D0
                     end if
                   end if
                      End Do
                    End Do
                  End Do
                End Do
              End Do
            End If
          End Do
        End Do
      End Do

!Construction of Fock matrix pieces
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
!I think I can generate the contributions from the potentials to the new
!FI and FA terms at this level.  Then there will be no need to store a
!large number of V_TE potentials.
!      CALL GETMEM('FI_V','ALLO','REAL',ifiv,ntot1)
!      CALL GETMEM('FA_V','ALLO','REAL',ifav,ntot1)
!      Call Get_dArray('FI_V',Work(ifiv),ntot1)
!      Call Get_dArray('FA_V',Work(ifav),ntot1)

      Do iSym = 1,nSym !sym for p
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iPUVX = off_PUVX(iSym)
        Do jSym = 1,nSym !sym for u
          jOrb = nOrb(jSym)
          jAsh = nAsh(jSym)
          jIsh = nIsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym !Sym for v
            kOrb = nOrb(kSym)
            kAsh = nAsh(kSym)
            kIsh = nIsh(kSym)
            Do lSym = 1,kSym !sym for x
              lOrb = nOrb(lSym)
              lAsh = nAsh(lSym)
              lIsh = nIsh(lSym)
              klSym = 1 + ieor(kSym-1,lSym-1)

*             find cases
              case = 4
              If ( iSym.eq.jSym ) case = case-2
              If ( iSym.eq.kSym ) case = case-1

              If ( ijSym.eq.klSym .and.
!     &             iAsh*jAsh*kAsh*lAsh.ne.0 ) then
     &             iOrb*jOrb*kOrb*lOrb.ne.0 ) then

                Goto (100,200,300,400) case

!Symmetry case (II|II)
100             Continue
                iFoff = off_Fmat(iSym)
                iDoff = off_Dmat(iSym)
                Do iV = 1,kIsh
                  jV = iV + off_basIsh(kSym)
                  iX=iV
                  jX=jV
                  DVX=2.0D0
                  If ( iX.eq.iV ) DVX = DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                 SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                      d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                    V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                       d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                        V_PUVX = V_PUVX + 0.0d0
                      end if
                    end if
                End Do!gridpt
           Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                       d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                       else
                             V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                   SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                    V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                    *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                    else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                    V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                       else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                       else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !V
!ACTIVE CONTRIBUTIONS
                Do iV = 1,kAsh
                  jV = iV + off_basAsh(ksym)
                do iX=1,iV
                    jX = iX + off_basAsh(ksym)
                    iVX = iTri(iV+off_Ash(ksym)) + iX + off_Ash(ksym)
                    DVX=D1MO(iVX)
                    !iVX = iTri(iV) + iX
                    !DVX=2.0D0*D1MO(iDoff+iVX)!Why *2?
                    If (iX.eq.iV) DVX=DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                       else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iU + off_basAsh(jSym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(iSym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                        else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                    V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                             V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                   !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                    V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !X
                end do !V
                Goto 500
!Symmetry case (II|KK)
200             Continue
                iFoff = off_Fmat(iSym)
                iDoff = off_Dmat(kSym)
                Do iV = 1,kIsh
                  jV = iV + off_basIsh(ksym)
                  iX=iV
                  jX=jV
                  DVX=2.0D0
                  If ( iX.eq.iV ) DVX = DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+nAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                   !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !V
!ACTIVE CONTRIBUTIONS
                Do iV = 1,kAsh
                  jV = iV + off_basAsh(ksym)
                  do iX=1,iV
                    jX = iX + off_basAsh(ksym)
                    iVX = iTri(iV+off_Ash(ksym)) + iX + off_Ash(ksym)
                    DVX=D1MO(iVX)
                    !iVX = iTri(iV) + iX
                    !DVX=2.0D0*D1MO(iDoff+iVX)
                    If (iX.eq.iV) DVX=DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+nAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP + iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                       d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                     else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                   !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                  SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                  *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(2,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(2,iGrid,jU)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(3,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(3,iGrid,jU)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(4,iGrid,jP)*TabMO(1,iGrid,jU)
     &                 +TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
     &                 *TabMO(1,iGrid,jP)*TabMO(4,iGrid,jU)
                d_ratio = 4.0d0/(dTot**2.0d0)*SQMO
                      if(((1.0d0-ratio).gt.thrsrho2)
     &                    .and.(ratio.lt.thrsrho3))  then
                      Zeta  = sqrt(1.0d0-ratio)
                        d_Zeta = -2.0D0/(Zeta*dTot**2.0D0)*SQMO
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                   SQMO/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                     +(0.5D0*RHOPx*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPx
     &                     +P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                     +(-0.5D0*RHOPx*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPx
     &                     -P2_ontop(2,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPx/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPx/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                     +(0.5D0*RHOPy*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPy
     &                     +P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                     +(-0.5D0*RHOPy*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPy
     &                     -P2_ontop(3,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPy/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPy/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(0.5D0*RHOPz*d_Zeta
     &                     -1.0D0/(Zeta*dTot)*SQMOPz
     &                     +P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     +(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     -(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(-0.5D0*RHOPz*d_Zeta
     &                     +1.0D0/(Zeta*dTot)*SQMOPz
     &                     -P2_ontop(4,iGrid)/(dTot*Zeta**2)*d_Zeta
     &                     -(RHOPz/(2.0d0*Zeta)*d_ratio)
     &                     +(Ratio*RHOPz/(2.0d0*Zeta**2.0d0)*d_Zeta))
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      elseif((ratio.ge.thrsrho3)
     &                  .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0
     &                          +4.0D0*Bb1*(ratio-1.15d0)**3.0d0
     &                          +3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                      d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0
     &                          +12.0D0*Bb1*(ratio-1.15d0)**2.0d0
     &                          +6.0D0*Cb1*(ratio-1.15d0)**1.0d0)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                    /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                   /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
                        d_Zetax =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPx
     &                          -2.0D0*RHOPx/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(2,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPx/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetay =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPy
     &                          -2.0D0*RHOPy/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(3,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPy/dTot))
     &                          *d_Deriv*d_ratio
                       d_Zetaz =Deriv*(4.0D0/(dTot**2.0D0)*SQMOPz
     &                          -2.0D0*RHOPz/dTot*d_ratio)
     &                          +((4.0D0*P2_ontop(4,iGrid)
     &                     /(dTot**2.0D0))-(2.0D0*ratio*RHOPz/dTot))
     &                          *d_Deriv*d_ratio
                   V_PUVX = V_PUVX +
     &                   Weights(iGrid)*(
     &                    (1.0D0/(2.0D0*dTot))
     &                   *SQMO*4.0D0*Deriv*
     &                   (dF_dRho(1,iGrid)-dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +(2.0D0*RHOPx/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetax)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +(2.0D0*RHOPy/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetay)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)
     &                     *(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +(2.0D0*RHOPz/(dTot**2.0D0)*SQMO*Deriv
     &                    +0.5D0*dTot*d_Zetaz)* (-1.0D0)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
              else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !X
                end do !V
                Goto 500
*               symmetry case (IJ!IJ)
300             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                   ! off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                   ! off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do
                Goto 500
*               symmetry case (IJ!KL)
400             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                  !  off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                  !  off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do
                Goto 500

500             Continue
              End If

            End Do
          End Do
        End Do
      End Do

!      Call Put_dArray('FI_V',Work(ifiv),ntot1)
!      Call Put_dArray('FA_V',Work(ifav),ntot1)
!      CALL GETMEM('FI_V','FREE','REAL',ifiv,ntot1)
!      CALL GETMEM('FA_V','FREE','REAL',ifav,ntot1)
      lsym=lsym_tmp
*
      Return
      End



      Subroutine Calc_OTOEGGA(OE,TabMO,mAO,nCoor,nTabMOs
     &                       ,P2_ontop,nP2_ontop,Rho,nRho,
     &                      dF_dRho,ndF_dRho,RhoI,RhoA,mRho,Weights,
     &                      nIrrep)
      Implicit Real*8 (A-H,O-Z)
      Dimension OE(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer nIrrep
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_Bas(mxSym)!,off_Fmat(mxSym),
!     &        off_Dmat(mxSym)
      Integer off_basIsh(mxSym),off_Fmat(mxSym)
      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),P2_ontop(nP2_ontop,nCoor),
     &       dF_dRho(ndF_dRho,nCoor),Rho(nRho,nCoor),
     &       RhoI(mRho,nCoor),RhoA(mRho,nCoor)
      Integer VX,ix,jx,iv,jv
*
      !iTri(i) = (i*i-i)/2
       iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      thrsrho=1.0d-15
      thrsrho2=1.0d-15
*
      lsym_tmp=lsym
      iStack  = 0
      iStack1 = 0
      Do iSym = 1,nSym
        off_Ash(iSym)    = iStack
        off_Bas(iSym)    = iStack1
        off_BasAsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)
        off_BasIsh(iSym) = iStack1+nFro(iSym)
        iStack1 = iStack1 + nBas(iSym)
        iStack  = iStack  + nAsh(iSym)
      End Do

      iStack = 0
      do iSym=1,nSym
        off_Fmat(iSym) = iStack
        iorb = nOrb(iSym)
        iStack = iStack + (iOrb*iOrb + iOrb)/2
      end do
*
      Call Unused_real_array(RhoI)
      Call Unused_real_array(RhoA)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!OE pieces - both orbs must be in the same irrep, eh?
      Do iSym = 1,nSym
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        Do iV = 1,iOrb!iIsh+iAsh
          jV = iV + off_BasIsh(iSym)
          do iX = 1,iV
            jX = iX + off_BasIsh(iSym)
            VX = off_Fmat(iSym) + iTri(iV,iX)
            fact=1.0d0
            Do iGrid = 1, nCoor
            dTot=Rho(1,iGrid)+Rho(2,iGrid)
              if (dTot.ge.thrsrho) then
                ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                      +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                      +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                      +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
                if((1.0d0-ratio).gt.thrsrho2) then
                  Zeta  = sqrt(1.0d0-ratio)
                       FTERMa=0.5D0*(1.0D0+Zeta)
                       FTERMb=0.5D0*(1.0D0-Zeta)
                       STERM=2.0D0*P2_ontop(1,iGrid)/(Zeta*dTot**2)
                       STERMG=2.0D0*P2_ontop(1,iGrid)
     &                       /(Zeta*dTot**3)
                       RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                       RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                       RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     OE(VX)=OE(VX)
     &                    +(SQMO
     &                    *(dF_dRho(1,iGrid)*(FTERMa
     &                    +STERM)
     &                    +dF_dRho(2,iGrid)*(FTERMb
     &                    -STERM))
                          !Now gradient part
                          !x alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
     &                     *(FTERMa*SQMOPx+
     &                     STERMG*RHOPx*SQMO)
                          !x beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
     &                     *(FTERMb*SQMOPx
     &                     -STERMG*RHOPx*SQMO)
                          !y alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
     &                     *(FTERMa*SQMOPy+
     &                     STERMG*RHOPy*SQMO)
                          !y beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
     &                     *(FTERMb*SQMOPy
     &                     -STERMG*RHOPy*SQMO)
                          !z alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
     &                     *(FTERMa*SQMOPz+
     &                     STERMG*RHOPz*SQMO)
                          !z beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                     *(FTERMb*SQMOPz
     &                     -STERMG*RHOPz*SQMO)
     &                    )*Weights(iGrid)*dble(nIrrep)
                  else
                         OE(VX)=OE(VX)
     &                          +(SQMO*
     &                     (dF_dRho(1,iGrid)
     &                     +dF_dRho(2,iGrid))
                            !gradient part
                          !x alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
     &                     *SQMOPx
                          !x beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
     &                     *SQMOPx
                          !y alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
     &                     *SQMOPy
                          !y beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
     &                     *SQMOPy
                          !z alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
     &                     *SQMOPz
                          !z beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                     *SQMOPz
     &                     )*Weights(iGrid)*0.5D0*dble(nIrrep)
                  end if
                 end if
            End Do     ! iGrid
          End do
        End Do
      End Do
           lsym=lsym_tmp
      Return
      End

      Subroutine Calc_OTOEGGA_ft(OE,TabMO,mAO,nCoor,nTabMOs
     &                       ,P2_ontop,nP2_ontop,Rho,nRho,
     &                      dF_dRho,ndF_dRho,RhoI,RhoA,mRho,Weights,
     &                      nIrrep)
      Implicit Real*8 (A-H,O-Z)
      Dimension OE(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer nIrrep
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_Bas(mxSym)!,off_Fmat(mxSym),
!     &        off_Dmat(mxSym)
      Integer off_basIsh(mxSym),off_Fmat(mxSym)
      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),P2_ontop(nP2_ontop,nCoor),
     &       dF_dRho(ndF_dRho,nCoor),Rho(nRho,nCoor),
     &       RhoI(mRho,nCoor),RhoA(mRho,nCoor)
      Integer VX,ix,jx,iv,jv
*
      !iTri(i) = (i*i-i)/2
       iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      thrsrho=1.0d-15
      thrsrho2=1.0d-15
      thrsrho3=0.9000000000d0
      thrsrho4=1.1500000000d0
      Ab1=-4.756065601d+2
      Bb1=-3.794733192d+2
      Cb1=-8.538149682d+1
*
      lsym_tmp=lsym
      iStack  = 0
      iStack1 = 0
      Do iSym = 1,nSym
        off_Ash(iSym)    = iStack
        off_Bas(iSym)    = iStack1
        off_BasAsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)
        off_BasIsh(iSym) = iStack1+nFro(iSym)
        iStack1 = iStack1 + nBas(iSym)
        iStack  = iStack  + nAsh(iSym)
      End Do

      Call Unused_real_array(RhoI)
      Call Unused_real_array(RhoA)

      iStack = 0
      do iSym=1,nSym
        off_Fmat(iSym) = iStack
        iorb = nOrb(iSym)
        iStack = iStack + (iOrb*iOrb + iOrb)/2
      end do
*

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!OE pieces - both orbs must be in the same irrep, eh?
      Do iSym = 1,nSym
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        Do iV = 1,iOrb!iIsh+iAsh
          jV = iV + off_BasIsh(iSym)
          do iX = 1,iV
            jX = iX + off_BasIsh(iSym)
            VX = off_Fmat(iSym) + iTri(iV,iX)
            fact=1.0d0
            Do iGrid = 1, nCoor
            dTot=Rho(1,iGrid)+Rho(2,iGrid)
              if (dTot.ge.thrsrho) then
                ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                SQMO=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
                SQMOPx=TabMO(2,iGrid,jV)*TabMO(1,iGrid,jX)
     &                      +TabMO(1,iGrid,jV)*TabMO(2,iGrid,jX)
                SQMOPy=TabMO(3,iGrid,jV)*TabMO(1,iGrid,jX)
     &                      +TabMO(1,iGrid,jV)*TabMO(3,iGrid,jX)
                SQMOPz=TabMO(4,iGrid,jV)*TabMO(1,iGrid,jX)
     &                      +TabMO(1,iGrid,jV)*TabMO(4,iGrid,jX)
                d_ratio = -8.0d0*P2_ontop(1,iGrid)/(dTot**3.0d0)*SQMO
                if(((1.0d0-ratio).gt.thrsrho2)
     &              .and.(ratio.lt.thrsrho3)) then
                  Zeta  = sqrt(1.0d0-ratio)
                  d_Zeta = -0.5D0/Zeta*d_ratio
                     FTERMa=0.5D0*(1.0D0+Zeta)
                       FTERMb=0.5D0*(1.0D0-Zeta)
                       STERM=2.0D0*P2_ontop(1,iGrid)/(Zeta*dTot**2)
                       STERMG=2.0D0*P2_ontop(1,iGrid)
     &                       /(Zeta*dTot**3)
                       RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                       RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                       RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                     OE(VX)=OE(VX)
     &                    +(SQMO
     &                    *(dF_dRho(1,iGrid)*(FTERMa
     &                    +STERM)
     &                    +dF_dRho(2,iGrid)*(FTERMb
     &                    -STERM))
                          !Now gradient part
                          !x alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
     &                     *(FTERMa*SQMOPx+
     &                     STERMG*RHOPx*SQMO
     &                    +(P2_ontop(2,iGrid)/((Zeta*dTot)**2.0D0)*
     &                      (SQMO*Zeta+d_Zeta*dTot))
     &                    +(ratio/(2*Zeta)*SQMOPx)
     &                    +(RHOPx/(2*Zeta)*d_ratio)
     &                    -(RHOPx*Ratio/(2.0D0*Zeta**2.0D0)*d_Zeta))
                          !x beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
     &                     *(FTERMb*SQMOPx
     &                     -STERMG*RHOPx*SQMO
     &                     -(P2_ontop(2,iGrid)/((Zeta*dTot)**2.0D0)*
     &                      (SQMO*Zeta+d_Zeta*dTot))
     &                    -(ratio/(2*Zeta)*SQMOPx)
     &                    -(RHOPx/(2*Zeta)*d_ratio)
     &                    +(RHOPx*Ratio/(2.0D0*Zeta**2.0D0)*d_Zeta))
                          !y alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
     &                     *(FTERMa*SQMOPy+
     &                     STERMG*RHOPy*SQMO
     &                    +(P2_ontop(3,iGrid)/((Zeta*dTot)**2.0D0)*
     &                      (SQMO*Zeta+d_Zeta*dTot))
     &                    +(ratio/(2*Zeta)*SQMOPy)
     &                    +(RHOPy/(2*Zeta)*d_ratio)
     &                    -(RHOPy*Ratio/(2.0D0*Zeta**2.0D0)*d_Zeta))
                          !y beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
     &                     *(FTERMb*SQMOPy
     &                     -STERMG*RHOPy*SQMO
     &                    -(P2_ontop(3,iGrid)/((Zeta*dTot)**2.0D0)*
     &                      (SQMO*Zeta+d_Zeta*dTot))
     &                    -(ratio/(2*Zeta)*SQMOPy)
     &                    -(RHOPy/(2*Zeta)*d_ratio)
     &                    +(RHOPy*Ratio/(2.0D0*Zeta**2.0D0)*d_Zeta))
                          !z alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
     &                     *(FTERMa*SQMOPz+
     &                     STERMG*RHOPz*SQMO
     &                    +(P2_ontop(4,iGrid)/((Zeta*dTot)**2.0D0)*
     &                      (SQMO*Zeta+d_Zeta*dTot))
     &                    +(ratio/(2*Zeta)*SQMOPz)
     &                    +(RHOPz/(2*Zeta)*d_ratio)
     &                    -(RHOPz*Ratio/(2.0D0*Zeta**2.0D0)*d_Zeta))
                          !z beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                     *(FTERMb*SQMOPz
     &                     -STERMG*RHOPz*SQMO
     &                    -(P2_ontop(4,iGrid)/((Zeta*dTot)**2.0D0)*
     &                      (SQMO*Zeta+d_Zeta*dTot))
     &                    -(ratio/(2*Zeta)*SQMOPz)
     &                    -(RHOPz/(2*Zeta)*d_ratio)
     &                    +(RHOPz*Ratio/(2.0D0*Zeta**2.0D0)*d_Zeta))
     &                    )*Weights(iGrid)*dble(nIrrep)
                 else if((ratio.ge.thrsrho3)
     &              .and.(ratio.le.thrsrho4)) then
                       Zeta = (Ab1*(ratio-1.15d0)**5.0d0) +
     &                        (Bb1*(ratio-1.15d0)**4.0d0) +
     &                        (Cb1*(ratio-1.15d0)**3.0d0)
                       Deriv = (5.0D0*Ab1*(ratio-1.15d0)**4.0d0) +
     &                        (4.0D0*Bb1*(ratio-1.15d0)**3.0d0) +
     &                        (3.0D0*Cb1*(ratio-1.15d0)**2.0d0)
                       d_Deriv = (20.0D0*Ab1*(ratio-1.15d0)**3.0d0) +
     &                        (12.0D0*Bb1*(ratio-1.15d0)**2.0d0) +
     &                        (6.0D0*Cb1*(ratio-1.15d0))
                       FTERMa=0.5D0*(1.0D0+Zeta)
                       FTERMb=0.5D0*(1.0D0-Zeta)
                       STERM=-4.0D0*P2_ontop(1,iGrid)/(dTot**2)*Deriv
                       STERMG=-4.0D0*P2_ontop(1,iGrid)
     &                       /(dTot**3)*Deriv
                       RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                       RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                       RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                       Zetax =Deriv*((4.0D0*P2_ontop(2,iGrid)
     &                          /(dTot**2.0D0))-(2.0D0*ratio*RHOPx
     &                          /dTot))
                       Zetay =Deriv*((4.0D0*P2_ontop(3,iGrid)
     &                         /(dTot**2.0D0))-(2.0D0*ratio*RHOPy
     &                         /dTot))
                       Zetaz =Deriv*((4.0D0*P2_ontop(4,iGrid)
     &                         /(dTot**2.0D0))-(2.0D0*ratio*RHOPz
     &                         /dTot))
                       d_Zetax =2.0D0*Deriv*((-4.0D0*P2_ontop(2,iGrid)
     &                         /(dTot**3.0d0)*SQMO)
     &                        -(ratio/(dTot)*SQMOPx)-(RHOPx*
     &                        (d_ratio/(dTot)
     &                         -ratio*SQMO/(dTot**2.0D0))))
     &                        +((4.0D0*P2_ontop(2,iGrid)
     &                          /(dTot**2.0D0))-(2.0D0*ratio*RHOPx
     &                          /dTot))*d_Deriv*d_ratio
                       d_Zetay =2.0D0*Deriv*((-4.0D0*P2_ontop(3,iGrid)
     &                          /(dTot**3.0d0)*SQMO)
     &                        -(ratio/dTot*SQMOPy)-(RHOPy*
     &                        (d_ratio/dTot
     &                         -ratio*SQMO/(dTot**2.0D0))))
     &                         +((4.0D0*P2_ontop(3,iGrid)
     &                         /(dTot**2.0D0))-(2.0D0*ratio*RHOPy
     &                         /dTot))*d_Deriv*d_ratio
                       d_Zetaz =2.0D0*Deriv*((-4.0D0*P2_ontop(4,iGrid)
     &                          /(dTot**3.0d0)*SQMO)
     &                        -(ratio/(dTot)*SQMOPz)
     &                        -(RHOPz*(d_ratio/(dTot)
     &                         -ratio*SQMO/(dTot**2.0D0))))
     &                        +((4.0D0*P2_ontop(4,iGrid)
     &                         /(dTot**2.0D0))-(2.0D0*ratio*RHOPz
     &                         /dTot))*d_Deriv*d_ratio
                    OE(VX)=OE(VX)
     &                    +(SQMO
     &                    *(dF_dRho(1,iGrid)*(FTERMa
     &                    +STERM)
     &                    +dF_dRho(2,iGrid)*(FTERMb
     &                    -STERM))
                          !Now gradient part
                          !x alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
     &                     *(FTERMa*SQMOPx+
     &                     STERMG*RHOPx*SQMO+0.5D0*dTot*d_Zetax
     &                      +0.5d0*Zetax*SQMO)
                          !x beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
     &                     *(FTERMb*SQMOPx
     &                     -STERMG*RHOPx*SQMO-0.5D0*dTot*d_Zetax
     &                      -0.5d0*Zetax*SQMO)
                          !y alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
     &                     *(FTERMa*SQMOPy+
     &                     STERMG*RHOPy*SQMO+0.5D0*dTot*d_Zetay
     &                      +0.5d0*Zetay*SQMO)
                          !y beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
     &                     *(FTERMb*SQMOPy
     &                     -STERMG*RHOPy*SQMO-0.5D0*dTot*d_Zetay
     &                      -0.5d0*Zetay*SQMO)
                          !z alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
     &                     *(FTERMa*SQMOPz+
     &                     STERMG*RHOPz*SQMO+0.5D0*dTot*d_Zetaz
     &                      +0.5d0*Zetaz*SQMO)
                          !z beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                     *(FTERMb*SQMOPz
     &                     -STERMG*RHOPz*SQMO-0.5D0*dTot*d_Zetaz
     &                      -0.5d0*Zetaz*SQMO)
     &                    )*Weights(iGrid)*dble(nIrrep)
                  else
                         OE(VX)=OE(VX)
     &                          +(SQMO*
     &                     (dF_dRho(1,iGrid)
     &                     +dF_dRho(2,iGrid))
                            !gradient part
                          !x alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
     &                     *SQMOPx
                          !x beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
     &                     *SQMOPx
                          !y alpha component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
     &                     *SQMOPy
                          !y beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
     &                     *SQMOPy
                          !z beta component
     &                    +(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
     &                     *SQMOPz
                          !z beta component
     &                    +(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                     *SQMOPz
     &                     )*Weights(iGrid)*0.5D0*dble(nIrrep)
                  end if
                 end if
            End Do     ! iGrid
          End do
        End Do
      End Do
           lsym=lsym_tmp
      Return
      End
      Subroutine Calc_OTPUVXGGA_2(PUVX,TabMO,mAO,nCoor,nTabMOs
     &                       ,P2_ontop,nP2_ontop,Rho,nRho,
     &                      dF_dRho,ndF_dRho,RhoI,RhoA,mRho,
     &                      Weights,D1MO,nD1MO,nIrrep)
      Implicit Real*8 (A-H,O-Z)
      Dimension PUVX(*)
#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "ksdft.fh"
      Integer nIrrep
      Integer off_Ash(mxSym), off_BasAsh(mxSym),
     &        off_PUVX(mxSym),off_Bas(mxSym),
     &        off_ish(mxSym),off_BasIsh(mxSym),
     &        off_BasVsh(mxSym)
      Integer   off_Dmat, off_Fmat
      Dimension off_Dmat(mxSym), off_Fmat(mxSym)

      Dimension TabMO(mAO,nCoor,nTabMOs),
     &       Weights(nCoor),P2_ontop(nP2_ontop,nCoor),
     &       dF_dRho(ndF_dRho,nCoor),Rho(nRho,nCoor),
     &       RhoI(mRho,nCoor),RhoA(mRho,nCoor)
      Real*8 D1MO(nD1MO)
      Integer count_tmp
      Real*8 DVX
      real*8 Fact
      Real*8 time1,time2
      Integer iftmpo,sztmp,ipq
      Real*8 junk_test,junk_test_t
      Integer case
      iTrii(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
      iTri(i)=(i*i-i)/2
*
      thrsrho=1.0d-15
      thrsrho2=1.0d-15
      thrspi=1.0d-30
*
      lsym_tmp=lsym

      Call Unused_real_array(RhoI)
      Call Unused_real_array(RhoA)

      Call CPU_Time(time1)

      iStack  = 0
      iStack1 = 0
      iStack2 = 0
      off_ish(:) = 0
      off_Ash(:) = 0
      off_Bas(:) = 0
      off_BasAsh(:) = 0
      off_BasIsh(:) = 0
      off_BasVsh(:) = 0
      ntot1 = 0

      Do iSym = 1,nSym
        off_ish(isym)    = iStack2
        off_Ash(iSym)    = iStack
        off_Bas(iSym)    = iStack1
        off_BasVsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)+nAsh(iSym)
        off_BasAsh(iSym) = iStack1+nIsh(iSym)+nFro(iSym)
        off_BasIsh(iSym) = iStack1+nFro(iSym)
        ntot1 = iTrii(nBas(iSym),nBas(iSym)) + ntot1
        iStack2 = iStack2 + nIsh(iSym)
        iStack1 = iStack1 + nBas(iSym)
        iStack  = iStack  + nAsh(iSym)
      End Do

!count the number of tmp_pot:
      count_tmp = 0
      do isym=1,nsym
        do jsym=1,nsym
        count_tmp = count_tmp + nIsh(isym)*(nIsh(jsym)+nAsh(jsym))**2
        end do
      end do

*
!Calculate PUVX offsets
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
*
!First build the "other stuff" - non-orbital grid-based info
!       Call GetMEM('other_st','ALLO','REAL',izet_stu,nCoor)
!       Call DCOPY_(1,0d0,0,Work(izet_stu),1)
!         Do iGrid=1,nCoor
!            dTot=Rho(1,iGrid)+Rho(2,iGrid)
!            ratio = 0.0d0
!            if(dTot.ge.thrsrho.and.
!     &           P2_ontop(1,iGrid).ge.thrspi) then
!               ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
!                if((1.0d0-ratio).gt.thrsrho2) then
!                    Zeta  = sqrt(1.0d0-ratio)
!                    RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
!                    RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
!                    RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
!                         izet_stu =
!     &                   Weights(iGrid)*dble(nIrrep)*(
!     &                   1.0D0/(Zeta*dTot)*
!     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
!                         !gradient part
!                         !x alpha part
!     &                    +RHOPx/(dTot**2*Zeta)
!     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
!     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
!                         !x beta part
!     &                    +RHOPx/(dTot**2*Zeta)
!     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
!     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
!                         !y alpha part
!     &                    +RHOPy/(dTot**2*Zeta)
!     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
!     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
!                         !y beta part
!     &                    +RHOPy/(dTot**2*Zeta)
!     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
!     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
!                         !z alpha part
!     &                    +RHOPz/(dTot**2*Zeta)
!     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
!     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
!     &                    +RHOPz/(dTot**2*Zeta)
!     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
!     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
!     &                    )
!                else
!                  PUVX(iPUVX) = PUVX(iPUVX) +0.0D0
!                end if
!              end if
!           End Do
!Now build pairs of active orbitals:
!      NASHT = 0
!      Do iSym = 1,nSym
!        NASHT = NASHT + nAsh(isym)
!      end do
!
!      counter =0
!      counter1 = 0
!      c_grid = 0
!       Call GetMEM('orb_pair','ALLO','REAL',iorb_pair,nCoor*nD1MO)
!      Do iSym = 1,nSym
!        iOrb = nOrb(iSym)
!        iAsh = nAsh(iSym)
!        iIsh = nIsh(iSym)
!        Do iV = 1,iAsh
!           jV = iV + off_basAsh(isym)
!           do iX=1,iV
!              jX = iX + off_basAsh(isym)
!              counter1 = counter1+1
!             do iGrid=1,nCoor
!               Work(iorb_pair+c_grid+counter1) =
!     &         TabMO(1,iGrid,jV)*TabMO(1,iGrid,JX)
!             end do
!           end do
!         end do
!       end do
!
!       Call GetMEM('orb_pair','FREE','REAL',iorb_pair,nCoor*nD1MO)
!       Call GetMEM('other_st','FREE','REAL',izet_stu,nCoor)


***********************************************************************
* BUILD PUVX potentials
!Note - I really only need the TUVX potentials.  Area for speedup.
***********************************************************************
      Fact = 1.00d0
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
!                        MO1=TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
                  Do iU = 1,jAsh
                  jU = iU + off_BasAsh(jSym)
                    Do iP = 1,iOrb
                      iT = iP - iIsh
                      iPUVX=iPUVX+1
                      jP = iP     +   off_Bas(iSym)
                    Do iGrid=1,nCoor
                       dTot=Rho(1,iGrid)+Rho(2,iGrid)
                     ratio = 0.0d0
                  if(dTot.ge.thrsrho.and.
     &               P2_ontop(1,iGrid).ge.thrspi) then
              ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                      if((1.0d0-ratio).gt.thrsrho2) then
                      Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
!                      interm1 = 1/(zeta*dTot)
!                      interm2 = interm1/dTot
                      PUVX(iPUVX) = PUVX(iPUVX) +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*dble(nIrrep)*(
     &                   1.0D0/(Zeta*dTot)*
!     &                   interm1*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
!     &                    +RHOPx/interm2
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
!     &                    +RHOPx/interm2
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
!     &                    +RHOPy/interm2
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
!     &                    +RHOPy/interm2
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
!     &                    +RHOPz/interm2
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
!     &                    +RHOPz/interm2
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                      else
                     PUVX(iPUVX) = PUVX(iPUVX) +0.0D0
                     end if
                   end if
                      End Do
                    End Do
                  End Do
                End Do
              End Do
            End If
          End Do
        End Do
      End Do

!Construction of Fock matrix pieces
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
!I think I can generate the contributions from the potentials to the new
!FI and FA terms at this level.  Then there will be no need to store a
!large number of V_TE potentials.
!*********************************
      Call CPU_Time(time3)
      PUVX_Time = PUVX_time + (time3-time1)

      iFoff=0
      ksym=1
      junk_test=0d0
      do iSym=1,nSym
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iPUVX = off_PUVX(iSym)
        sztmp = iTrii(nOrb(iSym),nOrb(iSym))
        if(iOrb.eq.0) cycle
        CALL GETMEM('TMP_other','ALLO','REAL',iftmpo,sztmp*nCoor)
        CALL DCOPY_(sztmp*nCoor,[0.0d0],0,Work(iftmpo),1)

        do iP=1,iOrb
          jP = iP + off_basIsh(iSym)
          do iU=1,iOrb
            jU = iU + off_basIsh(iSym)
            iPU = iTrii(iP,iU)
            do iGrid=1,nCoor
              dTot=Rho(1,iGrid)+Rho(2,iGrid)
              ratio = 0.0d0
              if(dTot.ge.thrsrho.and.
     &                P2_ontop(1,iGrid).ge.thrsrho) then
                      ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                if((1.0d0-ratio).gt.thrsrho2) then
                   Zeta  = sqrt(1.0d0-ratio)
                   Work(iftmpo + (iPU-1)*nCoor + iGrid - 1) =
     &             TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)
     &             *Weights(iGrid)
     &             *1.0D0/(Zeta*dTot)
                else
                  Work(iftmpo + (iPU-1)*nCoor + iGrid - 1) = 0d0
                end if
              end if
            end do!iGrid
          end do!iU
        end do!iP

        do iGrid=1,nCoor
          dTot=Rho(1,iGrid)+Rho(2,iGrid)
          if(dTot.ge.thrsrho.and.
     &     P2_ontop(1,iGrid).ge.thrsrho) then
            ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
            if((1.0d0-ratio).gt.thrsrho2) then
              RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
              RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
              RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
              junk_test = !Work(iftmpo+(ipq-1)*nCoor+iGrid)*
!     &               TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                ((-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !x alpha part
     &                   +(RHOPx/dTot)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
!            Work(ifav_n+iFoff+ipq-1) = Work(ifav_n+iFoff+ipq-1) +
!     &                               DVX*junk_test*dble(nIrrep)
            else
              junk_test=0
              cycle!skip this gridpt
            end if
          end if
          do iPQ=1,sztmp
            junk_test_t = Work(iftmpo+(iPQ-1)*nCoor+iGrid-1)
            do kSym = 1, nSym
              kOrb = nOrb(kSym)
              kAsh = nAsh(kSym)
              kIsh = nIsh(kSym)
              do ik = 1,nIsh(kSym)
                jK = iK + off_basIsh(kSym)
                Work(ifiv+iFoff+iPQ-1) = Work(ifiv+iFoff+iPQ-1) +
     &                               junk_test*dble(nIrrep)*junk_test_t*
     &                               TabMO(1,iGrid,jK)**2
              end do!iK
              do iV = 1,nAsh(kSym)
                jV = iV + off_BasAsh(kSym)
                do iX=1,iV
                  jX = iX + off_BasAsh(kSym)
                  iVX = iTri(iV+off_Ash(ksym)) + iX + off_Ash(ksym)
                  DVX=D1MO(iVX)
                  If (iX.eq.iV) DVX=DVX*0.5d0
                  Work(ifav+iFoff+iPQ-1) = Work(ifav+iFoff+iPQ-1) +
     &                DVX*junk_test*dble(nIrrep)*junk_test_t*
     &                TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)
                end do!1X
              end do!iV
            end do!kSym
          end do!ipq
        end do!igrid
!
        iFoff = iFoff + sztmp
        CALL GETMEM('TMP_other','FREE','REAL',iftmpo,sztmp*nCoor)
!
      end do!iSym
      Call CPU_Time(time2)
      SP_Time = SP_time + (time2-time3)

      Return
!*********************************
!      CALL GETMEM('FI_V','ALLO','REAL',ifiv,ntot1)
!      CALL GETMEM('FA_V','ALLO','REAL',ifav,ntot1)
!      Call Get_dArray('FI_V',Work(ifiv),ntot1)
!      Call Get_dArray('FA_V',Work(ifav),ntot1)

      Do iSym = 1,nSym !sym for p
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iPUVX = off_PUVX(iSym)
        Do jSym = 1,nSym !sym for u
          jOrb = nOrb(jSym)
          jAsh = nAsh(jSym)
          jIsh = nIsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym !Sym for v
            kOrb = nOrb(kSym)
            kAsh = nAsh(kSym)
            kIsh = nIsh(kSym)
            Do lSym = 1,kSym !sym for x
              lOrb = nOrb(lSym)
              lAsh = nAsh(lSym)
              lIsh = nIsh(lSym)
              klSym = 1 + ieor(kSym-1,lSym-1)

*             find cases
              case = 4
              If ( iSym.eq.jSym ) case = case-2
              If ( iSym.eq.kSym ) case = case-1

              If ( ijSym.eq.klSym .and.
!     &             iAsh*jAsh*kAsh*lAsh.ne.0 ) then
     &             iOrb*jOrb*kOrb*lOrb.ne.0 ) then

                Goto (100,200,300,400) case

!Symmetry case (II|II)
100             Continue
                iFoff = off_Fmat(iSym)
                iDoff = off_Dmat(iSym)
                Do iV = 1,kIsh
                  jV = iV + off_basIsh(kSym)
                  iX=iV
                  jX=jV
                  DVX=2.0D0
                  If ( iX.eq.iV ) DVX = DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !V

      Call CPU_Time(time1)
      FI_time = FI_Time + (time1-time2)
!ACTIVE CONTRIBUTIONS
                Do iV = 1,kAsh
                  jV = iV + off_basAsh(ksym)
                do iX=1,iV
                    jX = iX + off_basAsh(ksym)
                    iVX = iTri(iV+off_Ash(ksym)) + iX + off_Ash(ksym)
                    DVX=D1MO(iVX)
                    !iVX = iTri(iV) + iX
                    !DVX=2.0D0*D1MO(iDoff+iVX)!Why *2?
                    If (iX.eq.iV) DVX=DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iU + off_basAsh(jSym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(iSym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !X
                end do !V
      Call CPU_Time(time2)
      FA_time = FA_Time + (time2-time1)
                Goto 500
!Symmetry case (II|KK)
200             Continue
                iFoff = off_Fmat(iSym)
                iDoff = off_Dmat(kSym)
                Do iV = 1,kIsh
                  jV = iV + off_basIsh(ksym)
                  iX=iV
                  jX=jV
                  DVX=2.0D0
                  If ( iX.eq.iV ) DVX = DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+nAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifiv+iFoff+iPU-1) = Work(ifiv+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !V
      Call CPU_Time(time1)
      FI_time = FI_Time + (time1-time2)
!ACTIVE CONTRIBUTIONS
                Do iV = 1,kAsh
                  jV = iV + off_basAsh(ksym)
                  do iX=1,iV
                    jX = iX + off_basAsh(ksym)
                    iVX = iTri(iV+off_Ash(ksym)) + iX + off_Ash(ksym)
                    DVX=D1MO(iVX)
                    !iVX = iTri(iV) + iX
                    !DVX=2.0D0*D1MO(iDoff+iVX)
                    If (iX.eq.iV) DVX=DVX*0.5
                  !inact/inact case, p <= u
                  do iU=1,jIsh
                    jU = iu + off_basIsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !Inact/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+nAsh(jsym)
                    do iP = 1,iIsh
                      jP = iP + off_basIsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/act case
                  do iU=1,jAsh
                    jU = iu + off_basAsh(jsym)
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+iU) + iP+iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !act/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iAsh
                      jP = iP + off_basAsh(isym)
                      iPU  = iTri(jIsh+jAsh+iU) + iP + iIsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                  !virt/virt case
                  do iU=1,jOrb-jAsh-jIsh
                    jU = iu + off_basAsh(jsym)+jAsh
                    do iP = 1,iU
                      jP = iP + off_basAsh(isym)+iAsh
                      iPU  = iTri(jIsh+jAsh+iU) + iP+iIsh+iAsh
                      V_PUVX = 0.0d0
                      Do iGrid=1,nCoor
                        dTot=Rho(1,iGrid)+Rho(2,iGrid)
                        ratio = 0.0d0
                        if(dTot.ge.thrsrho.and.
     &                  P2_ontop(1,iGrid).ge.thrsrho) then
                          ratio = 4.0d0*P2_ontop(1,iGrid)/(dTot**2.0d0)
                          if((1.0d0-ratio).gt.thrsrho2) then
                            Zeta  = sqrt(1.0d0-ratio)
                      RHOPx=Rho(3,iGrid)+Rho(6,iGrid)
                      RHOPy=Rho(4,iGrid)+Rho(7,iGrid)
                      RHOPz=Rho(5,iGrid)+Rho(8,iGrid)
                      V_PUVX = V_PUVX +
     &                   TabMO(1,iGrid,jP)*TabMO(1,iGrid,jU)*
     &                   TabMO(1,iGrid,jV)*TabMO(1,iGrid,jX)*
     &                   Weights(iGrid)*(
     &                   1.0D0/(Zeta*dTot)*
     &                   (-dF_dRho(1,iGrid)+dF_dRho(2,iGrid))
                         !gradient part
                         !x alpha part
     &                    +RHOPx/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(3,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(6,iGrid))
                         !x beta part
     &                    +RHOPx/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(6,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(3,iGrid))
                         !y alpha part
     &                    +RHOPy/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(4,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(7,iGrid))
                         !y beta part
     &                    +RHOPy/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(7,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(4,iGrid))
                         !z alpha part
     &                    +RHOPz/(dTot**2*Zeta)
     &                   * (-1.0D0)*(2.0D0*dF_dRho(3,iGrid)*Rho(5,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(8,iGrid))
                         !z beta part
     &                    +RHOPz/(dTot**2*Zeta)
     &                     *(2.0D0*dF_dRho(5,iGrid)*Rho(8,iGrid)
     &                     +dF_dRho(4,iGrid)*Rho(5,iGrid))
     &                    )
                          else
                            V_PUVX = V_PUVX + 0.0d0
                          end if
                        end if
                      End Do!gridpt
            Work(ifav+iFoff+iPU-1) = Work(ifav+iFoff+iPU-1) +
     &                               DVX*V_PUVX*dble(nIrrep)
                    end do
                  end do
                end do !X
                end do !V
      Call CPU_Time(time2)
      FA_time = FA_Time + (time2-time1)
                Goto 500
*               symmetry case (IJ!IJ)
300             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                   ! off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                   ! off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do
                Goto 500
*               symmetry case (IJ!KL)
400             Continue
                Do iV = 1,kAsh
                  Do iX = 1,lAsh
                  !  off_PUVX(iSym) = off_PUVX(iSym) + jAsh*iOrb
                  !  off_PUVX(jSym) = off_PUVX(jSym) + iAsh*jOrb
                  End Do
                End Do
                Goto 500

500             Continue
              End If

            End Do
          End Do
        End Do
      End Do

!      Call Put_dArray('FI_V',Work(ifiv),ntot1)
!      Call Put_dArray('FA_V',Work(ifav),ntot1)
!      CALL GETMEM('FI_V','FREE','REAL',ifiv,ntot1)
!      CALL GETMEM('FA_V','FREE','REAL',ifav,ntot1)
      lsym=lsym_tmp
*
      Return
      End
