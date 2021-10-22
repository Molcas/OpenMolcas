!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE CHO_eval_waxy(irc,Scr,ChoV1,ChoV2,W_PWXY,nAorb,        &
     &                         JSYM,NUMV,DoTraInt,CMO)

      Use Data_structures, only: SBA_Type, twxy_Type, DSBA_Type
      Implicit Real*8 (a-h,o-z)
      Real*8 W_PWXY(*)
      Integer nAorb(*)
      Type (DSBA_Type) CMO
      Type (SBA_Type) ChoV1, ChoV2
      Type (twxy_type) Scr
      Integer off_PWXY(8,8,8)
      Logical DoTraInt

#include "real.fh"
#include "Molcas.fh"
#include "general.fh"
#include "wadr.fh"

! ************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
! ************************************************

      If (NumV .lt. 1) Return

! --- Computing the integrals
! ---------------------------------------------------------
! --- (wa|xy)  <-  (wa|xy)  +  sum_J  L(wa,#J) * L(xy,#J)
!==========================================================

      Do iSymy=1,nSym

         iSymx=MulD2h(iSymy,JSYM)

         Nxy=SIZE(ChoV2%SB(iSymx)%A2,1)

         If (iSymx.le.iSymy.and.Nxy.gt.0) then

            Do iSyma=1,nSym

               iSymw=MulD2h(iSyma,JSYM)

               Nwa  = SIZE(ChoV1%SB(iSymw)%A3,1)*                       &
     &                SIZE(ChoV1%SB(iSymw)%A3,2)

               If (Nwa<=0) Cycle

               CALL DGEMM_('N','T',Nwa,Nxy,NumV,                        &
     &                    ONE,ChoV1%SB(iSymw)%A3,Nwa,                   &
     &                        ChoV2%SB(iSymx)%A2,Nxy,ONE,               &
     &                    Scr%SB(iSymw,iSymx)%A,Nwa)


            End Do

         End If

      End Do


! --- MO-transformation to build the (pw|xy) integrals
! --- is performed (p is a general index).
! --- The storage of the result is the one required by
! --- the RASSCF program : LT-storage in xy, with sym(x)<=sym(y)
! ---
! ---       (pw|xy)  =  sum_a  C(a,p) * (wa|xy)
! ---
! ------------------------------------------------------------
      IF (DoTraInt) THEN

!
!     generate offsets to (pw|xy)
!
!     Important: the way GET_TUVX is written requires
!                a strict loop structure for the definition
!                of the offsets
!
         iStack = 0
         Do iSymp = 1,nSym
            iOrb = nOrb(iSymp)
            Do iSymw = 1,nSym
               jAsh = nAorb(iSymw)
               ijSym=MulD2h(iSymp,iSymw)
               Do iSymy = 1,nSym
                 kAsh = nAorb(iSymy)
                 iSymx=MulD2h(ijSym,iSymy)
                 If (iSymx.le.iSymy) Then
                  lAsh = nAorb(iSymx)
                   kl_Orb_pairs = kAsh*lAsh                             &
     &                          + Min(0,ijSym-2)*kAsh*(lAsh-1)/2
                   off_PWXY(iSymp,iSymw,iSymx)=iStack
                   iStack = iStack + iOrb*jAsh*kl_Orb_pairs
                 End If
               End Do
            End Do
         End Do
         nPWXY=iStack
!
!        Reordering and MO-transformation
!
         Do iSymy=1,nSym

            iSymx=MulD2h(iSymy,JSYM)

            If (iSymx.le.iSymy) then

               Nxy  = nAorb(iSymx)*nAorb(iSymy)                         &
     &              + Min(0,JSYM-2)*nAorb(iSymx)*(nAorb(iSymy)-1)/2


               Do iSymw=1,nSym

                  iSyma = MulD2h(iSymw,JSYM) ! =iSymp

                  Nwa = nAorb(iSymw)*nBas(iSyma)
                  Npw = nOrb(iSyma)*nAorb(iSymw)

                  iS = 1 + nBas(iSyma)*nFro(iSyma)

                  Do ixy = 1,Nxy

                        ipMpw = off_PWXY(iSyma,iSymw,iSymx)             &
     &                        +  1 + Npw*(ixy-1)

! --------------------------------------------------------
! ---       M(p,w)[xy]  =  sum_a  C(a,p) * M(w,a)[xy]
! --------------------------------------------------------
                        nBas_a = max(nBas(iSyma),1)
                        nAob_w = max(nAorb(iSymw),1)
                        nOrb_a = max(nOrb(iSyma),1)

                        CALL DGEMM_('T','T',                            &
     &                             nOrb(iSyma),nAorb(iSymw),nBas(iSyma),&
     &                             ONE,CMO%SB(iSyma)%A1(iS:),nBas_a,    &
     &                              Scr%SB(iSymw,iSymx)%A(:,ixy),nAob_w,&
     &                            ZERO,W_PWXY(ipMpw),nOrb_a)


                  End Do

               End Do

            End If

         End Do


      ENDIF


      irc=0

      Return
      END

!*************************************************************
