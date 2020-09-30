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
      Subroutine set_nm( exch, ncut, encut_definition, nk, mg,
     &                   nTempMagn, hmax, w, encut_rate, TempMagn,
     &                   nM, EM, dbg )

      Implicit None
#include "warnings.fh"
c input data:
      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)        :: exch, ncut, encut_definition, nk,
     &                              mg, nTempMagn
      Real(kind=8), intent(in)  :: hmax, W(exch), encut_rate,
     &                              TempMagn(nTempMagn)
      Logical, intent(in)        :: dbg
c output data:
      Integer,       intent(out) :: nM
      Real(kind=8), intent(out) :: EM
c local variables:
      Integer       :: i
      Real(kind=8) :: diff, T_High
      Real(kind=8) :: boltz_k, mu_bohr

      ! Constants:
      boltz_k=0.6950356_wp                    !   in cm^-1*K-1
      mu_bohr=0.466864374_wp                  !   in cm-1*T-1

      nM     = 1
      EM     = 0.0_wp
      diff   = 0.0_wp
      T_High = 0.0_wp
      If(nTempMagn>0) T_High = MAXVAL(TempMagn(1:nTempMagn))

      If (dbg) Write(6,*) 'exch             = ', exch
      If (dbg) Write(6,*) 'ncut             = ', ncut
      If (dbg) Write(6,*) 'encut_definition = ', encut_definition
      If (dbg) Write(6,*) 'nk               = ', nk
      If (dbg) Write(6,*) 'mg               = ', mg
      If (dbg) Write(6,*) 'nM               = ', nM
      If (dbg) Write(6,*) 'nTempMagn        = ', nTempMagn
      If (dbg) Write(6,*) 'hmax             = ', hmax
      If (dbg) Write(6,*) 'encut_rate       = ', encut_rate
      If (dbg) Write(6,*) 'EM               = ', EM
      If (dbg) Write(6,*) 'TempMagn()       = ', TempMagn(1:nTempMagn)
      If (dbg) Write(6,*) 'W()              = ', W(1:exch)


      If(encut_definition==1) Then

         If ( ncut>exch ) Then
            nm = exch
            em = w(exch)
         Else
            nm = ncut
            em = w(nm)
         End If


      Else If(encut_definition==2) Then

         nm = exch
         em = nk*boltz_k*T_High + mg*mu_bohr*abs(hmax)

         Do i=1,exch
            If ( i>1 ) diff = w(i) - w(i-1)
            If( ( w(i)>em ) .and. ( diff>1.0d-4 ) ) Then
               nm=i-1
               Go To 309
            End If
         End Do


      Else If(encut_definition==3) Then

         nm = exch
         em = w(exch)*encut_rate

         Do i=1,exch
            If(i>1) diff = w(i) - w(i-1)
            If( ( w(i)>em ) .and. ( diff>1.0d-4 ) ) Then
               nm=i-1
               Go To 309
            End If
         End Do


      Else

         Write(6,'(A)') 'something is wrong with "encut_definition" '
         Call quit(_RC_INPUT_ERROR_)

      End If !encut_definition

309   Continue

      Return
      End subroutine set_nm

