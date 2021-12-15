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
      Subroutine Get_D1I_RASSCF_m(CMO,D1I)
************************************************************************
*                                                                      *
*     Compute one-body density (AO-basis)                              *
*     from frozen+inactive orbitals                                    *
*                                                                      *
*     calling arguments:                                               *
*     CMO     : input, array of real*8                                 *
*               MO-coefficients                                        *
*     D1I     : input, array of real*8                                 *
*               one-electron density matrix (frozen+inactive)          *
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

      Dimension CMO(*), D1I(*)

#include "rasdim.fh"
#include "general.fh"

      Parameter ( Zero=0.0d0 , Two=2.0d0 )


* PAM Aug 2006, replace the following explicit code
* with the corresponding DGEMM calls.
*      iOff1 = 0
*      Do iSym = 1,nSym
*        iBas = nBas(iSym)
*        iOrb = nFro(iSym) + nIsh(iSym)
*        If ( iBas.ne.0 ) then
*          iOff2 = iOff1
*          Do i = 1,iBas
*            Do j = 1,iBas
*              Sum = Zero
*              Do k = 0,iOrb-1
*                Sum = Sum + Two * CMO(iOff1+k*iBas+i)
*     &                          * CMO(iOff1+k*iBas+j)
*              End Do
*              D1I(iOff2 + j) = Sum
*            End Do
*            iOff2 = iOff2 + iBas
*          End Do
*          iOff1 = iOff1 + iBas*iBas
*        End If
*      End Do
* Replaced by:
      ista=1
      do isym=1,nsym
        nb=nbas(isym)
        nbsq=nb**2
        nfi=nfro(isym)+nish(isym)
        if(nb.gt.0) then
          call dcopy_(nbsq,[0.0d0],0,d1i(ista),1)
          if(nfi.gt.0) then
            call DGEMM_('n','t',nb,nb,nfi,two,
     &           cmo(ista),nb,cmo(ista),nb,zero,d1i(ista),nb)
          end if
          ista=ista+nbsq
        end if
      end do


      Return
      End
