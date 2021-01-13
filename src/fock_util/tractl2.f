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
      Subroutine TraCtl2(CMO,PUVX,TUVX,D1I,FI,D1A,FA,IPR,lSquare,ExFac)

************************************************************************
*                                                                      *
*     main control section for Cholesky-based vs conventional          *
*     - transformation of ERIs from AO to MO basis                     *
*     - Fock matrix generation                                         *
*                                                                      *
************************************************************************
#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      Implicit Real*8 (A-H,O-Z)
      Parameter ( Zero=0.0d0 )

#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "wadr.fh"

*
      Dimension CMO(*), PUVX(*), TUVX(*)
      Dimension D1I(*), D1A(*), FI(*), FA(*)
      Logical lSquare,DoCholesky,TraOnly
      Integer ALGO
      Common /CHLCAS / DoCholesky,ALGO
*
*
*
*      Call DecideOnCholesky(DoCholesky)

      If (.not. DoCholesky) Then

        Call TRA_CTL2(CMO,PUVX,TUVX,D1I,FI,D1A,FA,IPR,lSquare,ExFac)

      ElseIf (ALGO.eq.1) Then

        lpwxy = ip_of_work(PUVX(1))

        TraOnly=.false.
       Call CHO_CAS_DRV(irc,CMO,D1I,FI,D1A,FA,WORK(LPMAT),TraOnly)

#if defined (_MOLCAS_MPP_)
c --------------------------------------------------
C  Synchronize Fock matrices if running parallel:
        If (nProcs .gt. 1 .and. Is_Real_Par()) then
           Call GADsum(FI,nTot1)
           Call GADsum(FA,nTot1)
C  Synchronize PUVX if running parallel:
           Call GADsum(PUVX,nPWXY)
        EndIf
#endif
c --------------------------------------------------
*     select integrals TUVX
        Call Get_TUVX(PUVX,TUVX)
*     save integrals on disk
*     nPWXY is computed in cho_eval_waxy
*           and stored in wadr.fh
        iDisk=0
        Call DDaFile(LUINTM,1,PUVX,nPWXY,iDisk)

      ElseIf (ALGO.eq.2) Then

        TraOnly=.false.
       Call CHO_CAS_DRV(irc,CMO,D1I,FI,D1A,FA,WORK(LPMAT),TraOnly)

        If (irc.ne.0) Then
         write(6,*)'TRACTL2: Cho_cas_drv non-Zero return code. rc= ',irc
         Call Abend
        EndIf

C  Synchronization for parallel runs is done in cho_cas_drv

      EndIf


      Return
      End
