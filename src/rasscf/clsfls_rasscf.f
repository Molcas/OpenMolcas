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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!***********************************************************************
      Subroutine ClsFls_RASSCF()
!***********************************************************************
!                                                                      *
!     Close files.                                                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************
#ifdef _HDF5_
      use mh5, only: mh5_close_file
      use RASWfn, only: wfn_fileid
#endif
      use general_data, only: JOBOLD,JOBIPH,ITERFILE,LUDAVID,LUINTM,    &
     &                        LUQUNE

      Implicit None
      Logical DoCholesky
      Integer iRC

!----------------------------------------------------------------------*
!     Start                                                            *
!-------------------------------------- -------------------------------*
! Local print level (if any)
!---  close the JOBOLD file -------------------------------------------*
      If(JOBOLD.gt.0.and.JOBOLD.ne.JOBIPH) Then
        Call DaClos(JOBOLD)
        JOBOLD=-1
      Else If (JOBOLD.gt.0) Then
        JOBOLD=-1
      End If
!---  close the JOBIPH file -------------------------------------------*
      If(JOBIPH.gt.0) Then
        Call DaClos(JOBIPH)
        JOBIPH=-1
      End If
#ifdef _HDF5_
      if (wfn_fileid.ne.0) then
        call mh5_close_file(wfn_fileid)
        wfn_fileid=0
      end if
#endif
!---  close the ORDINT file -------------------------------------------*
      CALL DecideOnCholesky(DoCholesky)
       If (.not.DoCholesky) then
         iRc=-1
         Call ClsOrd(iRc)
         If ( iRc.ne.0 ) Then
           Call WarningMessage(1,'Failed to close the ORDINT file.')
         End If
       End If
!---  close the file carrying the transformed two-electron integrals --*
      Call DaClos(LUINTM)
!---  close the DAVID file carrying temporary CI and sigma vectros ----*
      Call DaClos(LUDAVID)
!---  open the file carrying the hessian update vectors ---------------*
      Call DaClos(LuQune)
!---  close the file for storage of informations on CI-iterations
      Close(ITERFILE)
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
      End Subroutine ClsFls_RASSCF
