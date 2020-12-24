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
* Copyright (C) 1991, Markus P. Fuelscher                              *
*               1999, Roland Lindh                                     *
************************************************************************
      subroutine Motra(ireturn)
************************************************************************
*                                                                      *
*     Objective: AO to MO integral transformation                      *
*                                                                      *
*     Modify one-electron integrals to use dynamic memory allocation.  *
*     R. Lindh, March 1999.                                            *
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************

      !> module dependencies
#ifdef _HDF5_QCM_
      use hdf5_utils
#endif

#include "motra_global.fh"
#include "trafo_motra.fh"
#include "WrkSpc.fh"
      COMMON / CHO_Minp / iCTonly, iDoInt
      Character*3  tv2disk
      COMMON / CHOTRAW /tv2disk
      Logical DoCholesky, Do_int
*----------------------------------------------------------------------*
*     Start program and say Hello                                      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*     ( Dynamic work area has been allocated in Start() )              *
*----------------------------------------------------------------------*
*      Call IniMem
*----------------------------------------------------------------------*
*     Run through the input section                                    *
*----------------------------------------------------------------------*
      Call init_motra
      If (iPrintLevel(-1).LE.0) iPrint=-1
      Call InpCtl_Motra(ipOvlp,ipHOne,ipKine,ipCMO)

      !> initialize HDF5 interface
#ifdef _HDF5_QCM_
      if(ihdf5 == 1)then
        !> enable HDF5 support and open the file ijklname
        call hdf5_init()
        call hdf5_create(ijklname, file_id(1))
      end if
#endif

*----------------------------------------------------------------------*
* --- Cholesky check
      Call DecideOnCholesky(DoCholesky)
*----------------------------------------------------------------------*
* --- Use of MOTRA only for AO-->MO transf of the Cholesky vectors
      If (iCTonly.eq.1) Then
         If (.not.DoCholesky) Then
            write(6,*)'      Warning! This is not RI/CD calculation: '
            write(6,*)'                      keyword CTonly ignored! '
         Else
#ifdef _HDF5_QCM_
            If ((ihdf5.eq.1).and.(tv2disk.ne.'KPQ')) Then
              Write(6,*)' Transformed Cholesky vectors cannot be '//
     &          ' written as (pq,K) in HDF5 file as of now. Activate'//
     &          ' the KPQ option to store them or disable'//
     &          ' the HDF5 option.'
              Call Abend()
            End If
#endif
            write(6,*)
            write(6,*)'      ... Skipping MoTRA of ERIs ...'
            write(6,*)
            write(6,*)'      ... but Cholesky vectors will be MoTRA.'
            write(6,*)
!           Cholesky vectors in HDF5 must be stored as KPQ format (for now)
            Do_int=.false.
            If (iDoInt.eq.1) Do_int=.true.
            Call Cho_MOtra(Work(ipCMO),nTot2,Do_int,ihdf5)
            iOneOnly=666
            Go To 100  ! nothing else to be done except OneEl part
         EndIf
      EndIf
*----------------------------------------------------------------------*
* --- Preliminary step for integral transformation with CD
      If (DoCholesky) Then
         CALL CWTIME(TCR1,TWR1)
         Call Cho_X_init(irc,0.0)
         If (irc.ne.0) Then
           write(6,*) ' In MoTRA : Cho_X_Init returned non-zero'//
     &                ' rc = ',irc
           Call Abend()
         EndIf
         Call Cho_X_ReoVec(irc) ! get (if not there) CD vecs full stor
         If (irc.ne.0) Then
           write(6,*) ' In MoTRA : Cho_X_ReoVec returned non-zero'//
     &                ' rc = ',irc
           Call Abend()
         EndIf
         Call Cho_X_final(irc)
         CALL CWTIME(TCR2,TWR2)
         tcpu_reo=(TCR2-TCR1)
         write(6,*)
         write(6,*)'      Reordering Cholesky vectors to full storage.'
         write(6,*)'       Elapsed time for the reordering : ',tcpu_reo
         write(6,*)'       CPU time for the reordering     : ',tcpu_reo
         write(6,*)
      EndIf


*----------------------------------------------------------------------*
*     Transform the one-electron integrals                             *
*----------------------------------------------------------------------*
100   Continue
      Call Tr1Ctl(Work(ipOvlp),Work(ipHOne),Work(ipKine),Work(ipCMO))
      If ( iOneOnly.ne.0 ) Goto 900
*----------------------------------------------------------------------*
*     Transform the two-electron integrals                             *
*----------------------------------------------------------------------*
      Call Tr2Ctl(Work(ipCMO))
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
 900  Continue
      write(6,*)
      Call GetMem('CMO','Free','Real',ipCMO,nTot2)
      Call GetMem('Kine','Free','Real',ipKine,nTot1+4)
      Call GetMem('HOne','Free','Real',ipHOne,nTot1+4)
      Call GetMem('Ovlp','Free','Real',ipOvlp,nTot1+4)
*
#ifdef _HDF5_QCM_
      if(ihdf5 == 1)then
        !> close the file ijkl.h5 and turn off HDF5 support.
        call hdf5_close(file_id(1))
        call hdf5_exit()
      end if
#endif

      Call FastIO('STATUS')
      ireturn=0
      return
      End
