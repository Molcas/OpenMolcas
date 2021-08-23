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
! Copyright (C) 1991, Per-Olof Widmark                                 *
!               1993,1996,1997, Markus P. Fuelscher                    *
!               1996, Luis Serrano-Andres                              *
!               2012, Victor P. Vysotskiy                              *
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
!                                                                      *
!     This file defines the module Fast_IO                             *
!     including all entries required for fast I/O system.              *
!                                                                      *
!     Following the list of entries and their usage:                   *
!     LuName    : Name of the file                                     *
!     isOpen    : open/close flag                                      *
!     FSCB      : file descriptors (C language)                        *
!     Addr      : pointer to the current position                      *
!     Trace     : enable/disable debugging output                      *
!     Query     : enable/disable the traceback facilities              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, IBM Sweden, 1991                                   *
!     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
!     L. Serrano-Andres, University of Lund, Sweden, 1996              *
!                                                                      *
! History:                                                             *
!     New I/O Stat, V.P. Vysotskiy, University of Lund, Sweden, 2012   *
!     OpenMP,       V.P. Vysotskiy, University of Lund, Sweden, 2012   *
!     Cleanup and conversion to module,  I. Fdez. Galvan, 2021         *
!***********************************************************************

module Fast_IO

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: MxFile = 199, MaxSplitFile = 20
integer(kind=iwp) :: Addr(MxFile), FlsSize(MxFile), FSCB(MxFile), isFiM(MxFile), isOpen(MxFile), MaxFileSize, MBL(MxFile), &
                     MPUnit(0:MaxSplitFile-1,MxFile), NProfFiles
logical(kind=iwp) :: Multi_File(MxFile), Query, Trace
real(kind=wp) :: ProfData(8,MxFile)
character(len=8) :: LuName(MxFile), LuNameProf(MxFile)
!-----------------------------------------------
! Error codes for AIX/IO and FIM/IO routines.
integer(kind=iwp), parameter :: eEof   = 1024, &
                                eNtOpn = 1025, &
                                eInErr = 1026, &
                                eTmF   = 1027, &
                                eTlFn  = 1028, &
                                eBlNme = 1029, &
                                eNoMsg = 1030, &
                                eFiMFo = 1031
!-----------------------------------------------
! Control blocks for AIX I/O routines.
integer(kind=iwp), parameter :: pHndle = 1, &
                                pWhere = 2, &
                                pDesc  = 3, &
                                pStat  = 4
integer(kind=iwp) :: CtlBlk(4,MxFile)
character(len=80) :: FCtlBlk(MxFile)
!-----------------------------------------------
! Define the largest file that can be handled by the IO system (units=Bytes)
#ifdef _I8_
integer(kind=iwp), parameter :: Max_File_Length = 214748364800
#else
integer(kind=iwp), parameter :: Max_File_Length = 2147483647
#endif
! Define minimal block size handled by the IO system (units=Bytes)
integer(kind=iwp), parameter :: MBl_wa = 8, MBl_nwa = 512

public :: Addr, CtlBlk, eBlNme, eEof, eFiMFo, eInErr, eNoMsg, eNtOpn, eTlFn, eTmF, FCtlBlk, FlsSize, FSCB, isFiM, isOpen, LuName, &
          LuNameProf, Max_File_Length, MaxFileSize, MaxSplitFile, MBL, MBl_nwa, MBl_wa, MPUnit, Multi_File, MxFile, NProfFiles, &
          pDesc, pHndle, ProfData, pStat, pWhere, Query, Trace

end module Fast_IO
