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
* Copyright (C) 1991, Per-Olof Widmark                                 *
*               1993,1996,1997, Markus P. Fuelscher                    *
*               1996, Luis Serrano-Andres                              *
*               2012, Victor P. Vysotskiy                              *
************************************************************************
#if defined(_I8_) || defined(_OPENMP)
#define NO_SPLITTING
#endif
      Subroutine DaName_Main(Lu,String,mf,wa)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Open unit Lu for direct access I/O and link the data stream to   *
*     the logical file name Name.                                      *
*                                                                      *
*     calling arguments:                                               *
*     Lu      : integer, input                                         *
*               logical unit number                                    *
*     LuName  : character string, input                                *
*               logical file name                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, IBM Sweden, 1991                                   *
*     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
*     L. Serrano-Andres, University of Lund, Sweden, 1996              *
*----------------------------------------------------------------------*
*                                                                      *
* History:                                                             *
*          V.P. Vysotskiy, University of Lund, Sweden, 2012            *
*                                                                      *
************************************************************************

#ifndef _HAVE_EXTRA_
      Use Prgm
#endif
      Implicit Integer (A-Z)

#include "fio.fh"
#include "switch.fh"
#ifdef _OLD_IO_STAT_
#include "ofio.fh"
#else
#include "pfio.fh"
#endif
#include "blksize.fh"
#include "filesize.fh"
      Character*(*) String
      Logical mf, wa
      Character*8 StdNam
      Character*80 Text

      Character*16 TheName
      Data TheName/'DaName_Main'/
      If ( Query ) Call qEnter(TheName)

      If ( Trace ) then
        Write (6,*) ' >>> Enter DaName_Main <<<'
        Write (6,*) ' unit :',Lu
        Write (6,*) ' name :',String, mf, wa
      End If

      tmp=Lu
      Lu=isfreeunit(tmp)
*     Check calling arguments
      If ( (Lu.le.0) .or. (Lu.gt.MxFile) )
     * Call SysFileMsg(TheName,'MSG: unit', Lu,String)

*     Check for consistency
      If ( isOpen(Lu).ne.0 )
     * Call SysFileMsg(TheName,'MSG: used', Lu,String)

*     Reformat file name to standard notation
*     (capital letters, no leading blank)
      Call StdFmt(String,StdNam)

*     If no file name has been given link it to the default name
      If ( StdNam.eq.'        ' )
     &   Write(StdNam,'(A,I2.2,A)')'FT',Lu,'F001'

*     Check the storage scheme
#ifndef _GA_
      isFiM(Lu)=0
      isFiM(Lu)=isinmem(StdNam)
#ifdef _DEBUG_IO_
      if(isFiM(Lu)>0) write(6,*) "The file ",StdNam," will be kept in
     & memory"
#endif
*     Open file
      temp = isFiM(Lu)
#else
      temp=0
#endif
      iRc = AixOpn(temp,StdNam,.true.)
#ifndef _GA_
      if(iRc.eq.eFiMFo) Then
#ifdef _DEBUG_IO_
         write(6,*) "Failed to open file in memory"
#endif
         isFiM(Lu)=0
         iRc=0
      end if
#endif
      If ( iRc.ne.0 ) then
        iRc = AixErr(Text)
      Call SysFileMsg(TheName,'MSG: open', Lu,Text)
      End If
      isOpen(Lu)  = 1
      FSCB(Lu)    = temp
      LuName(Lu)    = StdNam
#ifndef _OLD_IO_STAT_
      inUse=0
      Do i=1,NProfFiles
         If(LuNameProf(i).eq.StdNam) Then
            inUse=1
         End If
      End Do

      If(inUse.eq.0) Then
        NProfFiles=NProfFiles+1
        LuNameProf(NProfFiles)=StdNam
      End If
#endif
      Call SetLuMark(Lu)
      Addr(Lu)    = 0
      MPUnit(0,Lu)=Lu
#ifdef _OLD_IO_STAT_
      MxAddr(Lu)  = 0
      Count(1,Lu) = 0
      Count(2,Lu) = 0
      Count(3,Lu) = 0
      Count(4,Lu) = 0
#endif
      Multi_File(Lu)=.False.
      MBL(Lu)=MBl_nwa
      if(wa) MBL(Lu)=MBl_wa
#ifndef NO_SPLITTING
      If(mf) then
        Multi_File(Lu)=.True.
        MaxFileSize = AllocDisk()
        If (MaxFileSize.gt.Max_File_Length/(1024**2)) Then
            Write (6,*)
            Write (6,*) 'DANAME_MF: Requested MaxFileSize is too large!'
            Write (6,*) ' Requested value of ',MaxFileSize
            MaxFileSize=Max_File_Length/(1024**2)
            Write (6,*) ' has been reset to  ',MaxFileSize
        Else If ( MaxFileSize.ne.0 ) then
          If ( Trace ) Write (6,*) ' This is a partitioned data set'
          lName = StrnLn(String)
          If ( (lName.eq.0).or.(lName.eq.8) )
     *   Call SysFileMsg(TheName,'Invalid file name. \n '//
     *   'File names used in '//
     *   'multiple unit files must be less than 8 characters.',
     *                  Lu,String)
          Do i = 1,MaxSplitFile-1
            MPUnit(i,Lu)=-99
          End Do
        End If
      End If
#endif
      If ( Trace ) then
        Write (6,*) ' >>> Exit DaName_Main <<<'
      End If

      If ( Query ) Call qExit(TheName)

      Return
      End
      Subroutine SetLuMark(Lu)
#include "fio.fh"
      LuMark(Lu)=0
      return
      end
