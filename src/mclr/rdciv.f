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
      Subroutine RdCIV
************************************************************************
*                                                                      *
*     Read the contents of the JOBIPH file.                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "negpre.fh"
#include "Files_mclr.fh"
      Call DaName(LuCIV,'ROOTS')
      iDisk=0
*----------------------------------------------------------------------*
*     Load the CI vector for all roots                            *
*----------------------------------------------------------------------*
      iDisk=iToc(4)
      idisk1=0
      Do i=1,lroots
       Call GetMem('OCIvec','Allo','Real',ipTmp,nConf)
       Call dDaFile(LuJob,2,Work(ipTmp),nConf,iDisk)
       Call Gugactl(ipTmp,1)
       Call dDafile(LuCIV,1,Work(ipTmp),nconf,iDisk1)
       Call GetMem('OCIVEC','Free','Real',ipTmp,nConf)
      End Do
*----------------------------------------------------------------------*
*     Load state energy                                                *
*----------------------------------------------------------------------*

      Call GetMem('Temp2','Allo','Real',ipTmp,mxRoot*mxIter)
      iDisk=iToc(6)
      Call dDaFile(LuJob,2,Work(ipTmp),mxRoot*mxIter,iDisk)
      Do i=0,lroots-1
      ERAS(i+1)=0.0d0
      Do  iter=0,mxIter-1
         Temp=Work(ipTmp+iter*mxRoot+i)
         If ( Temp.ne.0.0D0 ) ERAS(i+1)=Temp
      End Do
      End Do
      Call GetMem('Temp2','Free','Real',ipTmp,mxRoot*mxiter)
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
