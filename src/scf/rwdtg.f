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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
************************************************************************
      SubRoutine RWDTG(Num,DMat,lth,Option,DT,iDisk,MaxNum)
************************************************************************
*                                                                      *
* purpose: Read / write density matrix, two-electron hamiltonian       *
*          or gradient to the file if they are on disk / can't be      *
*          stored in core                                              *
*                                                                      *
* input:                                                               *
*   Num       - 'map number'                                           *
*   DMat(lth) - matrix to be written or read (Dens or TwoHam)          *
*   Option    - R- read, W- write                                      *
*   DT        - 'DENS  ', 'dVxcdR', 'TWOHAM', 'GRAD  '                 *
*   iDisk     - position in the file                                   *
*   MaxNum    - maximal number of density matrices that can be         *
*               stored on the disk                                     *
*                                                                      *
* output:                                                              *
*   DMat      - matrix read from the disk (if the option was read)     *
*                                                                      *
* called from: Diis, DMat, GetMat, GrdClc, MinDns, TraClc              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* written by:                                                          *
* P.O. Widmark, M.P. Fuelscher and P. Borowski                         *
* mod. by M. Schuetz                                                   *
* University of Lund, Sweden, 1992                                     *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Real*8 DMat(lth)
      Integer iDisk(MaxNum)
      Character Option,DT*6
#include "file.fh"

#include "SysDef.fh"
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
#ifdef _DEBUG_
#endif
*
* Check Num; Subroutine is called with Num = - MapDns(i)
*
      If (Num.le.0) Then
         Write (6,*) 'RWDTG: Num.le.0'
         Write (6,*) 'Num=',Num
         Write (6,*) 'Wrong density number supplied.'
         Call Abend()
      End If
      If (Num.gt.MaxNum) Then
         Write (6,*) 'RWDTG: Num.gt.MaxNum'
         Write (6,*) 'Num,MaxNum=',Num,MaxNum
         Write (6,*) 'Wrong density number supplied.'
         Call Abend()
      End If

      If (DT.ne.'DENS  '.and.DT.ne.'TWOHAM'.and.DT.ne.'GRAD  ' .AND.
     &    DT.ne.'dVxcdR') Then
         Write (6,*) 'RWDTG: invalid value of DT'
         Write (6,*) '->DT<-=->',DT,'<-'
         Write (6,*) 'Valid values: "DENS  "'
         Write (6,*) '              "dVxcdR"'
         Write (6,*) '              "TWOHAM"'
         Write (6,*) '              "GRAD  "'
         Call Abend()
      End If

      If (Option.ne.'W'.and.Option.ne.'R') Then
         Write (6,*) 'RWDTG: invalid Option'
         Write (6,*) '->Option<-=->',Option,'<-'
         Write (6,*) 'Valid Options: R'
         Write (6,*) '               W'
      End If
*
      If (DT.eq.'DENS  ') Then
         LU =  LuDSt
      Else If (DT.eq.'TWOHAM') Then
         LU = LuTSt
      Else If (DT.eq.'GRAD  ') Then
         LU = LuGrd
      Else
         LU = LuOSt
      End If
*
      If (Option.eq.'W') Then
*
*        Write density matrix to DNS
*
         If (Num.eq.1) iDisk(Num) = 0
         jDisk = iDisk(Num)
         If (jDisk.eq.-1) Then
            Write (6,*) 'RWDTG: jDisk.eq.-1'
            Write (6,*) 'Num,MaxNum=',Num,MaxNum
            Write (6,*) 'The preceeding block was not written.'
            Call Abend()
         End If
         Call dDaFile(LU,1,DMat,lth,jDisk)
         If (Num + 1.le.MaxNum) iDisk(Num + 1) = jDisk
*
      Else If (Option.eq.'R') Then
*
*        Read density matrix from DNS
*
         jDisk = iDisk(Num)
         Call dDaFile(LU,2,DMat,lth,jDisk)
      End If
#ifdef _DEBUG_
#endif
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
