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
* Copyright (C) 2020, Jie J. Bao                                       *
************************************************************************
      subroutine prrotmat(NGRP,U0,HEFF,NSTATE,Silent)

      INTEGER NGRP,NSTATE
      real(8) Heff(Nstate,Nstate)
      real(8) U0(Nstate,Nstate)
      Logical Silent

      LOGICAL FOUND
      INTEGER LUXMS,IsFreeUnit
      CHARACTER(len=128) filename,swapname
      CHARACTER(len=11)xmsfmt1
      CHARACTER(len=12)xmsfmt2
      External IsFreeUnit

      if(.not.silent) Then
      write(6,*) 'Writing Hamiltonian matrix for rotated states ',
     &'and the rotation matrix'
      End if
      write(FileName,'(a)') 'ROT_VEC'
      write(SwapName,'(a)') 'ROT_VEC0'
      Call F_Inquire(FileName,Found)
      If(Found)  Call RENAME(FileName,SwapName)
      LUXMS=233
      LUXMS=IsFreeUnit(LUXMS)
      Call Molcas_Open(LUXMS,FileName)
      if(NGRP.LT.10) then
       write(xmsfmt1,'(a4,I1,a6)') "(1x,",NGRP,"F16.8)"
       DO J=1,NGRP
        WRITE(LUXMS,xmsfmt1)(U0(I,J),I=1,NGRP)
       END DO
      else if(NGRP.LT.100) then
       write(xmsfmt2,'(a4,I2,a6)') "(1x,",NGRP,"F16.8)"
       DO J=1,NGRP
        WRITE(LUXMS,xmsfmt2)(U0(I,J),I=1,NGRP)
       END DO
      end if
      close (LUXMS)
      write(FileName,'(a)') 'ROT_HAM'
      write(SwapName,'(a)') 'ROT_HAM0'
      Call F_Inquire(FileName,Found)
      If(Found)  Call RENAME(FileName,SwapName)
      LUXMS=IsFreeUnit(LUXMS)
      Call Molcas_Open(LUXMS,FileName)
      if(NGrp.lt.10)then
       DO J1=1,NSTATE
        WRITE(LUXMS,xmsfmt1)(HEFF(J1,J2),J2=1,NSTATE)
       END DO
      else
       DO J1=1,NSTATE
        WRITE(LUXMS,xmsfmt2)(HEFF(J1,J2),J2=1,NSTATE)
       END DO
      endif
      Close(LUXMS)

      end subroutine
