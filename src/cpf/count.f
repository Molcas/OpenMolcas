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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE COUNT(NINTGR,NSYM,NORB,MUL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NORB(*),MUL(8,8)
C     COUNT TWO-ELECTRON INTEGRALS
      NINTGR=0
      DO 313 NSP=1,NSYM
      NOP=NORB(NSP)
      DO 312 NSQ=1,NSP
      NSPQ=MUL(NSP,NSQ)
      NOQ=NORB(NSQ)
      DO 311 NSR=1,NSP
      NSPQR=MUL(NSPQ,NSR)
      NOR=NORB(NSR)
      NSSM=NSR
      IF(NSR.EQ.NSP)NSSM=NSQ
      DO 310 NSS=1,NSSM
      IF(NSS.NE.NSPQR)GO TO 310
      NOS=NORB(NSS)
      NORBP=NOP*NOQ*NOR*NOS
      IF(NORBP.EQ.0)GO TO 310
      DO 309 NV=1,NOR
      NXM=NOS
      IF(NSR.EQ.NSS)NXM=NV
      DO 308 NX=1,NXM
      NTM=1
      IF(NSP.EQ.NSR)NTM=NV
      DO 307 NT=NTM,NOP
      NUMIN=1
      IF(NSP.EQ.NSR.AND.NT.EQ.NV)NUMIN=NX
      NUMAX=NOQ
      IF(NSP.EQ.NSQ)NUMAX=NT
      DO 306 NU=NUMIN,NUMAX
      NINTGR=NINTGR+1
306   CONTINUE
307   CONTINUE
308   CONTINUE
309   CONTINUE
310   CONTINUE
311   CONTINUE
312   CONTINUE
313   CONTINUE
      RETURN
      END
