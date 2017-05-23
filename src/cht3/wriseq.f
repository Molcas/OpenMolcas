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
      SUBROUTINE WRISEQ(G,IDIMENS,IFILE)
      implicit none
      integer IDIMENS,IFILE,IAS
      REAL*8 G(IDIMENS)
      integer INDI(IDIMENS)
C     write(6,*)'wrote',ifile,idimens
      WRITE(IFILE)G
      RETURN

      ENTRY REASEQ(G,IDIMENS,IFILE)
C     write(6,*)'reading',ifile,idimens
      READ(IFILE)G
      RETURN

      ENTRY WRIDIR(G,IDIMENS,IFILE,IAS)
      WRITE(IFILE,REC=IAS)G
      RETURN

      ENTRY READIR(G,IDIMENS,IFILE,IAS)
      READ(IFILE,REC=IAS)G
      RETURN


C     WRITED BY J. SIMUNEK JULY 2006
C     START PART

      ENTRY WRI_I(INDI,IDIMENS,IFILE,IAS)
      WRITE(IFILE,REC=IAS)INDI
      RETURN

      ENTRY REA_I(INDI,IDIMENS,IFILE,IAS)
      READ(IFILE,REC=IAS)INDI
      RETURN

      ENTRY WRI_IR(G,INDI,IDIMENS,IFILE,IAS)
      WRITE(IFILE,REC=IAS)INDI,G
      RETURN

      ENTRY REA_IR(G,INDI,IDIMENS,IFILE,IAS)
      READ(IFILE,REC=IAS)INDI,G
      RETURN
C     END PART
      END
