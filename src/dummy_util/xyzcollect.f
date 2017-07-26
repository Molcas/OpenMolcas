#ifndef _HAVE_EXTRA_

      Subroutine XYZcollect(iCoord, nCoord, OrigTrans, OrigRot,
     &                      nFragment)
      Integer :: iCoord, nCoord
      Real*8 :: OrigTrans(3,nFragment), OrigRot(3,3,nFragment)
      Call Unused_Integer(iCoord)
      Call Unused_Integer(nCoord)
      Call Unused_Real_Array(OrigTrans)
      Call Unused_Real_Array(OrigRot)
      Call Unused_Integer(nFragment)
      End Subroutine XYZcollect

#endif
