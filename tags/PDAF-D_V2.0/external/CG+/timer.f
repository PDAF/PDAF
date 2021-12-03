
      subroutine timer(ttime)
      double precision ttime
c     *********
c
c     Subroutine timer
c
c     This subroutine is used to determine user time. In a typical 
c     application, the user time for a code segment requires calls 
c     to subroutine timer to determine the initial and final time.
c
c     The subroutine statement is
c
c       subroutine timer(ttime)
c
c     where
c
c       ttime is an output variable which specifies the user time.
c
c     Argonne National Laboratory and University of Minnesota.
c     MINPACK-2 Project.
c
c     Modified October 1990 by Brett M. Averick.
c
c     **********
      real temp
      real tarray(2)
      real etime

c     The first element of the array tarray specifies user time

      temp = etime(tarray) 

      ttime = dble(tarray(1))
 
      return

      end
      
c====================== The end of timer ===============================
