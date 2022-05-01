      subroutine fcn( n, x, f, g )
      integer n
      double precision x(n), f, g(n)

c Rosenbrock 
      f = 100.*((x(2) - x(1)**2)**2) + (1. - x(1))**2
      g(1) = 200*(x(2) - x(1)**2)*(-2*x(1)) - 2*(1 - x(1))
      g(2) = 200*(x(2) - x(1)**2)

      return
      end

