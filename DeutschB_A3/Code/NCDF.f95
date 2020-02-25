program NCDF

	use numType
	implicit none
	real(dp) :: x, xmin, xmax, dx, a, F
	real(dp) :: mu, sig
	integer :: i, imax
	
	xmin = -4
	xmax = 8
	imax = 50
	dx =(xmax - xmin)/imax
	
	mu = 2._dp
	sig = 1.5_dp
	
	do i = 0, imax
		x = xmin + i*dx
		
		a = (x-mu)/(sig*sqrt(2._dp))
		F = 1/2._dp*(1+erf_cf(a))
		
		
	open(unit = 1, file='NCDF.dat', status='unknown')
	write(1,*) x, F	
		
	end do
	
	contains
	
		function w_cf(z) result(ss)
		
			implicit none
			complex(dp) :: z, ss
			integer :: nn, n
			
			nn = 10000
			ss = 0
			
			do n = nn, 1, -1
             ss = n/2._dp/(z-ss)
                
            end do    
            if (z == 0) then 
            
            	ss = 1._dp	
            else 
            	ss = iic/sqrt(pi) * 1/(z-ss)
            end if	
		
		end function w_cf
	
		function erfc_cf(x) result(ss)
			implicit none
			real(dp) :: x, ss, xx
			
			xx = abs(x)
			
			ss = exp(-xx**2) * w_cf(iic*xx)			
			
			if(x < 0) ss = 2 - ss
		
		end function erfc_cf
		
		function erf_cf(x) result(ss)
			implicit none
			real(dp) :: x, ss
			
			ss = 1 - erfc_cf(x)
			
		end function erf_cf

end program NCDF