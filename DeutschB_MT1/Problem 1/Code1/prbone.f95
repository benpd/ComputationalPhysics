program one
	use numtype
	implicit none
	real(dp) :: x, xmin, xmax, dx, f, g 
	integer :: i, imax
	
! step size selection in problem 	
	xmin = -3
	xmax = 3
	imax = 50
	dx =(xmax - xmin)/imax
	
	do i = 0, imax
		x = xmin + i*dx
		
! correction of function @ 0 
		if (x==0)  f = 1
		if (x > 0) f = g(x**2) / (sqrt(pi)) 
	 	if (x < 0) f = 2 - (g(x**2)) / (sqrt(pi)) 
	 	
! write 
		write (777,*) x, f	
		
	end do
end program one 
	
	function g(z) result(f) 
		use  numtype
		implicit none
		real(dp) :: z, f, a 
		integer :: n, nmax 			
		
		f = 0
		a = 0.5
		
		nmax = 1000 
		do n = nmax, 1, -1 
		
			f = z +((n-a)/(1 + (n/f)))
		end do 
			f = exp(-z)*(z**a)*1/f 
	
	end function g 