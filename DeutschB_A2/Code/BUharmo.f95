Program BUharmo
	
!setup
	use numtype
	implicit none
	
!variables
	real(dp) :: x, E
	real(dp), parameter :: w    = 2
	real(dp), parameter :: hbar = 1
	real(dp), parameter :: m    = 1
	real(dp), dimension (0:10) ::  psi
	integer :: j, i, n


!energy levels 
	do j = 0, 10 
	
!step size
		do i = -500, 500 
		
			x = float(i)*.01 

!hermitian of wave eq		
			do n = 0, 10 
			
			  psi(0) = (m*w/hbar/pi)**(.25) * exp(-m*w*x**2/2/hbar)
 			
 			psi(n+1) = (2*sqrt(m*w/hbar)*x)/((2.0*(n+1.0))**(.5))*psi(n)& 
 				 		   -(float(n)/(n+1.0))**(.5)*psi(n-1) 
 				
				   E = hbar*w*(float(j)+.5)
 			
	    	end do 

!write statement inside loop
!to plot each point 
	    write (10+j,*) x , psi(j) + E  
	    print *, psi(j) 
	    end do 
	     
	end do 
	
End program BUharmo