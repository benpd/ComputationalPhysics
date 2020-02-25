module setup
	use numtype
	implicit none
	
	integer, parameter :: n_eq = 12		  				!!!EDITED
	real(dp), parameter :: grav = 6.673e-11_dp
	real(dp), parameter :: mass_sun = 1.9891e+30_dp
	real(dp), parameter :: mass_earth = 5.9736e+24_dp
	
	real(dp), parameter :: mass_moon = 7.347e+22_dp		!!!EDITED
	real(dp), parameter :: etom = 3.844e+8_dp			!!!EDITED
	
end module setup

program SEM
	use setup
	implicit none
	real(dp) :: t, dt, tmax
	real(dp), dimension(n_eq) :: y
	

	 
	tmax = 60*60*24*365								    !!!EDITED
	
	t = 0._dp
	dt = 60*60*24										!!!EDITED
	
	
	
	!Earth pos
	y(1) = 1.496e+11_dp 			
	y(2) = 0._dp 
	y(3) = 0._dp
	!velo
	y(4) = 0._dp
	y(5) = 29.783e+3_dp 
	y(6) = 0._dp
	
	!Moon pos
	y(7) = y(1) + etom							!!!EDITED
	y(8) = 0._dp								!!!EDITED
	y(9) = 0._dp								!!!EDITED
	!velo
	y(10) =	0._dp								!!!EDITED
	y(11) =	y(5) + ((2*pi*etom)/(27*24*60*60))	!!!EDITED
	y(12) = 0._dp								!!!EDITED
	
	
	do while (t<tmax) 

		write(11,*) 0._dp, 0._dp
		
		write(12,*) y(1), y(2)
		
		write(13,*) y(7), y(8)							!!!EDITED
		
		if  (t<2628000) then							!!!EDITED
		
			write(14,*) y(7)-y(1), y(8)-y(2)  			!!!EDITED
		
		endif 
		
		call rkf45step(t, y, dt)
	
	
	end do


end program SEM