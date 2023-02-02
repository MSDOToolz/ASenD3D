module ASenD_constantVals

  ! Collection of constant values used throughout the structural modelling process, defined in either real or complex
    ! format.

    real*8 :: r_0 = 0d0
    real*8 :: r_1 = 1d0
    real*8 :: r_2 = 2d0
    real*8 :: r_p5 = 0.5d0
    real*8 :: r_p25 = 0.25d0
    real*8 :: r_p125 = 0.125d0
    real*8 :: r_1o6 = 1d0/6d0
    real*8 :: r_1o3 = 1d0/3d0
    real*8 :: r_1rt3 = 1d0/sqrt(3d0)
    real*8 :: r_1rt2 = 1d0/sqrt(2d0)
	real*8 :: r_pi = 3.141592653589793
    real*8 :: r_pi180 = 0.0174532925199433
    real*8 :: r_penFact = 1e+30_8
    real*8 :: r_cgTol = 1e-15_8
    real*8 :: r_0vec(3) = (/0d0,0d0,0d0/)
	real*8 :: compStepInv = 1e+14_8
	
	complex*16 :: c_0 = (0d0,0d0)
    complex*16 :: c_1 = (1d0,0d0)
    complex*16 :: c_2 = (2d0,0d0)
    complex*16 :: c_p5 = (0.5d0,0d0)
    complex*16 :: c_p25 = (0.25d0,0d0)
    complex*16 :: c_p125 = (0.125d0,0d0)
    complex*16 :: c_1o6 = (0.166666666666666d0,0d0)
    complex*16 :: c_1o3 = (0.333333333333333d0,0d0)
    complex*16 :: c_1rt3 = (0.577350269189626,0d0)
    complex*16 :: c_1rt2 = (0.707106781186547,0d0)
	complex*16 :: c_pi = (3.14159265358979323846,0d0)
    complex*16 :: c_pi180 = (0.0174532925199433,0d0)
    complex*16 :: c_penFact = (1e+30_8,0d0)
    complex*16 :: c_cgTol = (1e-30_8,0d0)
    complex*16 :: c_0vec(3) = (/(0d0,0d0),(0d0,0d0),(0d0,0d0)/)
	complex*16 :: compStep = (0d0,1e-14_8)
	
end module ASenD_constantVals