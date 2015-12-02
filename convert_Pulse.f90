PROGRAM bin_convert
IMPLICIT NONE

INTEGER :: t1,t2, clock_rate, clock_max
CHARACTER(LEN=1) :: junk
CHARACTER(LEN=10) :: str_c
character(LEN=10) :: x1
INTEGER:: MAX_REC = 1e9
INTEGER:: length1,count1,num_open,sum1,sum2
INTEGER:: ios, int_temp,int_check,int_pos1,int_min1,int_pos2,int_min2
DOUBLE PRECISION:: infinity = 1e30
INTEGER::yc,YMAX,I,J,K,IMAX,JMAX,KMAX,timesteps,RMAX,THMAX,ZMAX,rc,zc,tc,t,I_yp1,I_ym1,I_zp1,I_zm1,I_xp1,I_xm1,temp_rc
INTEGER::write_size,write_end,loop_open,fid_temp,fid_EP_P,fid_EP_G,fid_U, fid_ISO6,fid_GRAD4pt, fid_GRADsize,fid_GradV,fid_GradE  
INTEGER::fid_EP_G_t,fid_U_t,fid_ISO3,fid_Ri,fid_ROP1,fid_ROP2,fid_Dot,fid_c

INTEGER::Z_minus,Z_plus,X_minus,Y_plus,Z_total,I_local

DOUBLE PRECISION,ALLOCATABLE::EP_P(:,:,:),Iso_6(:,:,:)
DOUBLE PRECISION,ALLOCATABLE::Iso_3(:,:,:),four_point_Iso3(:,:,:),four_point(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: T_G(:,:,:),U_G(:,:,:), VEL_6(:,:,:),VEL_3(:,:,:),VEL_ALL(:,:,:), VEL_temp(:,:), VEL_temp2(:,:), Richardson(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: Ri_Dilute(:,:,:), Ri_Dense(:,:,:),char_dense(:), char_dilute(:),VEL_DENSE(:,:), VEL_DILUTE(:,:),ambient_height(:)


DOUBLE PRECISION, ALLOCATABLE :: topo2(:),topography(:),EP_G1(:,:),XXX(:,:),YYY(:,:),ZZZ(:,:),RO_G(:,:)
DOUBLE PRECISION, ALLOCATABLE :: T_G1(:,:), V_G1(:,:),U_G1(:,:),W_G1(:,:),T_S1(:,:),T_S2(:,:),ROP_S1(:,:),ROP_S2(:,:),ROP_S3(:,:)
DOUBLE PRECISION, ALLOCATABLE :: RO_C(:,:,:)

INTEGER, ALLOCATABLE::Location_I(:,:)

DOUBLE PRECISION, ALLOCATABLE :: XX(:,:,:),YY(:,:,:),ZZ(:,:,:),juive_topo(:,:)
DOUBLE PRECISION, ALLOCATABLE :: DY(:),y(:),x(:),theta(:),z(:),DX(:),DTH(:),DZ(:),GZ(:)
DOUBLE PRECISION:: y_boundary,sum_p1, sum_p2, sum_p3, sum_p4,sum_mass,sum_E1,sum_E2,sum_EG
DOUBLE PRECISION:: P_const, R_vapor, R_dryair,char_length, Reynolds_N,rho_g,mu_g,rho_dry, rho_p,T_amb
DOUBLE PRECISION:: max_mag,top,bottom, max_dilute, min_dilute, max_dense, min_dense
DOUBLE PRECISION:: initial_ep_g,ROP1,ROP2,initial_vel,initial_temp,Area_flux,Gas_flux,Solid_flux,gravity,temp_val,local_vel
DOUBLE PRECISION::c_pos1, c_pos2, c_min1, c_min2, delta_V1, delta_V3, shear_v,delta_rho, Ri
DOUBLE PRECISION::Gas_Volume(1,5),Volume_Unit,Energy_In_G, CP_Go, Energy_In_S1,CP_S1o, Energy_In_S2, CP_S2o
DOUBLE PRECISION::mag_grad,norm_grad(1,3),temp_dot,rho_p1,rho_p2,rho_p3,sum_gas

CALL SYSTEM_CLOCK(t1, clock_rate, clock_max )
print *, t1, clock_rate, clock_max

!------------OPEN files-----------!
OPEN(100,FILE='EP_G',form='unformatted')
OPEN(101,FILE='T_G',form='unformatted')
OPEN(102,FILE='U_G',form='unformatted')
OPEN(103,FILE='V_G',form='unformatted')
OPEN(104,FILE='W_G',form='unformatted')
OPEN(105,FILE='RO_G',form='unformatted')

OPEN(106,FILE='ROP_S1',form='unformatted')
OPEN(107,FILE='ROP_S2',form='unformatted')
OPEN(108,FILE='ROP_S3',form='unformatted')


!------------OPEN files-----------!
initial_vel = 152.61
!-------- Boundaries for Gradient Calculations -----------!
max_dilute  = 6.5
min_dilute  = 5.5
!-------- Boundaries for Gradient Calculations -----------!

!--------------------------- Constants ------------------------------------------!
Volume_Unit = 50.*50.*50.
gravity = 9.81 !m^2/s
P_const = 1.0e5 !Pa
R_vapor = 461.5 !J/kg K
R_dryair = 287.058 !J kg/K
T_amb    = 273.0 !K
rho_dry  = P_const/(R_dryair*T_amb)  !kg/m**3
char_length = 100.0
mu_g        = 2.0e-5 !Pa s
rho_p       = 2500.0 !kg/m**3
!--------------------------- Constants ------------------------------------------!

!-------- Set Size, Timesteps, and write size ------------!
timesteps=81  !0
RMAX=164
ZMAX=164
YMAX=204

length1 = RMAX*ZMAX*YMAX

rho_p1  = 2740.0
rho_p2  = 2740.0
rho_p3  = 1950.0

!------------------------ OPEN Gradient Location Files --------------------------!
  ! NOTE: If timesteps are greater than 100 this format will be an issue
print *, "Warning if timesteps ", timesteps, " is greater than 99 will have errors"
str_c = '(I2.2)'
DO t=1,timesteps
  write(x1,str_c) t

  ! ------------- DILUTE PORTION OF CURRENT-----------------!
  num_open = 200+t
  OPEN(num_open,FILE='GRAD_ISO6_t'//trim(x1)//'.txt')   !Used in *.general file for OPENDX Simluation

  !num_open = 300+t
  !OPEN(num_open,FILE='GRAD_4_point_t'//trim(x1)//'.txt')  

  !num_open = 400+t
  !OPEN(num_open,FILE='GRAD_Velocity_t'//trim(x1)//'.txt')

  !num_open = 500+t
  !OPEN(num_open,FILE='Surface_Richardson_t'//trim(x1)//'.txt')

  num_open = 600+t
  OPEN(num_open,FILE='EE_t'//trim(x1)//'.txt')
  

  ! ------------- DILUTE PORTION OF CURRENT-----------------!

  ! ---------------- ENTIRE CURRENT ---------------------!
  num_open = 1100+t
  OPEN(num_open,FILE='EP_t'//trim(x1)//'.txt')

  num_open = 1200+t
  OPEN(num_open,FILE='U_G_t'//trim(x1)//'.txt')

  num_open = 1300+t
  OPEN(num_open,FILE='T_G_t'//trim(x1)//'.txt')

  num_open = 1700+t
  OPEN(num_open,FILE='Richardson_t'//trim(x1)//'.txt')

  num_open = 1800+t
  OPEN(num_open,FILE='Current_Density_t'//trim(x1)//'.txt')
  ! ---------------- ENTIRE CURRENT ---------------------!
 
END DO

!------------------------ OPEN Gradient Location Files --------------------------!

OPEN(701,FILE='EP_G_sum1', form='formatted')
OPEN(702,FILE='EP_G_sum2', form='formatted')
OPEN(703,FILE='EP_G_sum3', form='formatted')

sum_p1 = 0.0 
sum_p2 = 0.0 
sum_p3 = 0.0 
sum_p4 = 0.0 
 
!----------------------Allocate the conversion variables-------------------------!
ALLOCATE( EP_G1(length1,timesteps))
ALLOCATE( T_G1(length1,timesteps))
ALLOCATE( V_G1(length1,timesteps))
ALLOCATE( U_G1(length1,timesteps))
ALLOCATE( W_G1(length1,timesteps))
ALLOCATE( RO_G(length1,timesteps))
ALLOCATE( ROP_S1(length1,timesteps))
ALLOCATE( ROP_S2(length1,timesteps))
ALLOCATE( ROP_S3(length1,timesteps))

ALLOCATE( XXX(length1,1))
ALLOCATE( ZZZ(length1,1))
ALLOCATE( YYY(length1,1))
ALLOCATE( XX(RMAX,YMAX,ZMAX))
ALLOCATE( YY(RMAX,YMAX,ZMAX))
ALLOCATE( ZZ(RMAX,YMAX,ZMAX))
ALLOCATE( x(RMAX))
ALLOCATE( y(YMAX))
ALLOCATE( z(ZMAX))
ALLOCATE( DX(RMAX))
ALLOCATE( DY(YMAX))
ALLOCATE( DZ(ZMAX))
ALLOCATE( ambient_height(YMAX))
!----------------------Allocate the conversion variables-------------------------!

!----------------------Allocate the Gradient variables-------------------------!
ALLOCATE(EP_P(length1,7,timesteps))   ! EP_P, EP_S1,EP_S2, EP_S3
ALLOCATE(T_G(length1,4,timesteps))
ALLOCATE(U_G(length1,6,timesteps))
ALLOCATE(Richardson(length1,5,timesteps))
ALLOCATE(RO_C(length1,7,timesteps))

ALLOCATE(Location_I(length1,timesteps))

ALLOCATE(Iso_6(length1,9,timesteps))

ALLOCATE(four_point(length1,11,timesteps))

ALLOCATE(VEL_temp(length1,6))
ALLOCATE(VEL_temp2(length1,6))

ALLOCATE(VEL_6(RMAX,7,timesteps))
ALLOCATE(VEL_ALL(RMAX,7,timesteps))

ALLOCATE(VEL_DILUTE(7,timesteps))


ALLOCATE(char_dense(timesteps))
ALLOCATE(char_dilute(timesteps))
!----------------------Allocate the Gradient variables-------------------------!

!------------------------ Create spatial deltas and distances -------------------!
DX(1)=50.
x(1)=DX(1)
DO rc=2,RMAX
  DX(rc)=DX(rc-1)
  x(rc)=DX(rc)+x(rc-1)
END DO

DY(1)=50.0!0.375
y(1)=DY(1)
DO zc=2,YMAX
  DY(zc)=DY(zc-1)
  y(zc)=DY(zc)+y(zc-1)
END DO

DZ(1)=50.0
z(1)=DZ(1) !Z is in reverse compared to X & Y
DO zc=2,ZMAX
 DZ(zc)=DZ(zc-1)
 z(zc)=z(zc-1)+DZ(zc)
END DO
!------------------------ Create spatial deltas and distances -------------------!

!------------------------------ Create grid matrices of length RMAX*ZMAX*YMAX ----------------------------!
I=1
DO zc=1,ZMAX
     DO rc=1,RMAX
          DO yc=1,YMAX
                XX(rc,yc,zc)=x(rc)
                XXX(I,1)=XX(rc,yc,zc)
                YY(rc,yc,zc)=y(yc)
                YYY(I,1)=YY(rc,yc,zc)
                ZZ(rc,yc,zc)=z(zc)
                ZZZ(I,1)=ZZ(rc,yc,zc)
                I=I+1
          END DO
     END DO
END DO
!------------------------------ Create grid matrices of length RMAX*ZMAX*YMAX ----------------------------!

!------------------------------ Read and Set to Variable ----------------------------!
REWIND(100)
DO I=1,timesteps
     READ(100) EP_G1(:,I) ! Format is RMAX*YMAX*ZMAX by timesteps, so looping
                          ! through timesteps
END DO

REWIND(101)
DO I=1,timesteps
     READ(101) T_G1(:,I)
END DO

REWIND(102)
DO I=1,timesteps
     READ(102) U_G1(:,I)
END DO

REWIND(103)
DO I=1,timesteps
     READ(103) V_G1(:,I)
END DO

REWIND(104)
DO I=1,timesteps
     READ(104) W_G1(:,I)
END DO

REWIND(105)
DO I=1,timesteps
     READ(105) RO_G(:,I)
END DO

REWIND(106)
DO I=1,timesteps
     READ(106) ROP_S1(:,I)
END DO

REWIND(107)
DO I=1,timesteps
     READ(107) ROP_S2(:,I)
END DO

REWIND(108)
DO I=1,timesteps
     READ(108) ROP_S3(:,I)
END DO

!------------------------------ Read and Set to Variable ----------------------------!


!------------------- Write Variable and set to 3D variable for gradient calculation ------------------!
print *, "Start writing 3D variables"

DO t=1,timesteps
        fid_EP_G_t  = 1100+t
        fid_U_t     = 1200+t
        fid_temp    = 1300+t
        fid_c       = 1800+t
       DO I=1,RMAX*ZMAX*YMAX
          !------------------ Volume Fraction of Gas or Particles ------------------!
          if (T_G1(I,t)==0.0) THEN !  .OR. EP_G1(I,t) > 0.9999999) Try to make sure not to calculate infinity densities or 
                                                               ! density for concentration of gas outside of the current
          RO_C(I,1,t) = 0.0
          else
          RO_C(I,1,t) = EP_G1(I,t)*RO_G(I,t)+ROP_S1(I,t)+ROP_S2(I,t)+ROP_S3(I,t)
          end if
          RO_C(I,2,t) = RO_G(I,t)
          RO_C(I,3,t) = 0.000
          RO_C(I,4,t) = 0.000
          
          RO_C(I,5,t) = XXX(I,1)
          RO_C(I,6,t) = YYY(I,1)
          RO_C(I,7,t) = ZZZ(I,1)

          !------------------------ Temperature of Gas ----------------------------!
          T_G(I,1,t) = T_G1(I,t) 
          T_G(I,2,t) = XXX(I,1)
          T_G(I,3,t) = YYY(I,1)
          T_G(I,4,t) = ZZZ(I,1)
 
          EP_P(I,1,t) = -LOG10(1-EP_G1(I,t)+1e-14)
          EP_P(I,2,t) = XXX(I,1)
          EP_P(I,3,t) = YYY(I,1)
          EP_P(I,4,t) = ZZZ(I,1)
          WRITE(fid_EP_G_t,407) EP_G1(I,t),ROP_S1(I,t)/rho_p1,ROP_S2(I,t)/rho_p2,ROP_S3(I,t)/rho_p3,XXX(I,1),YYY(I,1),ZZZ(I,1)
          !------------------------ Velocity of Gas --------------------------------!
          U_G(I,1,t) = U_G1(I,t)          
          U_G(I,2,t) = V_G1(I,t)          
          U_G(I,3,t) = W_G1(I,t)          
          U_G(I,4,t) = XXX(I,1)
          U_G(I,5,t) = YYY(I,1)
          U_G(I,6,t) = ZZZ(I,1)
          !WRITE(fid_U,300) U_G(I,1:6,t) WRITE LATER DOWN AFTER VELOCITY
          !CORRECTION

       END DO
!------------------- Write Variable and set to 3D variable for gradient calculation ------------------!
print *, "Done writing 3D variables"

!------------------- CALCULATE RICHARDSON #  ------------------!
fid_Ri     = 1700+t

! CORRECT THE LOCATION OF THE X AND Z DIRECTION VELOCITIES--------------------------!
  DO rc =4,RMAX-4           !X
   DO zc = 4,ZMAX-4         !Z
    DO yc = YMAX-4,4,-1     !Y
                ! FUNIJK_GL (LI, LJ, LK) = 1 + (LJ - jmin3) + (LI-imin3)*(jmax3-jmin3+1) &
                ! + (LK-kmin3)*(jmax3-jmin3+1)*(imax3-imin3+1)
             I       = 1 + (yc-1) +(rc-1)*(YMAX-1+1) + (zc-1)*(YMAX-1+1)*(RMAX-1+1) 
             
             I_xm1   = 1 + (yc-1) +(rc-1-1)*(YMAX-1+1) + (zc-1)*(YMAX-1+1)*(RMAX-1+1)
             I_zm1   = 1 + (yc-1) +(rc-1)*(YMAX-1+1) + (zc-1-1)*(YMAX-1+1)*(RMAX-1+1)
            
             U_G(I,1,t) = U_G1(I_xm1,t)
             U_G(I,2,t) = V_G1(I,t)
             U_G(I,3,t) = W_G1(I_zm1,t)
          !------------------------ Temperature of Gas ----------------------------!
          T_G(I,1,t) = T_G1(I,t)-.0098*50.0*REAL(yc)        

    END DO
   END DO
  END DO
! CORRECT THE LOCATION OF THE X AND Z DIRECTION VELOCITIES--------------------------!

     DO I=1,RMAX*ZMAX*YMAX
               Richardson(I,1,t) = 1.000e3
               Richardson(I,2,t) = XXX(I,1)
               Richardson(I,3,t) = YYY(I,1)
               Richardson(I,4,t) = ZZZ(I,1)
          WRITE(fid_temp,400) T_G(I,1:4,t)
     END DO

! --------CALCULATE AVERAGE AIR DENSITY AT EACH HEIGHT--------------------------!
DO yc = 1,YMAX
   ambient_height(yc) = 0.00
END DO


 DO yc = YMAX-2,2,-1     !Y
sum1 = 1
sum_gas = 0.0
   DO rc =2,RMAX-2           !X
     DO zc = 2,ZMAX-2         !Z
     I       = 1 + (yc-1) +(rc-1)*(YMAX-1+1) + (zc-1)*(YMAX-1+1)*(RMAX-1+1)
        IF (EP_G1(I,t)>0.9999999) THEN
                sum_gas = sum_gas + RO_G(I,t)
                sum1 = sum1+1
        END IF
     END DO
   END DO
   
     ambient_height(yc) = sum_gas/sum1 

 END DO
sum1 = 1

 DO yc = YMAX-2,2,-1     !Y
   DO rc =2,RMAX-2           !X
     DO zc = 2,ZMAX-2 
             I       = 1 + (yc-1) +(rc-1)*(YMAX-1+1) +(zc-1)*(YMAX-1+1)*(RMAX-1+1)
             RO_C(I,3,t) = RO_C(I,1,t)/ambient_height(yc)
             RO_C(I,4,t) = ambient_height(yc)
     END DO
   END DO
 END DO

! --------CALCULATE AVERAGE AIR DENSITY AT EACH HEIGHT--------------------------!

! NOW THE VELOCITIES ARE CORRECT FOR THE CALCULATION OF GRADIENT RICHARDSON
  DO rc =4,RMAX-4           !X
   DO zc = 4,ZMAX-4         !Z
    DO yc = YMAX-4,4,-1     !Y
                ! FUNIJK_GL (LI, LJ, LK) = 1 + (LJ - jmin3) + (LI-imin3)*(jmax3-jmin3+1) &
                ! + (LK-kmin3)*(jmax3-jmin3+1)*(imax3-imin3+1)
             I       = 1 + (yc-1) +(rc-1)*(YMAX-1+1) + (zc-1)*(YMAX-1+1)*(RMAX-1+1) 
             I_yp1   = 1 + (yc+1-1) +(rc-1)*(YMAX-1+1) + (zc-1)*(YMAX-1+1)*(RMAX-1+1) !Cell above
             I_ym1   = 1 + (yc-1-1) +(rc-1)*(YMAX-1+1) + (zc-1)*(YMAX-1+1)*(RMAX-1+1) !Cell above
             
 
             IF(T_G(I_ym1,1,t)>272.0 .AND. EP_P(I,1,t)<8.0) THEN
                                !------Y-Richardson-------!
                                   ! Current density at each position
                                int_pos1           = I_yp1
                                int_min1           = I_ym1
                                bottom             = 2*((abs(EP_P(int_pos1,3,t)-EP_P(I,3,t))))
                                c_pos1             = rho_p*(10**-(EP_P(int_pos1,1,t)))+(1-(10**-(EP_P(int_pos1,1,t))))*(P_const/(R_vapor*T_G(int_pos1,1,t))) 
                                c_min1             = rho_p*(10**-(EP_P(int_min1,1,t)))+(1-(10**-(EP_P(int_min1,1,t))))*(P_const/(R_vapor*T_G(int_min1,1,t))) 
                                   ! Shear Velocity components, X & Z
                                delta_V1           = (U_G(int_pos1,1,t)-U_G(int_min1,1,t))/bottom
                                !delta_V3           = (U_G(int_pos1,3,t)-U_G(int_min1,3,t))/bottom
                                !shear_v            = sqrt(delta_V1**2+delta_V3**2)
                                shear_v            = delta_V1
                                delta_rho          =(c_pos1-c_min1)/bottom
                                Ri                 =((-gravity/rho_dry)*delta_rho)/(shear_v**2)
                                
                                if (Ri>5) THEN 
                                Richardson(I,1,t) = 5.0
                                !print *, "Ri", I, t, Ri
                                elseif (Ri<-5) THEN
                                Richardson(I,1,t) = -5.0
                                !print *, "Ri", I, t, Ri
                                else
                                Richardson(I,1,t) = Ri
                                end if

                                Richardson(I,2,t) = XXX(I,1)
                                Richardson(I,3,t) = YYY(I,1)
                                Richardson(I,4,t) = ZZZ(I,1)
                                Richardson(I,5,t) = Ri
                                !------Y-Richardson-------!
           END IF
    END DO
   END DO
  END DO 
  

  DO I=1,RMAX*ZMAX*YMAX
    WRITE(fid_Ri,401) Richardson(I,1:5,t)
    WRITE(fid_U_t,300) U_G(I,1:6,t)
    WRITE(fid_c,407) RO_C(I,1:7,t)
  END DO

END DO
!------------------- CALCULATE RICHARDSON #  ------------------!






!--------------------------------------- DILUTE GRADIENT CALCULATION ---------------------------------!
print *, "Start Gradient Calculation for dilute case"
!------------------------------- Calculate Characteristic Velocity --------------------------------!
DO t=1,timesteps
 VEL_DILUTE(1:7,t)=0.0
END DO


Do t=1,timesteps
   sum1 = 1
   sum2 = 1
   DO I=1,RMAX*ZMAX*YMAX
     IF (EP_P(I,1,t)<min_dilute .AND. EP_P(I,1,t)>0.5) THEN
          VEL_temp(sum1,1:6) = U_G(I,1:6,t)
                VEL_DILUTE(1,t) = VEL_DILUTE(1,t) + VEL_temp(sum1,1)
                VEL_DILUTE(2,t) = VEL_DILUTE(2,t) + VEL_temp(sum1,2)
                VEL_DILUTE(3,t) = VEL_DILUTE(3,t) + VEL_temp(sum1,3)
                VEL_DILUTE(4,t) = VEL_DILUTE(4,t) + VEL_temp(sum1,4)
                VEL_DILUTE(5,t) = VEL_DILUTE(5,t) + VEL_temp(sum1,5)
                VEL_DILUTE(6,t) = VEL_DILUTE(6,t) + VEL_temp(sum1,6) 

          sum1 = sum1 + 1
     END IF


   END DO
     VEL_DILUTE(1,t) = VEL_DILUTE(1,t)/real(sum1-1)
     VEL_DILUTE(2,t) = VEL_DILUTE(2,t)/real(sum1-1)
     VEL_DILUTE(3,t) = VEL_DILUTE(3,t)/real(sum1-1)
     VEL_DILUTE(4,t) = VEL_DILUTE(4,t)/real(sum1-1)
     VEL_DILUTE(5,t) = VEL_DILUTE(5,t)/real(sum1-1)
     VEL_DILUTE(6,t) = VEL_DILUTE(6,t)/real(sum1-1)
     VEL_DILUTE(7,t) = sqrt(VEL_DILUTE(1,t)**2+VEL_DILUTE(2,t)**2+VEL_DILUTE(3,t)**2)
END DO
print *, "Done Velocity"
!------------------------------- Calculate Characteristic Velocity --------------------------------!

!----------- Search entire volume for crossing threshold of Volume Fraction of Particles  ------------!
DO t=1,timesteps
print *, "Done Velocity"
print *, t
fid_ISO6     = 200+t
!fid_GRAD4pt  = 300+t
!fid_GradV    = 400+t
!fid_Ri       = 500+t
fid_Dot      = 600+t

sum1=1
!--------------------- Searching entire volume but avoiding ghost cells and boundary cells -----------!
     DO rc=5,RMAX-4 !Looping through the x axis.
        VEL_6(rc,1:7,t) = 0.0
        sum2=1 
       DO zc=7,ZMAX-7
          DO yc=YMAX-4,4,-1
             I       = 1 + (yc-1) +(rc-1)*(YMAX-1+1) + (zc-1)*(YMAX-1+1)*(RMAX-1+1)
             I_yp1   = 1 + (yc+1-1) +(rc-1)*(YMAX-1+1) + (zc-1)*(YMAX-1+1)*(RMAX-1+1) !Cell above
             I_zp1   = 1 + (yc-1) +(rc-1)*(YMAX-1+1) + (zc+1-1)*(YMAX-1+1)*(RMAX-1+1) !Cell to the left in Z direction
             I_zm1   = 1 + (yc-1) +(rc-1)*(YMAX-1+1) + (zc-1-1)*(YMAX-1+1)*(RMAX-1+1) !Cell to the right in Z direction
             I_xm1   = 1 + (yc-1) +(rc-1-1)*(YMAX-1+1) + (zc-1)*(YMAX-1+1)*(RMAX-1+1) !Cell in front in the X direction
             I_xp1   = 1 + (yc-1) +(rc-1+1)*(YMAX-1+1) + (zc-1)*(YMAX-1+1)*(RMAX-1+1) !Cell in behind in the X direction
             IF (RO_C(I,1,t) < rho_dry .AND. EP_P(I,1,t)<max_dilute .AND. EP_P(I,1,t)>min_dilute) THEN  
             ! Making sure this portion is buoyant, and within a small range of volume fraction of particles
               IF (EP_P(I_yp1,1,t)>max_dilute .or. EP_P(I_zp1,1,t)>max_dilute .or. EP_P(I_zm1,1,t)>max_dilute .or. EP_P(I_xm1,1,t)>max_dilute .or. EP_P(I_xp1,1,t)>max_dilute) THEN
                        !Reset local position
                        Z_minus  = 0;
                        Z_plus   = 0;
                        X_minus  = 0;
                        Y_plus   = 0; 
                        Z_total  = 0; 
                        if (EP_P(I_yp1,1,t)>max_dilute) Y_plus = -1
                        if (EP_P(I_zp1,1,t)>max_dilute) Z_plus = -1
                        if (EP_P(I_zm1,1,t)>max_dilute) Z_minus = 1
                        if (EP_P(I_xm1,1,t)>max_dilute) X_minus = 1
                        if (Z_plus ==-1 .AND. Z_minus ==1) print *, "Error Z plus and minus",I,t
                         Z_total = Z_plus + Z_minus 
                        !Make sure the cells below are not below topography for gradient calculation
                        I_local = 1 + (yc+Y_plus-1) +(rc+X_minus-1)*(YMAX-1+1) + (zc+Z_total-1)*(YMAX-1+1)*(RMAX-1+1)
                  IF(T_G(I-1,1,t)>272.0 .AND. T_G(I-2,1,t) >272.0 .AND. EP_P(I_local,1,t)<max_dilute) THEN
                         
                 ! Making sure this is the point of change into the dilute region of interest
                        Iso_6(sum1,1:4,t) = EP_P(I,:,t)
                        Iso_6(sum1,5,t)   = DBLE(I)
                        Iso_6(sum1,6,t)   = T_G(I,1,t)
                        Iso_6(sum1,7,t)   = RO_C(I,1,t)
                        Iso_6(sum1,8,t)   = rho_dry
                        Iso_6(sum1,9,t)   = 1000.00
                        !WRITE(fid_GRADV,300) U_G(I,1:6,t)
                        VEL_6(rc,1,t) = VEL_6(rc,1,t) + U_G(I,1,t)
                        VEL_6(rc,2,t) = VEL_6(rc,2,t) + U_G(I,2,t)
                        VEL_6(rc,3,t) = VEL_6(rc,3,t) + U_G(I,3,t)
                        VEL_6(rc,4,t) = VEL_6(rc,4,t) + U_G(I,4,t)
                        VEL_6(rc,5,t) = VEL_6(rc,5,t) + U_G(I,5,t)
                        VEL_6(rc,6,t) = VEL_6(rc,6,t) + U_G(I,6,t)
                       ! ----------------------- 4-point central difference stencil -----------------------!
                                int_temp = I
                                tc       = sum1
                                four_point(tc,1:5,t) = Iso_6(tc,1:5,t)
                                !------X-direction-------!
                                int_pos1           = int_temp+YMAX
                                int_pos2           = int_temp+2*YMAX
                                int_min1           = int_temp-YMAX
                                int_min2           = int_temp-2*YMAX
                                top                = (EP_P(int_min2,1,t))-8*(EP_P(int_min1,1,t))+8*(EP_P(int_pos1,1,t))-(EP_P(int_pos2,1,t))
                                bottom             = 12*((abs(EP_P(int_pos1,2,t)-EP_P(int_temp,2,t))))
                                four_point(tc,6,t) = top/bottom
                                !------X-direction-------!

                                !------Y-direction-------!
                                int_pos1           = int_temp+1
                                int_pos2           = int_temp+2*1
                                int_min1           = int_temp-1
                                int_min2           = int_temp-2*1
                                local_vel          = sqrt(U_G(I_local,1,t)**2+U_G(I_local,2,t)**2+U_G(I_local,3,t)**2)
                                top                = (EP_P(int_min2,1,t))-8*(EP_P(int_min1,1,t))+8*(EP_P(int_pos1,1,t))-(EP_P(int_pos2,1,t))
                                bottom             = 12*((abs(EP_P(int_pos1,3,t)-EP_P(int_temp,3,t))))
                                four_point(tc,7,t) = top/bottom
                                !------Y-direction-------!

                                !------Y-Richardson-------!
                                   ! Current density at each position
                                int_pos1           = int_temp+1
                                int_pos2           = int_temp+2
                                int_min1           = int_temp-1
                                int_min2           = int_temp-2
                                if (T_G(int_min1,1,t) > 0.0) then
                                bottom             = 2*((abs(EP_P(int_pos1,3,t)-EP_P(int_temp,3,t))))
                                c_pos1             = rho_p*(10**-(EP_P(int_pos1,1,t)))+(1-(10**-(EP_P(int_pos1,1,t))))*(P_const/(R_vapor*T_G(int_pos1,1,t))) 
                                c_min1             = rho_p*(10**-(EP_P(int_min1,1,t)))+(1-(10**-(EP_P(int_min1,1,t))))*(P_const/(R_vapor*T_G(int_min1,1,t))) 
                                   ! Shear Velocity components, X & Z
                                delta_V1           = (U_G(int_pos1,1,t)-U_G(int_min1,1,t))/bottom
                                !delta_V3           = (U_G(int_pos1,3,t)-U_G(int_min1,3,t))/bottom
                                !shear_v            = sqrt(delta_V1**2+delta_V3**2)
                                shear_v            = delta_V1
                                delta_rho          =(c_pos1-c_min1)/bottom
                                Ri                 =((-gravity/rho_dry)*delta_rho)/(shear_v**2)

                                if (Ri > 5.0) then
                                Ri = 5.0
                                elseif (Ri< -5.0) then
                                Ri = -5.0
                                end if
                                ISO_6(tc,9,t) = Ri
                                !WRITE(fid_Ri,301) Iso_6(tc,1:4,t),c_pos1,c_min1,Ri,delta_V1,delta_V3,T_G(int_pos1,1,t),T_G(int_min1,1,t)
                                end if
                                !------Y-Richardson-------!

                                !------Z-direction-------!
                                int_pos1           = int_temp-RMAX*YMAX
                                int_pos2           = int_temp-2*RMAX*YMAX
                                int_min1           = int_temp+RMAX*YMAX
                                int_min2           = int_temp+2*RMAX*YMAX
                                top                = (EP_P(int_min2,1,t))-8*(EP_P(int_min1,1,t))+8*(EP_P(int_pos1,1,t))-(EP_P(int_pos2,1,t))
                                bottom             = 12*((abs(EP_P(int_pos1,4,t)-EP_P(int_temp,4,t))))
                                four_point(tc,8,t) = top/bottom
                                !------Z-direction-------!
                                Location_I(sum1,t)= rc
                                four_point(tc,9,t) = DBLE(rc)
                        ! ----------------------- 4-point central difference stencil -----------------------!
                        sum1 = sum1+1
                        sum2 = sum2+1
                        mag_grad  = sqrt(four_point(tc,6,t)**2 + four_point(tc,7,t)**2 + four_point(tc,8,t)**2)
                        norm_grad(1,1) = four_point(tc,6,t)/mag_grad
                        norm_grad(1,2) = four_point(tc,7,t)/mag_grad
                        norm_grad(1,3) = four_point(tc,8,t)/mag_grad
                        temp_dot  = dot_product(U_G(I,1:3,t),norm_grad(1,1:3))
                        WRITE(fid_Dot,302)four_point(tc,1:8,t),initial_vel,norm_grad(1,1:3),U_G(I,1:3,t),temp_dot/initial_vel 
                        WRITE(fid_ISO6,308) Iso_6(tc,1:9,t)
                        !WRITE(fid_Dot,302)four_point(tc,1:8,t),VEL_DILUTE(7,t),norm_grad(1,1:3),U_G(I,1:3,t),temp_dot/VEL_DILUTE(7,t)
                        !WRITE(fid_Dot,302)four_point(tc,1:8,t),mag_grad,norm_grad(1,1:3),U_G(I,1:3,t),temp_dot/local_vel 
                 END IF
                END IF
             END IF
          END DO
       END DO
         !WRITE(3825,806) VEL_ALL(rc,1:7,t), rc, t, sum2
     IF (sum2==1) THEN
           !WRITE(fid_temp,407) VEL_6(rc,1:7,t),rc,t
     ELSE
           VEL_6(rc,1,t) = VEL_6(rc,1,t)/real(sum2-1)
           VEL_6(rc,2,t) = VEL_6(rc,2,t)/real(sum2-1)
           VEL_6(rc,3,t) = VEL_6(rc,3,t)/real(sum2-1)
           VEL_6(rc,4,t) = VEL_6(rc,4,t)/real(sum2-1)
           VEL_6(rc,5,t) = VEL_6(rc,5,t)/real(sum2-1)
           VEL_6(rc,6,t) = VEL_6(rc,6,t)/real(sum2-1)
           VEL_6(rc,7,t) = sqrt(VEL_6(rc,1,t)**2+VEL_6(rc,2,t)**2+VEL_6(rc,3,t)**2)
           !WRITE(fid_temp,407) VEL_6(rc,1:7,t),rc,t
     END IF
    END DO

 DO K=1,sum1-1
   temp_rc = int(four_point(K,9,t))
   !J = Location_I(K,t)
   !print *, K, temp_rc, t
   four_point(K,10,t) = VEL_ALL(temp_rc,7,t)
   !print *, "Vel ALL (temp_rc)", VEL_ALL(temp_rc,1:7,t)
   !print *, "Vel ALL (J)", VEL_ALL(J,1:7,t)
   four_point(K,11,t) = VEL_6(temp_rc,7,t)
   !WRITE(fid_GRAD4pt,311) four_point(K,1:11,t)
 END DO 
END DO ! End timesteps
print *, "End Gradient Calculation"
!--------------------- Searching entire volume but avoiding ghost cells and boundary cells -----------!
!--------------------------------------- DILUTE GRADIENT CALCULATION ---------------------------------!



!-------------------------CALCULATE NUMBER OF GRIDS WITH SPECIFIC VOLUME FRACTION OF GAS -------------!
DO t=1,timesteps
sum_p1 = 0.0
  DO I=1,RMAX*ZMAX*YMAX
     IF (EP_G1(I,t) <0.999999) THEN
     sum_p1 = sum_p1 + 1.0
     END IF
  END DO
  WRITE(701,802) t, sum_p1
END DO

DO t=1,timesteps
sum_p2 = 0.0
  DO I=1,RMAX*ZMAX*YMAX
     IF (EP_G1(I,t) <0.99999) THEN
     sum_p2 = sum_p2 + 1.0
     END IF
  END DO
  WRITE(702,802) t, sum_p2
END DO

DO t=1,timesteps
sum_p3 = 0.0
  DO I=1,RMAX*ZMAX*YMAX
     IF (EP_G1(I,t) <0.9999) THEN
     sum_p3 = sum_p3 + 1.0
     END IF
  END DO
  WRITE(703,802) t, sum_p3
END DO
!-------------------------CALCULATE NUMBER OF GRIDS WITH SPECIFIC VOLUME FRACTION OF GAS -------------!
print *, "Program Complete"
CALL SYSTEM_CLOCK(t2, clock_rate, clock_max )
            print *, t2, clock_rate, clock_max
            print *, "Elapsed Timei = ", real(t2-t1)/real(clock_rate)/60.0, "minutes"

100 FORMAT(F22.12,3ES22.5,F22.12)
101 FORMAT(i7)
102 FORMAT(i7,1X(3ES22.5),1X,i7)
103 FORMAT(i7,1X(6ES22.5),1X,i7)
400 FORMAT(4F22.12)
401 FORMAT(4F22.12,1x,ES22.5)
300 FORMAT(6F22.12)
301 FORMAT(11F22.10)
302 FORMAT(8F22.10,ES22.7,7F22.10)
502 FORMAT(10F22.10)
308 FORMAT(9F22.10)
311 FORMAT(11F22.12)
802 FORMAT(i7,F22.5)
805 FORMAT(6F22.5,3i7)
806 FORMAT(7F22.5,3i7)
803 FORMAT(1X(5ES22.5),1X(3F22.12))
407 FORMAT(7F22.12)

END PROGRAM


