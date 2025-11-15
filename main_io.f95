MODULE POINT   ! Define the variables for each point   
! It is just the style of programming, but here in Fortran I prefer to use 'allocable' to 'pointer' since 
! if a pointer is not deallocated, there will be memory leaks, while an allocable variable will be automatically deallocated. 
 IMPLICIT NONE 
 TYPE :: POINT_COE
 	REAL (KIND = 8) :: CTH,CTHN  ! Concentration at old and new time steps
 	REAL (KIND = 8) :: KTH(2,2)   ! Dispersion tensor
 	REAL (KIND = 8) :: X,Y  ! X and Y coordinate
 	REAL (KIND = 8) :: XY_PXX1  ! Sum of the derivatives related to the dispersion 
 	REAL (KIND = 8), ALLOCATABLE :: A_POLY (:,:)  ! Collecting the derivatives related to DX, DY, DDX, DDY and DXDY
 	REAL (KIND = 8), ALLOCATABLE :: DIS(:)   ! Distance to neighboring points
 	INTEGER, ALLOCATABLE :: XYZ_LOCAL(:)  ! Neighboring point index (in global numbering system)
 	REAL (KIND = 8) :: VX,VY  ! Velocity at a point
 	REAL (KIND = 8) :: R_SMOOTH ! Length of smoothing kernel (to normalize the distance so that the matrix has a better condition number)
 	REAL (KIND = 8) :: BZ   ! The sum of terms from the old time step
 	INTEGER :: INT_FLAG,INT_NUM  ! They are both legacy parameters, where INT_FLAG is the boundary type, and INT_NUM is the index of BC
 END TYPE 

END MODULE POINT 


MODULE DEFINE
 USE kd_tree
 USE POINT
 INTEGER :: NCYC,NCYC_OUT,NCYC_CHECK  ! NCYC is the total number of iteration
 INTEGER :: XYZ_TNUM,XYZ_NUM   !XYZ_TNUM is total number of points; XYZ_NUM is the total number of neighboring points 
 INTEGER :: TIME_NUM
 INTEGER :: DT_INPUT_TYPE   ! Control how to determine the time step (= 0 manually input; = 1 compute from the formula)
 INTEGER :: T_TYPE   ! = 0, Uniform flow, so no need to search each time step; = 1, need to research each time step using k-d tree
 TYPE (POINT_COE), ALLOCATABLE :: POINT_XY(:)
!------For BC      This boundary I/O parameters are only legacy parts from the Eulerian program. 
! For the current code, the points at the boundary are all allowed to flow freely. 
! There used to be a version to use the periodical boundary condition, but it may cause some gaps when re-injecting points into the domain 
 INTEGER :: HEAD_NUM
 REAL (KIND = 8), ALLOCATABLE, DIMENSION (:) :: H_INPUT
 INTEGER :: FLUX_NUM
 REAL (KIND = 8), ALLOCATABLE, DIMENSION (:) :: QBX,QBY
!-------For basic parameters
 REAL (KIND = 8) :: T,DT,T_START,DT_AMP,T_LIM
 REAL (KIND = 8) :: ALPHA_T,ALPHA_L,DM ! Parameters for the dispersion 
 INTEGER :: DRAW_NUM ! Total number of points to output 
 REAL (KIND = 8), ALLOCATABLE, DIMENSION (:) :: DRAW_TIME ! When to output points
!-------Name for I/O files
 CHARACTER (LEN = 256) :: NAME1,NAME2,NAME3,NAME4,NAME5,NAME6
!-------For advection part  (Just need l+1/2 and l time step, but import into XYZ to do the genearl calculation)
 INTEGER :: VEL_NUM    ! Total number of input velocity nodes
 REAL (KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: VEL_XY
 REAL (KIND = 8), ALLOCATABLE, DIMENSION (:,:) :: VEL_INPUT
 INTEGER :: VEL_LOCAL_NUM
 REAL (KIND = 8) :: X_BOUND_L,X_BOUND_R,Y_BOUND_N,Y_BOUND_S  ! Legacy parameters to impose periodical BC
 type (tree_master_record), Pointer :: tree_v  ! Tree for velocity interpolation
 type (tree_master_record), Pointer :: tree_p  ! Tree for GMLS matrix construction
 CHARACTER (LEN = 256) :: NAME_FILE

END MODULE 

PROGRAM MAIN
 USE DEFINE
 IMPLICIT NONE 
 REAL (KIND = 8) :: XK
 INTEGER :: I,J,K,IK,JK,KK,N
 INTEGER :: IKKK
 CHARACTER (LEN = 256) :: PATH
 REAL (KIND = 8) :: START_T,END_T
 
!-----For the input file 
 CALL cpu_time(START_T)
 NCYC = 0
 CALL INITIAL_BASIC
 CALL ADV_HALF
 CALL INTERPOLATE
 
 CALL HEAD_READ 
 CALL FLUX_READ
 
 call getcwd(PATH)
 PATH = trim(PATH)//'/'
 NAME1 = trim(PATH)//'head_ini.dat'
 OPEN (1,FILE = trim(NAME1))
 DO IK = 1,XYZ_TNUM
 	READ (1,*) POINT_XY(IK)%CTH   ! read initial concentration
 END DO
 CLOSE (1)
 
 DO IK = 1,XYZ_TNUM
 	IF (POINT_XY(IK)%INT_FLAG.EQ.1) THEN ! Constant Concentration condition
 		N = POINT_XY(IK)%INT_NUM
 		POINT_XY(IK)%CTH = H_INPUT(N)
 	END IF
 END DO
 
 DO IK = 1,XYZ_TNUM
 	POINT_XY(IK)%CTHN = POINT_XY(IK)%CTH
 END DO

 CALL K_CAL
 CALL GAU_COE 
 
 IKKK = 0
 TIME_NUM = 1
 
 DO WHILE (T.LT.T_LIM)
 	IF (DT_INPUT_TYPE.NE.0) THEN 
 		DO IK = 1,XYZ_TNUM
 			XK = POINT_XY(IK)%KTH(1,1)+POINT_XY(IK)%KTH(2,2)
 			IF (XK.NE.0.0) THEN
 				DT = MIN(DT,0.1*(POINT_XY(IK)%R_SMOOTH)**2/XK)  
 			END IF
 		END DO
 	END IF
 	T = T+DT
 	NCYC = NCYC+1
 	IF (MOD(NCYC,NCYC_CHECK).EQ.0) THEN 
 		WRITE (*,*) 'NCYC = ',NCYC,T,DT
 	END IF
 	CALL INTERPOLATE
 	CALL HEAD_READ
 	CALL FLUX_READ
 	CALL K_CAL
 	CALL GAU_COE
 	CALL DIFF
 	CALL ADV_HALF
 	
 	CALL INTERPOLATE
 	CALL HEAD_READ
 	CALL FLUX_READ
 	CALL K_CAL
 	CALL GAU_COE
 	CALL DIFF
 	CALL ADV_HALF
 	
 	IF (TIME_NUM.LE.DRAW_NUM) THEN 
 		IF (T.LE.DRAW_TIME(TIME_NUM).AND.(T+DT).GT.DRAW_TIME(TIME_NUM)) THEN 
 			NAME6 = trim(PATH)//'output_con_'//CHAR(TIME_NUM+48)//'_'//trim(NAME_FILE)//'.dat'
 			OPEN (993,FILE = trim(NAME6))
 			WRITE (993,*) T
 			DO IK = 1,XYZ_TNUM
 				WRITE (993,'(5F30.20)') POINT_XY(IK)%X,POINT_XY(IK)%Y,POINT_XY(IK)%CTH
 		
 			END DO	 
 			CLOSE (993) 
 			
 			NAME6 = trim(PATH)//'KTH'//'_'//CHAR(TIME_NUM+48)//'_'//trim(NAME_FILE)//'.DAT'
 			OPEN (993,FILE = trim(NAME6))
 			DO IK = 1,XYZ_TNUM
 				WRITE (993,'(20F20.10)') POINT_XY(IK)%KTH(1,1),POINT_XY(IK)%KTH(1,2), &
 				POINT_XY(IK)%KTH(2,1),POINT_XY(IK)%KTH(2,2)
 			END DO
 			CLOSE (993)	
 			
 			TIME_NUM = TIME_NUM+1
 		END IF
 	END IF
 END DO
 CLOSE (993)
 CALL cpu_time(END_T)
 OPEN (96,FILE = trim(NAME_FILE)//'_CPU_TIME.DAT')
 WRITE (96,*) END_T-START_T
 CLOSE (96)
 
END PROGRAM 


SUBROUTINE INITIAL_BASIC
 USE DEFINE 
 CHARACTER (LEN = 256) :: PATH
 INTEGER :: IK,JK,KK
 
 CALL getcwd(PATH)
 PATH = trim(PATH)//'/'
 NAME1 = trim(PATH)//'basic_input.dat'
 OPEN (1,FILE = trim(NAME1))
 READ (1,*) XYZ_TNUM,XYZ_NUM  ! Total Points, Total Local Points
 READ (1,*) DT_INPUT_TYPE  ! Time step type
 READ (1,*) T_TYPE ! Whether it is uniform flow
 READ (1,*) T_LIM,DT,DT_AMP,T_START  ! Time Limits, Time Step, Amplification, Starting time
 READ (1,*) NCYC_OUT,NCYC_CHECK
 READ (1,*) HEAD_NUM,FLUX_NUM   ! HEAD_NUM is the total number of the essentail BC
 READ (1,*) DRAW_NUM   ! The total number of the output times
 READ (1,*) ALPHA_T,ALPHA_L,DM  ! Coefficients for the dispersion along the tranvers and longitudinual directions. DM is the diffusion coefficient
 READ (1,*) NAME_FILE  ! Name of the case
 CLOSE (1)
 
 T = T_START
 ALLOCATE (POINT_XY(XYZ_TNUM))
 DO IK = 1,XYZ_TNUM
 	ALLOCATE (POINT_XY(IK)%A_POLY(XYZ_NUM,5))  ! Because it is quadratic polynomial
 	ALLOCATE (POINT_XY(IK)%DIS(XYZ_NUM))
 	ALLOCATE (POINT_XY(IK)%XYZ_LOCAL(XYZ_NUM))
 END DO
 
!-------For boundary conditions 
 ALLOCATE (H_INPUT(HEAD_NUM),QBX(FLUX_NUM),QBY(FLUX_NUM))
!-------For the output
 ALLOCATE (DRAW_TIME(DRAW_NUM))

! Read velocity files
 NAME1 = trim(PATH)//'velocity_ini.dat'
 OPEN (1,FILE = trim(NAME1))
 READ (1,*) VEL_NUM, VEL_LOCAL_NUM  ! VEL_NUM is the total number of input velocity nodes; VEL_LOCAL_NUM is the numbers of local nodes
 READ (1,*) X_BOUND_L,X_BOUND_R,Y_BOUND_S,Y_BOUND_N  ! Boundary limits
 CLOSE (1)
 ALLOCATE (VEL_XY(VEL_NUM,2),VEL_INPUT(VEL_NUM,2))
 
 NAME1 = trim(PATH)//'velocity.dat'
 OPEN (1,FILE = trim(NAME1))
 DO IK = 1,XYZ_TNUM
 	READ (1,*) VEL_XY(IK,1),VEL_XY(IK,2),VEL_INPUT(IK,1),VEL_INPUT(IK,2)
 END DO
 CLOSE (1)
 
 
! Read Output Coordinates
 NAME1 = trim(PATH)//'output_time.dat'
 OPEN (1,FILE = trim(NAME1))
 DO IK = 1,DRAW_NUM
 	READ (1,*)DRAW_TIME(IK)
 END DO
 CLOSE (1)
 
! Read XYZ Coordinate 
 NAME1 = trim(PATH)//'xyz.geo'
 OPEN (1,FILE = trim(NAME1))
 DO IK = 1,XYZ_TNUM
 	READ (1,*) POINT_XY(IK)%X,POINT_XY(IK)%Y,POINT_XY(IK)%INT_FLAG,POINT_XY(IK)%INT_NUM ! Read the coordinates of the point 
 END DO
 CLOSE (1)
 
END SUBROUTINE

SUBROUTINE HEAD_READ   ! Reading essential boundary condition (legacy subroutine)
 USE DEFINE 
 IMPLICIT NONE 
 INTEGER :: I,J,K,N
 INTEGER :: ITYPE
 INTEGER :: N1,N2,N3,N4
 REAL (KIND = 8) :: H1,H2,H_TIME1,H_TIME2
 REAL (KIND = 8) :: H_SET
 CHARACTER (LEN = 256) :: PATH
 
 call getcwd(PATH)
 PATH = trim(PATH)//'/'
 DO N = 1,HEAD_NUM
 	IF (N.LT.10) THEN 
 		NAME3 = trim(PATH)//'head'//CHAR(N+48)//'.bound'
 	ELSE 
 		IF (N.LT.100) THEN 
 			N1 = FLOOR(N/10.0)
 			N2 = N-N1*10
 			NAME3 = trim(PATH)//'head'//CHAR(N1+48)//CHAR(N2+48)//'.bound'	
 		ELSE 
 			IF (N.LT.1000) THEN 
 				N1 = FLOOR(N/100.0)
 				N2 = FLOOR((N-N1*100.0)/10.0)
 				N3 = N-N1*100-N2*10
 				NAME3 = trim(PATH)//'head'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//'.bound'
 			ELSE 
 				N1 = FLOOR(N/1000.0)
 				N2 = FLOOR((N-N1*1000.0)/100.0)
 				N3 = FLOOR((N-N1*1000.0-N2*100.0)/10.0)
 				N4 = N-N1*1000-N2*100-N3*10
 				NAME3 = trim(PATH)//'head'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//CHAR(N4+48)//'.bound'
 			END IF
 		END IF
 	END IF
 	OPEN (3,FILE = trim(NAME3))
 	READ (3,*) ITYPE
 	IF (ITYPE.EQ.0) THEN ! Constant concentration
 		READ (3,*) H_SET
 		H_INPUT(N) = H_SET
 	ELSE  ! Time series data
 		READ (3,*) H_TIME1,H1
 		IF (NCYC.EQ.0) THEN 
 			H_SET = H1
 		END IF
 		DO WHILE (NCYC.NE.0.AND.T.GE.H_TIME1)
 			READ (3,*) H_TIME2,H2
 			IF (H_TIME2.GE.T) THEN 
 				H_SET = (H2*(T-H_TIME1)+H1*(H_TIME2-T))/(H_TIME2-H_TIME1)
 			END IF
 			H_TIME1 = H_TIME2
 			H1 = H2
 		END DO
 		H_INPUT(N) = H_SET
 	END IF
 	CLOSE (3)
 END DO

END SUBROUTINE 

SUBROUTINE FLUX_READ ! Reading flux boundary condition (legacy subroutine)
 USE DEFINE 
 IMPLICIT NONE 
 INTEGER :: I,J,K,N
 INTEGER :: ITYPE 
 INTEGER :: N1,N2,N3,N4
 REAL (KIND = 8) :: F_TIME1,F_TIME2
 REAL (KIND = 8) :: FX1,FX2,FY1,FY2
 REAL (KIND = 8) :: FX_SET,FY_SET
 CHARACTER (LEN = 256) :: PATH
 
 call getcwd(PATH)
 PATH = trim(PATH)//'/'
 DO N = 1,FLUX_NUM
 	IF (N.LT.10) THEN 
 		NAME3 = trim(PATH)//'flux'//CHAR(N+48)//'.bound'
 	ELSE 
 		IF (N.LT.100) THEN 
 			N1 = FLOOR(N/10.0)
 			N2 = N-N1*10
 			NAME3 = trim(PATH)//'flux'//CHAR(N1+48)//CHAR(N2+48)//'.bound'
 		ELSE 
 			IF (N.LT.1000) THEN 
 				N1 = FLOOR(N/100.0)
 				N2 = FLOOR((N-N1*100.0)/10.0)
 				N3 = N-N1*100-N2*10
 				NAME3 = trim(PATH)//'flux'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//'.bound'
 			ELSE 
 				N1 = FLOOR(N/1000.0)
 				N2 = FLOOR((N-N1*1000.0)/100.0)
 				N3 = FLOOR((N-N1*1000.0-N2*100.0)/10.0)
 				N4 = N-N1*1000-N2*100-N3*10
 				NAME3 = trim(PATH)//'flux'//CHAR(N1+48)//CHAR(N2+48)//CHAR(N3+48)//CHAR(N4+48)//'.bound'
 			END IF
 		END IF
 	END IF
 	OPEN (3,FILE = trim(NAME3))
 	READ (3,*) ITYPE
 	IF (ITYPE.EQ.0) THEN ! Constant flux
 		READ (3,*) FX_SET,FY_SET
 		QBX(N) = FX_SET; QBY(N) = FY_SET
 	ELSE ! Time-series data
 		READ (3,*) F_TIME1,FX1,FY1
 		IF (NCYC.EQ.0) THEN 
 			QBX(N) = FX1; QBY(N) = FY1
 		END IF
 		DO WHILE (NCYC.NE.0.AND.T.GE.F_TIME1)
 			READ (3,*) F_TIME2,FX2,FY2
 			IF (F_TIME2.GE.T) THEN 
 				FX_SET = (FX2*(T-F_TIME1)+FX1*(F_TIME2-T))/(F_TIME2 - F_TIME1)
 				FY_SET = (FY2*(T-F_TIME1)+FY1*(F_TIME2-T))/(F_TIME2 - F_TIME1) 
 			END IF
 			F_TIME1 = F_TIME2
 			FX1 = FX2; FY1 = FY2
 		END DO
 		QBX(N) = FX_SET; QBY(N) = FY_SET
 	END IF
 	CLOSE (3)
 END DO
 
END SUBROUTINE 
