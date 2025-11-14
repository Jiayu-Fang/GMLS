SUBROUTINE DIFF    ! Compute the dispersion term
 USE DEFINE 
 IMPLICIT NONE 
 INTEGER :: I,J,K
 INTEGER :: IK,JK,KK
 INTEGER :: NK,NI,NJ
 REAL (KIND = 8) :: KTH1,KS1
 REAL (KIND = 8) :: KTH(2,2)
 REAL (KIND = 8) :: C_COM
 
 DO IK = 1,XYZ_TNUM
 	POINT_XY(IK)%BZ = 0.0   ! It was written for L2 load term (though not needed for the cases in the paper)
 	POINT_XY(IK)%CTHN = POINT_XY(IK)%CTH
 END DO
 
 DO IK = 1,XYZ_TNUM
 	IF (POINT_XY(IK)%INT_FLAG.EQ.1) THEN ! Essential BC (subject to change) Useless here though (I wrote this part in case it is needed in the future)
 		I = POINT_XY(IK)%INT_NUM
 		POINT_XY(IK)%CTH = H_INPUT(I)
 	ELSE 
 		C_COM = -POINT_XY(IK)%BZ+POINT_XY(IK)%XY_PXX1
 		POINT_XY(IK)%CTH = C_COM*(DT/2.0)+POINT_XY(IK)%CTHN
 	END IF
 END DO 


END SUBROUTINE 

SUBROUTINE K_CAL ! Compute the dispersion coefficient
 USE DEFINE 
 IMPLICIT NONE 
 INTEGER :: I,J,K
 INTEGER :: IK,JK,KK
 INTEGER :: NK,NI,NJ
 REAL (KIND = 8) :: VXL,VYL,VL
 
 DO IK = 1,XYZ_TNUM
 	VXL = POINT_XY(IK)%VX; VYL = POINT_XY(IK)%VY
 	VL = SQRT(VXL**2+VYL**2)
 	POINT_XY(IK)%KTH(1,1) = ALPHA_T*VL+DM+(ALPHA_L-ALPHA_T)*VXL*VXL/VL
 	POINT_XY(IK)%KTH(1,2) = (ALPHA_L-ALPHA_T)*VXL*VYL/VL
 	POINT_XY(IK)%KTH(2,1) = POINT_XY(IK)%KTH(1,2)
 	POINT_XY(IK)%KTH(2,2) = ALPHA_T*VL+DM+(ALPHA_L-ALPHA_T)*VYL*VYL/VL
 END DO
 
END SUBROUTINE 

SUBROUTINE GAU_COE
 USE DEFINE 
 IMPLICIT NONE 
 REAL (KIND = 8) :: XYZ_L(XYZ_NUM,2)  ! Coordinates of local point
 REAL (KIND = 8) :: WL(XYZ_NUM),RHO1(XYZ_NUM)  ! Weight function and distance 
 INTEGER :: ORD_NUM,RF_NUM   ! Order of polynomial; Total number of polynomial terms 
 INTEGER :: IK,JK,KK,I,J,K
 INTEGER :: NK,NI,NJ
 INTEGER :: FLAG_NUM
 REAL (KIND = 8) :: KTH(2,2)
 REAL (KIND = 8) :: XK,XK1,XK2,RH_INT,RHO,X0,Y0,X1,Y1,X2,Y2,WLR,WLD
 REAL (KIND = 8) :: X_INI,Y_INI,RR
 REAL (KIND = 8) :: C_COE
 REAL (KIND = 8) :: KXXX(XYZ_TNUM),KXYX(XYZ_TNUM),KYXY(XYZ_TNUM),KYYY(XYZ_TNUM)
 REAL (KIND = 8) :: XY_PXX2(XYZ_TNUM,5)
 
 IF (NCYC.EQ.0.OR.T_TYPE.NE.0) THEN    ! The GMLS coefficient matrix does not need to be recomputed if the flow is uniform
 	DO IK = 1,XYZ_TNUM
 		POINT_XY(IK)%A_POLY = 0.0
 		RHO = 0.0; RHO1 = 0.0; WL = 0.0
 		X_INI = POINT_XY(IK)%X; Y_INI = POINT_XY(IK)%Y
 		DO JK = 1,XYZ_NUM
 			KK = POINT_XY(IK)%XYZ_LOCAL(JK)
 			XYZ_L(JK,1) = POINT_XY(KK)%X 
 			XYZ_L(JK,2) = POINT_XY(KK)%Y
 			RHO1(JK) = POINT_XY(IK)%DIS(JK)
 			RHO = MAX(RHO,RHO1(JK))
 		END DO
 		RHO = 1.05*RHO
 		POINT_XY(IK)%R_SMOOTH = RHO
 		DO JK = 1,XYZ_NUM
 			RR = RHO1(JK)
 			CALL WL_CAL1(RR,RHO,WLR,WLD)   ! WLD is the derivative of the kernel function, but it is not needed in GMLS 
 			                               ! It is included here only because the legacy from the MLS
 			WL(JK) = WLR
 		END DO
 		ORD_NUM = 2
 		RF_NUM = (ORD_NUM+1)*(ORD_NUM+2)/2
 		CALL MATRIX_INV_POLY(XYZ_NUM,ORD_NUM,RF_NUM,XYZ_L,WL,RHO,POINT_XY(IK)%A_POLY,X_INI,Y_INI)
 	END DO
 END IF

 XY_PXX2 = 0.0
 KXXX = 0.0; KXYX = 0.0; KYXY = 0.0; KYYY = 0.0
 DO IK = 1,XYZ_TNUM
 	DO JK = 1,XYZ_NUM
 		KK = POINT_XY(IK)%XYZ_LOCAL(JK)
 		XY_PXX2(IK,1) = XY_PXX2(IK,1)+POINT_XY(IK)%A_POLY(JK,1)*POINT_XY(KK)%CTH   ! DX
 		XY_PXX2(IK,2) = XY_PXX2(IK,2)+POINT_XY(IK)%A_POLY(JK,2)*POINT_XY(KK)%CTH   ! DY
 		XY_PXX2(IK,3) = XY_PXX2(IK,3)+POINT_XY(IK)%A_POLY(JK,3)*POINT_XY(KK)%CTH   ! DDX
 		XY_PXX2(IK,4) = XY_PXX2(IK,4)+POINT_XY(IK)%A_POLY(JK,4)*POINT_XY(KK)%CTH   ! DDY
 		XY_PXX2(IK,5) = XY_PXX2(IK,5)+POINT_XY(IK)%A_POLY(JK,5)*POINT_XY(KK)%CTH   ! DXDY
 		
 		KXXX(IK) = KXXX(IK)+POINT_XY(IK)%A_POLY(JK,1)*POINT_XY(KK)%KTH(1,1)    ! DKXXDX
 		KXYX(IK) = KXYX(IK)+POINT_XY(IK)%A_POLY(JK,1)*POINT_XY(KK)%KTH(1,2)    ! DKXYDX
 		KYXY(IK) = KYXY(IK)+POINT_XY(IK)%A_POLY(JK,2)*POINT_XY(KK)%KTH(2,1)    ! DKYXDY
 		KYYY(IK) = KYYY(IK)+POINT_XY(IK)%A_POLY(JK,2)*POINT_XY(KK)%KTH(2,2)    ! DKYYDY
 	END DO
 END DO
 
 DO IK = 1,XYZ_TNUM
 	KTH(1,1) = POINT_XY(IK)%KTH(1,1)
 	KTH(1,2) = POINT_XY(IK)%KTH(1,2)
 	KTH(2,1) = POINT_XY(IK)%KTH(2,1)
 	KTH(2,2) = POINT_XY(IK)%KTH(2,2)
 	POINT_XY(IK)%XY_PXX1 = KTH(1,1)*XY_PXX2(IK,3)+ &
 			         +KTH(2,2)*XY_PXX2(IK,4)+ &
 			         +(KTH(1,2)+KTH(2,1))*XY_PXX2(IK,5)+ &
 	              KXXX(IK)*XY_PXX2(IK,1)+KXYX(IK)*XY_PXX2(IK,2)+ &
 	              KYXY(IK)*XY_PXX2(IK,1)+KYYY(IK)*XY_PXX2(IK,2)     
 END DO

END SUBROUTINE  

SUBROUTINE WL_CAL1(RR,RM,WLR,WLD)
 IMPLICIT NONE
 REAL (KIND = 8) :: RM,WLR,WLD,RR
 REAL (KIND = 8) :: KV,CV,Q
 
 WLR = 0.0; WLD = 0.0
 KV = 9.0/3.1415926/RM/RM
 Q = RR/RM
 IF (Q.LE.1.0) THEN 
 	WLR = KV*((1.0-Q)**6)*(1.0+6.0*Q+35.0/3.0*Q*Q)
 	WLD = -6.0*KV*((1.0-Q)**5)*(1.0+6.0*Q+35.0/3.0*Q*Q)+  &
 	       KV*((1.0-Q)**6)*(6.0+70.0/3.0*Q)
 END IF
 
END SUBROUTINE

SUBROUTINE WL_CAL(RM,WLR,WLD)
 IMPLICIT NONE 
 REAL (KIND = 8) :: RM,WLR,WLD
 
 WLR = 0.0; WLD = 0.0
 
 IF (RM.GE.0.0.AND.RM.LE.1.0/2.0) THEN 
 	WLR = 2.0/3.0-4.0*RM*RM+4.0*RM*RM*RM
 	WLD = -8.0*RM+12.0*RM*RM
 END IF
 
 IF (RM.GT.1.0/2.0.AND.RM.LE.1.0) THEN 
 	WLR = 4.0/3.0-4.0*RM+4.0*RM*RM-4.0/3.0*RM*RM*RM
 	WLD = -4.0+8.0*RM-4.0*RM*RM
 END IF
 
END SUBROUTINE


SUBROUTINE LU_INVERSE(RF_NUM,A,A_INV)   ! Using LU decomposition to compute the inverse of matrix (I have a SVD version, but it is very slow)
 IMPLICIT NONE 
 INTEGER :: RF_NUM
 REAL (KIND = 8) :: A(RF_NUM,RF_NUM),A_INV(RF_NUM,RF_NUM),A_TEM(RF_NUM)
 REAL (KIND = 8) :: L(RF_NUM,RF_NUM),U(RF_NUM,RF_NUM),LI(RF_NUM,RF_NUM),UI(RF_NUM,RF_NUM)
 REAL (KIND = 8) :: XK
 INTEGER :: IK,JK,KK,NK_B
 INTEGER :: IND_B(RF_NUM,2)
 
 NK_B = 0; IND_B = 0
 DO IK = 1,RF_NUM
 	IF (ABS(A(IK,IK)).LT.1.0E-15) THEN 
 		JK = IK+1
 		DO WHILE (ABS(A(JK,IK)).LT.1.0E-15.AND.JK.LT.RF_NUM)
 			JK = JK+1
 		END DO			
 		DO KK = 1,RF_NUM
 			A_TEM(KK) = A(JK,KK)
 		END DO
 		DO KK = 1,RF_NUM
 			A(JK,KK) = A(IK,KK)
 		END DO
 		DO KK = 1,RF_NUM
 			A(IK,KK) = A_TEM(KK)
 		END DO
 		NK_B = NK_B+1
 		IND_B(NK_B,1) = JK; IND_B(NK_B,2) = IK
 	END IF
 END DO
 
 
 L = 0.0; U = 0.0
 DO IK = 1,RF_NUM
 	JK = IK
 	L(IK,JK) = 1.0
 END DO
 
! LU Decomposition
 DO JK = 1,RF_NUM
 	U(1,JK) = A(1,JK)
 END DO
 DO IK = 2,RF_NUM
 	L(IK,1) = A(IK,1)/U(1,1)
 	DO JK = 2,IK-1
 		XK = A(IK,JK)
 		DO KK = 1,JK-1
 			XK = XK-L(IK,KK)*U(KK,JK)
 		END DO
 		IF (ABS(U(JK,JK)).LT.1.0E-18) THEN 
 			WRITE (*,*) 'WARNING'
 			L(IK,JK) = 0.0
 		ELSE 
 			L(IK,JK) = XK/U(JK,JK)
 		END IF
 	END DO
 	DO JK = IK,RF_NUM
 		XK = A(IK,JK)
 		DO KK = 1,JK-1
 			XK = XK-L(IK,KK)*U(KK,JK)
 		END DO
 		U(IK,JK) = XK/1.0
 	END DO
 END DO
 
! Computing the inverse of L AND U MATRIX
! For U
 LI = 0.0; UI = 0.0
 DO IK = 1,RF_NUM
 	DO JK = 1,RF_NUM
 		IF (IK.EQ.JK) THEN 
 			XK = 1.0
 		ELSE 
 			XK = 0.0
 		END IF
 		DO KK = 1,JK-1
 			XK = XK-UI(IK,KK)*U(KK,JK)
 		END DO
 		IF (ABS(U(JK,JK)).LT.1.0E-18) THEN 
 			UI(IK,JK) = 0.0
 		ELSE  
 			UI(IK,JK) = XK/U(JK,JK)
 		END IF
 	END DO
 END DO

! For L
 DO IK = 1,RF_NUM
 	DO JK = RF_NUM,1,-1
 		IF (IK.EQ.JK) THEN 
 			XK = 1.0
 		ELSE 
 			XK = 0.0
 		END IF
 		DO KK = RF_NUM,JK+1,-1
 			XK = XK-LI(IK,KK)*L(KK,JK)
 		END DO
 		LI(IK,JK) = XK
 	END DO
 END DO

! For the final inverse of A
 A_INV = 0.0
 DO IK = 1,RF_NUM
 	DO JK = 1,RF_NUM
 		DO KK = 1,RF_NUM
 			A_INV(IK,JK) = A_INV(IK,JK)+UI(IK,KK)*LI(KK,JK)
 		END DO
 	END DO
 END DO

 DO IK = RF_NUM,1,-1
 	IF (IND_B(IK,1).NE.0) THEN 
 		DO KK = 1,RF_NUM
 			A_TEM(KK) = A_INV(KK,IND_B(IK,1))
 		END DO
 		DO KK = 1,RF_NUM
 			A_INV(KK,IND_B(IK,1)) = A_INV(KK,IND_B(IK,2)) 
 		END DO
 		DO KK = 1,RF_NUM
 			A_INV(KK,IND_B(IK,2)) = A_TEM(KK)
 		END DO
 	END IF
 END DO
 
 DO IK = RF_NUM,1,-1
 	IF (IND_B(IK,1).NE.0) THEN 
 		DO KK = 1,RF_NUM
 			A_TEM(KK) = A(KK,IND_B(IK,1))
 		END DO
 		DO KK = 1,RF_NUM
 			A(KK,IND_B(IK,1)) = A(KK,IND_B(IK,2)) 
 		END DO
 		DO KK = 1,RF_NUM
 			A(KK,IND_B(IK,2)) = A_TEM(KK)
 		END DO
 	END IF
 END DO
 
END SUBROUTINE 

 
SUBROUTINE MATRIX_INV_POLY(XYZ_NUM,ORD_NUM,RF_NUM,XYZ_L,WL,RM,A_POLY,X_INI,Y_INI)
! ORD_NUM is the order of polygonal (such as 2 for quadratic polynomial); RF_NUM is the total term of the polynomial
 IMPLICIT NONE 
 INTEGER :: RF_NUM,XYZ_NUM,ORD_NUM
 INTEGER :: IK,JK,KK,N,NK,I,J,K
 REAL (KIND = 8) :: XYZ_L(XYZ_NUM,2),AK(XYZ_NUM,RF_NUM)
 REAL (KIND = 8) :: A(RF_NUM,RF_NUM),L(RF_NUM,RF_NUM),U(RF_NUM,RF_NUM)
 REAL (KIND = 8) :: LI(RF_NUM,RF_NUM),UI(RF_NUM,RF_NUM),A_INV(RF_NUM,RF_NUM)
 REAL (KIND = 8) :: WL(XYZ_NUM)
 REAL (KIND = 8) :: A_POLY(XYZ_NUM,5)
 REAL (KIND = 8) :: A_POLY1(RF_NUM,XYZ_NUM)
 REAL (KIND = 8) :: DX,DY,XK,YK,RM
 REAL (KIND = 8) :: X_INI,Y_INI
 
 AK = 0.0
 DO IK = 1,XYZ_NUM
 	DX = (XYZ_L(IK,1)-X_INI)/RM
 	DY = (XYZ_L(IK,2)-Y_INI)/RM
 	K = 1
 	DO JK = 0,ORD_NUM
 		DO I = 0,JK   ! I is x order and J is y order
 			J = JK-I
 			AK(IK,K) = (DX**I)*(DY**J)
 			K = K+1
 		END DO
 	END DO
 END DO
 
 A = 0.0; 
 DO IK = 1,RF_NUM
 	DO JK = 1,RF_NUM
 		DO KK = 1,XYZ_NUM
 			A(IK,JK) = A(IK,JK)+AK(KK,IK)*AK(KK,JK)*WL(KK)
 		END DO
 	END DO
 END DO

 L = 0.0; U = 0.0
 DO IK = 1,RF_NUM
 	JK = IK
 	L(IK,JK) = 1.0
 END DO
 
! LU Decomposition
 DO JK = 1,RF_NUM
 	U(1,JK) = A(1,JK)
 END DO
 DO IK = 2,RF_NUM
 	L(IK,1) = A(IK,1)/U(1,1)
 	DO JK = 2,IK-1
 		XK = A(IK,JK)
 		DO KK = 1,JK-1
 			XK = XK-L(IK,KK)*U(KK,JK)
 		END DO
 		IF (ABS(U(JK,JK)).LT.1.0E-10) THEN 
 			L(IK,JK) = 0.0
 			WRITE (*,*) 'WARNING'   ! Meaning it is nearly singular
 		ELSE 
 			L(IK,JK) = XK/U(JK,JK)
 		END IF
 	END DO
 	DO JK = IK,RF_NUM
 		XK = A(IK,JK)
 		DO KK = 1,JK-1
 			XK = XK-L(IK,KK)*U(KK,JK)
 		END DO
 		U(IK,JK) = XK/1.0
 	END DO
 END DO
 
! Computing the inverse of L AND U MATRIX
! For U
 LI = 0.0; UI = 0.0
 DO IK = 1,RF_NUM
 	DO JK = 1,RF_NUM
 		IF (IK.EQ.JK) THEN 
 			XK = 1.0
 		ELSE 
 			XK = 0.0
 		END IF
 		DO KK = 1,JK-1
 			XK = XK-UI(IK,KK)*U(KK,JK)
 		END DO
 		IF (ABS(U(JK,JK)).LT.1.0E-10) THEN 
 			UI(IK,JK) = 0.0
 			WRITE (*,*) 'WARNING'
 		ELSE  
 			UI(IK,JK) = XK/U(JK,JK)
 		END IF
 	END DO
 END DO

! For L
 DO IK = 1,RF_NUM
 	DO JK = RF_NUM,1,-1
 		IF (IK.EQ.JK) THEN 
 			XK = 1.0
 		ELSE 
 			XK = 0.0
 		END IF
 		DO KK = RF_NUM,JK+1,-1
 			XK = XK-LI(IK,KK)*L(KK,JK)
 		END DO
 		LI(IK,JK) = XK
 	END DO
 END DO

! For the final inverse of A
 A_INV = 0.0
 DO IK = 1,RF_NUM
 	DO JK = 1,RF_NUM
 		DO KK = 1,RF_NUM
 			A_INV(IK,JK) = A_INV(IK,JK)+UI(IK,KK)*LI(KK,JK)
 		END DO
 	END DO
 END DO
 
 A_POLY = 0.0; A_POLY1 = 0.0
 DO IK = 1,RF_NUM
 	DO JK = 1,XYZ_NUM
 		DO KK = 1,RF_NUM
 			A_POLY1(IK,JK) = A_POLY1(IK,JK)+A_INV(IK,KK)*WL(JK)*AK(JK,KK)
 		END DO
 	END DO
 END DO
 
 DO JK = 1,XYZ_NUM  ! It is subject to change since it assumes the use of qudratic polynomial 
 	A_POLY(JK,1) = 1.0/RM*A_POLY1(3,JK)  ! dx
 	A_POLY(JK,2) = 1.0/RM*A_POLY1(2,JK)  ! dy
  	A_POLY(JK,3) = 2.0/RM/RM*A_POLY1(6,JK) ! dxx 
 	A_POLY(JK,4) = 2.0/RM/RM*A_POLY1(4,JK) ! dyy
 	A_POLY(JK,5) = 1.0/RM/RM*A_POLY1(5,JK) ! dxdy
 END DO
 
 
END SUBROUTINE 

 



















