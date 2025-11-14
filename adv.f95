SUBROUTINE ADV_HALF   ! Calculating advection process for the half time step
 USE DEFINE
 IMPLICIT NONE 
 INTEGER :: IK,JK,KK,NK
 INTEGER :: I,J,K
 INTEGER :: R_POINTS(VEL_LOCAL_NUM) 
 REAL :: R_DIST(VEL_LOCAL_NUM)
 REAL :: XYZ_LL(VEL_NUM,2)
 REAL :: query_vec(2)
 CHARACTER (LEN = 256) :: PATH
 REAL (KIND = 8) :: COE_SUM,VX_COE,VY_COE
 
 IF (NCYC.NE.0) THEN  
 	DO IK = 1,XYZ_TNUM
 		POINT_XY(IK)%X = POINT_XY(IK)%X+POINT_XY(IK)%VX*DT*0.5
 		POINT_XY(IK)%Y = POINT_XY(IK)%Y+POINT_XY(IK)%VY*DT*0.5
 	END DO
 	IF (T_TYPE.NE.0) THEN  ! With non-uniform flow, k-d tree needs to be rebuilt
 		CALL DESTROY_TREE(tree_v)
 	END IF
 END IF
 
 IF (NCYC.EQ.0.OR.T_TYPE.NE.0) THEN 
 	XYZ_LL = VEL_XY
 	tree_v => create_tree(XYZ_LL)
 	DO IK = 1,XYZ_TNUM
 		query_vec(1)=POINT_XY(IK)%X; query_vec(2)=POINT_XY(IK)%Y
 		CALL n_nearest_to(tp=tree_v, qv=query_vec, n=VEL_LOCAL_NUM, indexes=R_POINTS, distances=R_DIST)
 		COE_SUM = 0.0; VX_COE = 0.0; VY_COE = 0.0
 		DO JK = 1,VEL_LOCAL_NUM
 			COE_SUM = COE_SUM+1.0/MAX(R_DIST(JK),1.0E-2)
 			KK = R_POINTS(JK)
 			VX_COE = VX_COE+VEL_INPUT(KK,1)/MAX(R_DIST(JK),1.0E-2)
 			VY_COE = VY_COE+VEL_INPUT(KK,2)/MAX(R_DIST(JK),1.0E-2)
 		END DO
 		POINT_XY(IK)%VX = VX_COE/COE_SUM/1.0
 		POINT_XY(IK)%VY = VY_COE/COE_SUM/1.0
 	END DO
 END IF
 
END SUBROUTINE 


SUBROUTINE INTERPOLATE   ! Obtaining the K neighboring nodes with K-D tree  
 USE DEFINE
 IMPLICIT NONE
! type (tree_master_record), Pointer :: tree
 INTEGER :: IK,JK,KK,NK
 INTEGER :: I,J,K
 INTEGER :: R_POINTS(XYZ_NUM) 
 REAL :: R_DIST(XYZ_NUM),XK
 REAL :: XYZ_LL(XYZ_TNUM,2)
 REAL :: query_vec(2)
 CHARACTER (LEN = 256) :: PATH
 
 
 IF (NCYC.EQ.0.OR.T_TYPE.NE.0) THEN 
 	DO IK = 1,XYZ_TNUM
 		XYZ_LL(IK,1) = POINT_XY(IK)%X
 		XYZ_LL(IK,2) = POINT_XY(IK)%Y
 	END DO
 	tree_p => create_tree(XYZ_LL)
 	DO IK = 1,XYZ_TNUM
 		query_vec(1)=XYZ_LL(IK,1); query_vec(2)=XYZ_LL(IK,2)
 		CALL n_nearest_to(tp=tree_p, qv=query_vec, n=XYZ_NUM, indexes=R_POINTS, distances=R_DIST)
 		JK = 1
 		DO WHILE (R_POINTS(JK).NE.IK)
 			JK = JK+1
 		END DO
 	
 		KK = R_POINTS(JK)
	 	R_POINTS(JK) = R_POINTS(1)
 		R_POINTS(1) = KK    ! Make sure the first neighboring point is the point itself
 		
 		XK = R_DIST(JK)
 		R_DIST(JK) = R_DIST(1)
 		R_DIST(1) = XK
 	
 		DO JK = 1,XYZ_NUM
 			POINT_XY(IK)%XYZ_LOCAL(JK) = R_POINTS(JK)
 			POINT_XY(IK)%DIS(JK) = SQRT(R_DIST(JK))
 		END DO
	 END DO
	 IF (T_TYPE.NE.0) THEN 
	 	CALL DESTROY_TREE(tree_p)
	 END IF
 END IF
 
 ! Notice since the explicit scheme is adopted here, it is no need to rearrange the lcoal indices in an increasing order to reduce the bandwidth

  
END SUBROUTINE  


