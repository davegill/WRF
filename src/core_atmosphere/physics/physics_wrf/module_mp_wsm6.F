MODULE module_mp_wsm6
!
   USE mp_wsm6 
!
CONTAINS
!===================================================================
!
  SUBROUTINE wsm6(th, q, qc, qr, qi, qs, qg                        &
                 ,den, pii, p, delz                                &
                 ,delt,g, cpd, cpv, rd, rv, t0c                    &
                 ,ep1, ep2, qmin                                   &
                 ,XLS, XLV0, XLF0, den0, denr                      &
                 ,cliq,cice,psat                                   &
                 ,rain, rainncv                                    &
                 ,snow, snowncv                                    &
                 ,sr                                               &
                 ,refl_10cm, diagflag, do_radar_ref                &
                 ,graupel, graupelncv                              &
                 ,has_reqc, has_reqi, has_reqs                     &  ! for radiation
                 ,re_cloud, re_ice,   re_snow                      &  ! for radiation   
                 ,errmsg, errflg                                   &
                 ,ids,ide, jds,jde, kds,kde                        &
                 ,ims,ime, jms,jme, kms,kme                        &
                 ,its,ite, jts,jte, kts,kte                        &
#ifdef WRF_CHEM
                 ,evapprod, rainprod                               &
#endif
                                                                   )
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(INOUT) ::                                          &
                                                             th,  &
                                                              q,  &
                                                              qc, &
                                                              qi, &
                                                              qr, &
                                                              qs, &
                                                              qg
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(IN   ) ::                                          &
                                                             den, &
                                                             pii, &
                                                               p, &
                                                            delz
  REAL, INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                              rd, &
                                                              rv, &
                                                             t0c, &
                                                            den0, &
                                                             cpd, &
                                                             cpv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL, DIMENSION( ims:ime , jms:jme ),                           &
        INTENT(INOUT) ::                                    rain
  REAL, DIMENSION( ims:ime , jms:jme ),                           &
        INTENT(OUT) ::                                      rainncv, sr
! for radiation connecting
  INTEGER, INTENT(IN)::                                           &
                                                        has_reqc, &
                                                        has_reqi, &
                                                        has_reqs
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),                     &
        INTENT(INOUT)::                                           &
                                                        re_cloud, &
                                                          re_ice, &
                                                         re_snow
!+---+-----------------------------------------------------------------+
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme), OPTIONAL,           &   ! GT
        INTENT(INOUT) ::                               refl_10cm
!+---+-----------------------------------------------------------------+

  REAL, DIMENSION( ims:ime , jms:jme ), OPTIONAL,                 &
        INTENT(INOUT) ::                               snow, graupel
  REAL, DIMENSION( ims:ime , jms:jme ), OPTIONAL,                 &
        INTENT(OUT) ::                                 snowncv, graupelncv

  character(len=*), intent(out)    ::                  errmsg
  integer,          intent(out)    ::                  errflg

#ifdef WRF_CHEM
  REAL, DIMENSION( ims:ime , kms:kme, jms:jme ), INTENT(OUT) :: &
                                                      rainprod,   &
                                                      evapprod
! local variable
  REAL, DIMENSION( its:ite , kts:kte )                 ::         &
                                                      rainprod2d, &
                                                      evapprod2d
#endif

! LOCAL VAR
  REAL, DIMENSION( its:ite , kts:kte ) :: delz_hv, p_hv, q_hv, den_hv
  REAL, DIMENSION( its:ite ) ::  rain_hv, rainncv_hv, sr_hv
  REAL, DIMENSION( its:ite ) ::  snow_hv, snowncv_hv, graupel_hv, graupelncv_hv
  REAL, DIMENSION( its:ite , kts:kte ) ::   t_hv, th_hv, pii_hv
  REAL, DIMENSION( its:ite , kts:kte ) ::   qc_hv, qi_hv
  REAL, DIMENSION( its:ite , kts:kte ) ::   qr_hv, qs_hv, qg_hv
  INTEGER ::               i,j,k

!+---+-----------------------------------------------------------------+
      REAL, DIMENSION(kts:kte):: qv1d, t1d, p1d, qr1d, qs1d, qg1d, dBZ
      REAL, DIMENSION(kts:kte):: th_hv2, pii_hv2    
      LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
      INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref
!+---+-----------------------------------------------------------------+
! to calculate effective radius for radiation
  REAL, DIMENSION( kts:kte ) :: den1d
  REAL, DIMENSION( kts:kte ) :: qc1d
  REAL, DIMENSION( kts:kte ) :: qi1d
  REAL, DIMENSION( kts:kte ) :: re_qc, re_qi, re_qs

      DO j=jts,jte
         DO k=kts,kte
         DO i=its,ite
            den_hv(i,k)  = den(i,k,j)
            delz_hv(i,k) = delz(i,k,j)
            p_hv(i,k)    = p(i,k,j)
            q_hv(i,k)    = q(i,k,j)
            pii_hv(i,k)  = pii(i,k,j)
            th_hv(i,k)   = th(i,k,j)
            qc_hv(i,k)   = qc(i,k,j)
            qi_hv(i,k)   = qi(i,k,j)
            qr_hv(i,k)   = qr(i,k,j)
            qs_hv(i,k)   = qs(i,k,j)
            qg_hv(i,k)   = qg(i,k,j)
         ENDDO
         ENDDO
         DO i=its,ite
           rain_hv(i) = rain(i,j)
         ENDDO
         IF(PRESENT (snowncv) .AND. PRESENT (snow)) THEN
            DO i=its,ite
               snow_hv(i) = snow(i,j)
            ENDDO
         ENDIF
         IF(PRESENT (graupelncv) .AND. PRESENT (graupel)) THEN
            DO i=its,ite
               graupel_hv(i) = graupel(i,j)
            ENDDO
         ENDIF
         !
         CALL wsm62D(t_hv,th_hv,pii_hv                             &
                    ,q_hv,qc_hv,qi_hv,qr_hv,qs_hv,qg_hv            &
                    ,den_hv                                        &
                    ,p_hv, delz_hv                                 &
                    ,delt,g, cpd, cpv, rd, rv, t0c                 &
                    ,ep1, ep2, qmin                                &
                    ,XLS, XLV0, XLF0, den0, denr                   &
                    ,cliq,cice,psat                                &
                    ,rain_hv,rainncv_hv                            &
                    ,sr_hv                                         &
                    ,its,ite, kts,kte                              &
                    ,snow_hv,snowncv_hv                            &
                    ,graupel_hv,graupelncv_hv                      &
#ifdef WRF_CHEM
                    ,rainprod2d, evapprod2d                        &
#endif
                    ,errmsg, errflg                                &
                                                                   )
         DO k=kts,kte
         DO i=its,ite
            th(i,k,j) = t_hv(i,k)
            qc(i,k,j) = qc_hv(i,k)
            qi(i,k,j) = qi_hv(i,k)
            qr(i,k,j) = qr_hv(i,k)
            qs(i,k,j) = qs_hv(i,k)
            qg(i,k,j) = qg_hv(i,k)
            q(i,k,j)  = q_hv(i,k)
         ENDDO
         ENDDO
         DO i=its,ite
            rain(i,j)    = rain_hv(i)
            rainncv(i,j) = rainncv_hv(i)
            sr(i,j)      = sr_hv(i)
         ENDDO
         IF(PRESENT (snowncv) .AND. PRESENT (snow)) THEN
            DO i=its,ite
               snow(i,j) = snow_hv(i)
               snowncv(i,j) = snowncv_hv(i)
            ENDDO
         ENDIF
         IF(PRESENT (graupelncv) .AND. PRESENT (graupel)) THEN
            DO i=its,ite
               graupel(i,j) = graupel_hv(i)
               graupelncv(i,j) = graupelncv_hv(i)
            ENDDO
         ENDIF

         IF (has_reqc.ne.0 .and. has_reqi.ne.0 .and. has_reqs.ne.0) THEN
            do i=its,ite
               do k=kts,kte
                  re_qc(k)   = 2.51E-6
                  re_qi(k)   = 10.01E-6
                  re_qs(k)   = 25.E-6
                  th_hv2(k)  = th(i,k,j) 
                  pii_hv2(k) = pii(i,k,j) 
                  den1d(k)   = den(i,k,j)
                  qc1d(k)    = qc(i,k,j)
                  qi1d(k)    = qi(i,k,j)
                  qs1d(k)    = qs(i,k,j)
               enddo
            CALL effectRad_wsm6(t1d, pii_hv2,th_hv2,qc1d, qi1d, qs1d, den1d, &
                                qmin, t0c, re_qc, re_qi, re_qs,              &
                                kts, kte)
            do k=kts,kte
              re_cloud(i,k,j) = MAX(2.51E-6,  MIN(re_qc(k),  50.E-6))
              re_ice(i,k,j)   = MAX(10.01E-6, MIN(re_qi(k), 125.E-6))
              re_snow(i,k,j)  = MAX(25.E-6,   MIN(re_qs(k), 999.E-6))
            enddo 
          enddo   
        endif     ! has_reqc, etc...
!+---+-----------------------------------------------------------------+
#ifdef WRF_CHEM
        do i=its,ite
          do k=kts,kte
            rainprod(i,k,j) = rainprod2d(i,k)
            evapprod(i,k,j) = evapprod2d(i,k)
          enddo
        enddo
#endif
      ENDDO
  END SUBROUTINE wsm6
END MODULE module_mp_wsm6
