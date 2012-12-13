      subroutine wbfcuts(p,maxparts,passed)
c--- given the event momenta in p (with maxparts entries),
c--- performs some generic WBF cuts (a la Del Duca et al.)
c---  a point that fails the cuts returns passed=.false.
c--- Note: implements Eq. (3.2) of CEZ paper
c PG---      use sub_defs_io 
      implicit none
      include 'constants.f'
      logical passed, first 
      integer j,maxparts,found,j1,j2,j3
      double precision pt,etarap,p(mxpart,4),ptj,ptj1,ptj2,
     . etaj1,etaj2,mj1j2
      double precision ptj3,etaj3
      character*2 plabel(mxpart)
      common/plabel/plabel
      data first/.true./
      logical jetveto
      double precision JetPtMin,JetEtaCut,JetEtaMax,JetDEtaMin,
     . JetMjjMin,JetMjjMax
      

************************* Selected CUTS *******************************                                                                                                                                                                                                                                    
      JetPtMin   = 30d0
      JetEtaCut  = 4.7d0
      JetEtaMax  = 4.7d0
      JetDEtaMin = 3.5d0
      JetMjjMin  = 500d0
      JetMjjMax  = 10000d0
      jetveto    = .true.
      passed     = .false.      

***************************** START CUTS *******************************

      found=0
      ptj1=0d0
      ptj2=0d0

!      write(*,*) 'maxparts, plabel(3:)',maxparts, plabel(3:)
c--- check for 2 tagging jets (highest pt) 

      do j=3,maxparts
       if (  (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'qj')
     .  .or. (plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')) then
          ptj=pt(j,p)
!          write(*,*) 'j', j, ptj 
          if (ptj .gt. min(ptj1,ptj2)) then
            found=found+1
            if (found .eq. 1) then
              ptj1=ptj
              j1=j
            endif
            if (found .ge. 2) then
              if (ptj .gt. ptj1) then
                 ptj2=ptj1
                 ptj1=ptj
                 j2=j1
                 j1=j
              else
                 ptj2=ptj
                 j2=j
              endif
            endif
          endif
        endif
      enddo
      
c     At least two jet + basic pT and eta cut

      if (found .lt. 2) goto 999
      if (ptj1  .lt. JetPtMin) goto 999
      if (ptj2  .lt. JetPtMin) goto 999
c      write(*,*)'ptj1 ',ptj1,' ptj2',ptj2

      etaj1= 0d0
      etaj1= etarap(j1,p)
      etaj2= 0d0
      etaj2=etarap(j2,p)

      if(abs(etaj1)>JetEtaMax) goto 999
      if(abs(etaj2)>JetEtaMax) goto 999

c     In addition tighter selection on eta
      if(abs(etaj1)>JetEtaCut) goto 999
      if(abs(etaj2)>JetEtaCut) goto 999
c      write(*,*)'etaj1 ',etaj1,' etaj2',etaj2

c     Search third jet 
      ptj3 = 0d0 
      j3 = 0 
c      -- find, if it exists a third jet 
      do j=3,maxparts
        if (  (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'qj')
     .   .or. (plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')) then
           if ((j.ne. j1) .and. (j.ne.j2)) then  
              if (j3 .eq. 0) then 
                 ptj3=pt(j,p)
                 j3 = j 
              else
c                 write(*,*) 'ptj3 already defined?',j,j3 
              endif
           endif
        endif
      enddo

c     Eta acceptance also on the third jet
      etaj3 = 0d0
      if (j3 .ne. 0 .and. jetveto) then 
       etaj3=etarap(j3,p)
       if(abs(etaj3) .gt. JetEtaMax) goto 999
      endif
c      if(ptj3.ne.0) write(*,*)' ptj3',ptj3,' etaj3',etaj3
      if (jetveto) then 
         if ((ptj3 > JetPtMin) .and. 
     .        ((etaj1 < etaj3 .and. etaj3 < etaj2) .or. 
     .        (etaj2 < etaj3 .and. etaj3 < etaj1)))then
c                 write(*,*)'CJV Applied'
                 goto 999
              endif 
      endif

c--- ensure a rapidity gap of at least 4.2 between the tagging jets
       if (abs(etaj1-etaj2) .lt. JetDEtaMin) goto 999 
c       write(*,*) 'abs(etaj1-etaj2)',abs(abs(etaj1-etaj2))

c--- ensure the tagging jets lie in opposite hemispheres
c      if (etaj1*etaj2 .ge. 0d0) goto 999      
       mj1j2 = 0d0      
       mj1j2=dsqrt(max(0d0,(p(j1,4)+p(j2,4))**2
     . -(p(j1,1)+p(j2,1))**2-(p(j1,2)+p(j2,2))**2-(p(j1,3)+p(j2,3))**2))

c--- ensure the tagging jets have an invariant mass larger than mjjmin GeV
      if ((mj1j2 .lt. JetMjjMin) .or. (mj1j2 .gt. JetMjjMax)) goto 999      
c      write(*,*) 'mj1j2 ',mj1j2
********************** END OF CUT CROSS SECTIONS ***********************

      passed=.true.

c--- this is the alternate return if the cuts are failed    
  999 continue
      
      return
      
      end
      
     
