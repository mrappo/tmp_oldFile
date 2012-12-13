      Subroutine vhcuts(p,maxparts,passed)
c--- given the event momenta in p (with maxparts entries),
c--- performs some generic WBF cuts (a la Del Duca et al.)
c---  a point that fails the cuts returns passed=.false.
c--- Note: implements Eq. (3.2) of CEZ paper
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
      double precision JetPtMin, JetEtaCut, JetEtaMax,JetDEtaMin, 
     . JetMjjMin, JetMjjMax
      logical jetveto

************************* Selected CUTS *******************************
     
      JetPtMin   = 30d0
      JetEtaCut  = 2.5d0
      JetEtaMax  = 4.7d0
      JetDEtaMin = 2.1d0
      JetMjjMin  = 65d0
      JetMjjMax  = 105d0
      jetveto    = .true.
      passed     = .false.

***************************** START CUTS *******************************
      found=0
      ptj1=0d0
      ptj2=0d0     
      

c--- check for 2 tagging jets (highest pt) 
      do j=3,maxparts
       if (  (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'qj')
     .  .or. (plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')) then
          ptj=pt(j,p)
c          write(*,*) 'j', j, ptj 
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

c      write(*,*) 'found ',found 

      if (found .lt. 2) goto 999
      if (ptj1  .lt. JetPtMin) goto 999
      if (ptj2  .lt. JetPtMin) goto 999

      ptj3 = 0d0 
      j3 = 0 
      etaj3 = 0d0
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
c
c     Jet Veto over 30 GeV for the
      if (jetveto .and. (ptj3 .gt. JetPtMin) .and. (j3 .gt. 0)) goto 999
      if (jetveto .and. (j3 .gt. 0)) then
         etaj3=etarap(j3,p) 
         if(abs(etaj3) .gt. JetEtaMax) goto 999
      endif
 
      etaj1= 0d0
      etaj2= 0d0
      etaj1=etarap(j1,p)
      etaj2=etarap(j2,p)

c     eta acceptance selection

      if(abs(etaj1) .gt. JetEtaMax) goto 999 
      if(abs(etaj2) .gt. JetEtaMax) goto 999

      if(abs(etaj1) .gt. JetEtaCut) then 
c      write(*,*) 'abs(eta1) ',abs(etaj1)
       goto 999 
      endif
  
      if(abs(etaj2) .gt. JetEtaCut) goto 999


c--- ensure a rapidity gap  between the tagging jets
      if (abs(etaj1-etaj2) .gt. JetDEtaMin)then
c        write(*,*) 'abs(eta1-etaj2) ',abs(etaj1-etaj2)
        goto 999      
      endif
c--- Mjj Selection
      mj1j2=0d0      
      mj1j2=dsqrt(max(0d0,(p(j1,4)+p(j2,4))**2
     . -(p(j1,1)+p(j2,1))**2-(p(j1,2)+p(j2,2))**2-(p(j1,3)+p(j2,3))**2))

c--- ensure the tagging jets have an invariant mass larger than mjjmin GeV
      if ((mj1j2 .gt. JetMjjMax) .or. (mj1j2 .lt. JetMjjMin)) then
c        write(*,*) 'mj1j2 ',mj1j2     
        goto 999      
      endif
********************** END OF CUT CROSS SECTIONS ***********************

      passed=.true.

c--- this is the alternate return if the cuts are failed    
  999 continue
      
      return
      
      end
      
     
