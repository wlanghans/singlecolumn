load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

cases=1
xstr=new(cases,string)
colors=(/"blue","red","orange"/)
colorsp=(/"blue","red","orange","darkgreen"/)


rmsedir="/global/scratch/langhans/metamodel_large3/"
;f=addfile(rmsedir+"/RMSE_mmthvdqtdrmse4_paperreview.nc","r")
f=addfile(rmsedir+"/RMSE_paperreview.nc","r")


p1=f->alphathv
p2=f->alphaqt
p3=f->ctketau
p4=f->pwmin
p5=f->tauneggers
p6=f->Wc

pref=(/1.,550.,1.,1000.,2.,2./)
ind_def=new(6,integer)
ind_def(0)=2;ind(p6.eq.pref(0))
ind_def(1)=3;ind(p5.eq.pref(1))
ind_def(2)=3;ind(p4.eq.pref(2))
ind_def(3)=4;ind(p3.eq.pref(3))
ind_def(4)=1;ind(p2.eq.pref(4))
ind_def(5)=1;ind(p1.eq.pref(5))
ind_def(0)=ind(p6.eq.pref(0))
ind_def(1)=ind(p5.eq.pref(1))
ind_def(2)=ind(p4.eq.pref(2))
ind_def(3)=ind(p3.eq.pref(3))
ind_def(4)=ind(p2.eq.pref(4))
ind_def(5)=ind(p1.eq.pref(5))



p=new((/6,dimsizes(p1)/),float)
p(0,:)=(/p6/)
p(1,:)=(/p5/)
p(2,:)=(/p4/)
p(3,:)=(/p3/)
p(4,:)=(/p2/)
p(5,:)=(/p1/)


pnames=(/"a!BF8!q!F!v!N!","a!B!qt!N!","!F8!t!FB!l!NF!","w!B!min!N!","!F8!t!F!!BF8!e!NF!","b!BF8!e!NF!"/)

 rmse=f->RMSE_thvd_tot
 r=0
 xstr(r)="total"
 ;if (r.eq.1) then
 ;rmse=f->RMSE_CPBL2
 ;xstr(r)="CPBL 0.06"
 ;end if
 ;if (r.eq.2) then
 ;rmse=f->RMSE_CPBL4
 ;xstr(r)="CPBL 0.12"
 ;end if

rmse1D = ndtooned(rmse)
dims=dimsizes(rmse)
ii = new(dims,integer)
ij = new(dims,integer)
ik = new(dims,integer)
il = new(dims,integer)
im = new(dims,integer)
in = new(dims,integer)
do i=0,dims(0)-1
do j=0,dims(1)-1
do k=0,dims(2)-1
do l=0,dims(3)-1
do m=0,dims(4)-1
do n=0,dims(5)-1
  ii(i,j,k,l,m,n)=i
  ij(i,j,k,l,m,n)=j
  ik(i,j,k,l,m,n)=k
  il(i,j,k,l,m,n)=l
  im(i,j,k,l,m,n)=m
  in(i,j,k,l,m,n)=n
end do
end do
end do
end do
end do
end do

ii1D = ndtooned(ii)
ij1D = ndtooned(ij)
ik1D = ndtooned(ik)
il1D = ndtooned(il)
im1D = ndtooned(im)
in1D = ndtooned(in)
delete(ii)
delete(ij)
delete(ik)
delete(il)
delete(im)
delete(in)

opt   = False
; get several statistics
stat=stat_dispersion(rmse1D,opt)
print("Lower octile: " + stat(4))
print("Lower dectile: " + stat(3))
print("Lower 5%: " + stat(24))
print("Lower 1%: " + stat(23))
print("Lower 0.1%: " + stat(22))
best1=stat(23)
best5=stat(24)
best10=stat(3)
best01=stat(22)
print(best01)
loct = stat(23)
num_min=100
ind_min=new((/num_min,6/),integer)
ind_lmin=new(6,integer)
;ind_lmin=(/1,3,6,0/)

val_min=new((/num_min/),float)
ind_max=new((/num_min,6/),integer)
val_max=new((/num_min/),float)

;get six minima
tmp=rmse1D
tmp@_FillValue=-999.

do i=0,num_min-1

ind1D=ind(tmp.eq.min(tmp))
if (dimsizes(ind1D).gt.1) then
 print(i +" More than one min found, stopping")
 tmpi=ind1D(0)
 delete(ind1D)
 ind1D=tmpi
 delete(tmpi)
 ;exit
end if

ind_min(i,0)=ii1D(ind1D)
ind_min(i,1)=ij1D(ind1D)
ind_min(i,2)=ik1D(ind1D)
ind_min(i,3)=il1D(ind1D)
ind_min(i,4)=im1D(ind1D)
ind_min(i,5)=in1D(ind1D)
val_min(i)=rmse1D(ind1D)
tmp(ind1D) = -999.
delete(ind1D)

end do

;get six maxima
tmp=rmse1D
tmp@_FillValue=-999.
do i=0,num_min-1

ind1D=ind(tmp.eq.max(tmp))
if (dimsizes(ind1D).gt.1) then
 print("WARNING: More than one max found; using first one")
 tmpi=ind1D(0)
 delete(ind1D)
 ind1D=tmpi
 delete(tmpi)
end if
ind_max(i,0)=ii1D(ind1D)
ind_max(i,1)=ij1D(ind1D)
ind_max(i,2)=ik1D(ind1D)
ind_max(i,3)=il1D(ind1D)
ind_max(i,4)=im1D(ind1D)
ind_max(i,5)=in1D(ind1D)
val_max(i)=rmse1D(ind1D)
tmp(ind1D) = -999.
delete(ind1D)

end do

print("Case: "+xstr(r))

do i=0,num_min-1
   print("RMSE min "+ sprinti("%i",i+1) + " RMSE="+rmse(ind_min(i,0),ind_min(i,1),ind_min(i,2),ind_min(i,3),ind_min(i,4),ind_min(i,5)) +" for a_thv="+p1(ind_min(i,5))+ " a_qt="+p2(ind_min(i,4))+" ctketau="+p3(ind_min(i,3))  + " pwmin="+p4(ind_min(i,2)) + " tauneggers="+p5(ind_min(i,1)) + " Wc="+p6(ind_min(i,0)) )
end do


pdf=new((/6,dims(0)/),float)
dx=new(6,float)
pdf=0.0
ntot=1
do i=0,dimsizes(dims)-1
  ntot=ntot*dims(i)
end do
nc=0
do i=0,ntot-1

if (rmse1D(i).le.best01) then
nc=nc+1
pdf(0,ii1D(i))=pdf(0,ii1D(i))+1. ; Wc
pdf(1,ij1D(i))=pdf(1,ij1D(i))+1.
pdf(2,ik1D(i))=pdf(2,ik1D(i))+1.
pdf(3,il1D(i))=pdf(3,il1D(i))+1.
pdf(4,im1D(i))=pdf(4,im1D(i))+1.
pdf(5,in1D(i))=pdf(5,in1D(i))+1. ; alphathv
end if 


end do

pdf=pdf/int2flt(nc)

do i=0,5
  dx(i)=(p(i,dimsizes(p(i,:))-1) - p(i,0))/(dimsizes(p(i,:))-1.)
  scaleing=1./dx(i)/(sum(pdf(i,1:dims(0)-2)) + 0.5*(pdf(i,0)+pdf(i,dims(0)-1)) ) 
  pdf(i,:) = pdf(i,:) * scaleing
end do


plotpdf=new(6,graphic)

  wks = gsn_open_wks("eps","pdf01%_true")               ; send graphics to PNG file
;**********************************************
  res=True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@tiYAxisString=" "
  res@xyLineThicknessF=3.
  ;res@gsnYRefLineThicknessF=2.
  ;res@gsnYRefLineColors="gray"
  ;res@gsnYRefLine=(/best1,best5/)
  res@xyMarkLineMode="MarkLines"
  res@xyLineColor= "black";colors(0)
  res@xyMarkerColor="black";colors(0)
  res@xyMarker = 1
  res@xyMarkerSizeF=0.04
  res@tiXAxisFontHeightF=0.03
  res@trYMaxF = 1.0
  res@trYMinF = 0.
  
  opt1=(/0.2,0.6,0.0,0.21,449.,0.55/)
  opt2=(/0.6,0.88,1187.,0.72,407.,0.48/)
  res@gsnXRefLineColor="grey22"
  res@gsnXRefLineThicknessF=4.

  ymaxvec=(/2.,1.2,0.002,2.,0.0085,4.0/)
  do i=0,5

  ;res@gsnXRefLine=(/opt1(i)/)
  res@tiXAxisString=pnames(i)
  res@tiYAxisString="f("+pnames(i)+")"+"/!F8!D!F!"+pnames(i)
  res@tmXBLabels =  p(5-i,:) ; labels for each box
  res@trXMinF =p(5-i,0)
  res@trXMaxF = p(5-i,dims(0)-1)
  res@trYMaxF = ymaxvec(i)

  plotpdf(i)  = gsn_csm_xy (wks,p(5-i,:),pdf(5-i,:),res)
  end do
  
  resp=True
  resp@gsnPanelXWhiteSpacePercent=5.
  gsn_panel(wks,plotpdf,(/2,3/),resp)

; now plot RMSE stats
; box whisker plot

;if (r.eq.0) then
yval=new((/cases,5/),float)
;x=new((/cases/),float)
;x=fspan(-3.,1.,cases)
;end if
yval(r,0) = stat(2)
yval(r,1) = stat(6)
yval(r,2) = stat(8)
yval(r,3) = stat(10)
yval(r,4) = stat(14)

if (r.eq.0) then
;**********************************************
; create plot
;**********************************************
  plot=new(15,graphic)
  dum=new(4,graphic)
  wks = gsn_open_wks("eps","mmrmse_at_def")               ; send graphics to PNG file
;**********************************************
  delete(res)
  res=True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@tiYAxisString="Objective Function Value"
  res@xyLineThicknessF=2.
  ;res@gsnYRefLineThicknessF=2.
  ;res@gsnYRefLineColors="gray"
  ;res@gsnYRefLine=(/best1,best5/)
  res@xyMarkLineMode="Lines"
  res@xyLineColor= colors(0)
  res@xyMarker = 1
  res@xyMarkerSizeF=0.05
  res@trYMaxF = 2.5
  res@trYMinF = 0.
  pmres=True
  pmres@gsMarkerColor = colors(0)
  pmres@gsMarkerIndex = 1
  pmres@gsMarkerSizeF = 0.05

  res@tiXAxisString=pnames(0)
  res@tmXBLabels = p1 ; labels for each box
  range=p1(dimsizes(p1)-1)-p1(0)
  res@trXMinF = p1(0)-0.1*range
  res@trXMaxF = p1(dimsizes(p1)-1)+0.1*range
  
  ;ind_min(0,:)=ind_lmin(:)
  plot(0)  = gsn_csm_xy (wks,p1,rmse(ind_def(0),ind_def(1),ind_def(2),ind_def(3),ind_def(4),0:dimsizes(p(5,:))-1),res)
  res@tiXAxisString=pnames(1)
  delete(res@tmXBLabels)
  res@tmXBLabels = p2 ; labels for each box
  range=p2(dimsizes(p2)-1)-p2(0)
  res@trXMinF = p2(0)-0.1*range
  res@trXMaxF = p2(dimsizes(p2)-1)+0.1*range
  plot(1)  = gsn_csm_xy (wks,p2,rmse(ind_def(0),ind_def(1),ind_def(2),ind_def(3),0:dimsizes(p(4,:))-1,ind_def(5)),res)
  ;dum(1) = gsn_add_polymarker(wks, plot(1), p2(ind_min(0,1)),rmse(ind_min(0,0),ind_min(0,1),ind_min(0,2),ind_min(0,3)), pmres)

  delete(res@tmXBLabels)
  res@tiXAxisString=pnames(2)
  res@tmXBLabels = p3 ; labels for each box
  range=p3(dimsizes(p3)-1)-p3(0)
  res@trXMinF = p3(0)-0.1*range
  res@trXMaxF = p3(dimsizes(p3)-1)+0.1*range
   
  plot(2)  = gsn_csm_xy (wks,p3,rmse(ind_def(0),ind_def(1),ind_def(2),0:dimsizes(p(3,:))-1,ind_def(4),ind_def(5)),res)

  res@tiXAxisString="a";pnames(2)
  delete(res@tmXBLabels)
;  res@trXMaxF = 0.5*(1.-erf((p3(0)-0.1*range)/sqrt(2.)))
;  res@trXMinF = 0.5*(1.-erf((p3(dimsizes(p3)-1)+0.1*range)/sqrt(2.)))
  af=0.5*(1.-erf(p4/sqrt(2.)))
  range=af(0)-af(dimsizes(af)-1)
  res@trXMinF = 0.;af(dimsizes(af)-1)+0.1*range
  res@trXMaxF = 0.5;af(0)-0.1*range
  plot(3)  = gsn_csm_xy (wks,af,rmse(ind_def(0),ind_def(1),0:dimsizes(p(2,:))-1,ind_def(3),ind_def(4),ind_def(5)),res)
  afx=0.5*(1.-erf((p3(ind_min(0,2)))/sqrt(2.)))
 ; dum(3) = gsn_add_polymarker(wks, plot(2), afx,rmse(ind_min(0,0),ind_min(0,1),ind_min(0,2),ind_min(0,3)), pmres)

  res@tiXAxisString=pnames(4)
  delete(res@tmXBLabels)
  res@tmXBLabels = p5 ; labels for each box
  range=p5(dimsizes(p5)-1)-p5(0)
  res@trXMinF = p5(0)-0.1*range
  res@trXMaxF = p5(dimsizes(p5)-1)+0.1*range
  plot(4)  = gsn_csm_xy (wks,p5,rmse(ind_def(0),0:dimsizes(p(1,:))-1,ind_def(2),ind_def(3),ind_def(4),ind_def(5)),res)

  res@tiXAxisString=pnames(5)
  delete(res@tmXBLabels)
  res@tmXBLabels = p6 ; labels for each box
  range=p6(dimsizes(p6)-1)-p6(0)
  res@trXMinF = p6(0)-0.1*range
  res@trXMaxF = p6(dimsizes(p5)-1)+0.1*range
  plot(5)  = gsn_csm_xy (wks,p6,rmse(0:dimsizes(p(0,:))-1,ind_def(1),ind_def(2),ind_def(3),ind_def(4),ind_def(5)),res)

  gsn_panel(wks,plot(0:5),(/2,3/),False)

  delete(res@tmXBLabels)


  
  
  wks = gsn_open_wks("eps","cont_rmse")
  resc=True
  resc@gsnDraw     = False
  resc@gsnFrame    = False
  resc@gsnSpreadColors = True
  resc@gsnSpreadColorStart=45
  resc@cnInfoLabelOn = False
  resc@cnFillOn = True
  resc@cnLinesOn = True
  resc@cnLevelSelectionMode="ManualLevels"
  resc@cnMinLevelValF=0.4
  resc@cnMaxLevelValF=2.0
  resc@cnLevelSpacingF=0.05


  plotc=new(20,graphic) 

  coorname=(/"d0","d1","d2","d3","d4","d5"/)
  rmse!0="d0"
  rmse&d0=p6
  rmse!1="d1"
  rmse&d1=p5
  rmse!2="d2"
  rmse&d2=p4
  rmse!3="d3"
  rmse&d3=p3
  rmse!4="d4"
  rmse&d4=p2
  rmse!5="d5"
  rmse&d5=p1

  resc@lbOrientation = "vertical"
 
  numc=0
  do i=0,5
  do j=0,5
    if (j.gt.i) then
    

      resc@tiYAxisString          = pnames(5-j)
      resc@tiXAxisString          = pnames(5-i)
      ind_xname=ind(coorname.eq.("d"+sprinti("%0.1i",i)))
      ind_yname=ind(coorname.eq.("d"+sprinti("%0.1i",j)))
      ind_other=ind(coorname.ne.("d"+sprinti("%0.1i",i)).and.coorname.ne.("d"+sprinti("%0.1i",j)))
      data=rmse($coorname(ind_other(0))$|ind_def(ind_other(0)),$coorname(ind_other(1))$|ind_def(ind_other(1)),$coorname(ind_other(2))$|ind_def(ind_other(2)),\
      $coorname(ind_other(3))$|ind_def(ind_other(3)),$coorname(ind_yname)$|:,$coorname(ind_xname)$|:)

      data!0=coorname(ind_yname)
      data&$coorname(ind_yname)$=p(ind_yname,:)
      data!1=coorname(ind_xname)
      data&$coorname(ind_xname)$=p(ind_xname,:)
 
      plotc(numc) = gsn_csm_contour(wks,data,resc)
      numc=numc+1

      delete(data)
    end if 

  end do
  end do
 
  ;gsn_panel(wks,plotc(0:numc-1),(/3,5/),False)
  gsn_panel(wks,plotc(9),(/1,1/),False)

 

lskip=True

  if (.not.lskip) then

  statsd=new((/4,dimsizes(stat)/),float)
  drmsedx=new(dimsizes(rmse),float)
  xp=fspan(-3.,0.,4)

  yvalp=new((/4,5/),float)

  best20stat=new((/3,4/),float)


  do ri=0,3
    if (ri.eq.0) then 
    dx=p1(dimsizes(p1)-1)-p1(0)
    end if
    if (ri.eq.1) then 
    dx=p2(dimsizes(p2)-1)-p2(0)
    end if
    if (ri.eq.2) then 
    dx=p3(dimsizes(p3)-1)-p3(0)
    end if
    if (ri.eq.3) then 
    dx=p4(dimsizes(p4)-1)-p4(0)
    end if

    if (ri.eq.0) then
      num20=0
      best20stat=0.0
    end if
    if (ri.eq.1) then
      best20stat(0,:) = best20stat(0,:)/int2flt(num20)
    end if

    do i=0,dims(0)-1
    do j=0,dims(1)-1
    do k=0,dims(2)-1
    do l=0,dims(3)-1

     if (ri.eq.0) then
     if (rmse(i,j,k,l).le.loct) then
       best20stat(0,0) = best20stat(0,0) + p1(i)
       best20stat(0,1) = best20stat(0,1) + p2(j)
       best20stat(0,2) = best20stat(0,2) + p3(k)
       best20stat(0,3) = best20stat(0,3) + p4(l)
       num20=num20+1
     end if
     end if

     if (ri.eq.1) then
     if (rmse(i,j,k,l).le.loct) then
       best20stat(1,0) = best20stat(1,0) + (p1(i)-best20stat(0,0))^2
       best20stat(1,1) = best20stat(1,1) + (p2(j)-best20stat(0,1))^2
       best20stat(1,2) = best20stat(1,2) + (p3(k)-best20stat(0,2))^2
       best20stat(1,3) = best20stat(1,3) + (p4(l)-best20stat(0,3))^2
     end if
     end if

     
   
     if (ri.eq.0) then
       indu = min((/i+1,dimsizes(p1)-1/))
       indd = max((/i-1,0/))
       drmsedx(i,j,k,l) = (rmse(indu,j,k,l)-rmse(indd,j,k,l)) / (p1(indu)-p1(indd))
     end if
     if (ri.eq.1) then
       indu = min((/j+1,dimsizes(p2)-1/))
       indd = max((/j-1,0/))
       drmsedx(i,j,k,l) = (rmse(i,indu,k,l)-rmse(i,indd,k,l)) / (p2(indu)-p2(indd))
     end if
     if (ri.eq.2) then
       indu = min((/k+1,dimsizes(p3)-1/))
       indd = max((/k-1,0/))
       drmsedx(i,j,k,l) = (rmse(i,j,indu,l)-rmse(i,j,indd,l)) / (p3(indu)-p3(indd))
     end if
     if (ri.eq.3) then
       indu = min((/l+1,dimsizes(p4)-1/))
       indd = max((/l-1,0/))
       drmsedx(i,j,k,l) = (rmse(i,j,k,indu)-rmse(i,j,k,indd)) / (p4(indu)-p4(indd))
     end if
    
    end do
    end do
    end do
    end do
 
 
    drmsedx = drmsedx * dx


    statsd(ri,:)=stat_dispersion(drmsedx,opt)
  end do


  best20stat(1,:) =best20stat(1,:) / int2flt(num20) 
  best20avg=best20stat(0,:)
  do ri=0,3
    if (ri.eq.0) then 
    dx=p1(dimsizes(p1)-1)-p1(0)
    ps=p1(0)
    end if
    if (ri.eq.1) then 
    dx=p2(dimsizes(p2)-1)-p2(0)
    ps=p2(0)
    end if
    if (ri.eq.2) then 
    dx=p3(dimsizes(p3)-1)-p3(0)
    ps=p3(0)
    end if
    if (ri.eq.3) then 
    dx=p4(dimsizes(p4)-1)-p4(0)
    ps=p4(0)
    end if
    print("Average among top 1% for parameter "+ ri +"="+best20stat(0,ri))
    best20stat(0,ri) = (best20stat(0,ri)-ps)/ dx
    best20stat(1,ri) = sqrt(best20stat(1,ri))/dx
    delete(ps)
  end do

  yvalp(:,0) = statsd(:,2)
  ;yvalp(:,1) = statsd(:,6)
  yvalp(:,1) = statsd(:,3)
  yvalp(:,2) = statsd(:,8)
  ;yvalp(:,3) = statsd(:,10)
  yvalp(:,3) = statsd(:,13)
  yvalp(:,4) = statsd(:,14)

  ;**********************************************
  ; create plot for drmsedx
  ;**********************************************
  wks = gsn_open_wks("eps","drmsedx_box")               ; send graphics to PNG file
  ;**********************************************
  ; resources for plot background
  ;**********************************************
  delete(res)
  res            = True                         ; plot mods desired
 ; res@gsnDraw     = False
 ; res@gsnFrame    = False
  pnames2=pnames
  pnames2(2)="a"
  res@tmXBLabels = pnames2 ; labels for each box
  res@trYMinF = min(statsd(:,2))
  res@trYMaxF = max(statsd(:,14))
  res@tiMainString = "min/10%/median/90%/max"
  res@tiYAxisString="!F18!s!F!r!H-10V2F35!D!FV-2H3!/!F18!s!F!x!H-10V2F35!D!FV-2H3!"
  ;**********************************************
  ; resources for polylines that draws the boxes
  ;**********************************************
  llres                   = True
  llres@gsLineThicknessF  = 2.5                 ; line thickness
  ;**********************************************
  ; resources that control color and width of boxes
  ;**********************************************
  opti          = True
  opti@boxWidth = .25                           ; Width of box (x units)
  opti@boxColors = colorsp       ; Color of box(es)
  ;***********************************************
  plot(5) = boxplot(wks,xp,yvalp,opti,res,llres)  ; All 3 options used...
  draw(plot(5))                                     ; box plot does not call
  frame(wks)  

  delete(res@tmXBLabels)

end if


  wks = gsn_open_wks("eps","best20average")                         ; send graphics to PNG file
 
  delete(res)

  res                        = True                      ; plot mods desired
  res@tiYAxisString          = "x!H-10V2F35!D!FV-2H3! !B1!-!S1!+!N! 1/2!F8!s!F! !B!x!H-10V2F35!D!FV-2H3!!N!" ; add units title 
  res@gsnFrame               = False                     ; don't draw yet
  res@gsnDraw                = False                     ; don't advance frame
  res@tmXBMode = "Explicit"
  res@tmXBLabels = pnames ; labels for each box
  ;res@tmXBValues = xp ; labels for each box
  
  res@xyMarkLineMode         = "Markers"
  res@xyMarker               = 16
  res@xyMarkerSizeF          = 0.02
  res@trYMaxF=1.
  res@trYMinF=0.
  ;res@trXMaxF=xp(3)+1.
  ;res@trXMinF=xp(0)-1.
  ;best20stat!1="param"
  ;best20stat&param=xp
  res@tiMainString = "Mean and !F8!s!F! for best 1%"

  ;plot(7) = gsn_csm_xy(wks,xp,best20stat(0,0:3),res) 

  polyres                   = True                       ; marker resources
  polyres@gsMarkerIndex     = 1                          ; index for circle
  polyres@gsMarkerSizeF     = .02                        ; size 
  polyres@gsLineColor= "blue"
  polyres@gsLineThicknessF=3.
  
; gsn_add* templates are functions that we set to dummy values. Since
; we are going to draw numerous error bars, we create two arrays to
; hold the dummy values.
  error_bar = new(4,graphic)

; loop through the points  
;  do t=0,3
; the error bar is draw in two segements, above and below the true
; value of u(t). I made up these value. You would change this to 
; true error.      
;    error_bar(t) = gsn_add_polyline(wks,plot(7),(/xp(t),xp(t)/),(/best20stat(0,t)-0.5*best20stat(1,t),best20stat(0,t)+0.5*best20stat(1,t)/),\
;    polyres)
;  end do
  ;draw(plot(7))
  ;frame(wks)


;**********************************************
; create plot
;**********************************************
  wks = gsn_open_wks("eps","box")               ; send graphics to PNG file
;**********************************************
; resources for plot background
;**********************************************
  rmsepdf=pdfx(rmse, 30, False)

  delete(res)
  ;delete(opti)
  res            = True                         ; plot mods desired
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@xyLineThicknessF       = 2
  res@tiYAxisString          = "PDF (%)"
  res@tiXAxisString          = "RMSE r!H-10V2F35!D!FV-2H3!"
  res@xyMarkLineMode         = "MarkLines"
  res@xyMarker               = 16
 
  plot(4) = gsn_csm_xy (wks, rmsepdf@bin_center, rmsepdf, res)
  draw(plot(4))                                     ; box plot does not call
  frame(wks)                                    ; these for you


  end if  ; lskip

  dims=dimsizes(rmse)
  numall= num(rmse.le.val_min(99))

  frac=new((/dims(2),dims(3)/),float)
  frac!0=rmse!2
  frac!1=rmse!3
  frac&d2=p4
  frac&d3=p3
  ;frac&pwmin=0.5*(1.-erf(rmse&pwmin/sqrt(2.)))
  ;frac&ctketau=rmse&ctketau
  ;rmse20b=new((/dims(0),dims(3)/),float)
  ;rmse20b!0=rmse!0
  ;rmse20b!1=rmse!3
  ;rmse20b&Wc=rmse&Wc
  ;rmse20b&ctketau=rmse&ctketau
  ;fracb=rmse20b

  rmse20=frac
  do i=0,dims(2)-1
  do j=0,dims(3)-1
    frac(i,j)=int2flt(num( rmse(:,:,i,j,:,:) .le.val_min(99)  ))
    ;rmse20(i,j)=sum(where(rmse(:,:,i,j,:,:).le.val_min(99),rmse(:,:,i,j,:,:),0.0))
    ;rmse20(i,j)=min(where(rmse(:,:,i,j,:,:).le.val_min(99),rmse(:,:,i,j,:,:),10000000.0))
    rmse20(i,j)=min(rmse(:,:,i,j,:,:))
  end do
  end do

 ; do i=0,dims(0)-1
 ; do j=0,dims(3)-1
 ;   fracb(i,j)=int2flt(num( rmse(i,:,:,j) .le.loct  ))
 ;   rmse20b(i,j)=sum(where(rmse(i,:,:,j).le.loct,rmse(i,:,:,j),0.0))
 ; end do
 ; end do
 ; 
  ;rmse20=rmse20/(frac+1.e-10)
  ;rmse20=where(frac.gt.0.0,rmse20,-999.)
  rmse20@_FillValue=-999.

  frac= frac / int2flt(numall) /(p4(1)-p4(0))/(p3(1)-p3(0))
  ;rmse20b=rmse20b/(fracb+1.e-10)
  ;rmse20b=where(fracb.gt.0.0,rmse20b,-999.)
  ;rmse20b@_FillValue=-999.
  ;fracb= fracb / int2flt(numall) * 100.
 
;**********************************************
; create plot
;**********************************************
  wks = gsn_open_wks("eps","jpdf")               ; send graphics to PNG file
;**********************************************
  gsn_define_colormap(wks,"matlab_jet")
  
  delete(stat)
  ;delete(res)
  res      = True
  res@gsnSpreadColors = True
  res@gsnSpreadColorStart=3
  res@gsnDraw       = False
  res@gsnFrame      = False
  res@cnInfoLabelOn = False
  res@tiXAxisString          = "w!B!min!N!"
  res@tiYAxisString          = "!F8!t!F!!B!l!N! (s)"
  res@cnFillOn = True
  res@cnLinesOn = False
  res@cnFillMode = "RasterFill"
  res@cnMissingValFillPattern = 0 
  res@cnMissingValFillColor = "white"
  res@lbOrientation = "Vertical"
  res@tmYLMinorPerMajor = 4
  res@tmYRMinorPerMajor = 4
  frac=where(frac.eq.0.0,-999.,frac)
;  fracb=where(fracb.eq.0.0,-999.,fracb)
;  fracb@_FillValue=-999.
  frac@_FillValue=-999.
  
  plot(6) = gsn_csm_contour (wks,frac(d3|:,d2|:), res)

  res@cnLevelSelectionMode="ManualLevels"
  res@cnMinLevelValF=0.38
  res@cnMaxLevelValF=0.58
  res@cnLevelSpacingF=0.02
  
  plot(9) = gsn_csm_contour (wks,rmse20(d3|:,d2|:), res)

  draw(plot(9))
  frame(wks)

 ; res@tiYAxisString          = "b"
 ; res@tiXAxisString          = "mixing time scale !F8!t!F!!B!l!N! (s)"
 ; plot(9) = gsn_csm_contour (wks,fracb, res)
 ; resP                    = True                 ; modify the panel plot
 ; resP@gsnFrame           = False                ; don't advance panel plot
 ; resP@gsnPanelLabelBar   = False                 ; add common colorbar
 ; resP@gsnPanelMainString = "Distribution (%) of best 10%"     ; new way of setting main title (as of NCL V6.4.0)
 ; gsn_panel(wks,plot((/6,9/)),(/1,2/),False)


exit
;**********************************************
; create plot
;**********************************************
  wks = gsn_open_wks("eps","rmse20")               ; send graphics to PNG file
;**********************************************
  gsn_define_colormap(wks,"MPL_BuGn")
  delete(res)
  res      = True
  res@gsnDraw       = False
  res@gsnFrame      = False
  res@cnInfoLabelOn = False
  res@tiYAxisString          = "fractional area a"
  res@tiXAxisString          = "mixing time scale !F8!t!F!!B!l!N! (s)"
  res@cnFillOn = True
  res@cnLinesOn = False
  res@cnFillMode = "RasterFill"
  res@lbOrientation = "Vertical"
  res@tmYLMinorPerMajor = 4
  res@tmYRMinorPerMajor = 4
  res@cnMissingValFillColor="gray"
  
  plot(8) = gsn_csm_contour (wks,rmse20, res)

  res@tiYAxisString          = "b"
  plot(10) = gsn_csm_contour (wks,rmse20b, res)
  resP                    = True                 ; modify the panel plot
  resP@gsnFrame           = False                ; don't advance panel plot
  resP@gsnPanelLabelBar   = False                 ; add common colorbar
  resP@gsnPanelMainString = "Average RMSE r!H-10V2F35!D!FV-2H3! of best 10%"     ; new way of setting main title (as of NCL V6.4.0)
  gsn_panel(wks,plot((/8,10/)),(/1,2/),resP)
 
  

  




