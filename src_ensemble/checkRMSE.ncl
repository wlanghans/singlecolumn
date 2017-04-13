load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

cases=1
xstr=new(cases,string)
colors=(/"blue","red","orange"/)
colorsp=(/"blue","red","orange","darkgreen"/)


rmsedir="/global/scratch/langhans/opt_ensemble_test/"
f=addfile(rmsedir+"/RMSE.nc","r")

p1=f->Wc
p2=f->tauneggers
p3=f->pwmin
p4=f->ctketau

pnames=(/"b!BF8!e!FN!","!F8!t!F!!BF8!e!NF!","P!B!wmin!N!","!F8!t!F!!B!l!N!"/)

do r=0,cases-1
 if (r.eq.0) then
 rmse=f->RMSE_BOMEX 
 xstr(r)="BOMEX"
 end if
 if (r.eq.1) then
 rmse=f->RMSE_CPBL2
 xstr(r)="CPBL 0.06"
 end if
 if (r.eq.2) then
 rmse=f->RMSE_CPBL4
 xstr(r)="CPBL 0.12"
 end if

rmse1D = ndtooned(rmse)
dims=dimsizes(rmse)
ii = new(dims,integer)
ij = new(dims,integer)
ik = new(dims,integer)
il = new(dims,integer)
do i=0,dims(0)-1
do j=0,dims(1)-1
do k=0,dims(2)-1
do l=0,dims(3)-1
  ii(i,j,k,l)=i
  ij(i,j,k,l)=j
  ik(i,j,k,l)=k
  il(i,j,k,l)=l
end do
end do
end do
end do

ii1D = ndtooned(ii)
ij1D = ndtooned(ij)
ik1D = ndtooned(ik)
il1D = ndtooned(il)
delete(ii)
delete(ij)
delete(ik)
delete(il)

opt   = False
; get several statistics
stat=stat_dispersion(rmse1D,opt)
ind_min=new((/3,4/),integer)
val_min=new((/3/),float)
ind_max=new((/3,4/),integer)
val_max=new((/3/),float)

;get three minima
ind1D=ind(rmse1D.eq.min(rmse1D))
if (dimsizes(ind1D)gt.1) then
 print("More than one min/max found, stopping")
 exit
end if

ind_min(0,0)=ii1D(ind1D)
ind_min(0,1)=ij1D(ind1D)
ind_min(0,2)=ik1D(ind1D)
ind_min(0,3)=il1D(ind1D)
val_min(0)=rmse1D(ind1D)
tmp=rmse1D
tmp(ind1D) = -999.
tmp@_FillValue=-999.
delete(ind1D)
ind1D=ind(tmp.eq.min(tmp))
if (dimsizes(ind1D)gt.1) then
 print("More than one min/max found, stopping")
 exit
end if
ind_min(1,0)=ii1D(ind1D)
ind_min(1,1)=ij1D(ind1D)
ind_min(1,2)=ik1D(ind1D)
ind_min(1,3)=il1D(ind1D)
val_min(1)=tmp(ind1D)
tmp(ind1D) = -999.
tmp@_FillValue=-999.
delete(ind1D)
ind1D=ind(tmp.eq.min(tmp))
if (dimsizes(ind1D)gt.1) then
 print("More than one min/max found, stopping")
 exit
end if
ind_min(2,0)=ii1D(ind1D)
ind_min(2,1)=ij1D(ind1D)
ind_min(2,2)=ik1D(ind1D)
ind_min(2,3)=il1D(ind1D)
val_min(2)=tmp(ind1D)
;get three maxima
delete(ind1D)
ind1D=ind(rmse1D.eq.max(rmse1D))
if (dimsizes(ind1D)gt.1) then
 print("More than one min/max found, stopping")
 exit
end if
ind_max(0,0)=ii1D(ind1D)
ind_max(0,1)=ij1D(ind1D)
ind_max(0,2)=ik1D(ind1D)
ind_max(0,3)=il1D(ind1D)
val_max(0)=rmse1D(ind1D)
tmp=rmse1D
tmp(ind1D) = -999.
tmp@_FillValue=-999.
delete(ind1D)
ind1D=ind(tmp.eq.max(tmp))
if (dimsizes(ind1D)gt.1) then
 print("More than one min/max found, stopping")
 exit
end if
ind_max(1,0)=ii1D(ind1D)
ind_max(1,1)=ij1D(ind1D)
ind_max(1,2)=ik1D(ind1D)
ind_max(1,3)=il1D(ind1D)
val_max(1)=tmp(ind1D(0))
tmp(ind1D) = -999.
tmp@_FillValue=-999.
delete(ind1D)
ind1D=ind(tmp.eq.max(tmp))
if (dimsizes(ind1D)gt.1) then
 print("More than one min/max found, stopping")
 exit
end if
ind_max(2,0)=ii1D(ind1D)
ind_max(2,1)=ij1D(ind1D)
ind_max(2,2)=ik1D(ind1D)
ind_max(2,3)=il1D(ind1D)
val_max(2)=tmp(ind1D)

print("Case: "+xstr(r))

do i=0,2
 print("RMSE min "+ i + " RMSE="+rmse(ind_min(i,0),ind_min(i,1),ind_min(i,2),ind_min(i,3)) +" for Wc="+p1(ind_min(i,0)) + " tauneggers="+p2(ind_min(i,1)) + " pwmin="+p3(ind_min(i,2)) +" ctketau="+p4(ind_min(i,3)) )
 print("RMSE max "+ i + " RMSE="+rmse(ind_max(i,0),ind_max(i,1),ind_max(i,2),ind_max(i,3)) +" for Wc="+p1(ind_max(i,0)) + " tauneggers="+p2(ind_max(i,1)) + " pwmin="+p3(ind_max(i,2)) +" ctketau="+p4(ind_max(i,3)) )
end do

; now plot RMSE stats
; box whisker plot

if (r.eq.0) then
yval=new((/cases,5/),float)
;x=new((/cases/),float)
;x=fspan(-3.,1.,cases)
end if
yval(r,0) = stat(2)
yval(r,1) = stat(6)
yval(r,2) = stat(8)
yval(r,3) = stat(10)
yval(r,4) = stat(14)

if (r.eq.0) then
;**********************************************
; create plot
;**********************************************
  plot=new(6,graphic)
  dum=new(4,graphic)
  wks = gsn_open_wks("eps","rmse_at_opt")               ; send graphics to PNG file
;**********************************************
  res=True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@tiYAxisString="rmse"
  res@xyLineThicknessF=3.
  ;res@gsnXRefLineThicknessF=7.
  ;res@gsnXRefLineColor="blue"
  ;res@gsnXRefLine=p1(ind_min(0,0))
  res@xyMarkLineMode="MarkLines"
  res@xyMarker = 1
  res@xyMarkerSizeF=0.05
  res@trYMaxF = max(rmse)
  res@trYMinF = 0.
  pmres=True
  pmres@gsMarkerColor = "blue"
  pmres@gsMarkerIndex = 1
  pmres@gsMarkerSizeF = 0.05

  res@tiXAxisString=pnames(0)
  res@tmXBLabels = p1 ; labels for each box
  range=p1(dimsizes(p1)-1)-p1(0)
  res@trXMinF = p1(0)-0.1*range
  res@trXMaxF = p1(dimsizes(p1)-1)+0.1*range
  plot(0)  = gsn_csm_xy (wks,p1,rmse(0:dimsizes(p1)-1,ind_min(0,1),ind_min(0,2),ind_min(0,3)),res)
  dum(0) = gsn_add_polymarker(wks, plot(0), p1(ind_min(0,0)),rmse(ind_min(0,0),ind_min(0,1),ind_min(0,2),ind_min(0,3)), pmres)
  res@tiXAxisString=pnames(1)
  delete(res@tmXBLabels)
  res@tmXBLabels = p2 ; labels for each box
  range=p2(dimsizes(p2)-1)-p2(0)
  res@trXMinF = p2(0)-0.1*range
  res@trXMaxF = p2(dimsizes(p2)-1)+0.1*range
  plot(1)  = gsn_csm_xy (wks,p2,rmse(ind_min(0,0),0:dimsizes(p2)-1,ind_min(0,2),ind_min(0,3)),res)
  dum(1) = gsn_add_polymarker(wks, plot(1), p2(ind_min(0,1)),rmse(ind_min(0,0),ind_min(0,1),ind_min(0,2),ind_min(0,3)), pmres)
  res@tiXAxisString=pnames(2)
  delete(res@tmXBLabels)
  res@tmXBLabels = p3 ; labels for each box
  range=p3(dimsizes(p3)-1)-p3(0)
  res@trXMinF = p3(0)-0.1*range
  res@trXMaxF = p3(dimsizes(p3)-1)+0.1*range
  plot(2)  = gsn_csm_xy (wks,p3,rmse(ind_min(0,0),ind_min(0,1),0:dimsizes(p3)-1,ind_min(0,3)),res)
  dum(2) = gsn_add_polymarker(wks, plot(2), p3(ind_min(0,2)),rmse(ind_min(0,0),ind_min(0,1),ind_min(0,2),ind_min(0,3)), pmres)
  res@tiXAxisString=pnames(3)
  delete(res@tmXBLabels)
  res@tmXBLabels = p4 ; labels for each box
  range=p4(dimsizes(p4)-1)-p4(0)
  res@trXMinF = p4(0)-0.1*range
  res@trXMaxF = p4(dimsizes(p4)-1)+0.1*range
  plot(3)  = gsn_csm_xy (wks,p4,rmse(ind_min(0,0),ind_min(0,1),ind_min(0,2),0:dimsizes(p4)-1),res)
  dum(3) = gsn_add_polymarker(wks, plot(3), p4(ind_min(0,3)),rmse(ind_min(0,0),ind_min(0,1),ind_min(0,2),ind_min(0,3)), pmres)

  gsn_panel(wks,plot(0:3),(/1,4/),False)

  delete(res@tmXBLabels)

  statsd=new((/4,dimsizes(stat)/),float)
  drmsedx=new(dimsizes(rmse),float)
  xp=fspan(-3.,0.,4)

  yvalp=new((/4,5/),float)


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

    do i=0,dims(0)-1
    do j=0,dims(1)-1
    do k=0,dims(2)-1
    do l=0,dims(3)-1
   
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
  res@tmXBLabels = pnames ; labels for each box
  res@trYMinF = min(statsd(:,2))
  res@trYMaxF = max(statsd(:,14))
  res@tiMainString = "dRMSE/dx!B!i!N!; min/10%/median/90%/max"
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

; end loop over cases
end do


;**********************************************
; create plot
;**********************************************
;  wks = gsn_open_wks("eps","box")               ; send graphics to PNG file
;**********************************************
; resources for plot background
;**********************************************
;  delete(res)
;  delete(opti)
;  res            = True                         ; plot mods desired
; ; res@gsnDraw     = False
; ; res@gsnFrame    = False
;  res@tmXBLabels = xstr ; labels for each box
;  res@trYMinF = 0.
;  res@trYMaxF = 1.
;  res@tiMainString = "RMSE; min/25%/median/75%/max"
;;**********************************************
;; resources for polylines that draws the boxes
;;**********************************************  
;  llres                   = True			
;  llres@gsLineThicknessF  = 2.5                 ; line thickness 
;;**********************************************
;; resources that control color and width of boxes
;;**********************************************  
;  opti          = True			
;  opti@boxWidth = .25				; Width of box (x units)
;  opti@boxColors = colors  	; Color of box(es)
;;***********************************************
;  plot(4) = boxplot(wks,x,yval,opti,res,llres)	; All 3 options used...
;  draw(plot(4))                                     ; box plot does not call
;  frame(wks)                                    ; these for you
;  delete(res@tmXBLabels)

  




