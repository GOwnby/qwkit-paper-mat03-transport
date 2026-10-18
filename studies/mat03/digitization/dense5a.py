"""Dense digitization of Lee et al. Fig. 5a (5 V and 7 V, 2 K): line centre per pixel column for |B| = 0.5-7 T."""
from PIL import Image; import numpy as np, json
a=np.asarray(Image.open('fig5a.png').convert('RGB')).astype(float)
X0,X1,Y0,Y1=68.5,1319.5,33.5,1911.5
C={'7V':(48,128,192),'5V':(32,32,128),'3V':(240,144,64)}
out={}
for v,c in C.items():
    m=np.sqrt(((a-np.array(c)-8)**2).sum(2))<30; m[int(0.75*a.shape[0]):,:]=False
    Bs=[];MR=[]
    for x in range(int(X0)+2,int(X1)-2):
        B=(x-X0)/(X1-X0)*20-10
        if not (0.5<=abs(B)<=7.0): continue
        ys=np.nonzero(m[:,x])[0]
        if len(ys)<3: continue
        # keep the run nearest the median (the curve), drop stray pixels
        y=np.median(ys)
        if np.ptp(ys)>40: continue
        Bs.append(B); MR.append(5-(y-Y0)/(Y1-Y0)*35)
    Bs=np.array(Bs); MR=np.array(MR)
    sel=np.abs(Bs)>=2
    A=np.c_[np.ones(sel.sum()),Bs[sel]**2]
    coef,res,_,_=np.linalg.lstsq(A,MR[sel],rcond=None)
    resid=MR[sel]-A@coef
    # parabola fit uncertainty
    cov=np.linalg.inv(A.T@A)*np.var(resid,ddof=2)
    d27=coef[1]*45; e27=np.sqrt(cov[1,1])*45
    asym=np.mean(MR[(Bs>=2)])-np.mean(MR[(Bs<=-2)])
    out[v]=dict(n=int(sel.sum()),b=float(coef[1]),change_2_7=float(d27),err=float(e27),rms=float(resid.std()),asym=float(asym),
               mean_2_7=float(MR[sel].mean()))
    print(v,{k:round(val,4) if isinstance(val,float) else val for k,val in out[v].items()})
json.dump(out,open('dense5a.json','w'),indent=1)
