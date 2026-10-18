"""Independent re-digitization of Lee et al. Fig. 3b: marker centroids per colour (600 dpi crop fig3b.png)."""
from PIL import Image; import numpy as np, json
from scipy import ndimage
a=np.asarray(Image.open('fig3b.png').convert('RGB')).astype(float)
TX={2:None,10:None,50:None,100:None,150:None}
x10,x100=488.5,853.0; dec=x100-x10
xT={T:x10+dec*np.log10(T/10) for T in TX}
top=dict(y4=21.0,y0=732.0)       # mobility 10^4 .. 10^0
bot=dict(y14=826.0,y12=1444.0)   # density 10^14 .. 10^12
# sample colours at known markers (x at 2 K, y read from the image)
def col_at(x,y,r=4):
    return a[int(y)-r:int(y)+r,int(x)-r:int(x)+r].reshape(-1,3).mean(0)
# reference colours: take from Fig 3a legend digitization classes, refined at markers
C={'0V':(30,50,62),'1V':(48,160,128),'2V':(224,184,72),'3V':(240,150,72),'4V':(230,80,64),'5V':(40,48,140),'6V':(40,100,160),'7V':(60,140,210),'8V':(40,150,70)}
out={}
for v,c in C.items():
    d=np.sqrt(((a-np.array(c))**2).sum(2)); m=d<45
    lab,nl=ndimage.label(m)
    res={'mu':{},'n':{}}
    for T,x in xT.items():
        for panel,(ylo,yhi) in (('mu',(25,728)),('n',(740,1440))):
            sub=lab[ylo:yhi,int(x-22):int(x+23)]
            ids,cnt=np.unique(sub[sub>0],return_counts=True)
            best=None
            for i,cn in zip(ids,cnt):
                if cn<120: continue
                ys,xs=np.nonzero(sub==i)
                # restrict to pixels within the marker window
                yc=ys.mean()+ylo
                if best is None or cn>best[1]: best=(yc,cn)
            if best is None: res[panel][T]=None; continue
            yc,cn=best
            if panel=='mu': val=10**(4-(yc-top['y4'])/(top['y0']-top['y4'])*4)
            else: val=10**(14-(yc-bot['y14'])/(bot['y12']-bot['y14'])*2)
            res[panel][T]=(round(float(val),3 if panel=='n' else 1) if panel=='mu' else float('%.4g'%val), int(cn))
    out[v]=res
    print(v,'mu',{T:r for T,r in res['mu'].items()})
    print(v,' n',{T:r for T,r in res['n'].items()})
json.dump(out,open('dig3b_v2.json','w'),indent=1)
