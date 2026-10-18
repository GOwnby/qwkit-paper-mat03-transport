"""Digitize Lee et al. Fig. 4 (mu vs n_2D at 10 K blue and 150 K green) at 600 dpi: marker centroids."""
from PIL import Image; import numpy as np, json
from scipy import ndimage
a=np.asarray(Image.open(__import__('sys').argv[1] if len(__import__('sys').argv)>1 else 'fig4.png').convert('RGB')).astype(float)
dark=a.sum(2)<200
xt=[x for x in range(30,1440) if dark[1440:1468,x].sum()>22]
yt=[y for y in range(30,1460) if dark[y,1420:1448].sum()>22]
print('x ticks',xt); print('y ticks',yt)
cb=a[260:275,890:905].reshape(-1,3).mean(0); cg=a[340:352,980:994].reshape(-1,3).mean(0)
print('blue',cb.round(),'green',cg.round())
out={}
for name,c in (('10K',cb),('150K',cg)):
    m=np.sqrt(((a-c)**2).sum(2))<45
    m[1180:,:]=False                       # legend
    lab,n=ndimage.label(m)
    pts=[]
    for i in range(1,n+1):
        ys,xs=np.nonzero(lab==i)
        if len(ys)<300: continue
        groups=[(xs,ys)]
        if len(ys)>2200:                     # two overlapping markers (4 V and 5 V at 10 K): 2-means split
            P=np.c_[xs,ys].astype(float); c=P[[P[:,0].argmin(),P[:,0].argmax()]]
            for _ in range(50):
                l=((P[:,None,:]-c[None])**2).sum(2).argmin(1); c=np.array([P[l==j].mean(0) for j in (0,1)])
            groups=[(xs[l==j],ys[l==j]) for j in (0,1)]
        for gx,gy in groups:
            x=gx.mean(); y=gy.mean()
            nn=10**(13+(x-486.5)/484); mu=10**(3-(y-384)/361.0)
            pts.append((float('%.4g'%nn),round(float(mu),1),len(gy)))
    pts.sort()
    out[name]=pts
    for p in pts: print(name,p)
json.dump(out,open('dig4.json','w'),indent=1)
